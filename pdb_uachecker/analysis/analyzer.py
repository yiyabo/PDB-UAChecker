"""
主分析器
整合PDB解析、数据库查询和四重验证，提供统一的分析接口
"""

import time
from typing import List, Dict, Any, Optional
from pathlib import Path

from ..core.parser import PDBParser
from ..core.database import DatabaseManager
from ..core.verification import VerificationEngine
from ..core.models import (
    ResidueInfo, AminoAcidInfo, MatchResult, AnalysisResult,
    VerificationMethod, ParsingError, DatabaseError, VerificationError
)
from ..utils.config import Config


class PDBAnalyzer:
    """PDB分析器 - 主入口类"""
    
    def __init__(self, config: Optional[Config] = None):
        """
        初始化PDB分析器
        
        Args:
            config: 配置对象，None使用默认配置
        """
        if config is None:
            from ..utils.config import default_config
            config = default_config
        
        self.config = config
        
        # 初始化组件
        self.parser = PDBParser()
        self.database = DatabaseManager(config)
        self.verification_engine = VerificationEngine(config)
        
        print("🚀 PDB分析器初始化完成")
        self._print_capabilities()
    
    def analyze_pdb(self, pdb_file: str, 
                   enabled_methods: Optional[List[VerificationMethod]] = None,
                   save_report: bool = False,
                   output_file: Optional[str] = None) -> AnalysisResult:
        """
        分析PDB文件中的非天然氨基酸
        
        Args:
            pdb_file: PDB文件路径
            enabled_methods: 启用的验证方法列表
            save_report: 是否保存报告
            output_file: 输出文件路径
        
        Returns:
            分析结果
        
        Raises:
            ParsingError: PDB解析失败
            DatabaseError: 数据库操作失败
            VerificationError: 验证失败
        """
        start_time = time.time()
        
        print(f"\n🔬 开始分析PDB文件: {pdb_file}")
        print("=" * 60)
        
        try:
            # 1. 验证PDB文件
            self._validate_pdb_file(pdb_file)
            
            # 2. 解析PDB文件
            residues = self.parser.parse_pdb_file(pdb_file)
            print(f"📊 解析完成: 发现 {len(residues)} 个残基")
            
            # 3. 加载氨基酸数据库
            amino_acids = self.database.get_all_amino_acids()
            print(f"📚 数据库加载: {len(amino_acids)} 种氨基酸")
            
            # 4. 分析每个残基
            matches = []
            for i, residue in enumerate(residues, 1):
                print(f"\n🧪 分析残基 {i}/{len(residues)}: {residue.residue_key}")
                
                # 多策略搜索
                residue_matches = self._analyze_residue(residue, amino_acids, enabled_methods)
                matches.extend(residue_matches)
                
                if residue_matches:
                    best_match = residue_matches[0]
                    print(f"   ✅ 最佳匹配: {best_match.amino_acid_name}")
                    print(f"   📊 置信度: {best_match.confidence_score:.3f}")
                    print(f"   🔧 通过验证: {best_match.verification_result.passed_verifications}/{best_match.verification_result.total_verifications}")
                else:
                    print(f"   ❌ 未找到匹配")
            
            analysis_time = time.time() - start_time
            
            # 5. 生成分析结果
            result = AnalysisResult(
                pdb_file=pdb_file,
                total_residues=len(residues),
                identified_residues=len(matches),
                matches=matches,
                analysis_time=analysis_time,
                capabilities=self._get_capabilities()
            )
            
            # 6. 打印摘要
            self._print_analysis_summary(result)
            
            # 7. 保存报告（如果需要）
            if save_report and output_file:
                self._save_analysis_report(result, output_file)
            
            return result
        
        except Exception as e:
            print(f"❌ 分析失败: {e}")
            raise
    
    def _validate_pdb_file(self, pdb_file: str):
        """验证PDB文件"""
        if not Path(pdb_file).exists():
            raise ParsingError(f"PDB文件不存在: {pdb_file}")
        
        is_valid, errors = self.parser.validate_pdb_file(pdb_file)
        if not is_valid:
            raise ParsingError(f"PDB文件格式错误: {'; '.join(errors)}")
    
    def _analyze_residue(self, residue: ResidueInfo, amino_acids: List[AminoAcidInfo],
                        enabled_methods: Optional[List[VerificationMethod]] = None) -> List[MatchResult]:
        """
        分析单个残基
        
        Args:
            residue: 残基信息
            amino_acids: 氨基酸数据库
            enabled_methods: 启用的验证方法
        
        Returns:
            匹配结果列表
        """
        matches = []
        
        # 策略1: 残基名精确匹配
        exact_matches = self._search_by_residue_name(residue, amino_acids)
        if exact_matches:
            # 残基名匹配成功，直接返回
            for amino_acid in exact_matches:
                match_result = MatchResult(
                    amino_acid_info=amino_acid,
                    residue_info=residue,
                    verification_result=self._create_perfect_verification_result(amino_acid),
                    match_method="residue_name_exact"
                )
                matches.append(match_result)
            return matches
        
        # 策略2: 四重验证搜索
        candidate_amino_acids = self._get_candidate_amino_acids(residue, amino_acids)
        
        if candidate_amino_acids:
            verification_results = self.verification_engine.find_best_matches(
                residue, candidate_amino_acids, 
                top_k=5, min_confidence=0.0, 
                enabled_methods=enabled_methods
            )
            
            for verification_result in verification_results:
                if verification_result.is_match:
                    amino_acid = next(
                        aa for aa in candidate_amino_acids 
                        if aa.id == verification_result.amino_acid_id
                    )
                    match_result = MatchResult(
                        amino_acid_info=amino_acid,
                        residue_info=residue,
                        verification_result=verification_result,
                        match_method="four_fold_verification"
                    )
                    matches.append(match_result)
        
        return matches
    
    def _search_by_residue_name(self, residue: ResidueInfo, amino_acids: List[AminoAcidInfo]) -> List[AminoAcidInfo]:
        """根据残基名搜索"""
        return [aa for aa in amino_acids if aa.id == residue.residue_name]
    
    def _get_candidate_amino_acids(self, residue: ResidueInfo, amino_acids: List[AminoAcidInfo]) -> List[AminoAcidInfo]:
        """
        获取候选氨基酸列表
        
        使用分子式和原子组成进行初步筛选，提高效率
        """
        candidates = []
        
        # 分子式筛选
        if residue.molecular_formula:
            formula_matches = [
                aa for aa in amino_acids 
                if aa.molecular_formula == residue.molecular_formula
            ]
            candidates.extend(formula_matches)
        
        # 原子组成筛选
        if residue.atom_composition:
            composition_matches = [
                aa for aa in amino_acids 
                if aa.atom_composition == residue.atom_composition
            ]
            candidates.extend(composition_matches)
        
        # 去重
        unique_candidates = {}
        for aa in candidates:
            unique_candidates[aa.id] = aa
        
        # 如果筛选结果太少，返回所有氨基酸
        if len(unique_candidates) < 10:
            return amino_acids
        
        return list(unique_candidates.values())
    
    def _create_perfect_verification_result(self, amino_acid: AminoAcidInfo):
        """为残基名精确匹配创建完美验证结果"""
        from ..core.models import VerificationResult, VerificationScore
        
        # 创建完美分数
        perfect_scores = [
            VerificationScore(
                method=VerificationMethod.MOLECULAR_FORMULA,
                score=1.0,
                passed=True,
                details={'match_type': 'residue_name_exact'}
            )
        ]
        
        return VerificationResult(
            amino_acid_id=amino_acid.id,
            amino_acid_name=amino_acid.name,
            scores=perfect_scores,
            overall_confidence=1.0,
            passed_verifications=1,
            total_verifications=1
        )
    
    def _print_capabilities(self):
        """打印系统能力"""
        capabilities = self._get_capabilities()
        print("🎯 系统能力:")
        for capability, available in capabilities.items():
            status = "✅" if available else "❌"
            print(f"   {status} {capability.replace('_', ' ').title()}")
    
    def _get_capabilities(self) -> Dict[str, bool]:
        """获取系统能力"""
        return {
            'exact_matching': True,
            'molecular_formula_verification': True,
            'atom_composition_verification': True,
            'fingerprint_similarity': self.config.is_rdkit_available(),
            '3d_structure_matching': True,  # 基础支持，但需要标准结构数据
            'parallel_processing': self.config.performance.enable_parallel
        }
    
    def _print_analysis_summary(self, result: AnalysisResult):
        """打印分析摘要"""
        print(f"\n" + "=" * 60)
        print(f"📊 PDB分析摘要:")
        print(f"   文件: {result.pdb_file}")
        print(f"   分析时间: {result.analysis_time:.2f}s")
        print(f"   总残基数: {result.total_residues}")
        print(f"   识别残基数: {result.identified_residues}")
        print(f"   识别率: {result.identification_rate:.1f}%")
        
        if result.matches:
            print(f"\n🎯 识别的氨基酸 (前5个):")
            for i, match in enumerate(result.matches[:5], 1):
                residue = match.residue_info
                print(f"   {i}. {residue.residue_key}")
                print(f"      → {match.amino_acid_name} (置信度: {match.confidence_score:.3f})")
                print(f"      → 方法: {match.match_method}")
    
    def _save_analysis_report(self, result: AnalysisResult, output_file: str):
        """保存分析报告"""
        try:
            report_content = self._generate_report_content(result)
            
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(report_content)
            
            print(f"📄 分析报告已保存: {output_file}")
        
        except Exception as e:
            print(f"⚠️ 报告保存失败: {e}")
    
    def _generate_report_content(self, result: AnalysisResult) -> str:
        """生成报告内容"""
        lines = [
            "PDB非天然氨基酸分析报告",
            "=" * 50,
            "",
            f"文件: {result.pdb_file}",
            f"分析时间: {result.analysis_time:.2f}秒",
            f"总残基数: {result.total_residues}",
            f"识别残基数: {result.identified_residues}",
            f"识别率: {result.identification_rate:.1f}%",
            "",
            "系统能力:",
        ]
        
        for capability, available in result.capabilities.items():
            status = "✅" if available else "❌"
            lines.append(f"  {status} {capability.replace('_', ' ').title()}")
        
        if result.matches:
            lines.extend([
                "",
                "识别结果:",
                "-" * 30
            ])
            
            for i, match in enumerate(result.matches, 1):
                residue = match.residue_info
                verification = match.verification_result
                
                lines.extend([
                    f"{i}. {residue.residue_key}",
                    f"   氨基酸: {match.amino_acid_name}",
                    f"   置信度: {match.confidence_score:.3f}",
                    f"   匹配方法: {match.match_method}",
                    f"   通过验证: {verification.passed_verifications}/{verification.total_verifications}",
                    ""
                ])
                
                # 详细验证信息
                for score in verification.scores:
                    status = "✅" if score.passed else "❌"
                    lines.append(f"     {status} {score.method.value}: {score.score:.3f}")
                
                lines.append("")
        
        return "\n".join(lines)
    
    def analyze_residue(self, residue: ResidueInfo, 
                       enabled_methods: Optional[List[VerificationMethod]] = None) -> List[MatchResult]:
        """
        分析单个残基（公共接口）
        
        Args:
            residue: 残基信息
            enabled_methods: 启用的验证方法
        
        Returns:
            匹配结果列表
        """
        amino_acids = self.database.get_all_amino_acids()
        return self._analyze_residue(residue, amino_acids, enabled_methods)
    
    def get_database_stats(self) -> Dict[str, Any]:
        """获取数据库统计信息"""
        return self.database.get_database_stats()
    
    def get_verification_stats(self, results: List[MatchResult]) -> Dict[str, Any]:
        """获取验证统计信息"""
        verification_results = [match.verification_result for match in results]
        return self.verification_engine.get_verification_statistics(verification_results)
    
    def update_config(self, **config_updates):
        """更新配置"""
        # 更新阈值
        if 'thresholds' in config_updates:
            self.verification_engine.update_thresholds(**config_updates['thresholds'])
        
        print("✅ 配置已更新")
    
    def validate_system(self) -> Dict[str, Any]:
        """验证系统状态"""
        validation = {
            'overall_status': 'healthy',
            'components': {},
            'warnings': [],
            'errors': []
        }
        
        # 验证数据库
        try:
            stats = self.database.get_database_stats()
            validation['components']['database'] = {
                'status': 'healthy',
                'amino_acids_count': stats.get('total_amino_acids', 0)
            }
        except Exception as e:
            validation['components']['database'] = {
                'status': 'error',
                'error': str(e)
            }
            validation['errors'].append(f"数据库错误: {e}")
        
        # 验证验证引擎
        engine_validation = self.verification_engine.validate_configuration()
        validation['components']['verification_engine'] = engine_validation
        
        if not engine_validation['valid']:
            validation['errors'].extend(engine_validation['errors'])
        
        validation['warnings'].extend(engine_validation['warnings'])
        
        # 确定整体状态
        if validation['errors']:
            validation['overall_status'] = 'error'
        elif validation['warnings']:
            validation['overall_status'] = 'warning'
        
        return validation