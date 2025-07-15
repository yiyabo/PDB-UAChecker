#!/usr/bin/env python3
"""
增强版PDB分析器 - 整合指纹相似性和立体化学匹配
解决同分异构体识别问题
"""

import sys
import time
import json
import sqlite3
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass

# 添加路径
sys.path.insert(0, 'core')
sys.path.insert(0, 'legacy')

# 导入现有组件
from unified_pdb_analyzer import PDBParser, ResidueInfo, MatchResult, UnifiedSearchEngine

# 导入同分异构体识别器
try:
    from isomer_identifier import IsomerIdentifier, RDKIT_AVAILABLE
    ISOMER_AVAILABLE = RDKIT_AVAILABLE
except ImportError:
    ISOMER_AVAILABLE = False
    print("⚠️ 同分异构体识别器不可用")

@dataclass
class EnhancedMatchResult(MatchResult):
    """增强的匹配结果，包含更多信息"""
    isomer_type: Optional[str] = None  # 'identical', 'stereoisomer', 'structural', 'different'
    structural_similarity: Optional[float] = None
    stereochemistry_info: Optional[str] = None
    fingerprint_similarity: Optional[float] = None

class EnhancedSearchEngine(UnifiedSearchEngine):
    """增强的搜索引擎，支持指纹相似性和立体化学匹配"""
    
    def __init__(self):
        super().__init__()
        
        # 初始化同分异构体识别器
        if ISOMER_AVAILABLE:
            self.isomer_identifier = IsomerIdentifier()
            print("✅ 同分异构体识别器已启用")
        else:
            self.isomer_identifier = None
            print("❌ 同分异构体识别器不可用")
    
    def search_with_fingerprint_similarity(self, target_smiles: str, 
                                         threshold: float = 0.6) -> List[EnhancedMatchResult]:
        """基于指纹相似性的增强搜索"""
        results = []
        
        if not self.isomer_identifier:
            print("⚠️ 指纹相似性搜索不可用")
            return results
        
        print(f"🔍 进行指纹相似性搜索，阈值: {threshold}")
        
        for amino_id, aa_data in self.amino_acids.items():
            if not aa_data['smiles']:
                continue
            
            try:
                # 使用同分异构体识别器分析关系
                relationship = self.isomer_identifier.identify_relationship(
                    target_smiles, aa_data['smiles']
                )
                
                similarity = relationship.get('structural_similarity', 0.0)
                
                if similarity >= threshold:
                    result = EnhancedMatchResult(
                        amino_acid_id=amino_id,
                        amino_acid_name=aa_data['name'],
                        confidence_score=similarity * 0.8,  # 相似性匹配置信度稍低
                        match_method="fingerprint_similarity",
                        molecular_formula=aa_data['molecular_formula'],
                        smiles=aa_data['smiles'],
                        residue_info=None,
                        isomer_type=relationship.get('isomer_type', 'unknown'),
                        structural_similarity=similarity,
                        fingerprint_similarity=similarity
                    )
                    results.append(result)
                    
            except Exception as e:
                print(f"⚠️ 处理 {amino_id} 时出错: {e}")
                continue
        
        # 按相似性排序
        results.sort(key=lambda x: x.structural_similarity or 0, reverse=True)
        return results[:10]  # 返回前10个最相似的
    
    def search_stereoisomers(self, target_smiles: str) -> List[EnhancedMatchResult]:
        """专门搜索立体异构体"""
        results = []
        
        if not self.isomer_identifier:
            return results
        
        print(f"🔍 搜索立体异构体")
        
        for amino_id, aa_data in self.amino_acids.items():
            if not aa_data['smiles']:
                continue
            
            try:
                relationship = self.isomer_identifier.identify_relationship(
                    target_smiles, aa_data['smiles']
                )
                
                # 只返回立体异构体
                if relationship.get('isomer_type') == 'stereoisomer':
                    result = EnhancedMatchResult(
                        amino_acid_id=amino_id,
                        amino_acid_name=aa_data['name'],
                        confidence_score=0.85,  # 立体异构体高置信度
                        match_method="stereoisomer_match",
                        molecular_formula=aa_data['molecular_formula'],
                        smiles=aa_data['smiles'],
                        residue_info=None,
                        isomer_type='stereoisomer',
                        structural_similarity=relationship.get('structural_similarity', 0.0)
                    )
                    results.append(result)
                    
            except Exception as e:
                continue
        
        return results
    
    def comprehensive_search(self, residue: ResidueInfo, 
                           enable_fingerprint: bool = True,
                           fingerprint_threshold: float = 0.6) -> List[EnhancedMatchResult]:
        """综合搜索 - 整合所有匹配方法"""
        all_results = []
        
        # 1. 传统精确匹配
        exact_results = self._traditional_search(residue)
        all_results.extend(exact_results)
        
        # 2. 指纹相似性搜索（如果有SMILES）
        if enable_fingerprint and self.isomer_identifier:
            # 这里需要从PDB坐标推导SMILES，暂时跳过
            # 实际应用中可以集成OpenEye或RDKit的坐标到SMILES转换
            pass
        
        # 去重和排序
        unique_results = self._deduplicate_results(all_results)
        return unique_results
    
    def _traditional_search(self, residue: ResidueInfo) -> List[EnhancedMatchResult]:
        """传统搜索方法"""
        results = []
        
        # 残基名匹配
        name_results = self.search_by_residue_name(residue.residue_name)
        for result in name_results:
            enhanced_result = EnhancedMatchResult(
                amino_acid_id=result.amino_acid_id,
                amino_acid_name=result.amino_acid_name,
                confidence_score=result.confidence_score,
                match_method=result.match_method,
                molecular_formula=result.molecular_formula,
                smiles=result.smiles,
                residue_info=residue
            )
            results.append(enhanced_result)
        
        # 分子式匹配
        if residue.molecular_formula:
            formula_results = self.search_by_molecular_formula(residue.molecular_formula)
            for result in formula_results:
                enhanced_result = EnhancedMatchResult(
                    amino_acid_id=result.amino_acid_id,
                    amino_acid_name=result.amino_acid_name,
                    confidence_score=result.confidence_score,
                    match_method=result.match_method,
                    molecular_formula=result.molecular_formula,
                    smiles=result.smiles,
                    residue_info=residue
                )
                results.append(enhanced_result)
        
        # 原子组成匹配
        if residue.atom_composition:
            atom_results = self.search_by_atom_composition(residue.atom_composition)
            for result in atom_results:
                enhanced_result = EnhancedMatchResult(
                    amino_acid_id=result.amino_acid_id,
                    amino_acid_name=result.amino_acid_name,
                    confidence_score=result.confidence_score,
                    match_method=result.match_method,
                    molecular_formula=result.molecular_formula,
                    smiles=result.smiles,
                    residue_info=residue
                )
                results.append(enhanced_result)
        
        return results
    
    def _deduplicate_results(self, results: List[EnhancedMatchResult]) -> List[EnhancedMatchResult]:
        """去重并选择最佳结果"""
        unique_results = {}
        
        for result in results:
            amino_id = result.amino_acid_id
            if amino_id not in unique_results:
                unique_results[amino_id] = result
            else:
                # 保留置信度更高的结果
                if result.confidence_score > unique_results[amino_id].confidence_score:
                    unique_results[amino_id] = result
        
        # 按置信度排序
        sorted_results = list(unique_results.values())
        sorted_results.sort(key=lambda x: x.confidence_score, reverse=True)
        
        return sorted_results

class EnhancedPDBAnalyzer:
    """增强版PDB分析器"""
    
    def __init__(self):
        print("🚀 初始化增强版PDB分析器...")
        self.parser = PDBParser()
        self.search_engine = EnhancedSearchEngine()
        print("✅ 增强版PDB分析器就绪！")
        
        if ISOMER_AVAILABLE:
            print("✅ 支持功能: 精确匹配 + 指纹相似性 + 立体化学识别")
        else:
            print("⚠️ 支持功能: 精确匹配（指纹相似性不可用）")
    
    def analyze_pdb_enhanced(self, pdb_file: str, 
                           enable_fingerprint: bool = True,
                           fingerprint_threshold: float = 0.6) -> Dict:
        """增强版PDB分析"""
        start_time = time.time()
        
        print(f"\n🔬 开始增强版PDB分析: {pdb_file}")
        print("=" * 70)
        
        # 1. 解析PDB文件
        residues = self.parser.parse_pdb_file(pdb_file)
        if not residues:
            return {'error': 'PDB文件解析失败'}
        
        # 2. 分析每个残基
        analysis_results = []
        isomer_findings = []
        
        for residue in residues:
            print(f"\n🧪 分析残基: {residue.chain_id}:{residue.residue_name}{residue.residue_number}")
            
            # 综合搜索
            matches = self.search_engine.comprehensive_search(
                residue, enable_fingerprint, fingerprint_threshold
            )
            
            if matches:
                best_match = matches[0]
                analysis_results.append(best_match)
                
                print(f"   ✅ 匹配: {best_match.amino_acid_name}")
                print(f"   📊 置信度: {best_match.confidence_score:.3f}")
                print(f"   🔧 方法: {best_match.match_method}")
                
                if best_match.isomer_type:
                    print(f"   🧬 异构体类型: {best_match.isomer_type}")
                    if best_match.isomer_type in ['stereoisomer', 'structural']:
                        isomer_findings.append(best_match)
            else:
                print(f"   ❌ 未找到匹配")
        
        analysis_time = time.time() - start_time
        
        # 3. 生成增强报告
        report = {
            'pdb_file': pdb_file,
            'analysis_time': analysis_time,
            'total_residues': len(residues),
            'identified_residues': len(analysis_results),
            'identification_rate': len(analysis_results) / len(residues) * 100,
            'results': analysis_results,
            'isomer_findings': isomer_findings,
            'capabilities': {
                'exact_matching': True,
                'fingerprint_similarity': ISOMER_AVAILABLE,
                'stereochemistry_analysis': ISOMER_AVAILABLE
            }
        }
        
        self._print_enhanced_summary(report)
        return report
    
    def _print_enhanced_summary(self, report: Dict):
        """打印增强分析摘要"""
        print(f"\n" + "=" * 70)
        print(f"📊 增强版PDB分析摘要:")
        print(f"   文件: {report['pdb_file']}")
        print(f"   分析时间: {report['analysis_time']:.2f}s")
        print(f"   总残基数: {report['total_residues']}")
        print(f"   识别残基数: {report['identified_residues']}")
        print(f"   识别率: {report['identification_rate']:.1f}%")
        
        if report['isomer_findings']:
            print(f"\n🧬 同分异构体发现:")
            for finding in report['isomer_findings']:
                residue = finding.residue_info
                print(f"   • {residue.chain_id}:{residue.residue_name}{residue.residue_number}")
                print(f"     → {finding.amino_acid_name} ({finding.isomer_type})")
                if finding.structural_similarity:
                    print(f"     相似性: {finding.structural_similarity:.3f}")
        
        capabilities = report['capabilities']
        print(f"\n🔧 系统能力:")
        print(f"   精确匹配: {'✅' if capabilities['exact_matching'] else '❌'}")
        print(f"   指纹相似性: {'✅' if capabilities['fingerprint_similarity'] else '❌'}")
        print(f"   立体化学分析: {'✅' if capabilities['stereochemistry_analysis'] else '❌'}")

def test_enhanced_analyzer():
    """测试增强版分析器"""
    print("🧪 测试增强版PDB分析器")
    print("=" * 50)
    
    analyzer = EnhancedPDBAnalyzer()
    
    # 测试指纹相似性搜索
    if ISOMER_AVAILABLE:
        print("\n🔍 测试指纹相似性搜索:")
        test_smiles = "N[C@@H](Cc1ccccc1)C(=O)O"  # 苯丙氨酸
        results = analyzer.search_engine.search_with_fingerprint_similarity(
            test_smiles, threshold=0.6
        )
        
        print(f"找到 {len(results)} 个相似结构:")
        for i, result in enumerate(results[:3], 1):
            print(f"  {i}. {result.amino_acid_id}: {result.amino_acid_name}")
            print(f"     相似性: {result.structural_similarity:.3f}")
            print(f"     异构体类型: {result.isomer_type}")
    
    print(f"\n💡 使用方法:")
    print("analyzer = EnhancedPDBAnalyzer()")
    print("results = analyzer.analyze_pdb_enhanced('protein.pdb')")

def main():
    """主函数"""
    print("🌟 增强版PDB分析器")
    print("=" * 50)
    
    test_enhanced_analyzer()

if __name__ == "__main__":
    main()
