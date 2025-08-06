"""
增强的数据库管理器
支持自动扫描data/structures目录，加载所有非天然氨基酸数据
"""

import os
import json
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass

from .connection import DatabaseManager
from ..models import AminoAcidInfo, DatabaseError
from ...utils.chemistry import ChemistryUtils, MolecularFingerprint
from ...utils.config import Config


@dataclass
class DataQualityReport:
    """数据质量报告"""
    aa_id: str
    smiles_valid: bool = False
    formula_consistent: bool = False
    structure_reasonable: bool = False
    files_complete: bool = False
    overall_quality: float = 0.0
    issues: List[str] = None
    warnings: List[str] = None
    
    def __post_init__(self):
        if self.issues is None:
            self.issues = []
        if self.warnings is None:
            self.warnings = []


class EnhancedDatabaseManager(DatabaseManager):
    """增强的数据库管理器"""
    
    def __init__(self, config: Optional[Config] = None, structures_path: str = "data/structures"):
        super().__init__(config)
        self.structures_path = Path(structures_path)
        self.chemistry_utils = ChemistryUtils()
        self.fingerprint_utils = MolecularFingerprint()
        self.quality_reports = {}
        
    def load_amino_acid_database(self) -> Dict[str, AminoAcidInfo]:
        """
        扫描并加载所有氨基酸数据
        
        Returns:
            加载的氨基酸数据字典 {id: AminoAcidInfo}
        """
        print(f"🔍 扫描氨基酸数据目录: {self.structures_path}")
        
        if not self.structures_path.exists():
            raise DatabaseError(f"数据目录不存在: {self.structures_path}")
        
        amino_acids = {}
        loaded_count = 0
        error_count = 0
        
        # 扫描所有氨基酸目录
        for aa_dir in self.structures_path.iterdir():
            if aa_dir.is_dir():
                aa_id = aa_dir.name
                try:
                    amino_acid = self._load_single_amino_acid(aa_id, aa_dir)
                    if amino_acid:
                        amino_acids[aa_id] = amino_acid
                        # 保存到数据库
                        if self.add_amino_acid(amino_acid):
                            loaded_count += 1
                        else:
                            error_count += 1
                            print(f"⚠️ 保存到数据库失败: {aa_id}")
                    else:
                        error_count += 1
                        
                except Exception as e:
                    error_count += 1
                    print(f"❌ 加载氨基酸失败 {aa_id}: {e}")
        
        print(f"✅ 数据库加载完成: 成功{loaded_count}个, 失败{error_count}个")
        return amino_acids
    
    def _load_single_amino_acid(self, aa_id: str, aa_dir: Path) -> Optional[AminoAcidInfo]:
        """
        加载单个氨基酸的数据
        
        Args:
            aa_id: 氨基酸ID
            aa_dir: 氨基酸数据目录
            
        Returns:
            氨基酸信息对象
        """
        try:
            # 检查必要文件
            smi_file = aa_dir / f"{aa_id}.smi"
            pdb_file = aa_dir / f"{aa_id}.pdb"
            mol2_file = aa_dir / f"{aa_id}.mol2"
            
            # 读取SMILES
            smiles = ""
            if smi_file.exists():
                with open(smi_file, 'r') as f:
                    smiles = f.read().strip().split()[0]  # 取第一个字段
            
            # 从SMILES计算基本信息
            molecular_formula = ""
            molecular_weight = 0.0
            atom_composition = {}
            
            if smiles:
                try:
                    molecular_formula = self.fingerprint_utils.generate_molecular_formula_from_smiles(smiles)
                    if not molecular_formula:
                        # 如果RDKit不可用，尝试从PDB文件计算
                        if pdb_file.exists():
                            formula, composition = self._extract_formula_from_pdb(pdb_file)
                            molecular_formula = formula
                            atom_composition = composition
                    else:
                        # 从分子式解析原子组成
                        atom_composition = self.chemistry_utils.parse_molecular_formula(molecular_formula)
                        molecular_weight = self._calculate_molecular_weight(atom_composition)
                        
                except Exception as e:
                    print(f"⚠️ 从SMILES计算分子信息失败 {aa_id}: {e}")
            
            # 如果SMILES处理失败，尝试从PDB文件获取信息
            if not molecular_formula and pdb_file.exists():
                formula, composition = self._extract_formula_from_pdb(pdb_file)
                molecular_formula = formula
                atom_composition = composition
                molecular_weight = self._calculate_molecular_weight(atom_composition)
            
            # 计算分子指纹
            fingerprints = {}
            if smiles:
                try:
                    morgan_fp = self.fingerprint_utils.calculate_morgan_fingerprint(smiles)
                    if morgan_fp:
                        fingerprints['ecfp2'] = morgan_fp
                except Exception as e:
                    print(f"⚠️ 计算分子指纹失败 {aa_id}: {e}")
            
            # 确定氨基酸名称
            name = self._get_amino_acid_name(aa_id)
            
            # 创建氨基酸信息对象
            amino_acid = AminoAcidInfo(
                id=aa_id,
                name=name,
                molecular_formula=molecular_formula,
                molecular_weight=molecular_weight,
                smiles=smiles,
                atom_composition=atom_composition,
                key_features=self._extract_key_features(aa_id, smiles),
                fingerprints=fingerprints
            )
            
            # 生成质量报告
            quality_report = self._validate_amino_acid_data(amino_acid, aa_dir)
            self.quality_reports[aa_id] = quality_report
            
            return amino_acid
            
        except Exception as e:
            print(f"❌ 加载氨基酸数据失败 {aa_id}: {e}")
            return None
    
    def _extract_formula_from_pdb(self, pdb_file: Path) -> Tuple[str, Dict[str, int]]:
        """
        从PDB文件提取分子式和原子组成
        
        Args:
            pdb_file: PDB文件路径
            
        Returns:
            (分子式, 原子组成)
        """
        try:
            from ..parser import PDBParser
            
            parser = PDBParser()
            residues = parser.parse_pdb_file(str(pdb_file))
            
            if residues:
                residue = residues[0]
                return residue.molecular_formula, residue.atom_composition
            
        except Exception as e:
            print(f"⚠️ 从PDB提取分子式失败: {e}")
        
        return "", {}
    
    def _calculate_molecular_weight(self, atom_composition: Dict[str, int]) -> float:
        """
        计算分子量
        
        Args:
            atom_composition: 原子组成
            
        Returns:
            分子量
        """
        # 原子量表（简化版）
        atomic_weights = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
            'S': 32.065, 'P': 30.974, 'F': 18.998, 'Cl': 35.453,
            'Br': 79.904, 'I': 126.904, 'Se': 78.971
        }
        
        molecular_weight = 0.0
        for element, count in atom_composition.items():
            if element in atomic_weights:
                molecular_weight += atomic_weights[element] * count
            else:
                print(f"⚠️ 未知元素: {element}")
        
        return round(molecular_weight, 2)
    
    def _get_amino_acid_name(self, aa_id: str) -> str:
        """
        获取氨基酸名称
        
        Args:
            aa_id: 氨基酸ID
            
        Returns:
            氨基酸名称
        """
        # 标准氨基酸名称映射
        standard_names = {
            'ALA': '丙氨酸', 'ARG': '精氨酸', 'ASN': '天冬酰胺', 'ASP': '天冬氨酸',
            'CYS': '半胱氨酸', 'GLN': '谷氨酰胺', 'GLU': '谷氨酸', 'GLY': '甘氨酸',
            'HIS': '组氨酸', 'ILE': '异亮氨酸', 'LEU': '亮氨酸', 'LYS': '赖氨酸',
            'MET': '蛋氨酸', 'PHE': '苯丙氨酸', 'PRO': '脯氨酸', 'SER': '丝氨酸',
            'THR': '苏氨酸', 'TRP': '色氨酸', 'TYR': '酪氨酸', 'VAL': '缬氨酸'
        }
        
        if aa_id in standard_names:
            return standard_names[aa_id]
        else:
            return f"氨基酸-{aa_id}"
    
    def _extract_key_features(self, aa_id: str, smiles: str) -> List[str]:
        """
        提取氨基酸的关键特征
        
        Args:
            aa_id: 氨基酸ID
            smiles: SMILES字符串
            
        Returns:
            关键特征列表
        """
        features = []
        
        # 基于ID的特征
        if aa_id in ['ALA', 'GLY', 'VAL', 'LEU', 'ILE']:
            features.append('疏水性')
        elif aa_id in ['SER', 'THR', 'ASN', 'GLN']:
            features.append('极性')
        elif aa_id in ['ASP', 'GLU']:
            features.append('酸性')
        elif aa_id in ['LYS', 'ARG', 'HIS']:
            features.append('碱性')
        elif aa_id in ['PHE', 'TYR', 'TRP']:
            features.append('芳香性')
        elif aa_id in ['CYS', 'MET']:
            features.append('含硫')
        
        # 基于SMILES的特征
        if smiles:
            if 'c1ccccc1' in smiles:
                features.append('苯环')
            if 'S' in smiles:
                features.append('含硫')
            if 'N' in smiles and smiles.count('N') > 1:
                features.append('多氮')
            if 'F' in smiles:
                features.append('含氟')
            if 'Cl' in smiles:
                features.append('含氯')
        
        # 非天然氨基酸标记
        standard_aas = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                       'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                       'THR', 'TRP', 'TYR', 'VAL'}
        
        if aa_id not in standard_aas:
            features.append('非天然氨基酸')
        
        return features
    
    def _validate_amino_acid_data(self, amino_acid: AminoAcidInfo, aa_dir: Path) -> DataQualityReport:
        """
        验证氨基酸数据质量
        
        Args:
            amino_acid: 氨基酸信息
            aa_dir: 氨基酸数据目录
            
        Returns:
            数据质量报告
        """
        report = DataQualityReport(aa_id=amino_acid.id)
        
        # 检查SMILES有效性
        if amino_acid.smiles:
            report.smiles_valid = len(amino_acid.smiles) > 0 and not amino_acid.smiles.isspace()
            if not report.smiles_valid:
                report.issues.append("SMILES无效")
        else:
            report.warnings.append("缺少SMILES")
        
        # 检查分子式一致性
        if amino_acid.molecular_formula and amino_acid.smiles:
            try:
                smiles_formula = self.fingerprint_utils.generate_molecular_formula_from_smiles(amino_acid.smiles)
                if smiles_formula:
                    report.formula_consistent = (smiles_formula == amino_acid.molecular_formula)
                    if not report.formula_consistent:
                        report.issues.append(f"分子式不一致: SMILES={smiles_formula}, 数据库={amino_acid.molecular_formula}")
                else:
                    report.warnings.append("无法从SMILES计算分子式")
            except Exception:
                report.warnings.append("SMILES分子式计算失败")
        
        # 检查文件完整性
        required_files = [f"{amino_acid.id}.smi", f"{amino_acid.id}.pdb"]
        optional_files = [f"{amino_acid.id}.mol2", f"{amino_acid.id}.hdb", f"{amino_acid.id}.top"]
        
        missing_required = []
        missing_optional = []
        
        for filename in required_files:
            if not (aa_dir / filename).exists():
                missing_required.append(filename)
        
        for filename in optional_files:
            if not (aa_dir / filename).exists():
                missing_optional.append(filename)
        
        report.files_complete = len(missing_required) == 0
        if missing_required:
            report.issues.append(f"缺少必要文件: {', '.join(missing_required)}")
        if missing_optional:
            report.warnings.append(f"缺少可选文件: {', '.join(missing_optional)}")
        
        # 计算总体质量分数
        quality_score = 0.0
        if report.smiles_valid:
            quality_score += 0.3
        if report.formula_consistent:
            quality_score += 0.3
        if report.files_complete:
            quality_score += 0.2
        if len(report.issues) == 0:
            quality_score += 0.2
        
        report.overall_quality = quality_score
        
        return report
    
    def build_search_indices(self):
        """构建搜索索引以提高查询效率"""
        print("🔍 构建搜索索引...")
        
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                
                # 创建额外的索引
                indices = [
                    'CREATE INDEX IF NOT EXISTS idx_molecular_weight ON amino_acids(molecular_weight)',
                    'CREATE INDEX IF NOT EXISTS idx_smiles ON amino_acids(smiles)',
                    'CREATE INDEX IF NOT EXISTS idx_key_features ON amino_acids(key_features)',
                ]
                
                for index_sql in indices:
                    cursor.execute(index_sql)
                
                conn.commit()
                print("✅ 搜索索引构建完成")
                
        except Exception as e:
            print(f"⚠️ 搜索索引构建失败: {e}")
    
    def get_quality_report(self, aa_id: str = None) -> Dict[str, DataQualityReport]:
        """
        获取数据质量报告
        
        Args:
            aa_id: 特定氨基酸ID，None表示获取所有报告
            
        Returns:
            质量报告字典
        """
        if aa_id:
            return {aa_id: self.quality_reports.get(aa_id)}
        else:
            return self.quality_reports.copy()
    
    def get_enhanced_database_stats(self) -> Dict[str, any]:
        """获取增强的数据库统计信息"""
        basic_stats = self.get_database_stats()
        
        # 添加质量统计
        quality_stats = {
            'high_quality': 0,    # 质量分数 >= 0.8
            'medium_quality': 0,  # 质量分数 0.5-0.8
            'low_quality': 0,     # 质量分数 < 0.5
            'with_issues': 0,     # 有问题的数据
            'with_warnings': 0,   # 有警告的数据
        }
        
        for report in self.quality_reports.values():
            if report.overall_quality >= 0.8:
                quality_stats['high_quality'] += 1
            elif report.overall_quality >= 0.5:
                quality_stats['medium_quality'] += 1
            else:
                quality_stats['low_quality'] += 1
            
            if report.issues:
                quality_stats['with_issues'] += 1
            if report.warnings:
                quality_stats['with_warnings'] += 1
        
        return {**basic_stats, **quality_stats}