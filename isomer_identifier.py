#!/usr/bin/env python3
"""
同分异构体识别模块 - 基于RDKit的ECFP分子指纹
解决同分异构体识别问题，提升搜索准确性
"""

import hashlib
import numpy as np
from typing import Dict, List, Tuple, Optional, Set, Any
from dataclasses import dataclass
from collections import defaultdict
import json

# 尝试导入RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, Descriptors, AllChem
    from rdkit import DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("警告: RDKit未安装，将使用简化的同分异构体识别")

@dataclass
class StructuralFingerprint:
    """结构指纹数据类"""
    ecfp2: Optional[List[int]] = None      # ECFP2指纹
    ecfp4: Optional[List[int]] = None      # ECFP4指纹
    ecfp6: Optional[List[int]] = None      # ECFP6指纹
    morgan2: Optional[List[int]] = None    # Morgan2指纹
    morgan3: Optional[List[int]] = None    # Morgan3指纹
    topological_hash: Optional[str] = None # 拓扑哈希
    structural_keys: Optional[List[str]] = None  # 结构关键词
    stereochemistry: Optional[str] = None  # 立体化学信息
    
    def to_dict(self) -> Dict:
        """转换为字典格式"""
        return {
            'ecfp2': self.ecfp2,
            'ecfp4': self.ecfp4,
            'ecfp6': self.ecfp6,
            'morgan2': self.morgan2,
            'morgan3': self.morgan3,
            'topological_hash': self.topological_hash,
            'structural_keys': self.structural_keys,
            'stereochemistry': self.stereochemistry
        }
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'StructuralFingerprint':
        """从字典创建对象"""
        return cls(**data)

@dataclass
class IsomerAnalysisResult:
    """同分异构体分析结果"""
    is_isomer: bool
    isomer_type: str  # 'structural', 'stereoisomer', 'identical', 'different'
    similarity_score: float
    structural_differences: List[str]
    confidence: float
    details: Dict[str, Any]

class ECFPGenerator:
    """基于RDKit的ECFP指纹生成器"""
    
    def __init__(self, fingerprint_size: int = 2048):
        self.fingerprint_size = fingerprint_size
        self.rdkit_available = RDKIT_AVAILABLE
        
        # 缓存分子对象
        self.mol_cache = {}
        self.cache_size = 1000
        
        # 预定义的氨基酸特征模式
        self.amino_acid_patterns = {
            'alpha_amino': '[NX3,NX4+][CX4H1,CX4H2]([*])[CX3](=[OX1])[OX2H,OX1-]',
            'beta_amino': '[NX3,NX4+][CX4H1,CX4H2]([*])[CX4H1,CX4H2]([*])[CX3](=[OX1])[OX2H,OX1-]',
            'aromatic_side': '[cX3]1[cX3H][cX3H][cX3H][cX3H][cX3H]1',
            'hydroxyl': '[OX2H]',
            'amino_group': '[NX3,NX4+]',
            'carboxyl': '[CX3](=[OX1])[OX2H,OX1-]',
            'sulfur': '[SX2H,SX1-]',
            'imidazole': '[nX3H][cX3H][nX3H]',
            'guanidino': '[NX3H2][CX3](=[NX3H])[NX3H2]',
            'indole': '[nX3H][cX3]1[cX3H][cX3H][cX3H][cX3H][cX3]1'
        }
        
        if self.rdkit_available:
            self.pattern_mols = {
                name: Chem.MolFromSmarts(pattern) 
                for name, pattern in self.amino_acid_patterns.items()
            }
    
    def generate_structural_fingerprint(self, smiles: str) -> StructuralFingerprint:
        """生成完整的结构指纹"""
        if not self.rdkit_available:
            return self._generate_simplified_fingerprint(smiles)
        
        try:
            mol = self._get_mol(smiles)
            if mol is None:
                return StructuralFingerprint()
            
            # 生成多种类型的指纹
            ecfp2 = self._generate_ecfp(mol, radius=1)
            ecfp4 = self._generate_ecfp(mol, radius=2)
            ecfp6 = self._generate_ecfp(mol, radius=3)
            
            morgan2 = self._generate_morgan(mol, radius=2)
            morgan3 = self._generate_morgan(mol, radius=3)
            
            topological_hash = self._generate_topological_hash(mol)
            structural_keys = self._extract_structural_keys(mol)
            stereochemistry = self._analyze_stereochemistry(mol)
            
            return StructuralFingerprint(
                ecfp2=ecfp2,
                ecfp4=ecfp4,
                ecfp6=ecfp6,
                morgan2=morgan2,
                morgan3=morgan3,
                topological_hash=topological_hash,
                structural_keys=structural_keys,
                stereochemistry=stereochemistry
            )
            
        except Exception as e:
            print(f"生成结构指纹失败: {e}")
            return StructuralFingerprint()
    
    def _get_mol(self, smiles: str) -> Optional[Chem.Mol]:
        """获取分子对象（带缓存）"""
        if smiles in self.mol_cache:
            return self.mol_cache[smiles]
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # 标准化分子
            mol = Chem.AddHs(mol)
            
            # 缓存管理
            if len(self.mol_cache) >= self.cache_size:
                # 移除最旧的缓存
                oldest_key = next(iter(self.mol_cache))
                del self.mol_cache[oldest_key]
            
            self.mol_cache[smiles] = mol
        
        return mol
    
    def _generate_ecfp(self, mol: Chem.Mol, radius: int = 2) -> List[int]:
        """生成ECFP指纹"""
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, 
            radius=radius, 
            nBits=self.fingerprint_size,
            useChirality=True,  # 重要：考虑手性
            useBondTypes=True   # 重要：考虑键类型
        )
        return [int(bit) for bit in fp.ToBitString()]
    
    def _generate_morgan(self, mol: Chem.Mol, radius: int = 2) -> List[int]:
        """生成Morgan指纹"""
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, 
            radius=radius, 
            nBits=self.fingerprint_size,
            useFeatures=True,
            useChirality=True  # 重要：考虑手性
        )
        return [int(bit) for bit in fp.ToBitString()]
    
    def _generate_topological_hash(self, mol: Chem.Mol) -> str:
        """生成拓扑哈希"""
        # 使用分子的拓扑不变量
        invariants = []
        
        # 原子不变量
        for atom in mol.GetAtoms():
            atom_invariant = (
                atom.GetAtomicNum(),
                atom.GetDegree(),
                atom.GetTotalNumHs(),
                atom.GetFormalCharge(),
                int(atom.GetIsAromatic())
            )
            invariants.append(atom_invariant)
        
        # 键不变量
        for bond in mol.GetBonds():
            bond_invariant = (
                bond.GetBondTypeAsDouble(),
                int(bond.GetIsAromatic()),
                int(bond.IsInRing())
            )
            invariants.append(bond_invariant)
        
        # 生成哈希
        invariants_str = str(sorted(invariants))
        return hashlib.md5(invariants_str.encode()).hexdigest()
    
    def _extract_structural_keys(self, mol: Chem.Mol) -> List[str]:
        """提取结构关键词"""
        keys = []
        
        # 检查预定义的氨基酸模式
        for pattern_name, pattern_mol in self.pattern_mols.items():
            if mol.HasSubstructMatch(pattern_mol):
                keys.append(pattern_name)
        
        # 环信息
        ring_info = mol.GetRingInfo()
        if ring_info.NumRings() > 0:
            keys.append(f"rings_{ring_info.NumRings()}")
            for ring in ring_info.AtomRings():
                keys.append(f"ring_size_{len(ring)}")
        
        # 手性中心
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if chiral_centers:
            keys.append(f"chiral_centers_{len(chiral_centers)}")
        
        # 双键数量
        double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
        if double_bonds > 0:
            keys.append(f"double_bonds_{double_bonds}")
        
        return sorted(keys)
    
    def _analyze_stereochemistry(self, mol: Chem.Mol) -> str:
        """分析立体化学"""
        stereo_info = []
        
        # 手性中心
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        for center in chiral_centers:
            atom_idx, chiral_tag = center
            stereo_info.append(f"chiral_{atom_idx}_{chiral_tag}")
        
        # 双键立体化学
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                if bond.GetStereo() != Chem.BondStereo.STEREONONE:
                    stereo_info.append(f"double_bond_{bond.GetIdx()}_{bond.GetStereo()}")
        
        return "|".join(sorted(stereo_info))
    
    def _generate_simplified_fingerprint(self, smiles: str) -> StructuralFingerprint:
        """生成简化的结构指纹（不依赖RDKit）"""
        # 基于SMILES字符串的简化分析
        structural_keys = []
        
        # 检查常见的氨基酸特征
        if 'c1ccccc1' in smiles or 'C1=CC=CC=C1' in smiles:
            structural_keys.append('aromatic_ring')
        if 'N' in smiles:
            structural_keys.append('amino_group')
        if 'O' in smiles:
            structural_keys.append('oxygen_containing')
        if 'S' in smiles:
            structural_keys.append('sulfur_containing')
        if '=' in smiles:
            structural_keys.append('double_bond')
        if '#' in smiles:
            structural_keys.append('triple_bond')
        
        # 简化的拓扑哈希
        char_counts = defaultdict(int)
        for char in smiles:
            char_counts[char] += 1
        
        topological_hash = hashlib.md5(str(sorted(char_counts.items())).encode()).hexdigest()
        
        return StructuralFingerprint(
            structural_keys=structural_keys,
            topological_hash=topological_hash
        )
    
    def calculate_structural_similarity(self, fp1: StructuralFingerprint, 
                                      fp2: StructuralFingerprint) -> float:
        """计算结构相似性"""
        similarities = []
        
        # ECFP相似性
        if fp1.ecfp4 and fp2.ecfp4:
            ecfp_sim = self._calculate_tanimoto_similarity(fp1.ecfp4, fp2.ecfp4)
            similarities.append(('ecfp4', ecfp_sim, 0.4))
        
        if fp1.ecfp2 and fp2.ecfp2:
            ecfp2_sim = self._calculate_tanimoto_similarity(fp1.ecfp2, fp2.ecfp2)
            similarities.append(('ecfp2', ecfp2_sim, 0.3))
        
        # Morgan相似性
        if fp1.morgan2 and fp2.morgan2:
            morgan_sim = self._calculate_tanimoto_similarity(fp1.morgan2, fp2.morgan2)
            similarities.append(('morgan2', morgan_sim, 0.2))
        
        # 拓扑哈希相似性
        if fp1.topological_hash and fp2.topological_hash:
            topo_sim = 1.0 if fp1.topological_hash == fp2.topological_hash else 0.0
            similarities.append(('topological', topo_sim, 0.1))
        
        # 加权平均
        if not similarities:
            return 0.0
        
        total_weight = sum(weight for _, _, weight in similarities)
        weighted_sum = sum(sim * weight for _, sim, weight in similarities)
        
        return weighted_sum / total_weight if total_weight > 0 else 0.0
    
    def _calculate_tanimoto_similarity(self, fp1: List[int], fp2: List[int]) -> float:
        """计算Tanimoto相似性"""
        if len(fp1) != len(fp2):
            return 0.0
        
        fp1_int = [int(x) for x in fp1]
        fp2_int = [int(x) for x in fp2]
        
        intersection = sum(a & b for a, b in zip(fp1_int, fp2_int))
        union = sum(a | b for a, b in zip(fp1_int, fp2_int))
        
        return intersection / union if union > 0 else 0.0

class IsomerIdentifier:
    """同分异构体识别器"""
    
    def __init__(self, fingerprint_size: int = 2048):
        self.ecfp_generator = ECFPGenerator(fingerprint_size)
        self.fingerprint_cache = {}
        self.analysis_cache = {}
        
        # 相似性阈值 - 调整为更合理的值
        self.thresholds = {
            'identical': 0.95,        # 认为是同一分子
            'stereoisomer': 0.70,     # 立体异构体（降低阈值以提高检测）
            'structural_isomer': 0.60, # 结构异构体（降低阈值）
            'different': 0.40         # 不同分子
        }
    
    def analyze_isomers(self, smiles1: str, smiles2: str) -> IsomerAnalysisResult:
        """分析两个分子的同分异构体关系"""
        # 生成结构指纹
        fp1 = self._get_fingerprint(smiles1)
        fp2 = self._get_fingerprint(smiles2)
        
        # 计算结构相似性
        similarity = self.ecfp_generator.calculate_structural_similarity(fp1, fp2)
        
        # 判断异构体类型
        isomer_type = self._classify_isomer_type(fp1, fp2, similarity)
        
        # 分析结构差异
        differences = self._analyze_structural_differences(fp1, fp2)
        
        # 计算置信度
        confidence = self._calculate_confidence(fp1, fp2, similarity)
        
        return IsomerAnalysisResult(
            is_isomer=isomer_type in ['structural', 'stereoisomer'],
            isomer_type=isomer_type,
            similarity_score=similarity,
            structural_differences=differences,
            confidence=confidence,
            details={
                'fp1_keys': fp1.structural_keys or [],
                'fp2_keys': fp2.structural_keys or [],
                'stereochemistry_match': fp1.stereochemistry == fp2.stereochemistry,
                'topological_match': fp1.topological_hash == fp2.topological_hash
            }
        )
    
    def _get_fingerprint(self, smiles: str) -> StructuralFingerprint:
        """获取分子指纹（带缓存）"""
        if smiles in self.fingerprint_cache:
            return self.fingerprint_cache[smiles]
        
        fp = self.ecfp_generator.generate_structural_fingerprint(smiles)
        self.fingerprint_cache[smiles] = fp
        return fp
    
    def _classify_isomer_type(self, fp1: StructuralFingerprint, 
                            fp2: StructuralFingerprint, 
                            similarity: float) -> str:
        """分类异构体类型"""
        
        # 检查立体化学差异
        has_stereo_diff = (fp1.stereochemistry != fp2.stereochemistry and 
                          fp1.stereochemistry and fp2.stereochemistry)
        
        # 检查拓扑结构差异
        has_topo_diff = (fp1.topological_hash != fp2.topological_hash and 
                        fp1.topological_hash and fp2.topological_hash)
        
        # 检查结构关键词差异
        keys1 = set(fp1.structural_keys) if fp1.structural_keys else set()
        keys2 = set(fp2.structural_keys) if fp2.structural_keys else set()
        has_structural_diff = keys1 != keys2
        
        # 完全相同（高相似性且无结构差异）
        if (similarity >= self.thresholds['identical'] and 
            not has_stereo_diff and not has_topo_diff and not has_structural_diff):
            return 'identical'
        
        # 立体异构体（拓扑相同但立体化学不同，无其他重大结构差异）
        if (not has_topo_diff and has_stereo_diff and 
            similarity >= self.thresholds['stereoisomer'] and
            not has_structural_diff):
            return 'stereoisomer'
        
        # 结构异构体（有结构差异但相似性较高，且分子式可能相同）
        if ((has_topo_diff or has_structural_diff) and 
            similarity >= self.thresholds['structural_isomer']):
            return 'structural'
        
        # 如果有立体化学差异但相似性很高，且没有其他重大差异，可能是立体异构体
        if (has_stereo_diff and similarity >= 0.75 and not has_structural_diff):
            return 'stereoisomer'
        
        # 不同分子
        return 'different'
    
    def _analyze_structural_differences(self, fp1: StructuralFingerprint, 
                                      fp2: StructuralFingerprint) -> List[str]:
        """分析结构差异"""
        differences = []
        
        # 比较结构关键词
        if fp1.structural_keys and fp2.structural_keys:
            keys1 = set(fp1.structural_keys)
            keys2 = set(fp2.structural_keys)
            
            unique_to_1 = keys1 - keys2
            unique_to_2 = keys2 - keys1
            
            if unique_to_1:
                differences.append(f"unique_to_mol1: {list(unique_to_1)}")
            if unique_to_2:
                differences.append(f"unique_to_mol2: {list(unique_to_2)}")
        
        # 立体化学差异
        if fp1.stereochemistry != fp2.stereochemistry:
            differences.append("stereochemistry_different")
        
        # 拓扑差异
        if fp1.topological_hash != fp2.topological_hash:
            differences.append("topological_different")
        
        return differences
    
    def _calculate_confidence(self, fp1: StructuralFingerprint, 
                            fp2: StructuralFingerprint, 
                            similarity: float) -> float:
        """计算分析置信度"""
        confidence_factors = []
        
        # 指纹完整性
        if fp1.ecfp4 and fp2.ecfp4:
            confidence_factors.append(0.4)
        if fp1.ecfp2 and fp2.ecfp2:
            confidence_factors.append(0.3)
        if fp1.morgan2 and fp2.morgan2:
            confidence_factors.append(0.2)
        if fp1.topological_hash and fp2.topological_hash:
            confidence_factors.append(0.1)
        
        # 结构关键词一致性
        if fp1.structural_keys and fp2.structural_keys:
            key_overlap = len(set(fp1.structural_keys) & set(fp2.structural_keys))
            key_total = len(set(fp1.structural_keys) | set(fp2.structural_keys))
            if key_total > 0:
                confidence_factors.append(key_overlap / key_total * 0.2)
        
        # 基于相似性的置信度调整
        if similarity > 0.9:
            confidence_factors.append(0.1)
        elif similarity > 0.7:
            confidence_factors.append(0.05)
        
        return min(1.0, sum(confidence_factors))
    
    def enhanced_composition_similarity(self, smiles1: str, smiles2: str,
                                      basic_similarity: float) -> float:
        """增强的组成相似性计算"""
        # 结构相似性
        fp1 = self._get_fingerprint(smiles1)
        fp2 = self._get_fingerprint(smiles2)
        structural_similarity = self.ecfp_generator.calculate_structural_similarity(fp1, fp2)
        
        # 加权组合
        # 如果基础相似性很高（同分异构体），结构相似性权重更大
        if basic_similarity > 0.95:
            return 0.3 * basic_similarity + 0.7 * structural_similarity
        else:
            return 0.6 * basic_similarity + 0.4 * structural_similarity

# 工具函数
def detect_potential_isomers(smiles_list: List[str]) -> Dict[str, List[str]]:
    """检测潜在的同分异构体组"""
    identifier = IsomerIdentifier()
    groups = defaultdict(list)
    
    for i, smiles1 in enumerate(smiles_list):
        for j, smiles2 in enumerate(smiles_list[i+1:], i+1):
            result = identifier.analyze_isomers(smiles1, smiles2)
            if result.is_isomer:
                group_key = f"group_{min(i, j)}"
                groups[group_key].extend([smiles1, smiles2])
    
    # 去重
    for group_key in groups:
        groups[group_key] = list(set(groups[group_key]))
    
    return dict(groups)

def main():
    """测试模块"""
    # 测试同分异构体识别
    test_pairs = [
        # 同分异构体：异亮氨酸 vs 亮氨酸
        ("CC(C)C[C@@H](C(=O)O)N", "C[C@H](C)C[C@@H](C(=O)O)N"),
        
        # 立体异构体：L-丙氨酸 vs D-丙氨酸
        ("C[C@@H](C(=O)O)N", "C[C@H](C(=O)O)N"),
        
        # 相同分子
        ("C[C@@H](C(=O)O)N", "C[C@@H](C(=O)O)N"),
        
        # 不同分子
        ("C[C@@H](C(=O)O)N", "CC(C)C[C@@H](C(=O)O)N")
    ]
    
    identifier = IsomerIdentifier()
    
    for i, (smiles1, smiles2) in enumerate(test_pairs):
        print(f"\n测试对 {i+1}:")
        print(f"分子1: {smiles1}")
        print(f"分子2: {smiles2}")
        
        result = identifier.analyze_isomers(smiles1, smiles2)
        print(f"是否为异构体: {result.is_isomer}")
        print(f"异构体类型: {result.isomer_type}")
        print(f"相似性分数: {result.similarity_score:.3f}")
        print(f"置信度: {result.confidence:.3f}")
        if result.structural_differences:
            print(f"结构差异: {result.structural_differences}")

if __name__ == "__main__":
    main()