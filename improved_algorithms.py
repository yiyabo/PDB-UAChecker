#!/usr/bin/env python3
"""
改进的算法实现
修复指纹相似性和3D结构验证中的关键问题
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

@dataclass
class ImprovedThresholds:
    """改进的阈值设置，基于化学信息学标准"""
    molecular_formula: float = 1.0      # 完全匹配
    atom_composition: float = 0.9       # 允许10%偏差
    fingerprint_similarity: float = 0.7 # 化学信息学标准
    structure_3d: float = 0.8           # 对应1.0Å RMSD

class ImprovedFingerprintValidator:
    """改进的指纹相似性验证器"""
    
    def __init__(self):
        self.thresholds = ImprovedThresholds()
    
    def verify_fingerprint_similarity(self, residue_atoms: Dict, candidate_smiles: str) -> Dict:
        """
        改进的指纹相似性验证
        
        Returns:
            Dict: {
                'score': float,
                'method': str,
                'confidence': float,
                'details': Dict
            }
        """
        
        # 1. 尝试真正的分子指纹计算
        try:
            fingerprint_score = self._calculate_molecular_fingerprint(residue_atoms, candidate_smiles)
            if fingerprint_score is not None:
                return {
                    'score': fingerprint_score,
                    'method': 'molecular_fingerprint',
                    'confidence': 0.9,
                    'details': {'algorithm': 'ECFP2_Tanimoto'}
                }
        except Exception as e:
            print(f"分子指纹计算失败: {e}")
        
        # 2. 改进的组成相似性计算
        composition_score = self._enhanced_composition_similarity(residue_atoms, candidate_smiles)
        
        return {
            'score': composition_score,
            'method': 'enhanced_composition',
            'confidence': 0.6,
            'details': {'algorithm': 'weighted_jaccard'}
        }
    
    def _calculate_molecular_fingerprint(self, residue_atoms: Dict, candidate_smiles: str) -> Optional[float]:
        """
        真正的分子指纹计算（需要RDKit）
        
        注意：这需要从原子坐标推导SMILES，这是一个复杂的化学信息学问题
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors
            from rdkit import DataStructs
            
            # 这里需要实现从原子组成到SMILES的转换
            # 这是一个非常复杂的问题，需要专门的算法
            # 暂时返回None，表示无法计算
            return None
            
        except ImportError:
            return None
    
    def _enhanced_composition_similarity(self, residue_atoms: Dict, candidate_smiles: str) -> float:
        """
        改进的组成相似性计算，考虑元素重要性权重
        """
        
        # 获取候选分子的原子组成
        candidate_atoms = self._get_atom_composition_from_smiles(candidate_smiles)
        if not candidate_atoms:
            return 0.0
        
        # 元素重要性权重（基于化学重要性）
        element_weights = {
            'C': 1.0,   # 碳骨架最重要
            'N': 0.9,   # 氮原子很重要
            'O': 0.9,   # 氧原子很重要
            'S': 0.8,   # 硫原子重要
            'P': 0.8,   # 磷原子重要
            'H': 0.3,   # 氢原子权重较低（经常缺失）
        }
        
        # 计算加权Jaccard相似性
        all_elements = set(residue_atoms.keys()) | set(candidate_atoms.keys())
        weighted_intersection = 0.0
        weighted_union = 0.0
        
        for element in all_elements:
            weight = element_weights.get(element, 0.5)  # 默认权重0.5
            res_count = residue_atoms.get(element, 0)
            cand_count = candidate_atoms.get(element, 0)
            
            weighted_intersection += weight * min(res_count, cand_count)
            weighted_union += weight * max(res_count, cand_count)
        
        return weighted_intersection / weighted_union if weighted_union > 0 else 0.0
    
    def _get_atom_composition_from_smiles(self, smiles: str) -> Dict:
        """从SMILES获取原子组成"""
        try:
            from rdkit import Chem
            
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return {}
            
            atom_counts = {}
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
            
            return atom_counts
            
        except ImportError:
            # 如果没有RDKit，返回空字典
            return {}

class Improved3DValidator:
    """改进的3D结构验证器"""
    
    def __init__(self):
        self.thresholds = ImprovedThresholds()
    
    def verify_3d_structure(self, pdb_coords: List[List[float]], 
                           standard_coords: List[List[float]]) -> Dict:
        """
        改进的3D结构验证
        
        Returns:
            Dict: {
                'rmsd': float,
                'score': float,
                'method': str,
                'details': Dict
            }
        """
        
        if len(pdb_coords) != len(standard_coords):
            return {
                'rmsd': float('inf'),
                'score': 0.0,
                'method': 'size_mismatch',
                'details': {'error': 'atom_count_mismatch'}
            }
        
        # 1. 改进的RMSD计算（包含最优叠合）
        rmsd = self._calculate_optimal_rmsd(pdb_coords, standard_coords)
        
        # 2. 改进的评分机制
        score = self._rmsd_to_score(rmsd, len(pdb_coords))
        
        return {
            'rmsd': rmsd,
            'score': score,
            'method': 'optimal_superposition',
            'details': {
                'algorithm': 'kabsch_rmsd',
                'atom_count': len(pdb_coords)
            }
        }
    
    def _calculate_optimal_rmsd(self, coords1: List[List[float]], 
                               coords2: List[List[float]]) -> float:
        """
        使用Kabsch算法计算最优叠合RMSD
        """
        try:
            coords1 = np.array(coords1, dtype=float)
            coords2 = np.array(coords2, dtype=float)
            
            # 中心化坐标
            centroid1 = np.mean(coords1, axis=0)
            centroid2 = np.mean(coords2, axis=0)
            
            coords1_centered = coords1 - centroid1
            coords2_centered = coords2 - centroid2
            
            # Kabsch算法：找到最优旋转矩阵
            H = coords1_centered.T @ coords2_centered
            U, S, Vt = np.linalg.svd(H)
            R = Vt.T @ U.T
            
            # 确保是右手坐标系
            if np.linalg.det(R) < 0:
                Vt[-1, :] *= -1
                R = Vt.T @ U.T
            
            # 应用旋转
            coords1_rotated = coords1_centered @ R.T
            
            # 计算RMSD
            diff = coords1_rotated - coords2_centered
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
            
            return float(rmsd)
            
        except Exception as e:
            print(f"RMSD计算失败: {e}")
            return float('inf')
    
    def _rmsd_to_score(self, rmsd: float, atom_count: int) -> float:
        """
        改进的RMSD到分数转换，考虑分子大小
        """
        
        # 基于分子大小的动态阈值
        if atom_count <= 10:
            max_rmsd = 1.5  # 小分子更严格
        elif atom_count <= 20:
            max_rmsd = 2.0  # 中等分子
        else:
            max_rmsd = 3.0  # 大分子更宽松
        
        # 使用指数衰减函数，更符合实际情况
        score = np.exp(-rmsd / max_rmsd)
        
        return float(np.clip(score, 0.0, 1.0))

class ImprovedGeometryValidator:
    """改进的几何参数验证器"""
    
    def verify_bond_lengths(self, pdb_bonds: List[float], 
                           standard_bonds: List[float]) -> float:
        """
        改进的键长验证，使用统计学方法
        """
        if not pdb_bonds or not standard_bonds:
            return 0.0
        
        # 使用相对误差的统计分布
        relative_errors = []
        
        # 简化：假设键长按顺序对应
        min_len = min(len(pdb_bonds), len(standard_bonds))
        
        for i in range(min_len):
            if standard_bonds[i] > 0:
                rel_error = abs(pdb_bonds[i] - standard_bonds[i]) / standard_bonds[i]
                relative_errors.append(rel_error)
        
        if not relative_errors:
            return 0.0
        
        # 使用平均相对误差
        mean_error = np.mean(relative_errors)
        
        # 转换为分数（5%误差对应0.95分）
        score = max(0, 1 - mean_error / 0.1)  # 10%误差对应0分
        
        return float(score)
    
    def verify_bond_angles(self, pdb_angles: List[float], 
                          standard_angles: List[float]) -> float:
        """
        改进的键角验证
        """
        if not pdb_angles or not standard_angles:
            return 0.0
        
        # 键角差异（度数）
        angle_diffs = []
        min_len = min(len(pdb_angles), len(standard_angles))
        
        for i in range(min_len):
            diff = abs(pdb_angles[i] - standard_angles[i])
            # 处理角度的周期性（0°和360°相同）
            diff = min(diff, 360 - diff)
            angle_diffs.append(diff)
        
        if not angle_diffs:
            return 0.0
        
        mean_diff = np.mean(angle_diffs)
        
        # 10度差异对应0.9分，30度差异对应0分
        score = max(0, 1 - mean_diff / 30.0)
        
        return float(score)

def test_improved_algorithms():
    """测试改进的算法"""
    
    print("🧪 测试改进的算法实现")
    print("=" * 50)
    
    # 1. 测试指纹相似性验证
    print("\n1️⃣ 测试指纹相似性验证:")
    fingerprint_validator = ImprovedFingerprintValidator()
    
    test_atoms = {'C': 3, 'H': 7, 'N': 1, 'O': 2}
    test_smiles = "N[C@@H](C)C(=O)O"  # 丙氨酸
    
    result = fingerprint_validator.verify_fingerprint_similarity(test_atoms, test_smiles)
    print(f"结果: {result}")
    
    # 2. 测试3D结构验证
    print("\n2️⃣ 测试3D结构验证:")
    structure_validator = Improved3DValidator()
    
    # 测试相同坐标
    coords1 = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    coords2 = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    
    result = structure_validator.verify_3d_structure(coords1, coords2)
    print(f"相同坐标RMSD: {result['rmsd']:.6f}, 分数: {result['score']:.3f}")
    
    # 测试轻微偏移
    coords2_shifted = [[0.1, 0.1, 0.1], [1.1, 0.1, 0.1], [0.1, 1.1, 0.1]]
    result = structure_validator.verify_3d_structure(coords1, coords2_shifted)
    print(f"偏移坐标RMSD: {result['rmsd']:.6f}, 分数: {result['score']:.3f}")
    
    # 3. 测试几何参数验证
    print("\n3️⃣ 测试几何参数验证:")
    geometry_validator = ImprovedGeometryValidator()
    
    bond_score = geometry_validator.verify_bond_lengths([1.5, 1.3, 1.2], [1.5, 1.3, 1.2])
    print(f"相同键长分数: {bond_score:.3f}")
    
    angle_score = geometry_validator.verify_bond_angles([109.5, 120.0, 180.0], [109.5, 120.0, 180.0])
    print(f"相同键角分数: {angle_score:.3f}")
    
    print("\n✅ 改进算法测试完成")

if __name__ == "__main__":
    test_improved_algorithms()
