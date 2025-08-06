"""
3D结构验证器
基于Kabsch算法进行3D结构匹配验证
"""

import numpy as np
from typing import Dict, Any, List, Optional, Tuple

from .base import BaseVerifier
from ..models import ResidueInfo, AminoAcidInfo, VerificationMethod, Structure3DInfo


class Structure3DVerifier(BaseVerifier):
    """3D结构验证器"""
    
    def __init__(self, threshold: float = 0.5):
        super().__init__(threshold)
    
    def get_method(self) -> VerificationMethod:
        return VerificationMethod.STRUCTURE_3D
    
    def calculate_score(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> float:
        """
        计算3D结构验证分数
        
        验证逻辑：
        1. 提取残基和标准结构的坐标
        2. 使用Kabsch算法计算RMSD
        3. 转换为动态评分
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            验证分数 (0.0-1.0)
        """
        # 获取残基坐标
        residue_coords = self._extract_residue_coordinates(residue)
        if residue_coords is None or len(residue_coords) == 0:
            return 0.0
        
        # 获取标准结构坐标（这里需要从数据库或其他来源获取）
        standard_coords = self._get_standard_structure_coordinates(amino_acid)
        if standard_coords is None or len(standard_coords) == 0:
            return 0.0
        
        # 计算RMSD
        rmsd = self._calculate_kabsch_rmsd(residue_coords, standard_coords)
        if rmsd is None:
            return 0.0
        
        # 转换为动态评分
        atom_count = len(residue_coords)
        score = self._rmsd_to_score_dynamic(rmsd, atom_count)
        
        return score
    
    def _extract_residue_coordinates(self, residue: ResidueInfo) -> Optional[np.ndarray]:
        """
        提取残基坐标
        
        Args:
            residue: 残基信息
        
        Returns:
            坐标数组 (N, 3)，失败返回None
        """
        if not residue.atoms:
            return None
        
        # 过滤氢原子（与数据库保持一致）
        heavy_atoms = [atom for atom in residue.atoms if atom.element != 'H']
        
        if not heavy_atoms:
            return None
        
        coords = np.array([[atom.x, atom.y, atom.z] for atom in heavy_atoms])
        return coords
    
    def _get_standard_structure_coordinates(self, amino_acid: AminoAcidInfo) -> Optional[np.ndarray]:
        """
        获取标准结构坐标
        
        注意：这里需要实现从SMILES或其他来源生成3D坐标的功能
        目前返回None表示不支持
        
        Args:
            amino_acid: 氨基酸信息
        
        Returns:
            标准结构坐标数组，失败返回None
        """
        # 这里可以实现：
        # 1. 从SMILES生成3D坐标
        # 2. 从数据库中获取预计算的3D结构
        # 3. 使用分子力学优化生成标准构象
        
        # 目前返回None，表示不支持3D结构验证
        return None
    
    def _calculate_kabsch_rmsd(self, coords1: np.ndarray, coords2: np.ndarray) -> Optional[float]:
        """
        使用Kabsch算法计算RMSD
        
        Args:
            coords1: 第一组坐标 (N, 3)
            coords2: 第二组坐标 (M, 3)
        
        Returns:
            RMSD值，失败返回None
        """
        if coords1.shape[0] != coords2.shape[0]:
            # 原子数量不匹配，无法直接比较
            return None
        
        if coords1.shape[0] < 3:
            # 原子数量太少，无法进行有意义的结构比较
            return None
        
        try:
            # 1. 中心化坐标
            centroid1 = np.mean(coords1, axis=0)
            centroid2 = np.mean(coords2, axis=0)
            coords1_centered = coords1 - centroid1
            coords2_centered = coords2 - centroid2
            
            # 2. Kabsch算法：计算最优旋转矩阵
            H = coords1_centered.T @ coords2_centered
            U, S, Vt = np.linalg.svd(H)
            R = Vt.T @ U.T
            
            # 3. 确保是右手坐标系
            if np.linalg.det(R) < 0:
                Vt[-1, :] *= -1
                R = Vt.T @ U.T
            
            # 4. 应用最优旋转计算RMSD
            coords1_rotated = coords1_centered @ R.T
            diff = coords1_rotated - coords2_centered
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
            
            return float(rmsd)
        
        except Exception as e:
            print(f"⚠️ Kabsch算法计算失败: {e}")
            return None
    
    def _rmsd_to_score_dynamic(self, rmsd: float, atom_count: int) -> float:
        """
        将RMSD转换为动态评分
        
        Args:
            rmsd: RMSD值
            atom_count: 原子数量
        
        Returns:
            评分 (0.0-1.0)
        """
        # 基于原子数的动态阈值
        if atom_count <= 10:
            max_rmsd = 1.0  # 小分子更严格
        elif atom_count <= 20:
            max_rmsd = 1.5  # 中等分子
        else:
            max_rmsd = 2.0  # 大分子稍微宽松
        
        # 使用指数衰减函数
        score = np.exp(-rmsd / max_rmsd)
        return np.clip(score, 0.0, 1.0)
    
    def _get_verification_details(self, residue: ResidueInfo, amino_acid: AminoAcidInfo, score: float) -> Dict[str, Any]:
        """获取3D结构验证详情"""
        details = super()._get_verification_details(residue, amino_acid, score)
        
        # 提取坐标信息
        residue_coords = self._extract_residue_coordinates(residue)
        standard_coords = self._get_standard_structure_coordinates(amino_acid)
        
        # 计算RMSD
        rmsd = None
        if residue_coords is not None and standard_coords is not None:
            rmsd = self._calculate_kabsch_rmsd(residue_coords, standard_coords)
        
        atom_count = len(residue_coords) if residue_coords is not None else 0
        
        details.update({
            'residue_atom_count': atom_count,
            'rmsd': rmsd,
            'max_rmsd_threshold': self._get_max_rmsd_threshold(atom_count),
            'coordinates_available': {
                'residue': residue_coords is not None,
                'standard': standard_coords is not None
            }
        })
        
        return details
    
    def _get_max_rmsd_threshold(self, atom_count: int) -> float:
        """获取最大RMSD阈值"""
        if atom_count <= 10:
            return 1.0
        elif atom_count <= 20:
            return 1.5
        else:
            return 2.0
    
    def calculate_detailed_structure_analysis(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """
        计算详细的结构分析
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            详细结构分析结果
        """
        result = {
            'residue_id': residue.residue_key,
            'amino_acid_id': amino_acid.id
        }
        
        # 提取坐标
        residue_coords = self._extract_residue_coordinates(residue)
        standard_coords = self._get_standard_structure_coordinates(amino_acid)
        
        result['coordinates'] = {
            'residue_available': residue_coords is not None,
            'standard_available': standard_coords is not None,
            'residue_atom_count': len(residue_coords) if residue_coords is not None else 0,
            'standard_atom_count': len(standard_coords) if standard_coords is not None else 0
        }
        
        if residue_coords is not None and standard_coords is not None:
            # 计算RMSD
            rmsd = self._calculate_kabsch_rmsd(residue_coords, standard_coords)
            result['rmsd'] = rmsd
            
            if rmsd is not None:
                atom_count = len(residue_coords)
                score = self._rmsd_to_score_dynamic(rmsd, atom_count)
                result['kabsch_score'] = score
                result['max_rmsd_threshold'] = self._get_max_rmsd_threshold(atom_count)
                
                # 几何分析
                result['geometry_analysis'] = self._analyze_geometry(residue_coords, standard_coords)
        
        return result
    
    def _analyze_geometry(self, coords1: np.ndarray, coords2: np.ndarray) -> Dict[str, Any]:
        """
        分析几何特征
        
        Args:
            coords1: 第一组坐标
            coords2: 第二组坐标
        
        Returns:
            几何分析结果
        """
        analysis = {}
        
        try:
            # 计算质心
            centroid1 = np.mean(coords1, axis=0)
            centroid2 = np.mean(coords2, axis=0)
            analysis['centroids'] = {
                'residue': centroid1.tolist(),
                'standard': centroid2.tolist(),
                'distance': float(np.linalg.norm(centroid1 - centroid2))
            }
            
            # 计算惯性矩
            inertia1 = self._calculate_inertia_tensor(coords1 - centroid1)
            inertia2 = self._calculate_inertia_tensor(coords2 - centroid2)
            
            # 计算主轴
            eigenvals1, eigenvecs1 = np.linalg.eigh(inertia1)
            eigenvals2, eigenvecs2 = np.linalg.eigh(inertia2)
            
            analysis['inertia'] = {
                'eigenvalues_residue': eigenvals1.tolist(),
                'eigenvalues_standard': eigenvals2.tolist(),
                'shape_similarity': float(np.corrcoef(eigenvals1, eigenvals2)[0, 1])
            }
            
        except Exception as e:
            analysis['error'] = str(e)
        
        return analysis
    
    def _calculate_inertia_tensor(self, coords: np.ndarray) -> np.ndarray:
        """
        计算惯性张量
        
        Args:
            coords: 中心化后的坐标
        
        Returns:
            惯性张量矩阵 (3, 3)
        """
        n = coords.shape[0]
        I = np.zeros((3, 3))
        
        for i in range(n):
            x, y, z = coords[i]
            I[0, 0] += y*y + z*z
            I[1, 1] += x*x + z*z
            I[2, 2] += x*x + y*y
            I[0, 1] -= x*y
            I[0, 2] -= x*z
            I[1, 2] -= y*z
        
        # 对称化
        I[1, 0] = I[0, 1]
        I[2, 0] = I[0, 2]
        I[2, 1] = I[1, 2]
        
        return I / n
    
    def create_structure_3d_info(self, residue: ResidueInfo) -> Structure3DInfo:
        """
        创建3D结构信息对象
        
        Args:
            residue: 残基信息
        
        Returns:
            3D结构信息
        """
        coords = self._extract_residue_coordinates(residue)
        elements = [atom.element for atom in residue.atoms if atom.element != 'H']
        
        return Structure3DInfo(
            coordinates=coords.tolist() if coords is not None else [],
            elements=elements
        )