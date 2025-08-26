"""
3D结构匹配器
专注于高精度的3D结构比对和异构体区分
"""

import numpy as np
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import logging

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Geometry import Point3D

from .models_simple import ResidueInfo, AminoAcidInfo


@dataclass
class StructureAlignment:
    """结构对齐结果"""
    rmsd: float
    transformation_matrix: np.ndarray
    aligned_coords: np.ndarray
    atom_mapping: List[int]
    confidence: float


class StructureMatcher:
    """
    3D结构匹配器
    专注于准确的三维结构比对，用于区分立体异构体
    """
    
    def __init__(self, rmsd_threshold: float = 0.5):
        """
        初始化结构匹配器
        
        Args:
            rmsd_threshold: RMSD阈值，用于判断结构相似性
        """
        self.rmsd_threshold = rmsd_threshold
    
    def match_structures(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> StructureAlignment:
        """
        匹配两个结构
        
        Args:
            residue: PDB残基
            amino_acid: 参考氨基酸
            
        Returns:
            结构对齐结果
        """
        try:
            # 生成参考分子的3D构象
            ref_mol = self._generate_reference_conformer(amino_acid.smiles)
            if ref_mol is None:
                return self._create_failed_alignment()
            
            # 提取坐标
            residue_coords = self._extract_coordinates(residue)
            ref_coords = self._extract_mol_coordinates(ref_mol)
            
            # 原子匹配
            atom_mapping = self._match_atoms(residue, ref_mol)
            
            # 结构对齐
            if len(atom_mapping) >= 3:  # 至少需要3个原子进行对齐
                rmsd, transformation, aligned_coords = self._align_structures(
                    residue_coords, ref_coords, atom_mapping
                )
                
                confidence = self._calculate_alignment_confidence(rmsd, len(atom_mapping))
                
                return StructureAlignment(
                    rmsd=rmsd,
                    transformation_matrix=transformation,
                    aligned_coords=aligned_coords,
                    atom_mapping=atom_mapping,
                    confidence=confidence
                )
            else:
                return self._create_failed_alignment("原子匹配不足")
                
        except Exception as e:
            logging.error(f"结构匹配失败: {e}")
            return self._create_failed_alignment(str(e))
    
    def _generate_reference_conformer(self, smiles: str) -> Optional[Chem.Mol]:
        """从SMILES生成参考分子的最优构象"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            # 添加氢原子
            mol = Chem.AddHs(mol)
            
            # 生成3D构象
            AllChem.EmbedMolecule(mol, randomSeed=42)
            
            # 优化构象（使用MMFF94力场）
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except:
                # 如果MMFF94失败，尝试UFF
                try:
                    AllChem.UFFOptimizeMolecule(mol)
                except:
                    # 如果都失败，跳过优化
                    pass
            
            return mol
            
        except Exception as e:
            logging.warning(f"参考构象生成失败: {e}")
            return None
    
    def _extract_coordinates(self, residue: ResidueInfo) -> np.ndarray:
        """提取残基坐标"""
        coords = []
        for atom in residue.atoms:
            coords.append([atom.x, atom.y, atom.z])
        return np.array(coords)
    
    def _extract_mol_coordinates(self, mol: Chem.Mol) -> np.ndarray:
        """提取RDKit分子坐标"""
        conf = mol.GetConformer()
        coords = []
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append([pos.x, pos.y, pos.z])
        return np.array(coords)
    
    def _match_atoms(self, residue: ResidueInfo, ref_mol: Chem.Mol) -> List[int]:
        """
        匹配原子对应关系
        
        Returns:
            原子映射列表，索引为残基原子，值为参考分子原子索引
        """
        atom_mapping = []
        
        # 简单的元素匹配策略
        residue_elements = [atom.element for atom in residue.atoms]
        ref_elements = [ref_mol.GetAtomWithIdx(i).GetSymbol() for i in range(ref_mol.GetNumAtoms())]
        
        for i, residue_element in enumerate(residue_elements):
            # 找到第一个匹配的元素
            for j, ref_element in enumerate(ref_elements):
                if residue_element == ref_element and j not in atom_mapping:
                    atom_mapping.append(j)
                    break
            else:
                # 没有找到匹配的原子
                atom_mapping.append(-1)
        
        # 过滤掉无效映射
        valid_mapping = [idx for idx in atom_mapping if idx != -1]
        
        return valid_mapping
    
    def _align_structures(self, coords1: np.ndarray, coords2: np.ndarray, 
                         atom_mapping: List[int]) -> Tuple[float, np.ndarray, np.ndarray]:
        """
        使用Kabsch算法对齐结构
        
        Returns:
            (RMSD, 变换矩阵, 对齐后的坐标)
        """
        # 选择匹配的原子坐标
        mapped_coords1 = coords1[:len(atom_mapping)]
        mapped_coords2 = coords2[atom_mapping]
        
        # 中心化坐标
        center1 = np.mean(mapped_coords1, axis=0)
        center2 = np.mean(mapped_coords2, axis=0)
        
        centered_coords1 = mapped_coords1 - center1
        centered_coords2 = mapped_coords2 - center2
        
        # 计算协方差矩阵
        H = centered_coords1.T @ centered_coords2
        
        # SVD分解
        U, S, Vt = np.linalg.svd(H)
        
        # 计算旋转矩阵
        R = Vt.T @ U.T
        
        # 确保是右手坐标系
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        
        # 计算平移向量
        t = center2 - center1 @ R.T
        
        # 构建变换矩阵
        transformation = np.eye(4)
        transformation[:3, :3] = R
        transformation[:3, 3] = t
        
        # 应用变换
        aligned_coords = centered_coords1 @ R.T + center2
        
        # 计算RMSD
        diff = aligned_coords - centered_coords2
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        
        return rmsd, transformation, aligned_coords
    
    def _calculate_alignment_confidence(self, rmsd: float, num_atoms: int) -> float:
        """计算对齐置信度"""
        # 基于RMSD和匹配原子数计算置信度
        rmsd_score = max(0, 1 - rmsd / self.rmsd_threshold)
        
        # 原子数量奖励（更多匹配原子 = 更高置信度）
        atom_score = min(1.0, num_atoms / 20.0)  # 假设20个原子为满分
        
        # 综合评分
        confidence = 0.7 * rmsd_score + 0.3 * atom_score
        
        return max(0.0, min(1.0, confidence))
    
    def _create_failed_alignment(self, error_msg: str = "对齐失败") -> StructureAlignment:
        """创建失败的对齐结果"""
        return StructureAlignment(
            rmsd=float('inf'),
            transformation_matrix=np.eye(4),
            aligned_coords=np.array([]),
            atom_mapping=[],
            confidence=0.0
        )
    
    def detect_chirality_difference(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> Dict[str, bool]:
        """
        检测手性差异
        
        Returns:
            手性分析结果
        """
        try:
            # 从坐标构建分子
            residue_mol = self._build_molecule_from_coords(residue)
            ref_mol = Chem.MolFromSmiles(amino_acid.smiles)
            
            if residue_mol is None or ref_mol is None:
                return {"has_chirality": False, "chirality_matches": False}
            
            # 检测手性中心
            residue_chiral = Chem.FindMolChiralCenters(residue_mol, includeUnassigned=True)
            ref_chiral = Chem.FindMolChiralCenters(ref_mol, includeUnassigned=True)
            
            has_chirality = len(residue_chiral) > 0 or len(ref_chiral) > 0
            
            # 比较手性配置
            chirality_matches = True
            if has_chirality:
                # 简化比较：检查手性中心数量是否匹配
                if len(residue_chiral) != len(ref_chiral):
                    chirality_matches = False
                else:
                    # 详细的手性比较需要更复杂的算法
                    # 这里使用简化的方法
                    residue_stereo = Chem.MolToSmiles(residue_mol, isomericSmiles=True)
                    ref_stereo = Chem.MolToSmiles(ref_mol, isomericSmiles=True)
                    chirality_matches = residue_stereo == ref_stereo
            
            return {
                "has_chirality": has_chirality,
                "chirality_matches": chirality_matches,
                "residue_chiral_centers": len(residue_chiral),
                "reference_chiral_centers": len(ref_chiral)
            }
            
        except Exception as e:
            logging.warning(f"手性差异检测失败: {e}")
            return {"has_chirality": False, "chirality_matches": False}
    
    def _build_molecule_from_coords(self, residue: ResidueInfo) -> Optional[Chem.Mol]:
        """从坐标构建分子对象"""
        try:
            mol = Chem.RWMol()
            
            # 添加原子
            for atom in residue.atoms:
                rd_atom = Chem.Atom(atom.element)
                mol.AddAtom(rd_atom)
            
            # 设置坐标
            conf = Chem.Conformer(len(residue.atoms))
            for i, atom in enumerate(residue.atoms):
                point = Point3D(atom.x, atom.y, atom.z)
                conf.SetAtomPosition(i, point)
            mol.AddConformer(conf)
            
            # 推断化学键
            from rdkit.Chem import rdDetermineBonds
            rdDetermineBonds.DetermineBonds(mol)
            
            # 清理分子
            Chem.SanitizeMol(mol)
            
            return mol
            
        except Exception as e:
            logging.warning(f"从坐标构建分子失败: {e}")
            return None
    
    def calculate_shape_similarity(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> float:
        """
        计算分子形状相似度
        
        Returns:
            形状相似度分数 (0-1)
        """
        try:
            alignment = self.match_structures(residue, amino_acid)
            
            if alignment.confidence > 0:
                # 基于RMSD计算形状相似度
                max_rmsd = 5.0  # 最大可接受RMSD
                shape_similarity = max(0, 1 - alignment.rmsd / max_rmsd)
                return shape_similarity
            
            return 0.0
            
        except Exception as e:
            logging.warning(f"形状相似度计算失败: {e}")
            return 0.0


class StereochemistryAnalyzer:
    """
    立体化学分析器
    专门用于分析和比较分子的立体化学特征
    """
    
    def __init__(self):
        self.structure_matcher = StructureMatcher()
    
    def analyze_stereochemistry(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> Dict[str, any]:
        """
        完整的立体化学分析
        
        Returns:
            立体化学分析结果
        """
        result = {
            "is_stereoisomer": False,
            "chirality_analysis": {},
            "conformation_analysis": {},
            "confidence": 0.0
        }
        
        try:
            # 手性分析
            chirality_info = self.structure_matcher.detect_chirality_difference(residue, amino_acid)
            result["chirality_analysis"] = chirality_info
            
            # 构象分析
            alignment = self.structure_matcher.match_structures(residue, amino_acid)
            result["conformation_analysis"] = {
                "rmsd": alignment.rmsd,
                "structural_similarity": alignment.confidence
            }
            
            # 判断是否为立体异构体
            if chirality_info.get("has_chirality", False):
                if not chirality_info.get("chirality_matches", True):
                    result["is_stereoisomer"] = True
                    result["confidence"] = 0.9
                elif alignment.rmsd < 2.0:  # 结构相似但可能是构象异构体
                    result["is_stereoisomer"] = True
                    result["confidence"] = 0.7
            
            return result
            
        except Exception as e:
            logging.error(f"立体化学分析失败: {e}")
            return result