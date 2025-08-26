"""
异构体检测器
专注于从PDB坐标识别非天然氨基酸异构体的核心模块
"""

import numpy as np
from typing import List, Dict, Optional, Tuple, Set
from dataclasses import dataclass
import logging

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdDetermineBonds, Descriptors
from rdkit.Geometry import Point3D

from .models_simple import ResidueInfo, AminoAcidInfo, VerificationResult, VerificationScore, VerificationMethod


@dataclass
class IsomerMatch:
    """异构体匹配结果"""
    amino_acid_id: str
    amino_acid_name: str
    smiles_similarity: float
    structure_3d_rmsd: float
    fingerprint_similarity: float
    overall_confidence: float
    isomer_type: str  # "identical", "structural", "stereoisomer", "conformational"


class IsomerDetector:
    """
    异构体检测器
    专注于高精度的异构体识别，假设RDKit可用
    """
    
    def __init__(self, rmsd_threshold: float = 0.5, fingerprint_threshold: float = 0.8):
        """
        初始化检测器
        
        Args:
            rmsd_threshold: 3D结构RMSD阈值
            fingerprint_threshold: 分子指纹相似度阈值
        """
        self.rmsd_threshold = rmsd_threshold
        self.fingerprint_threshold = fingerprint_threshold
        
        # 原子半径用于键推断
        self.covalent_radii = {
            'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05,
            'P': 1.07, 'F': 0.57, 'Cl': 0.99, 'Br': 1.14, 'I': 1.33
        }
    
    def detect_isomer(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> IsomerMatch:
        """
        检测异构体类型和相似度
        
        Args:
            residue: 残基信息
            amino_acid: 候选氨基酸
            
        Returns:
            异构体匹配结果
        """
        try:
            # 1. 从坐标生成SMILES
            residue_mol = self._coordinates_to_molecule(residue)
            if residue_mol is None:
                return self._create_failed_match(amino_acid, "分子构建失败")
            
            residue_smiles = Chem.MolToSmiles(residue_mol, canonical=True, isomericSmiles=True)
            
            # 2. 计算SMILES相似度
            smiles_similarity = self._calculate_smiles_similarity(residue_smiles, amino_acid.smiles)
            
            # 3. 计算分子指纹相似度  
            fingerprint_similarity = self._calculate_fingerprint_similarity(residue_mol, amino_acid.smiles)
            
            # 4. 计算3D结构相似度
            structure_3d_rmsd = self._calculate_3d_similarity(residue, amino_acid)
            
            # 5. 确定异构体类型
            isomer_type = self._classify_isomer_type(residue_smiles, amino_acid.smiles, 
                                                   smiles_similarity, structure_3d_rmsd)
            
            # 6. 计算综合置信度
            overall_confidence = self._calculate_confidence(
                smiles_similarity, fingerprint_similarity, structure_3d_rmsd, isomer_type
            )
            
            return IsomerMatch(
                amino_acid_id=amino_acid.id,
                amino_acid_name=amino_acid.name,
                smiles_similarity=smiles_similarity,
                structure_3d_rmsd=structure_3d_rmsd,
                fingerprint_similarity=fingerprint_similarity,
                overall_confidence=overall_confidence,
                isomer_type=isomer_type
            )
            
        except Exception as e:
            logging.error(f"异构体检测失败: {e}")
            return self._create_failed_match(amino_acid, str(e))
    
    def _coordinates_to_molecule(self, residue: ResidueInfo) -> Optional[Chem.Mol]:
        """从坐标生成RDKit分子对象"""
        try:
            # 创建分子
            mol = Chem.RWMol()
            
            # 添加原子
            for atom in residue.atoms:
                rd_atom = Chem.Atom(atom.element)
                mol.AddAtom(rd_atom)
            
            # 设置3D坐标
            conf = Chem.Conformer(len(residue.atoms))
            for i, atom in enumerate(residue.atoms):
                point = Point3D(atom.x, atom.y, atom.z)
                conf.SetAtomPosition(i, point)
            mol.AddConformer(conf)
            
            # 推断化学键（更robust的方式）
            try:
                rdDetermineBonds.DetermineBonds(mol, charge=0)
                Chem.SanitizeMol(mol)
            except:
                # 如果自动推断失败，使用距离推断
                try:
                    mol_with_bonds = self._add_bonds_by_distance(mol)
                    if mol_with_bonds is not None:
                        mol = mol_with_bonds
                    else:
                        return None
                except Exception as e2:
                    logging.warning(f"距离推断也失败: {e2}")
                    return None
            
            return mol
            
        except Exception as e:
            logging.warning(f"分子构建失败: {e}")
            return None
    
    def _add_bonds_by_distance(self, mol):
        """基于距离添加化学键"""
        try:
            rwmol = Chem.RWMol(mol)
            conf = mol.GetConformer()
            num_atoms = mol.GetNumAtoms()
            
            # 清除现有键
            for i in range(rwmol.GetNumBonds() - 1, -1, -1):
                rwmol.RemoveBond(rwmol.GetBondWithIdx(i).GetBeginAtomIdx(),
                                rwmol.GetBondWithIdx(i).GetEndAtomIdx())
            
            # 基于距离添加键
            for i in range(num_atoms):
                for j in range(i + 1, num_atoms):
                    atom1 = mol.GetAtomWithIdx(i)
                    atom2 = mol.GetAtomWithIdx(j)
                    
                    pos1 = conf.GetAtomPosition(i)
                    pos2 = conf.GetAtomPosition(j)
                    distance = pos1.Distance(pos2)
                    
                    # 检查是否应该连接
                    if self._should_bond(atom1.GetSymbol(), atom2.GetSymbol(), distance):
                        rwmol.AddBond(i, j, Chem.BondType.SINGLE)
            
            # 尝试清理分子
            try:
                Chem.SanitizeMol(rwmol)
                return rwmol
            except:
                # 如果清理失败，返回原分子
                return None
                
        except Exception as e:
            logging.warning(f"距离推断失败: {e}")
            return None
    
    def _should_bond(self, elem1: str, elem2: str, distance: float) -> bool:
        """判断两个原子是否应该连接"""
        # 获取共价半径
        radius1 = self.covalent_radii.get(elem1, 0.8)
        radius2 = self.covalent_radii.get(elem2, 0.8)
        
        # 最大键长
        max_distance = (radius1 + radius2) * 1.3
        min_distance = (radius1 + radius2) * 0.7
        
        return min_distance <= distance <= max_distance
    
    def _calculate_smiles_similarity(self, smiles1: str, smiles2: str) -> float:
        """计算SMILES相似度"""
        try:
            mol1 = Chem.MolFromSmiles(smiles1)
            mol2 = Chem.MolFromSmiles(smiles2)
            
            if mol1 is None or mol2 is None:
                return 0.0
            
            # 计算Tanimoto相似度
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)
            
            return DataStructs.TanimotoSimilarity(fp1, fp2)
            
        except Exception:
            return 0.0
    
    def _calculate_fingerprint_similarity(self, mol1: Chem.Mol, smiles2: str) -> float:
        """计算分子指纹相似度"""
        try:
            mol2 = Chem.MolFromSmiles(smiles2)
            if mol2 is None:
                return 0.0
            
            # Morgan指纹
            fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
            fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)
            
            return DataStructs.TanimotoSimilarity(fp1, fp2)
            
        except Exception:
            return 0.0
    
    def _calculate_3d_similarity(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> float:
        """计算3D结构相似度(RMSD)"""
        try:
            # 从SMILES生成参考构象
            ref_mol = Chem.MolFromSmiles(amino_acid.smiles)
            if ref_mol is None:
                return float('inf')
            
            # 生成3D构象
            ref_mol = Chem.AddHs(ref_mol)
            AllChem.EmbedMolecule(ref_mol)
            AllChem.OptimizeMoleculeConfs(ref_mol)
            
            # 获取坐标
            residue_coords = np.array([[atom.x, atom.y, atom.z] for atom in residue.atoms])
            
            if ref_mol.GetNumConformers() > 0:
                ref_coords = ref_mol.GetConformer().GetPositions()
                
                # 只比较重原子
                if len(residue_coords) == len(ref_coords):
                    # 使用Kabsch算法对齐
                    rmsd = self._kabsch_rmsd(residue_coords, ref_coords)
                    return rmsd
            
            return float('inf')
            
        except Exception as e:
            logging.warning(f"3D相似度计算失败: {e}")
            return float('inf')
    
    def _kabsch_rmsd(self, coords1: np.ndarray, coords2: np.ndarray) -> float:
        """使用Kabsch算法计算RMSD"""
        # 中心化坐标
        center1 = np.mean(coords1, axis=0)
        center2 = np.mean(coords2, axis=0)
        
        coords1_centered = coords1 - center1
        coords2_centered = coords2 - center2
        
        # 计算协方差矩阵
        H = coords1_centered.T @ coords2_centered
        
        # SVD分解
        U, S, Vt = np.linalg.svd(H)
        
        # 计算旋转矩阵
        R = Vt.T @ U.T
        
        # 确保是右手坐标系
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        
        # 应用旋转
        coords1_aligned = coords1_centered @ R.T
        
        # 计算RMSD
        diff = coords1_aligned - coords2_centered
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        
        return rmsd
    
    def _classify_isomer_type(self, smiles1: str, smiles2: str, 
                             smiles_similarity: float, rmsd: float) -> str:
        """分类异构体类型"""
        # 完全相同
        if smiles1 == smiles2:
            return "identical"
        
        # 立体异构体（分子式相同但立体化学不同）
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        
        if mol1 and mol2:
            # 去除立体化学信息比较
            smiles1_no_stereo = Chem.MolToSmiles(mol1, isomericSmiles=False)
            smiles2_no_stereo = Chem.MolToSmiles(mol2, isomericSmiles=False) 
            
            if smiles1_no_stereo == smiles2_no_stereo:
                return "stereoisomer"
            
            # 检查是否为结构异构体
            formula1 = Chem.rdMolDescriptors.CalcMolFormula(mol1)
            formula2 = Chem.rdMolDescriptors.CalcMolFormula(mol2)
            
            if formula1 == formula2:
                return "structural"
        
        # 构象异构体（3D结构相似但SMILES不同）
        if smiles_similarity > 0.9 and rmsd < self.rmsd_threshold:
            return "conformational"
        
        return "different"
    
    def _calculate_confidence(self, smiles_sim: float, fp_sim: float, 
                             rmsd: float, isomer_type: str) -> float:
        """计算综合置信度"""
        # 根据异构体类型调整权重
        if isomer_type == "identical":
            return 1.0
        elif isomer_type == "stereoisomer":
            # 立体异构体主要看3D结构
            structure_score = max(0, 1 - rmsd / self.rmsd_threshold)
            return 0.4 * smiles_sim + 0.3 * fp_sim + 0.3 * structure_score
        elif isomer_type == "structural":
            # 结构异构体主要看指纹相似度
            return 0.3 * smiles_sim + 0.5 * fp_sim + 0.2 * max(0, 1 - rmsd / 10)
        else:
            # 其他情况综合考虑
            structure_score = max(0, 1 - rmsd / self.rmsd_threshold) if rmsd != float('inf') else 0
            return 0.4 * smiles_sim + 0.4 * fp_sim + 0.2 * structure_score
    
    def _create_failed_match(self, amino_acid: AminoAcidInfo, error_msg: str) -> IsomerMatch:
        """创建失败的匹配结果"""
        return IsomerMatch(
            amino_acid_id=amino_acid.id,
            amino_acid_name=amino_acid.name,
            smiles_similarity=0.0,
            structure_3d_rmsd=float('inf'),
            fingerprint_similarity=0.0,
            overall_confidence=0.0,
            isomer_type="error"
        )
    
    def batch_detect(self, residue: ResidueInfo, 
                    amino_acids: List[AminoAcidInfo], 
                    min_confidence: float = 0.6) -> List[IsomerMatch]:
        """
        批量检测异构体
        
        Args:
            residue: 残基信息
            amino_acids: 候选氨基酸列表
            min_confidence: 最小置信度阈值
            
        Returns:
            匹配结果列表，按置信度排序
        """
        matches = []
        
        for amino_acid in amino_acids:
            match = self.detect_isomer(residue, amino_acid)
            if match.overall_confidence >= min_confidence:
                matches.append(match)
        
        # 按置信度排序
        matches.sort(key=lambda x: x.overall_confidence, reverse=True)
        
        return matches
    
    def find_best_match(self, residue: ResidueInfo, 
                       amino_acids: List[AminoAcidInfo]) -> Optional[IsomerMatch]:
        """
        找到最佳匹配
        
        Args:
            residue: 残基信息
            amino_acids: 候选氨基酸列表
            
        Returns:
            最佳匹配结果
        """
        matches = self.batch_detect(residue, amino_acids, min_confidence=0.0)
        
        if matches and matches[0].overall_confidence > 0.0:
            return matches[0]
        
        return None


class EnhancedIsomerDetector(IsomerDetector):
    """
    增强的异构体检测器
    添加更多的检测策略用于提高准确率
    """
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
    def detect_isomer(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> IsomerMatch:
        """增强的异构体检测"""
        # 调用基础检测
        match = super().detect_isomer(residue, amino_acid)
        
        # 添加额外验证
        if match.overall_confidence > 0.8:
            # 高置信度结果，进行额外验证
            match = self._verify_high_confidence_match(match, residue, amino_acid)
        
        return match
    
    def _verify_high_confidence_match(self, match: IsomerMatch, 
                                     residue: ResidueInfo, 
                                     amino_acid: AminoAcidInfo) -> IsomerMatch:
        """验证高置信度匹配"""
        try:
            # 检查原子数量是否匹配
            if len(residue.atoms) != len(amino_acid.atom_composition):
                # 原子数不匹配，降低置信度
                match.overall_confidence *= 0.8
                
            # 检查分子式是否匹配
            if residue.molecular_formula != amino_acid.molecular_formula:
                match.overall_confidence *= 0.7
                
            return match
            
        except Exception as e:
            logging.warning(f"高置信度匹配验证失败: {e}")
            return match