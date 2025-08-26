"""
增强的3D结构验证模块
支持从SMILES生成标准构象、药效团匹配、立体化学验证
"""

import numpy as np
from typing import Dict, Any, List, Optional, Tuple, Set
from dataclasses import dataclass
import logging
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment

try:
    from rdkit import Chem
    from rdkit.Chem import rdDistGeom, AllChem, rdMolAlign
    from rdkit.Geometry import Point3D
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit未安装，3D结构验证功能受限")

from .structure_3d import Structure3DVerifier
from ..models import ResidueInfo, AminoAcidInfo, VerificationMethod


@dataclass
class PharmacophorePoint:
    """药效团点"""
    point_type: str  # "donor", "acceptor", "hydrophobic", "aromatic", "positive", "negative"
    coordinates: Tuple[float, float, float]
    radius: float = 1.5
    importance: float = 1.0


@dataclass 
class ConformerInfo:
    """构象信息"""
    coordinates: np.ndarray  # 原子坐标 (N, 3)
    energy: float = 0.0
    confidence: float = 1.0
    atom_types: List[str] = None


class Enhanced3DVerifier(Structure3DVerifier):
    """增强的3D结构验证器"""
    
    def __init__(self, threshold: float = 0.5):
        super().__init__(threshold)
        self.conformer_cache = {}  # 缓存生成的构象
        
    def calculate_score(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> float:
        """
        增强的3D结构验证分数计算
        
        验证策略：
        1. 从SMILES生成标准构象
        2. Kabsch RMSD比较
        3. 药效团匹配
        4. 立体化学验证
        5. 综合评分
        """
        try:
            # 获取残基坐标
            residue_coords = self._extract_residue_coordinates(residue)
            if residue_coords is None or len(residue_coords) < 3:
                return super().calculate_score(residue, amino_acid)
            
            # 生成标准构象
            standard_conformers = self._generate_standard_conformers(amino_acid)
            if not standard_conformers:
                return super().calculate_score(residue, amino_acid)
            
            # 计算多种3D相似性指标
            scores = []
            
            for conformer in standard_conformers:
                score_dict = self._calculate_comprehensive_3d_similarity(
                    residue_coords, conformer, residue, amino_acid
                )
                scores.append(score_dict)
            
            # 选择最佳匹配
            best_score = max(scores, key=lambda x: x['total_score'])
            
            return best_score['total_score']
            
        except Exception as e:
            logging.error(f"增强3D验证失败: {e}")
            return super().calculate_score(residue, amino_acid)
    
    def _generate_standard_conformers(self, amino_acid: AminoAcidInfo, 
                                    num_conformers: int = 5) -> List[ConformerInfo]:
        """
        从SMILES生成标准构象
        
        Args:
            amino_acid: 氨基酸信息
            num_conformers: 生成构象数量
            
        Returns:
            构象信息列表
        """
        if not amino_acid.smiles or not RDKIT_AVAILABLE:
            return []
        
        cache_key = f"{amino_acid.id}_{num_conformers}"
        if cache_key in self.conformer_cache:
            return self.conformer_cache[cache_key]
        
        try:
            conformers = []
            mol = Chem.MolFromSmiles(amino_acid.smiles)
            
            if mol is None:
                return []
            
            # 添加氢原子
            mol = Chem.AddHs(mol)
            
            # 生成多个构象
            conf_ids = rdDistGeom.EmbedMultipleConfs(
                mol, 
                numConfs=num_conformers,
                randomSeed=42,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True
            )
            
            if not conf_ids:
                # 如果失败，尝试更简单的方法
                conf_id = rdDistGeom.EmbedMolecule(mol, randomSeed=42)
                if conf_id >= 0:
                    conf_ids = [conf_id]
                else:
                    return []
            
            # 优化构象
            for conf_id in conf_ids:
                try:
                    # MMFF优化
                    AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
                    
                    # 提取坐标
                    conf = mol.GetConformer(conf_id)
                    coords = []
                    atom_types = []
                    
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() != 'H':  # 只考虑重原子
                            pos = conf.GetAtomPosition(atom.GetIdx())
                            coords.append([pos.x, pos.y, pos.z])
                            atom_types.append(atom.GetSymbol())
                    
                    if coords:
                        conformer_info = ConformerInfo(
                            coordinates=np.array(coords),
                            energy=0.0,  # 可以计算能量
                            confidence=1.0,
                            atom_types=atom_types
                        )
                        conformers.append(conformer_info)
                        
                except Exception as e:
                    logging.warning(f"构象优化失败 (conf_id={conf_id}): {e}")
                    continue
            
            # 缓存结果
            self.conformer_cache[cache_key] = conformers
            return conformers
            
        except Exception as e:
            logging.error(f"构象生成失败 {amino_acid.smiles}: {e}")
            return []
    
    def _calculate_comprehensive_3d_similarity(self, residue_coords: np.ndarray, 
                                             conformer: ConformerInfo,
                                             residue: ResidueInfo, 
                                             amino_acid: AminoAcidInfo) -> Dict[str, float]:
        """
        综合3D相似性计算
        
        Args:
            residue_coords: 残基坐标
            conformer: 标准构象
            residue: 残基信息
            amino_acid: 氨基酸信息
            
        Returns:
            相似性分数字典
        """
        scores = {
            'kabsch_rmsd_score': 0.0,
            'pharmacophore_score': 0.0,
            'shape_similarity_score': 0.0,
            'atom_mapping_score': 0.0,
            'total_score': 0.0
        }
        
        try:
            standard_coords = conformer.coordinates
            
            # 1. Kabsch RMSD比较（核心）
            rmsd = self._calculate_kabsch_rmsd(residue_coords, standard_coords)
            if rmsd is not None:
                atom_count = len(residue_coords)
                scores['kabsch_rmsd_score'] = self._rmsd_to_score_dynamic(rmsd, atom_count)
            
            # 2. 药效团匹配
            pharmacophore_score = self._calculate_pharmacophore_similarity(
                residue_coords, standard_coords, residue, amino_acid
            )
            scores['pharmacophore_score'] = pharmacophore_score
            
            # 3. 分子形状相似性
            shape_score = self._calculate_shape_similarity(residue_coords, standard_coords)
            scores['shape_similarity_score'] = shape_score
            
            # 4. 原子匹配评分
            mapping_score = self._calculate_atom_mapping_score(
                residue_coords, standard_coords, residue, conformer
            )
            scores['atom_mapping_score'] = mapping_score
            
            # 综合评分（加权平均）
            weights = {
                'kabsch_rmsd_score': 0.4,
                'pharmacophore_score': 0.25,
                'shape_similarity_score': 0.2,
                'atom_mapping_score': 0.15
            }
            
            total_score = sum(scores[key] * weights[key] for key in weights.keys())
            scores['total_score'] = max(0.0, min(1.0, total_score))
            
        except Exception as e:
            logging.error(f"综合3D相似性计算失败: {e}")
            scores['total_score'] = 0.0
        
        return scores
    
    def _calculate_pharmacophore_similarity(self, residue_coords: np.ndarray,
                                          standard_coords: np.ndarray,
                                          residue: ResidueInfo,
                                          amino_acid: AminoAcidInfo) -> float:
        """
        计算药效团相似性
        
        Args:
            residue_coords: 残基坐标
            standard_coords: 标准坐标
            residue: 残基信息
            amino_acid: 氨基酸信息
            
        Returns:
            药效团相似性分数
        """
        try:
            # 提取药效团特征点
            residue_pharmacophores = self._extract_pharmacophores_from_residue(residue)
            standard_pharmacophores = self._extract_pharmacophores_from_smiles(amino_acid.smiles)
            
            if not residue_pharmacophores or not standard_pharmacophores:
                return 0.5  # 无法提取特征时返回中性分数
            
            # 计算药效团点之间的最佳匹配
            return self._match_pharmacophore_points(residue_pharmacophores, standard_pharmacophores)
            
        except Exception as e:
            logging.error(f"药效团相似性计算失败: {e}")
            return 0.5
    
    def _extract_pharmacophores_from_residue(self, residue: ResidueInfo) -> List[PharmacophorePoint]:
        """从残基提取药效团特征点"""
        pharmacophores = []
        
        try:
            for atom in residue.atoms:
                element = atom.get('element', '').strip()
                coords = (
                    float(atom.get('x', 0)),
                    float(atom.get('y', 0)), 
                    float(atom.get('z', 0))
                )
                
                # 简化的药效团分类
                if element == 'O':
                    # 氧原子可能是氢键受体
                    pharmacophores.append(PharmacophorePoint(
                        point_type='acceptor',
                        coordinates=coords,
                        importance=0.9
                    ))
                elif element == 'N':
                    # 氮原子可能是氢键供体或受体
                    pharmacophores.append(PharmacophorePoint(
                        point_type='donor',
                        coordinates=coords,
                        importance=0.8
                    ))
                elif element == 'C':
                    # 碳原子可能是疏水中心
                    pharmacophores.append(PharmacophorePoint(
                        point_type='hydrophobic',
                        coordinates=coords,
                        importance=0.3
                    ))
            
            return pharmacophores
            
        except Exception as e:
            logging.error(f"残基药效团提取失败: {e}")
            return []
    
    def _extract_pharmacophores_from_smiles(self, smiles: str) -> List[PharmacophorePoint]:
        """从SMILES提取药效团特征点"""
        if not RDKIT_AVAILABLE or not smiles:
            return []
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return []
            
            # 生成3D构象
            mol = Chem.AddHs(mol)
            conf_id = rdDistGeom.EmbedMolecule(mol, randomSeed=42)
            if conf_id < 0:
                return []
            
            AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)
            conf = mol.GetConformer(conf_id)
            
            pharmacophores = []
            
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == 'H':
                    continue
                    
                pos = conf.GetAtomPosition(atom.GetIdx())
                coords = (pos.x, pos.y, pos.z)
                element = atom.GetSymbol()
                
                if element == 'O':
                    pharmacophores.append(PharmacophorePoint(
                        point_type='acceptor',
                        coordinates=coords,
                        importance=0.9
                    ))
                elif element == 'N':
                    pharmacophores.append(PharmacophorePoint(
                        point_type='donor',
                        coordinates=coords,
                        importance=0.8
                    ))
                elif element == 'C':
                    pharmacophores.append(PharmacophorePoint(
                        point_type='hydrophobic',
                        coordinates=coords,
                        importance=0.3
                    ))
            
            return pharmacophores
            
        except Exception as e:
            logging.error(f"SMILES药效团提取失败: {e}")
            return []
    
    def _match_pharmacophore_points(self, points1: List[PharmacophorePoint],
                                  points2: List[PharmacophorePoint]) -> float:
        """匹配药效团特征点"""
        if not points1 or not points2:
            return 0.5
        
        try:
            # 按类型分组
            groups1 = {}
            groups2 = {}
            
            for point in points1:
                if point.point_type not in groups1:
                    groups1[point.point_type] = []
                groups1[point.point_type].append(point)
            
            for point in points2:
                if point.point_type not in groups2:
                    groups2[point.point_type] = []
                groups2[point.point_type].append(point)
            
            # 计算各类型的匹配分数
            type_scores = []
            
            for point_type in set(groups1.keys()) | set(groups2.keys()):
                if point_type in groups1 and point_type in groups2:
                    # 计算同类型点之间的距离匹配
                    coords1 = np.array([p.coordinates for p in groups1[point_type]])
                    coords2 = np.array([p.coordinates for p in groups2[point_type]])
                    
                    # 使用匈牙利算法找最优匹配
                    distances = cdist(coords1, coords2)
                    if min(distances.shape) > 0:
                        row_ind, col_ind = linear_sum_assignment(distances)
                        matched_distances = distances[row_ind, col_ind]
                        
                        # 转换距离为相似性分数
                        similarity_scores = np.exp(-matched_distances / 2.0)  # 高斯相似性
                        type_score = np.mean(similarity_scores)
                        type_scores.append(type_score)
            
            return np.mean(type_scores) if type_scores else 0.5
            
        except Exception as e:
            logging.error(f"药效团匹配失败: {e}")
            return 0.5
    
    def _calculate_shape_similarity(self, coords1: np.ndarray, coords2: np.ndarray) -> float:
        """
        计算分子形状相似性
        
        Args:
            coords1: 第一组坐标
            coords2: 第二组坐标
            
        Returns:
            形状相似性分数
        """
        try:
            if len(coords1) != len(coords2):
                # 原子数不同时的处理
                return self._calculate_partial_shape_similarity(coords1, coords2)
            
            # 使用惯性矩阵比较分子形状
            inertia1 = self._calculate_inertia_tensor(coords1)
            inertia2 = self._calculate_inertia_tensor(coords2)
            
            # 计算惯性矩的特征值
            eigenvals1 = np.sort(np.linalg.eigvals(inertia1))
            eigenvals2 = np.sort(np.linalg.eigvals(inertia2))
            
            # 归一化特征值
            eigenvals1 = eigenvals1 / np.sum(eigenvals1) if np.sum(eigenvals1) > 0 else eigenvals1
            eigenvals2 = eigenvals2 / np.sum(eigenvals2) if np.sum(eigenvals2) > 0 else eigenvals2
            
            # 计算特征值相似性
            similarity = 1.0 - np.linalg.norm(eigenvals1 - eigenvals2) / 2.0
            
            return max(0.0, min(1.0, similarity))
            
        except Exception as e:
            logging.error(f"形状相似性计算失败: {e}")
            return 0.5
    
    def _calculate_inertia_tensor(self, coords: np.ndarray) -> np.ndarray:
        """计算惯性张量"""
        # 将坐标中心化
        centered_coords = coords - np.mean(coords, axis=0)
        
        # 计算惯性张量
        I = np.zeros((3, 3))
        for coord in centered_coords:
            x, y, z = coord
            I[0, 0] += y*y + z*z
            I[1, 1] += x*x + z*z
            I[2, 2] += x*x + y*y
            I[0, 1] = I[1, 0] = I[0, 1] - x*y
            I[0, 2] = I[2, 0] = I[0, 2] - x*z
            I[1, 2] = I[2, 1] = I[1, 2] - y*z
        
        return I
    
    def _calculate_partial_shape_similarity(self, coords1: np.ndarray, coords2: np.ndarray) -> float:
        """计算不同原子数分子的形状相似性"""
        try:
            # 使用分子的几何特征比较
            # 1. 计算质心距离
            center1 = np.mean(coords1, axis=0)
            center2 = np.mean(coords2, axis=0)
            
            # 2. 计算平均原子距离
            avg_dist1 = np.mean(np.linalg.norm(coords1 - center1, axis=1))
            avg_dist2 = np.mean(np.linalg.norm(coords2 - center2, axis=1))
            
            # 3. 距离相似性
            if avg_dist1 > 0 and avg_dist2 > 0:
                dist_ratio = min(avg_dist1, avg_dist2) / max(avg_dist1, avg_dist2)
            else:
                dist_ratio = 0.5
            
            # 4. 原子数相似性
            size_ratio = min(len(coords1), len(coords2)) / max(len(coords1), len(coords2))
            
            # 综合评分
            similarity = (dist_ratio * 0.6 + size_ratio * 0.4)
            
            return max(0.0, min(1.0, similarity))
            
        except Exception as e:
            logging.error(f"部分形状相似性计算失败: {e}")
            return 0.5
    
    def _calculate_atom_mapping_score(self, residue_coords: np.ndarray,
                                    standard_coords: np.ndarray,
                                    residue: ResidueInfo,
                                    conformer: ConformerInfo) -> float:
        """
        计算原子匹配评分
        
        Args:
            residue_coords: 残基坐标
            standard_coords: 标准坐标
            residue: 残基信息
            conformer: 构象信息
            
        Returns:
            原子匹配分数
        """
        try:
            # 提取原子类型
            residue_elements = []
            for atom in residue.atoms:
                element = atom.get('element', '').strip()
                if element != 'H':  # 只考虑重原子
                    residue_elements.append(element)
            
            conformer_elements = conformer.atom_types if conformer.atom_types else []
            
            if len(residue_elements) != len(conformer_elements):
                # 原子数不匹配时的处理
                common_elements = min(len(residue_elements), len(conformer_elements))
                if common_elements == 0:
                    return 0.0
                
                # 计算元素类型匹配度
                residue_composition = {}
                conformer_composition = {}
                
                for elem in residue_elements:
                    residue_composition[elem] = residue_composition.get(elem, 0) + 1
                
                for elem in conformer_elements:
                    conformer_composition[elem] = conformer_composition.get(elem, 0) + 1
                
                # 计算组成相似性
                all_elements = set(residue_composition.keys()) | set(conformer_composition.keys())
                similarity_sum = 0
                for elem in all_elements:
                    count1 = residue_composition.get(elem, 0)
                    count2 = conformer_composition.get(elem, 0)
                    if count1 + count2 > 0:
                        similarity_sum += min(count1, count2) / max(count1, count2)
                
                return similarity_sum / len(all_elements) if all_elements else 0.0
            
            else:
                # 原子数匹配时，计算类型匹配度
                matches = sum(1 for e1, e2 in zip(residue_elements, conformer_elements) if e1 == e2)
                return matches / len(residue_elements)
            
        except Exception as e:
            logging.error(f"原子匹配评分计算失败: {e}")
            return 0.0
    
    def _get_standard_structure_coordinates(self, amino_acid: AminoAcidInfo) -> Optional[np.ndarray]:
        """
        获取氨基酸的标准结构坐标（重写父类方法）
        
        Args:
            amino_acid: 氨基酸信息
            
        Returns:
            标准结构坐标，失败返回None
        """
        conformers = self._generate_standard_conformers(amino_acid, num_conformers=1)
        
        if conformers:
            return conformers[0].coordinates
        
        return None
    
    def get_detailed_analysis(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """
        获取详细的3D结构分析
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
            
        Returns:
            详细分析结果
        """
        analysis = {
            'conformers_generated': 0,
            'best_conformer_score': 0.0,
            'all_conformer_scores': [],
            'pharmacophore_match': False,
            'shape_similarity': 0.0,
            'rmsd_values': [],
            'error_messages': []
        }
        
        try:
            residue_coords = self._extract_residue_coordinates(residue)
            if residue_coords is None:
                analysis['error_messages'].append("无法提取残基坐标")
                return analysis
            
            # 生成构象
            conformers = self._generate_standard_conformers(amino_acid)
            analysis['conformers_generated'] = len(conformers)
            
            if not conformers:
                analysis['error_messages'].append("无法生成标准构象")
                return analysis
            
            # 分析每个构象
            for i, conformer in enumerate(conformers):
                try:
                    score_dict = self._calculate_comprehensive_3d_similarity(
                        residue_coords, conformer, residue, amino_acid
                    )
                    
                    analysis['all_conformer_scores'].append({
                        'conformer_id': i,
                        'scores': score_dict
                    })
                    
                    # 记录RMSD
                    rmsd = self._calculate_kabsch_rmsd(residue_coords, conformer.coordinates)
                    if rmsd is not None:
                        analysis['rmsd_values'].append(rmsd)
                    
                except Exception as e:
                    analysis['error_messages'].append(f"构象{i}分析失败: {e}")
            
            # 找到最佳匹配
            if analysis['all_conformer_scores']:
                best_match = max(analysis['all_conformer_scores'], 
                               key=lambda x: x['scores']['total_score'])
                analysis['best_conformer_score'] = best_match['scores']['total_score']
                analysis['pharmacophore_match'] = best_match['scores']['pharmacophore_score'] > 0.7
                analysis['shape_similarity'] = best_match['scores']['shape_similarity_score']
            
        except Exception as e:
            analysis['error_messages'].append(f"详细分析失败: {e}")
        
        return analysis