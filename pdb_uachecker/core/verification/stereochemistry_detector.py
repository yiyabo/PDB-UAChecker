"""
立体化学检测器
识别手性中心、确定R/S构型、比较立体异构体
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Set, Any
from dataclasses import dataclass
import logging

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, AllChem, rdDistGeom
    from rdkit.Geometry import Point3D
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit未安装，立体化学检测功能受限")

from ..models import ResidueInfo, AminoAcidInfo


@dataclass
class ChiralCenter:
    """手性中心信息"""
    atom_idx: int
    atom_symbol: str
    configuration: Optional[str] = None  # 'R', 'S', 或 None
    coordinates: Optional[Tuple[float, float, float]] = None
    neighbors: List[int] = None
    cip_priority: List[int] = None


@dataclass
class StereochemistryInfo:
    """立体化学信息"""
    chiral_centers: List[ChiralCenter]
    total_chiral_centers: int
    stereoisomer_type: str = "unknown"  # "enantiomer", "diastereomer", "meso", "achiral"
    confidence: float = 0.0


class StereochemistryDetector:
    """立体化学检测器"""
    
    def __init__(self):
        """初始化立体化学检测器"""
        self.atomic_numbers = {
            'H': 1, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'P': 15, 'S': 16,
            'Cl': 17, 'Br': 35, 'I': 53, 'Se': 34
        }
        
    def detect_stereochemistry_from_residue(self, residue: ResidueInfo) -> StereochemistryInfo:
        """
        从残基信息检测立体化学
        
        Args:
            residue: 残基信息
            
        Returns:
            立体化学信息
        """
        try:
            if RDKIT_AVAILABLE:
                return self._detect_with_rdkit_from_coords(residue)
            else:
                return self._detect_with_basic_algorithm(residue)
                
        except Exception as e:
            logging.error(f"残基立体化学检测失败: {e}")
            return StereochemistryInfo(chiral_centers=[], total_chiral_centers=0)
    
    def detect_stereochemistry_from_smiles(self, smiles: str) -> StereochemistryInfo:
        """
        从SMILES检测立体化学
        
        Args:
            smiles: SMILES字符串
            
        Returns:
            立体化学信息
        """
        if not RDKIT_AVAILABLE or not smiles:
            return StereochemistryInfo(chiral_centers=[], total_chiral_centers=0)
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return StereochemistryInfo(chiral_centers=[], total_chiral_centers=0)
            
            return self._analyze_molecule_stereochemistry(mol)
            
        except Exception as e:
            logging.error(f"SMILES立体化学检测失败: {e}")
            return StereochemistryInfo(chiral_centers=[], total_chiral_centers=0)
    
    def _detect_with_rdkit_from_coords(self, residue: ResidueInfo) -> StereochemistryInfo:
        """
        使用RDKit从坐标检测立体化学
        
        Args:
            residue: 残基信息
            
        Returns:
            立体化学信息
        """
        try:
            # 创建分子对象
            mol = self._create_molecule_from_residue(residue)
            if mol is None:
                return StereochemistryInfo(chiral_centers=[], total_chiral_centers=0)
            
            return self._analyze_molecule_stereochemistry(mol)
            
        except Exception as e:
            logging.error(f"RDKit坐标立体化学检测失败: {e}")
            return StereochemistryInfo(chiral_centers=[], total_chiral_centers=0)
    
    def _create_molecule_from_residue(self, residue: ResidueInfo):
        """从残基信息创建RDKit分子对象"""
        try:
            mol = Chem.RWMol()
            atom_map = {}
            
            # 添加原子
            for i, atom in enumerate(residue.atoms):
                element = atom.get('element', '').strip()
                if not element:
                    continue
                
                rd_atom = Chem.Atom(element)
                atom_idx = mol.AddAtom(rd_atom)
                atom_map[i] = atom_idx
            
            # 设置坐标
            if len(atom_map) > 0:
                conf = Chem.Conformer(len(atom_map))
                for i, atom in enumerate(residue.atoms):
                    if i in atom_map:
                        x = float(atom.get('x', 0))
                        y = float(atom.get('y', 0))
                        z = float(atom.get('z', 0))
                        point = Point3D(x, y, z)
                        conf.SetAtomPosition(atom_map[i], point)
                
                mol.AddConformer(conf)
            
            # 推断化学键（简化版）
            mol = self._infer_bonds_simple(mol)
            
            return mol
            
        except Exception as e:
            logging.error(f"从残基创建分子失败: {e}")
            return None
    
    def _infer_bonds_simple(self, mol):
        """简单的化学键推断"""
        try:
            # 使用RDKit的自动键推断
            from rdkit.Chem import rdDetermineBonds
            rdDetermineBonds.DetermineBonds(mol)
            
            # 清理和验证分子
            Chem.SanitizeMol(mol)
            
            return mol
            
        except Exception as e:
            logging.warning(f"简单键推断失败: {e}")
            
            # 如果失败，返回只有原子没有键的分子
            try:
                mol_copy = Chem.RWMol(mol)
                return mol_copy
            except:
                return None
    
    def _analyze_molecule_stereochemistry(self, mol) -> StereochemistryInfo:
        """
        分析分子的立体化学
        
        Args:
            mol: RDKit分子对象
            
        Returns:
            立体化学信息
        """
        try:
            chiral_centers = []
            
            # 检测手性中心
            chiral_centers_info = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            
            for center_idx, chirality_code in chiral_centers_info:
                atom = mol.GetAtomWithIdx(center_idx)
                
                # 获取坐标
                coords = None
                if mol.GetNumConformers() > 0:
                    conf = mol.GetConformer()
                    pos = conf.GetAtomPosition(center_idx)
                    coords = (pos.x, pos.y, pos.z)
                
                # 解析手性代码
                configuration = None
                if chirality_code == 'R':
                    configuration = 'R'
                elif chirality_code == 'S':
                    configuration = 'S'
                
                # 获取邻近原子
                neighbors = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()]
                
                # 计算CIP优先级
                cip_priority = self._calculate_cip_priority(mol, center_idx, neighbors)
                
                chiral_center = ChiralCenter(
                    atom_idx=center_idx,
                    atom_symbol=atom.GetSymbol(),
                    configuration=configuration,
                    coordinates=coords,
                    neighbors=neighbors,
                    cip_priority=cip_priority
                )
                
                chiral_centers.append(chiral_center)
            
            # 确定立体异构体类型
            stereoisomer_type = self._determine_stereoisomer_type(chiral_centers)
            
            # 计算置信度
            confidence = self._calculate_confidence(mol, chiral_centers)
            
            return StereochemistryInfo(
                chiral_centers=chiral_centers,
                total_chiral_centers=len(chiral_centers),
                stereoisomer_type=stereoisomer_type,
                confidence=confidence
            )
            
        except Exception as e:
            logging.error(f"分子立体化学分析失败: {e}")
            return StereochemistryInfo(chiral_centers=[], total_chiral_centers=0)
    
    def _calculate_cip_priority(self, mol, center_idx: int, neighbors: List[int]) -> List[int]:
        """
        计算CIP优先级规则的原子优先级
        
        Args:
            mol: 分子对象
            center_idx: 手性中心索引
            neighbors: 邻近原子索引列表
            
        Returns:
            按CIP优先级排序的原子索引列表
        """
        try:
            priority_list = []
            
            for neighbor_idx in neighbors:
                neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
                
                # 计算优先级（简化版CIP规则）
                priority = self._calculate_atom_priority(mol, neighbor_idx, visited=set([center_idx]))
                priority_list.append((neighbor_idx, priority))
            
            # 按优先级排序（高优先级在前）
            priority_list.sort(key=lambda x: x[1], reverse=True)
            
            return [idx for idx, _ in priority_list]
            
        except Exception as e:
            logging.error(f"CIP优先级计算失败: {e}")
            return neighbors
    
    def _calculate_atom_priority(self, mol, atom_idx: int, visited: Set[int], depth: int = 0) -> float:
        """
        计算原子的CIP优先级（递归）
        
        Args:
            mol: 分子对象
            atom_idx: 原子索引
            visited: 已访问的原子集合
            depth: 递归深度
            
        Returns:
            优先级分数
        """
        if depth > 3 or atom_idx in visited:  # 防止无限递归
            return 0.0
        
        atom = mol.GetAtomWithIdx(atom_idx)
        visited_copy = visited.copy()
        visited_copy.add(atom_idx)
        
        # 基础优先级：原子序数
        base_priority = self.atomic_numbers.get(atom.GetSymbol(), 0)
        
        # 考虑邻近原子的影响（递归）
        neighbor_priorities = []
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited:
                neighbor_priority = self._calculate_atom_priority(mol, neighbor_idx, visited_copy, depth + 1)
                neighbor_priorities.append(neighbor_priority)
        
        # 计算综合优先级
        if neighbor_priorities:
            avg_neighbor_priority = sum(neighbor_priorities) / len(neighbor_priorities)
            total_priority = base_priority + avg_neighbor_priority * 0.1
        else:
            total_priority = base_priority
        
        return total_priority
    
    def _determine_stereoisomer_type(self, chiral_centers: List[ChiralCenter]) -> str:
        """
        确定立体异构体类型
        
        Args:
            chiral_centers: 手性中心列表
            
        Returns:
            立体异构体类型
        """
        if not chiral_centers:
            return "achiral"
        
        if len(chiral_centers) == 1:
            return "enantiomer"
        
        # 多个手性中心的情况
        configured_centers = [c for c in chiral_centers if c.configuration is not None]
        
        if len(configured_centers) == len(chiral_centers):
            # 所有手性中心都有明确构型
            return "diastereomer"
        else:
            # 部分手性中心构型不明
            return "unknown"
    
    def _calculate_confidence(self, mol, chiral_centers: List[ChiralCenter]) -> float:
        """
        计算立体化学检测的置信度
        
        Args:
            mol: 分子对象
            chiral_centers: 手性中心列表
            
        Returns:
            置信度分数
        """
        try:
            if not chiral_centers:
                return 1.0  # 没有手性中心时置信度为1
            
            # 基于已确定构型的手性中心比例
            configured_centers = [c for c in chiral_centers if c.configuration is not None]
            configuration_ratio = len(configured_centers) / len(chiral_centers)
            
            # 基于分子的3D坐标完整性
            coord_completeness = 1.0
            if mol.GetNumConformers() > 0:
                conf = mol.GetConformer()
                valid_coords = 0
                for center in chiral_centers:
                    pos = conf.GetAtomPosition(center.atom_idx)
                    if abs(pos.x) + abs(pos.y) + abs(pos.z) > 1e-6:  # 非零坐标
                        valid_coords += 1
                coord_completeness = valid_coords / len(chiral_centers)
            else:
                coord_completeness = 0.5  # 没有3D坐标时降低置信度
            
            # 综合置信度
            confidence = (configuration_ratio * 0.7 + coord_completeness * 0.3)
            
            return max(0.0, min(1.0, confidence))
            
        except Exception as e:
            logging.error(f"置信度计算失败: {e}")
            return 0.0
    
    def _detect_with_basic_algorithm(self, residue: ResidueInfo) -> StereochemistryInfo:
        """
        使用基础算法检测立体化学（当RDKit不可用时）
        
        Args:
            residue: 残基信息
            
        Returns:
            立体化学信息
        """
        try:
            # 基础算法：基于几何分析识别可能的手性中心
            potential_chiral_centers = []
            
            for i, atom in enumerate(residue.atoms):
                element = atom.get('element', '').strip()
                
                # 只考虑碳原子作为潜在手性中心
                if element == 'C':
                    # 计算该原子与其他原子的距离
                    neighbors = self._find_neighbors_basic(residue, i)
                    
                    if len(neighbors) == 4:  # 四面体碳原子
                        coords = (
                            float(atom.get('x', 0)),
                            float(atom.get('y', 0)),
                            float(atom.get('z', 0))
                        )
                        
                        chiral_center = ChiralCenter(
                            atom_idx=i,
                            atom_symbol=element,
                            configuration=None,  # 无法确定R/S
                            coordinates=coords,
                            neighbors=neighbors
                        )
                        
                        potential_chiral_centers.append(chiral_center)
            
            return StereochemistryInfo(
                chiral_centers=potential_chiral_centers,
                total_chiral_centers=len(potential_chiral_centers),
                stereoisomer_type="unknown" if potential_chiral_centers else "achiral",
                confidence=0.3  # 基础算法置信度较低
            )
            
        except Exception as e:
            logging.error(f"基础立体化学检测失败: {e}")
            return StereochemistryInfo(chiral_centers=[], total_chiral_centers=0)
    
    def _find_neighbors_basic(self, residue: ResidueInfo, atom_idx: int, 
                            max_distance: float = 2.0) -> List[int]:
        """
        基础算法寻找邻近原子
        
        Args:
            residue: 残基信息
            atom_idx: 原子索引
            max_distance: 最大键长距离
            
        Returns:
            邻近原子索引列表
        """
        neighbors = []
        
        if atom_idx >= len(residue.atoms):
            return neighbors
        
        center_atom = residue.atoms[atom_idx]
        center_coords = np.array([
            float(center_atom.get('x', 0)),
            float(center_atom.get('y', 0)),
            float(center_atom.get('z', 0))
        ])
        
        for i, other_atom in enumerate(residue.atoms):
            if i == atom_idx:
                continue
            
            other_coords = np.array([
                float(other_atom.get('x', 0)),
                float(other_atom.get('y', 0)),
                float(other_atom.get('z', 0))
            ])
            
            distance = np.linalg.norm(center_coords - other_coords)
            
            if distance <= max_distance:
                neighbors.append(i)
        
        return neighbors
    
    def compare_stereochemistry(self, stereo1: StereochemistryInfo, 
                              stereo2: StereochemistryInfo) -> Dict[str, Any]:
        """
        比较两个立体化学信息
        
        Args:
            stereo1: 第一个立体化学信息
            stereo2: 第二个立体化学信息
            
        Returns:
            比较结果
        """
        comparison = {
            'match_type': 'unknown',
            'similarity_score': 0.0,
            'chiral_center_match': False,
            'configuration_match': False,
            'details': {}
        }
        
        try:
            # 比较手性中心数量
            if stereo1.total_chiral_centers != stereo2.total_chiral_centers:
                comparison['match_type'] = 'different_structure'
                comparison['similarity_score'] = 0.0
                comparison['details']['chiral_center_count_diff'] = abs(
                    stereo1.total_chiral_centers - stereo2.total_chiral_centers
                )
                return comparison
            
            if stereo1.total_chiral_centers == 0 and stereo2.total_chiral_centers == 0:
                comparison['match_type'] = 'both_achiral'
                comparison['similarity_score'] = 1.0
                comparison['chiral_center_match'] = True
                return comparison
            
            # 比较手性中心配置
            config_matches = 0
            total_comparisons = 0
            
            for center1 in stereo1.chiral_centers:
                for center2 in stereo2.chiral_centers:
                    if (center1.atom_symbol == center2.atom_symbol and
                        center1.configuration is not None and
                        center2.configuration is not None):
                        
                        total_comparisons += 1
                        if center1.configuration == center2.configuration:
                            config_matches += 1
            
            if total_comparisons > 0:
                config_similarity = config_matches / total_comparisons
                comparison['similarity_score'] = config_similarity
                comparison['configuration_match'] = config_matches == total_comparisons
                
                if config_matches == total_comparisons:
                    comparison['match_type'] = 'same_configuration'
                elif config_matches == 0:
                    comparison['match_type'] = 'opposite_configuration'
                else:
                    comparison['match_type'] = 'mixed_configuration'
            else:
                comparison['similarity_score'] = 0.5
                comparison['match_type'] = 'insufficient_data'
            
            comparison['chiral_center_match'] = True
            comparison['details']['config_matches'] = config_matches
            comparison['details']['total_comparisons'] = total_comparisons
            
        except Exception as e:
            logging.error(f"立体化学比较失败: {e}")
            comparison['similarity_score'] = 0.0
            comparison['details']['error'] = str(e)
        
        return comparison
    
    def generate_stereochemistry_report(self, stereo_info: StereochemistryInfo) -> Dict[str, Any]:
        """
        生成立体化学分析报告
        
        Args:
            stereo_info: 立体化学信息
            
        Returns:
            分析报告
        """
        report = {
            'total_chiral_centers': stereo_info.total_chiral_centers,
            'stereoisomer_type': stereo_info.stereoisomer_type,
            'confidence': stereo_info.confidence,
            'chiral_centers_detail': [],
            'summary': ''
        }
        
        try:
            # 详细手性中心信息
            for i, center in enumerate(stereo_info.chiral_centers):
                center_info = {
                    'index': i,
                    'atom_idx': center.atom_idx,
                    'atom_symbol': center.atom_symbol,
                    'configuration': center.configuration,
                    'coordinates': center.coordinates,
                    'neighbor_count': len(center.neighbors) if center.neighbors else 0
                }
                report['chiral_centers_detail'].append(center_info)
            
            # 生成摘要
            if stereo_info.total_chiral_centers == 0:
                report['summary'] = "分子不含手性中心，为非手性分子"
            elif stereo_info.total_chiral_centers == 1:
                config = stereo_info.chiral_centers[0].configuration
                if config:
                    report['summary'] = f"分子含1个手性中心，构型为{config}"
                else:
                    report['summary'] = "分子含1个手性中心，构型未确定"
            else:
                configured = len([c for c in stereo_info.chiral_centers if c.configuration])
                report['summary'] = f"分子含{stereo_info.total_chiral_centers}个手性中心，其中{configured}个构型已确定"
            
        except Exception as e:
            logging.error(f"立体化学报告生成失败: {e}")
            report['summary'] = f"报告生成失败: {e}"
        
        return report