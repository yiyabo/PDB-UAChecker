#!/usr/bin/env python3
"""
可扩展非天然氨基酸PDB搜索引擎 - 第三阶段高级功能
实现3D结构匹配、机器学习分类、异构体识别和高级分析工具

高级功能模块：
- Structure3DAnalyzer: 三维结构分析和匹配
- MLClassifier: 机器学习分类器
- StereochemistryAnalyzer: 立体化学和异构体分析
- MolecularVisualizer: 分子可视化
- AdvancedAnalytics: 高级分析工具
"""

import os
import json
import time
import numpy as np
import pandas as pd
from typing import Dict, List, Set, Optional, Tuple, Any, Union
from collections import defaultdict, Counter
from dataclasses import dataclass, asdict
import pickle
import warnings
warnings.filterwarnings('ignore')

# 导入前两阶段的核心组件
from scalable_search_engine import AminoAcidRecord, SearchResult
from performance_optimized_engine import (
    PerformanceOptimizedSearchEngine, 
    PerformanceConfig
)

# 科学计算库
try:
    from scipy.spatial.distance import cdist, euclidean
    from scipy.optimize import minimize
    from scipy.spatial.transform import Rotation
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import SVC
    from sklearn.neural_network import MLPClassifier
    from sklearn.model_selection import train_test_split, cross_val_score
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import classification_report, confusion_matrix
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    print("警告: scikit-learn未安装，机器学习功能将不可用")

@dataclass
class Structure3D:
    """三维结构数据类"""
    amino_acid_id: str
    atoms: List[Dict[str, Any]]  # 原子信息：{name, element, x, y, z}
    bonds: List[Tuple[int, int, str]]  # 键信息：(atom1_idx, atom2_idx, bond_type)
    center_of_mass: Tuple[float, float, float]
    molecular_volume: float
    surface_area: float
    
    def to_dict(self) -> Dict:
        return asdict(self)

@dataclass
class GeometricFeatures:
    """几何特征数据类"""
    bond_lengths: List[float]
    bond_angles: List[float]
    dihedral_angles: List[float]
    ring_systems: List[Dict[str, Any]]
    chiral_centers: List[Dict[str, Any]]
    planarity_deviations: List[float]

class Structure3DAnalyzer:
    """三维结构分析器"""
    
    def __init__(self):
        # 标准键长数据（Å）
        self.standard_bond_lengths = {
            ('C', 'C'): 1.54,  # C-C单键
            ('C', 'N'): 1.47,  # C-N单键
            ('C', 'O'): 1.43,  # C-O单键
            ('C', 'S'): 1.82,  # C-S单键
            ('N', 'H'): 1.01,  # N-H键
            ('O', 'H'): 0.96,  # O-H键
            ('C', 'H'): 1.09,  # C-H键
            # 双键
            ('C', '=C'): 1.34,  # C=C双键
            ('C', '=N'): 1.28,  # C=N双键
            ('C', '=O'): 1.23,  # C=O双键
            # 芳香键
            ('C', ':C'): 1.39,  # 芳香C-C键
        }
        
        # 标准键角数据（度）
        self.standard_bond_angles = {
            'sp3': 109.5,  # 四面体角
            'sp2': 120.0,  # 平面三角形
            'sp': 180.0,   # 线性
        }
        
        # 容差范围
        self.bond_length_tolerance = 0.1  # Å
        self.bond_angle_tolerance = 10.0   # 度
        self.planarity_tolerance = 0.1     # Å
    
    def parse_pdb_structure(self, pdb_content: str, residue_name: str) -> Optional[Structure3D]:
        """从PDB内容解析三维结构"""
        atoms = []
        bonds = []
        
        lines = pdb_content.strip().split('\n')
        
        # 解析原子坐标
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                if line[17:20].strip() == residue_name:
                    atom_info = {
                        'name': line[12:16].strip(),
                        'element': line[76:78].strip() or self._infer_element(line[12:16].strip()),
                        'x': float(line[30:38]),
                        'y': float(line[38:46]),
                        'z': float(line[46:54]),
                        'occupancy': float(line[54:60]) if line[54:60].strip() else 1.0,
                        'b_factor': float(line[60:66]) if line[60:66].strip() else 0.0
                    }
                    atoms.append(atom_info)
        
        if not atoms:
            return None
        
        # 解析键连接（CONECT记录）
        atom_name_to_idx = {atom['name']: i for i, atom in enumerate(atoms)}
        
        for line in lines:
            if line.startswith('CONECT'):
                try:
                    atom_idx = int(line[6:11]) - 1  # PDB索引从1开始
                    connected_atoms = []
                    
                    # 解析连接的原子
                    for i in range(11, len(line), 5):
                        if i + 5 <= len(line):
                            connected_idx = line[i:i+5].strip()
                            if connected_idx:
                                connected_atoms.append(int(connected_idx) - 1)
                    
                    # 添加键信息
                    for connected_idx in connected_atoms:
                        if 0 <= atom_idx < len(atoms) and 0 <= connected_idx < len(atoms):
                            bonds.append((atom_idx, connected_idx, 'single'))
                
                except (ValueError, IndexError):
                    continue
        
        # 如果没有CONECT记录，基于距离推断键
        if not bonds:
            bonds = self._infer_bonds_from_distance(atoms)
        
        # 计算分子属性
        center_of_mass = self._calculate_center_of_mass(atoms)
        molecular_volume = self._estimate_molecular_volume(atoms)
        surface_area = self._estimate_surface_area(atoms)
        
        return Structure3D(
            amino_acid_id=residue_name,
            atoms=atoms,
            bonds=bonds,
            center_of_mass=center_of_mass,
            molecular_volume=molecular_volume,
            surface_area=surface_area
        )
    
    def _infer_element(self, atom_name: str) -> str:
        """从原子名推断元素符号"""
        atom_name = atom_name.strip()
        if atom_name.startswith('C'):
            return 'C'
        elif atom_name.startswith('N'):
            return 'N'
        elif atom_name.startswith('O'):
            return 'O'
        elif atom_name.startswith('S'):
            return 'S'
        elif atom_name.startswith('P'):
            return 'P'
        elif atom_name.startswith('H'):
            return 'H'
        else:
            return atom_name[0] if atom_name else 'C'
    
    def _infer_bonds_from_distance(self, atoms: List[Dict]) -> List[Tuple[int, int, str]]:
        """基于原子间距离推断化学键"""
        bonds = []
        
        # 典型键长范围
        max_bond_distances = {
            ('C', 'C'): 1.8,
            ('C', 'N'): 1.7,
            ('C', 'O'): 1.6,
            ('C', 'S'): 2.1,
            ('N', 'N'): 1.6,
            ('N', 'O'): 1.5,
            ('O', 'O'): 1.5,
            ('C', 'H'): 1.3,
            ('N', 'H'): 1.2,
            ('O', 'H'): 1.1,
            ('S', 'H'): 1.5
        }
        
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                atom1 = atoms[i]
                atom2 = atoms[j]
                
                # 计算距离
                distance = np.sqrt(
                    (atom1['x'] - atom2['x'])**2 +
                    (atom1['y'] - atom2['y'])**2 +
                    (atom1['z'] - atom2['z'])**2
                )
                
                # 判断是否形成键
                element_pair = tuple(sorted([atom1['element'], atom2['element']]))
                max_distance = max_bond_distances.get(element_pair, 2.0)
                
                if distance <= max_distance:
                    bonds.append((i, j, 'single'))
        
        return bonds
    
    def _calculate_center_of_mass(self, atoms: List[Dict]) -> Tuple[float, float, float]:
        """计算质心"""
        # 原子质量（简化）
        atomic_masses = {
            'C': 12.01, 'N': 14.01, 'O': 16.00, 'S': 32.07,
            'P': 30.97, 'H': 1.008, 'F': 19.00, 'Cl': 35.45
        }
        
        total_mass = 0.0
        weighted_x = weighted_y = weighted_z = 0.0
        
        for atom in atoms:
            mass = atomic_masses.get(atom['element'], 12.0)
            total_mass += mass
            weighted_x += atom['x'] * mass
            weighted_y += atom['y'] * mass
            weighted_z += atom['z'] * mass
        
        if total_mass > 0:
            return (weighted_x / total_mass, weighted_y / total_mass, weighted_z / total_mass)
        else:
            return (0.0, 0.0, 0.0)
    
    def _estimate_molecular_volume(self, atoms: List[Dict]) -> float:
        """估算分子体积"""
        # 使用原子半径估算
        atomic_radii = {
            'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8,
            'P': 1.8, 'H': 1.2, 'F': 1.47, 'Cl': 1.75
        }
        
        total_volume = 0.0
        for atom in atoms:
            radius = atomic_radii.get(atom['element'], 1.7)
            volume = (4/3) * np.pi * radius**3
            total_volume += volume
        
        # 考虑重叠，乘以经验因子
        return total_volume * 0.7
    
    def _estimate_surface_area(self, atoms: List[Dict]) -> float:
        """估算分子表面积"""
        # 简化的表面积计算
        atomic_radii = {
            'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8,
            'P': 1.8, 'H': 1.2, 'F': 1.47, 'Cl': 1.75
        }
        
        total_surface = 0.0
        for atom in atoms:
            radius = atomic_radii.get(atom['element'], 1.7)
            surface = 4 * np.pi * radius**2
            total_surface += surface
        
        # 考虑重叠，乘以经验因子
        return total_surface * 0.6
    
    def extract_geometric_features(self, structure: Structure3D) -> GeometricFeatures:
        """提取几何特征"""
        bond_lengths = []
        bond_angles = []
        dihedral_angles = []
        
        # 计算键长
        for bond in structure.bonds:
            atom1_idx, atom2_idx, bond_type = bond
            atom1 = structure.atoms[atom1_idx]
            atom2 = structure.atoms[atom2_idx]
            
            distance = np.sqrt(
                (atom1['x'] - atom2['x'])**2 +
                (atom1['y'] - atom2['y'])**2 +
                (atom1['z'] - atom2['z'])**2
            )
            bond_lengths.append(distance)
        
        # 计算键角
        bond_angles = self._calculate_bond_angles(structure)
        
        # 计算二面角
        dihedral_angles = self._calculate_dihedral_angles(structure)
        
        # 检测环系统
        ring_systems = self._detect_ring_systems(structure)
        
        # 检测手性中心
        chiral_centers = self._detect_chiral_centers(structure)
        
        # 计算平面性偏差
        planarity_deviations = self._calculate_planarity_deviations(structure)
        
        return GeometricFeatures(
            bond_lengths=bond_lengths,
            bond_angles=bond_angles,
            dihedral_angles=dihedral_angles,
            ring_systems=ring_systems,
            chiral_centers=chiral_centers,
            planarity_deviations=planarity_deviations
        )
    
    def _calculate_bond_angles(self, structure: Structure3D) -> List[float]:
        """计算键角"""
        angles = []
        
        # 构建邻接表
        adjacency = defaultdict(list)
        for bond in structure.bonds:
            atom1_idx, atom2_idx, _ = bond
            adjacency[atom1_idx].append(atom2_idx)
            adjacency[atom2_idx].append(atom1_idx)
        
        # 计算每个原子的键角
        for center_idx, neighbors in adjacency.items():
            if len(neighbors) >= 2:
                center_atom = structure.atoms[center_idx]
                center_pos = np.array([center_atom['x'], center_atom['y'], center_atom['z']])
                
                # 计算所有邻居对的角度
                for i in range(len(neighbors)):
                    for j in range(i + 1, len(neighbors)):
                        atom1 = structure.atoms[neighbors[i]]
                        atom2 = structure.atoms[neighbors[j]]
                        
                        pos1 = np.array([atom1['x'], atom1['y'], atom1['z']])
                        pos2 = np.array([atom2['x'], atom2['y'], atom2['z']])
                        
                        # 计算向量
                        vec1 = pos1 - center_pos
                        vec2 = pos2 - center_pos
                        
                        # 计算角度
                        cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
                        cos_angle = np.clip(cos_angle, -1.0, 1.0)
                        angle = np.arccos(cos_angle) * 180.0 / np.pi
                        angles.append(angle)
        
        return angles
    
    def _calculate_dihedral_angles(self, structure: Structure3D) -> List[float]:
        """计算二面角"""
        dihedrals = []
        
        # 构建邻接表
        adjacency = defaultdict(list)
        for bond in structure.bonds:
            atom1_idx, atom2_idx, _ = bond
            adjacency[atom1_idx].append(atom2_idx)
            adjacency[atom2_idx].append(atom1_idx)
        
        # 查找四原子序列计算二面角
        for bond in structure.bonds:
            atom2_idx, atom3_idx, _ = bond
            
            # 找到atom2和atom3的其他邻居
            atom2_neighbors = [n for n in adjacency[atom2_idx] if n != atom3_idx]
            atom3_neighbors = [n for n in adjacency[atom3_idx] if n != atom2_idx]
            
            for atom1_idx in atom2_neighbors:
                for atom4_idx in atom3_neighbors:
                    # 计算二面角 atom1-atom2-atom3-atom4
                    dihedral = self._calculate_dihedral_angle_four_atoms(
                        structure.atoms[atom1_idx],
                        structure.atoms[atom2_idx],
                        structure.atoms[atom3_idx],
                        structure.atoms[atom4_idx]
                    )
                    dihedrals.append(dihedral)
        
        return dihedrals
    
    def _calculate_dihedral_angle_four_atoms(self, atom1: Dict, atom2: Dict, 
                                           atom3: Dict, atom4: Dict) -> float:
        """计算四个原子的二面角"""
        # 获取坐标
        p1 = np.array([atom1['x'], atom1['y'], atom1['z']])
        p2 = np.array([atom2['x'], atom2['y'], atom2['z']])
        p3 = np.array([atom3['x'], atom3['y'], atom3['z']])
        p4 = np.array([atom4['x'], atom4['y'], atom4['z']])
        
        # 计算向量
        v1 = p2 - p1
        v2 = p3 - p2
        v3 = p4 - p3
        
        # 计算法向量
        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)
        
        # 归一化
        n1 = n1 / np.linalg.norm(n1)
        n2 = n2 / np.linalg.norm(n2)
        
        # 计算二面角
        cos_angle = np.dot(n1, n2)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle = np.arccos(cos_angle) * 180.0 / np.pi
        
        # 确定符号
        if np.dot(np.cross(n1, n2), v2) < 0:
            angle = -angle
        
        return angle

    def _detect_ring_systems(self, structure: Structure3D) -> List[Dict[str, Any]]:
        """检测环系统"""
        rings = []

        # 构建邻接表
        adjacency = defaultdict(list)
        for bond in structure.bonds:
            atom1_idx, atom2_idx, _ = bond
            adjacency[atom1_idx].append(atom2_idx)
            adjacency[atom2_idx].append(atom1_idx)

        # 使用DFS查找环
        visited = set()

        def dfs_find_rings(start, current, path, depth):
            if depth > 8:  # 限制搜索深度
                return

            if current in path and len(path) >= 3:
                # 找到环
                ring_start = path.index(current)
                ring = path[ring_start:]
                if len(ring) >= 3 and len(ring) <= 8:  # 3-8元环
                    ring_info = {
                        'atoms': ring,
                        'size': len(ring),
                        'type': self._classify_ring_type(ring, structure)
                    }
                    rings.append(ring_info)
                return

            if current in visited:
                return

            path.append(current)

            for neighbor in adjacency[current]:
                if neighbor != path[-2] if len(path) > 1 else True:  # 避免回到上一个原子
                    dfs_find_rings(start, neighbor, path.copy(), depth + 1)

        # 从每个原子开始搜索
        for atom_idx in range(len(structure.atoms)):
            if atom_idx not in visited:
                dfs_find_rings(atom_idx, atom_idx, [], 0)
                visited.add(atom_idx)

        return rings

    def _classify_ring_type(self, ring_atoms: List[int], structure: Structure3D) -> str:
        """分类环的类型"""
        # 检查是否为芳香环
        carbon_count = sum(1 for idx in ring_atoms
                          if structure.atoms[idx]['element'] == 'C')

        if len(ring_atoms) == 6 and carbon_count == 6:
            return 'benzene'
        elif len(ring_atoms) == 5 and carbon_count == 4:
            return 'pyrrole'
        elif len(ring_atoms) == 6 and carbon_count == 5:
            return 'pyridine'
        elif len(ring_atoms) == 5:
            return 'five_membered'
        elif len(ring_atoms) == 6:
            return 'six_membered'
        else:
            return f'{len(ring_atoms)}_membered'

    def _detect_chiral_centers(self, structure: Structure3D) -> List[Dict[str, Any]]:
        """检测手性中心"""
        chiral_centers = []

        # 构建邻接表
        adjacency = defaultdict(list)
        for bond in structure.bonds:
            atom1_idx, atom2_idx, _ = bond
            adjacency[atom1_idx].append(atom2_idx)
            adjacency[atom2_idx].append(atom1_idx)

        # 检查每个碳原子
        for atom_idx, atom in enumerate(structure.atoms):
            if atom['element'] == 'C':
                neighbors = adjacency[atom_idx]

                # 手性中心需要4个不同的取代基
                if len(neighbors) == 4:
                    # 简化的手性检测：检查邻居的元素类型
                    neighbor_elements = [structure.atoms[n]['element'] for n in neighbors]

                    # 如果有4个不同的邻居类型，可能是手性中心
                    if len(set(neighbor_elements)) >= 3:
                        chiral_info = {
                            'atom_index': atom_idx,
                            'atom_name': atom['name'],
                            'neighbors': neighbors,
                            'neighbor_elements': neighbor_elements,
                            'chirality': self._determine_chirality(atom_idx, neighbors, structure)
                        }
                        chiral_centers.append(chiral_info)

        return chiral_centers

    def _determine_chirality(self, center_idx: int, neighbors: List[int],
                           structure: Structure3D) -> str:
        """确定手性（R/S）"""
        # 简化的手性判断
        # 实际应用中需要更复杂的CIP规则

        center_atom = structure.atoms[center_idx]
        center_pos = np.array([center_atom['x'], center_atom['y'], center_atom['z']])

        # 计算邻居原子的优先级（基于原子序数）
        atomic_numbers = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'S': 16, 'P': 15}

        neighbor_priorities = []
        for neighbor_idx in neighbors:
            neighbor_atom = structure.atoms[neighbor_idx]
            priority = atomic_numbers.get(neighbor_atom['element'], 6)
            neighbor_priorities.append((priority, neighbor_idx))

        # 按优先级排序
        neighbor_priorities.sort(reverse=True)

        # 简化判断：基于几何排列
        if len(neighbor_priorities) == 4:
            # 计算四面体的手性
            positions = []
            for _, neighbor_idx in neighbor_priorities:
                neighbor_atom = structure.atoms[neighbor_idx]
                pos = np.array([neighbor_atom['x'], neighbor_atom['y'], neighbor_atom['z']])
                positions.append(pos - center_pos)

            # 计算混合积判断手性
            if len(positions) >= 3:
                mixed_product = np.dot(positions[0], np.cross(positions[1], positions[2]))
                return 'R' if mixed_product > 0 else 'S'

        return 'unknown'

    def _calculate_planarity_deviations(self, structure: Structure3D) -> List[float]:
        """计算平面性偏差"""
        deviations = []

        # 检测可能的平面结构（如芳香环）
        rings = self._detect_ring_systems(structure)

        for ring in rings:
            if ring['size'] >= 4:  # 至少4个原子才能定义平面
                ring_atoms = ring['atoms']

                # 获取原子坐标
                coords = []
                for atom_idx in ring_atoms:
                    atom = structure.atoms[atom_idx]
                    coords.append([atom['x'], atom['y'], atom['z']])

                coords = np.array(coords)

                # 计算最佳拟合平面
                if len(coords) >= 3:
                    # 使用SVD计算平面
                    centroid = np.mean(coords, axis=0)
                    centered_coords = coords - centroid

                    _, _, vh = np.linalg.svd(centered_coords)
                    normal = vh[-1]  # 最小奇异值对应的向量

                    # 计算每个原子到平面的距离
                    distances = []
                    for coord in coords:
                        distance = abs(np.dot(coord - centroid, normal))
                        distances.append(distance)

                    # 平面性偏差为最大距离
                    max_deviation = max(distances)
                    deviations.append(max_deviation)

        return deviations

    def calculate_rmsd(self, structure1: Structure3D, structure2: Structure3D) -> float:
        """计算两个结构的RMSD"""
        if len(structure1.atoms) != len(structure2.atoms):
            return float('inf')

        # 提取坐标
        coords1 = np.array([[atom['x'], atom['y'], atom['z']] for atom in structure1.atoms])
        coords2 = np.array([[atom['x'], atom['y'], atom['z']] for atom in structure2.atoms])

        # 计算RMSD
        diff = coords1 - coords2
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

        return rmsd

    def superpose_structures(self, reference: Structure3D, mobile: Structure3D) -> Tuple[Structure3D, float]:
        """分子叠合算法"""
        if len(reference.atoms) != len(mobile.atoms):
            raise ValueError("结构的原子数量不匹配")

        # 提取坐标
        ref_coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in reference.atoms])
        mob_coords = np.array([[atom['x'], atom['y'], atom['z']] for atom in mobile.atoms])

        # 中心化坐标
        ref_centroid = np.mean(ref_coords, axis=0)
        mob_centroid = np.mean(mob_coords, axis=0)

        ref_centered = ref_coords - ref_centroid
        mob_centered = mob_coords - mob_centroid

        # 使用Kabsch算法计算最优旋转矩阵
        H = np.dot(mob_centered.T, ref_centered)
        U, S, Vt = np.linalg.svd(H)
        R = np.dot(Vt.T, U.T)

        # 确保旋转矩阵的行列式为正
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = np.dot(Vt.T, U.T)

        # 应用变换
        transformed_coords = np.dot(mob_centered, R.T) + ref_centroid

        # 创建变换后的结构
        transformed_atoms = []
        for i, atom in enumerate(mobile.atoms):
            transformed_atom = atom.copy()
            transformed_atom['x'] = transformed_coords[i, 0]
            transformed_atom['y'] = transformed_coords[i, 1]
            transformed_atom['z'] = transformed_coords[i, 2]
            transformed_atoms.append(transformed_atom)

        transformed_structure = Structure3D(
            amino_acid_id=mobile.amino_acid_id + '_aligned',
            atoms=transformed_atoms,
            bonds=mobile.bonds,
            center_of_mass=tuple(np.mean(transformed_coords, axis=0)),
            molecular_volume=mobile.molecular_volume,
            surface_area=mobile.surface_area
        )

        # 计算叠合后的RMSD
        rmsd = self.calculate_rmsd(reference, transformed_structure)

        return transformed_structure, rmsd

    def validate_geometry(self, structure: Structure3D) -> Dict[str, Any]:
        """验证分子几何"""
        validation_results = {
            'valid_bond_lengths': [],
            'invalid_bond_lengths': [],
            'valid_bond_angles': [],
            'invalid_bond_angles': [],
            'planarity_violations': [],
            'overall_valid': True
        }

        # 验证键长
        for bond in structure.bonds:
            atom1_idx, atom2_idx, bond_type = bond
            atom1 = structure.atoms[atom1_idx]
            atom2 = structure.atoms[atom2_idx]

            # 计算实际键长
            distance = np.sqrt(
                (atom1['x'] - atom2['x'])**2 +
                (atom1['y'] - atom2['y'])**2 +
                (atom1['z'] - atom2['z'])**2
            )

            # 获取标准键长
            element_pair = tuple(sorted([atom1['element'], atom2['element']]))
            standard_length = self.standard_bond_lengths.get(element_pair, 1.5)

            # 检查是否在容差范围内
            if abs(distance - standard_length) <= self.bond_length_tolerance:
                validation_results['valid_bond_lengths'].append({
                    'bond': bond,
                    'actual_length': distance,
                    'standard_length': standard_length
                })
            else:
                validation_results['invalid_bond_lengths'].append({
                    'bond': bond,
                    'actual_length': distance,
                    'standard_length': standard_length,
                    'deviation': abs(distance - standard_length)
                })
                validation_results['overall_valid'] = False

        # 验证键角
        geometric_features = self.extract_geometric_features(structure)
        for angle in geometric_features.bond_angles:
            # 简化的键角验证
            if 90 <= angle <= 130:  # 合理的键角范围
                validation_results['valid_bond_angles'].append(angle)
            else:
                validation_results['invalid_bond_angles'].append(angle)
                if angle < 60 or angle > 150:  # 严重偏差
                    validation_results['overall_valid'] = False

        # 验证平面性
        for deviation in geometric_features.planarity_deviations:
            if deviation > self.planarity_tolerance:
                validation_results['planarity_violations'].append(deviation)
                if deviation > 0.3:  # 严重偏差
                    validation_results['overall_valid'] = False

        return validation_results

class MLClassifier:
    """机器学习分类器"""

    def __init__(self):
        self.models = {}
        self.scalers = {}
        self.feature_names = []
        self.is_trained = False

        # 支持的模型类型
        if SKLEARN_AVAILABLE:
            self.model_types = {
                'random_forest': RandomForestClassifier(n_estimators=100, random_state=42),
                'svm': SVC(kernel='rbf', probability=True, random_state=42),
                'neural_network': MLPClassifier(hidden_layer_sizes=(100, 50), max_iter=1000, random_state=42)
            }
        else:
            self.model_types = {}
            print("警告: scikit-learn未安装，机器学习功能不可用")

    def extract_molecular_features(self, record: AminoAcidRecord,
                                 structure: Optional[Structure3D] = None) -> np.ndarray:
        """提取分子描述符特征"""
        features = []

        # 基本分子特征
        features.extend([
            record.molecular_weight,
            len(record.atom_composition),
            sum(record.atom_composition.values()),  # 总原子数
            record.atom_composition.get('C', 0),
            record.atom_composition.get('N', 0),
            record.atom_composition.get('O', 0),
            record.atom_composition.get('S', 0),
            record.atom_composition.get('P', 0),
            record.atom_composition.get('H', 0)
        ])

        # 原子比例特征
        total_atoms = sum(record.atom_composition.values())
        if total_atoms > 0:
            features.extend([
                record.atom_composition.get('C', 0) / total_atoms,
                record.atom_composition.get('N', 0) / total_atoms,
                record.atom_composition.get('O', 0) / total_atoms,
                (record.atom_composition.get('N', 0) + record.atom_composition.get('O', 0)) / total_atoms,  # 杂原子比例
            ])
        else:
            features.extend([0, 0, 0, 0])

        # 关键特征编码
        feature_encoding = {
            'aromatic_ring': 0,
            'carboxyl_group': 0,
            'amino_group': 0,
            'hydroxyl_group': 0,
            'methoxy_group': 0,
            'double_bond': 0,
            'triple_bond': 0,
            'sulfur_group': 0
        }

        for feature in record.key_features:
            if feature in feature_encoding:
                feature_encoding[feature] = 1

        features.extend(list(feature_encoding.values()))

        # 3D结构特征（如果可用）
        if structure:
            features.extend([
                structure.molecular_volume,
                structure.surface_area,
                len(structure.bonds),
                structure.center_of_mass[0],  # X坐标
                structure.center_of_mass[1],  # Y坐标
                structure.center_of_mass[2],  # Z坐标
            ])

            # 几何特征
            geometric_features = Structure3DAnalyzer().extract_geometric_features(structure)

            # 统计特征
            if geometric_features.bond_lengths:
                features.extend([
                    np.mean(geometric_features.bond_lengths),
                    np.std(geometric_features.bond_lengths),
                    min(geometric_features.bond_lengths),
                    max(geometric_features.bond_lengths)
                ])
            else:
                features.extend([0, 0, 0, 0])

            if geometric_features.bond_angles:
                features.extend([
                    np.mean(geometric_features.bond_angles),
                    np.std(geometric_features.bond_angles)
                ])
            else:
                features.extend([0, 0])

            # 环系统特征
            features.extend([
                len(geometric_features.ring_systems),
                len(geometric_features.chiral_centers)
            ])
        else:
            # 如果没有3D结构，用零填充
            features.extend([0] * 14)

        return np.array(features)

    def prepare_training_data(self, amino_acid_records: List[AminoAcidRecord],
                            structures: Optional[Dict[str, Structure3D]] = None,
                            labels: Optional[List[str]] = None) -> Tuple[np.ndarray, np.ndarray]:
        """准备训练数据"""
        if not SKLEARN_AVAILABLE:
            raise ImportError("scikit-learn未安装，无法进行机器学习")

        features_list = []
        labels_list = []

        for i, record in enumerate(amino_acid_records):
            structure = structures.get(record.id) if structures else None
            features = self.extract_molecular_features(record, structure)
            features_list.append(features)

            # 如果没有提供标签，使用氨基酸ID作为标签
            if labels:
                labels_list.append(labels[i])
            else:
                labels_list.append(record.id)

        X = np.array(features_list)
        y = np.array(labels_list)

        # 保存特征名称
        if not self.feature_names:
            self.feature_names = [
                'molecular_weight', 'num_atom_types', 'total_atoms',
                'C_count', 'N_count', 'O_count', 'S_count', 'P_count', 'H_count',
                'C_ratio', 'N_ratio', 'O_ratio', 'heteroatom_ratio',
                'aromatic_ring', 'carboxyl_group', 'amino_group', 'hydroxyl_group',
                'methoxy_group', 'double_bond', 'triple_bond', 'sulfur_group',
                'molecular_volume', 'surface_area', 'num_bonds',
                'center_x', 'center_y', 'center_z',
                'mean_bond_length', 'std_bond_length', 'min_bond_length', 'max_bond_length',
                'mean_bond_angle', 'std_bond_angle',
                'num_rings', 'num_chiral_centers'
            ]

        return X, y

    def train_models(self, X: np.ndarray, y: np.ndarray,
                    model_types: List[str] = None) -> Dict[str, Dict[str, Any]]:
        """训练多个模型"""
        if not SKLEARN_AVAILABLE:
            raise ImportError("scikit-learn未安装，无法进行机器学习")

        if model_types is None:
            model_types = list(self.model_types.keys())

        results = {}

        # 数据预处理
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # 分割训练和测试数据
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y, test_size=0.2, random_state=42, stratify=y
        )

        for model_name in model_types:
            if model_name not in self.model_types:
                continue

            print(f"训练 {model_name} 模型...")

            # 训练模型
            model = self.model_types[model_name]
            model.fit(X_train, y_train)

            # 预测
            y_pred = model.predict(X_test)
            y_pred_proba = model.predict_proba(X_test) if hasattr(model, 'predict_proba') else None

            # 交叉验证
            cv_scores = cross_val_score(model, X_scaled, y, cv=5)

            # 保存模型和缩放器
            self.models[model_name] = model
            self.scalers[model_name] = scaler

            # 计算性能指标
            results[model_name] = {
                'accuracy': np.mean(y_pred == y_test),
                'cv_mean': np.mean(cv_scores),
                'cv_std': np.std(cv_scores),
                'classification_report': classification_report(y_test, y_pred),
                'confusion_matrix': confusion_matrix(y_test, y_pred).tolist(),
                'feature_importance': self._get_feature_importance(model, model_name)
            }

            print(f"{model_name} 准确率: {results[model_name]['accuracy']:.3f}")
            print(f"{model_name} 交叉验证: {results[model_name]['cv_mean']:.3f} ± {results[model_name]['cv_std']:.3f}")

        self.is_trained = True
        return results

    def _get_feature_importance(self, model, model_name: str) -> List[Tuple[str, float]]:
        """获取特征重要性"""
        if model_name == 'random_forest' and hasattr(model, 'feature_importances_'):
            importances = model.feature_importances_
        elif model_name == 'svm' and hasattr(model, 'coef_'):
            # SVM的系数绝对值作为重要性
            importances = np.abs(model.coef_[0]) if len(model.coef_.shape) > 1 else np.abs(model.coef_)
        else:
            # 其他模型返回空列表
            return []

        # 组合特征名称和重要性
        feature_importance = list(zip(self.feature_names, importances))
        feature_importance.sort(key=lambda x: x[1], reverse=True)

        return feature_importance[:10]  # 返回前10个重要特征

    def predict(self, record: AminoAcidRecord, structure: Optional[Structure3D] = None,
               model_name: str = 'random_forest') -> Dict[str, Any]:
        """预测氨基酸类别"""
        if not self.is_trained or model_name not in self.models:
            raise ValueError(f"模型 {model_name} 未训练或不存在")

        # 提取特征
        features = self.extract_molecular_features(record, structure)
        features = features.reshape(1, -1)

        # 标准化
        scaler = self.scalers[model_name]
        features_scaled = scaler.transform(features)

        # 预测
        model = self.models[model_name]
        prediction = model.predict(features_scaled)[0]

        # 预测概率
        if hasattr(model, 'predict_proba'):
            probabilities = model.predict_proba(features_scaled)[0]
            class_probabilities = dict(zip(model.classes_, probabilities))
        else:
            class_probabilities = {}

        return {
            'predicted_class': prediction,
            'confidence': max(class_probabilities.values()) if class_probabilities else 0.0,
            'class_probabilities': class_probabilities,
            'model_used': model_name
        }

    def batch_predict(self, records: List[AminoAcidRecord],
                     structures: Optional[Dict[str, Structure3D]] = None,
                     model_name: str = 'random_forest') -> List[Dict[str, Any]]:
        """批量预测"""
        results = []

        for record in records:
            structure = structures.get(record.id) if structures else None
            prediction = self.predict(record, structure, model_name)
            prediction['amino_acid_id'] = record.id
            results.append(prediction)

        return results

    def incremental_learning(self, new_records: List[AminoAcidRecord],
                           new_structures: Optional[Dict[str, Structure3D]] = None,
                           new_labels: Optional[List[str]] = None,
                           model_name: str = 'random_forest'):
        """增量学习"""
        if not self.is_trained:
            raise ValueError("模型尚未进行初始训练")

        # 准备新数据
        X_new, y_new = self.prepare_training_data(new_records, new_structures, new_labels)

        # 标准化新数据
        scaler = self.scalers[model_name]
        X_new_scaled = scaler.transform(X_new)

        # 对于支持增量学习的模型
        model = self.models[model_name]

        if hasattr(model, 'partial_fit'):
            # 支持增量学习
            model.partial_fit(X_new_scaled, y_new)
            print(f"模型 {model_name} 增量学习完成，新增 {len(new_records)} 个样本")
        else:
            # 不支持增量学习，需要重新训练
            print(f"模型 {model_name} 不支持增量学习，建议重新训练")

    def save_models(self, filepath: str):
        """保存训练好的模型"""
        model_data = {
            'models': self.models,
            'scalers': self.scalers,
            'feature_names': self.feature_names,
            'is_trained': self.is_trained
        }

        with open(filepath, 'wb') as f:
            pickle.dump(model_data, f)

        print(f"模型已保存到: {filepath}")

    def load_models(self, filepath: str):
        """加载训练好的模型"""
        with open(filepath, 'rb') as f:
            model_data = pickle.load(f)

        self.models = model_data['models']
        self.scalers = model_data['scalers']
        self.feature_names = model_data['feature_names']
        self.is_trained = model_data['is_trained']

        print(f"模型已从 {filepath} 加载")

    def get_model_statistics(self) -> Dict[str, Any]:
        """获取模型统计信息"""
        return {
            'trained_models': list(self.models.keys()),
            'num_features': len(self.feature_names),
            'feature_names': self.feature_names,
            'is_trained': self.is_trained,
            'sklearn_available': SKLEARN_AVAILABLE
        }

class StereochemistryAnalyzer:
    """立体化学和异构体分析器"""

    def __init__(self):
        # CIP优先级规则的原子序数
        self.atomic_numbers = {
            'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
            'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
            'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Br': 35, 'I': 53
        }

    def detect_stereoisomers(self, structures: List[Structure3D]) -> Dict[str, Any]:
        """检测立体异构体"""
        stereoisomer_groups = []
        processed = set()

        for i, structure1 in enumerate(structures):
            if structure1.amino_acid_id in processed:
                continue

            group = [structure1]
            processed.add(structure1.amino_acid_id)

            for j, structure2 in enumerate(structures[i+1:], i+1):
                if structure2.amino_acid_id in processed:
                    continue

                # 检查是否为立体异构体
                if self._are_stereoisomers(structure1, structure2):
                    group.append(structure2)
                    processed.add(structure2.amino_acid_id)

            if len(group) > 1:
                stereoisomer_groups.append({
                    'group_id': f"stereoisomer_group_{len(stereoisomer_groups)}",
                    'structures': group,
                    'stereoisomer_type': self._classify_stereoisomerism(group)
                })

        return {
            'stereoisomer_groups': stereoisomer_groups,
            'total_groups': len(stereoisomer_groups),
            'total_stereoisomers': sum(len(group['structures']) for group in stereoisomer_groups)
        }

    def _are_stereoisomers(self, structure1: Structure3D, structure2: Structure3D) -> bool:
        """判断两个结构是否为立体异构体"""
        # 检查分子式是否相同
        formula1 = self._get_molecular_formula(structure1)
        formula2 = self._get_molecular_formula(structure2)

        if formula1 != formula2:
            return False

        # 检查连接性是否相同
        connectivity1 = self._get_connectivity_signature(structure1)
        connectivity2 = self._get_connectivity_signature(structure2)

        if connectivity1 != connectivity2:
            return False

        # 检查3D结构是否不同
        rmsd = Structure3DAnalyzer().calculate_rmsd(structure1, structure2)

        # 如果RMSD较大，可能是立体异构体
        return rmsd > 0.5  # 阈值可调

    def _get_molecular_formula(self, structure: Structure3D) -> str:
        """获取分子式"""
        element_count = defaultdict(int)
        for atom in structure.atoms:
            element_count[atom['element']] += 1

        # 按标准顺序排列
        formula_parts = []
        for element in ['C', 'H', 'N', 'O', 'S', 'P']:
            if element in element_count:
                count = element_count[element]
                formula_parts.append(f"{element}{count}" if count > 1 else element)

        # 其他元素按字母顺序
        for element in sorted(element_count.keys()):
            if element not in ['C', 'H', 'N', 'O', 'S', 'P']:
                count = element_count[element]
                formula_parts.append(f"{element}{count}" if count > 1 else element)

        return ''.join(formula_parts)

    def _get_connectivity_signature(self, structure: Structure3D) -> str:
        """获取连接性签名"""
        # 构建邻接表
        adjacency = defaultdict(list)
        for bond in structure.bonds:
            atom1_idx, atom2_idx, _ = bond
            adjacency[atom1_idx].append(atom2_idx)
            adjacency[atom2_idx].append(atom1_idx)

        # 生成每个原子的连接性描述
        atom_signatures = []
        for i, atom in enumerate(structure.atoms):
            neighbors = adjacency[i]
            neighbor_elements = sorted([structure.atoms[j]['element'] for j in neighbors])
            signature = f"{atom['element']}({','.join(neighbor_elements)})"
            atom_signatures.append(signature)

        # 排序后连接
        return '|'.join(sorted(atom_signatures))

    def _classify_stereoisomerism(self, structures: List[Structure3D]) -> str:
        """分类立体异构体类型"""
        if len(structures) < 2:
            return 'none'

        # 检查手性中心
        analyzer = Structure3DAnalyzer()
        chiral_centers_list = []

        for structure in structures:
            geometric_features = analyzer.extract_geometric_features(structure)
            chiral_centers_list.append(len(geometric_features.chiral_centers))

        # 如果有手性中心，可能是对映异构体或非对映异构体
        if any(count > 0 for count in chiral_centers_list):
            if len(structures) == 2:
                return 'enantiomers'  # 对映异构体
            else:
                return 'diastereomers'  # 非对映异构体

        # 检查双键
        double_bonds_list = []
        for structure in structures:
            double_bond_count = sum(1 for bond in structure.bonds if bond[2] == 'double')
            double_bonds_list.append(double_bond_count)

        if any(count > 0 for count in double_bonds_list):
            return 'geometric_isomers'  # 几何异构体（E/Z）

        return 'conformational_isomers'  # 构象异构体

    def analyze_chirality(self, structure: Structure3D) -> Dict[str, Any]:
        """分析手性"""
        analyzer = Structure3DAnalyzer()
        geometric_features = analyzer.extract_geometric_features(structure)

        chirality_analysis = {
            'chiral_centers': geometric_features.chiral_centers,
            'num_chiral_centers': len(geometric_features.chiral_centers),
            'is_chiral': len(geometric_features.chiral_centers) > 0,
            'stereochemical_descriptors': []
        }

        # 为每个手性中心生成立体化学描述符
        for chiral_center in geometric_features.chiral_centers:
            descriptor = {
                'atom_index': chiral_center['atom_index'],
                'atom_name': chiral_center['atom_name'],
                'configuration': chiral_center['chirality'],
                'priority_order': self._get_cip_priority_order(
                    chiral_center, structure
                )
            }
            chirality_analysis['stereochemical_descriptors'].append(descriptor)

        return chirality_analysis

    def _get_cip_priority_order(self, chiral_center: Dict, structure: Structure3D) -> List[int]:
        """获取CIP优先级顺序"""
        neighbors = chiral_center['neighbors']

        # 计算每个邻居的优先级
        priorities = []
        for neighbor_idx in neighbors:
            neighbor_atom = structure.atoms[neighbor_idx]
            priority = self.atomic_numbers.get(neighbor_atom['element'], 0)
            priorities.append((priority, neighbor_idx))

        # 按优先级排序（高优先级在前）
        priorities.sort(reverse=True)

        return [idx for _, idx in priorities]

    def detect_ez_isomers(self, structure: Structure3D) -> List[Dict[str, Any]]:
        """检测E/Z异构体"""
        ez_isomers = []

        # 查找双键
        double_bonds = [bond for bond in structure.bonds if bond[2] == 'double']

        for bond in double_bonds:
            atom1_idx, atom2_idx, _ = bond

            # 获取双键两端原子的取代基
            substituents1 = self._get_substituents(atom1_idx, atom2_idx, structure)
            substituents2 = self._get_substituents(atom2_idx, atom1_idx, structure)

            if len(substituents1) >= 2 and len(substituents2) >= 2:
                # 可能存在E/Z异构
                ez_info = {
                    'double_bond': bond,
                    'atom1_substituents': substituents1,
                    'atom2_substituents': substituents2,
                    'configuration': self._determine_ez_configuration(
                        atom1_idx, atom2_idx, substituents1, substituents2, structure
                    )
                }
                ez_isomers.append(ez_info)

        return ez_isomers

    def _get_substituents(self, center_idx: int, exclude_idx: int,
                         structure: Structure3D) -> List[int]:
        """获取原子的取代基（排除指定原子）"""
        substituents = []

        for bond in structure.bonds:
            atom1_idx, atom2_idx, _ = bond

            if atom1_idx == center_idx and atom2_idx != exclude_idx:
                substituents.append(atom2_idx)
            elif atom2_idx == center_idx and atom1_idx != exclude_idx:
                substituents.append(atom1_idx)

        return substituents

    def _determine_ez_configuration(self, atom1_idx: int, atom2_idx: int,
                                  substituents1: List[int], substituents2: List[int],
                                  structure: Structure3D) -> str:
        """确定E/Z构型"""
        # 简化的E/Z判断
        # 实际应用中需要更复杂的CIP规则

        if len(substituents1) < 2 or len(substituents2) < 2:
            return 'unknown'

        # 获取最高优先级的取代基
        priority1 = self._get_highest_priority_substituent(substituents1, structure)
        priority2 = self._get_highest_priority_substituent(substituents2, structure)

        # 计算几何关系
        atom1 = structure.atoms[atom1_idx]
        atom2 = structure.atoms[atom2_idx]
        sub1 = structure.atoms[priority1]
        sub2 = structure.atoms[priority2]

        # 计算向量
        bond_vector = np.array([atom2['x'] - atom1['x'],
                               atom2['y'] - atom1['y'],
                               atom2['z'] - atom1['z']])

        sub1_vector = np.array([sub1['x'] - atom1['x'],
                               sub1['y'] - atom1['y'],
                               sub1['z'] - atom1['z']])

        sub2_vector = np.array([sub2['x'] - atom2['x'],
                               sub2['y'] - atom2['y'],
                               sub2['z'] - atom2['z']])

        # 简化判断：基于二面角
        # 实际应用中需要更精确的几何计算
        cross_product = np.cross(sub1_vector, sub2_vector)
        dot_product = np.dot(cross_product, bond_vector)

        return 'E' if dot_product > 0 else 'Z'

    def _get_highest_priority_substituent(self, substituents: List[int],
                                        structure: Structure3D) -> int:
        """获取最高优先级的取代基"""
        max_priority = -1
        highest_priority_idx = substituents[0]

        for sub_idx in substituents:
            atom = structure.atoms[sub_idx]
            priority = self.atomic_numbers.get(atom['element'], 0)

            if priority > max_priority:
                max_priority = priority
                highest_priority_idx = sub_idx

        return highest_priority_idx
