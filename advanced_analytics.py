#!/usr/bin/env python3
"""
高级分析工具模块
包含分子可视化、药物相似性分析、ADMET预测和化学空间分析
"""

import os
import json
import numpy as np
import pandas as pd
from typing import Dict, List, Set, Optional, Tuple, Any, Union
from collections import defaultdict, Counter
from dataclasses import dataclass, asdict
import matplotlib.pyplot as plt
import seaborn as sns

# 导入核心组件
from scalable_search_engine import AminoAcidRecord
from advanced_features_engine import Structure3D, Structure3DAnalyzer

# 科学计算库
try:
    from sklearn.decomposition import PCA
    from sklearn.cluster import KMeans, DBSCAN
    from sklearn.manifold import TSNE
    from sklearn.metrics import silhouette_score
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

@dataclass
class DrugLikenessScore:
    """药物相似性评分"""
    lipinski_violations: int
    molecular_weight: float
    logp: float
    hbd: int  # 氢键供体
    hba: int  # 氢键受体
    tpsa: float  # 拓扑极性表面积
    rotatable_bonds: int
    overall_score: float
    drug_like: bool

@dataclass
class ADMETProperties:
    """ADMET性质预测"""
    absorption: Dict[str, float]
    distribution: Dict[str, float]
    metabolism: Dict[str, float]
    excretion: Dict[str, float]
    toxicity: Dict[str, float]
    overall_assessment: str

class MolecularVisualizer:
    """分子可视化器"""
    
    def __init__(self):
        self.color_scheme = {
            'C': '#909090',
            'N': '#3050F8',
            'O': '#FF0D0D',
            'S': '#FFFF30',
            'P': '#FF8000',
            'H': '#FFFFFF',
            'F': '#90E050',
            'Cl': '#1FF01F',
            'Br': '#A62929',
            'I': '#940094'
        }
    
    def create_3d_structure_plot(self, structure: Structure3D, 
                               title: str = None) -> Optional[Any]:
        """创建3D结构图"""
        if not PLOTLY_AVAILABLE:
            print("警告: Plotly未安装，无法创建3D可视化")
            return None
        
        # 提取原子坐标和元素
        x_coords = [atom['x'] for atom in structure.atoms]
        y_coords = [atom['y'] for atom in structure.atoms]
        z_coords = [atom['z'] for atom in structure.atoms]
        elements = [atom['element'] for atom in structure.atoms]
        atom_names = [atom['name'] for atom in structure.atoms]
        
        # 创建原子散点图
        colors = [self.color_scheme.get(element, '#808080') for element in elements]
        
        fig = go.Figure()
        
        # 添加原子
        fig.add_trace(go.Scatter3d(
            x=x_coords,
            y=y_coords,
            z=z_coords,
            mode='markers',
            marker=dict(
                size=8,
                color=colors,
                opacity=0.8
            ),
            text=[f"{name} ({element})" for name, element in zip(atom_names, elements)],
            hovertemplate='%{text}<br>坐标: (%{x:.2f}, %{y:.2f}, %{z:.2f})<extra></extra>',
            name='原子'
        ))
        
        # 添加化学键
        for bond in structure.bonds:
            atom1_idx, atom2_idx, bond_type = bond
            atom1 = structure.atoms[atom1_idx]
            atom2 = structure.atoms[atom2_idx]
            
            fig.add_trace(go.Scatter3d(
                x=[atom1['x'], atom2['x']],
                y=[atom1['y'], atom2['y']],
                z=[atom1['z'], atom2['z']],
                mode='lines',
                line=dict(
                    color='gray',
                    width=3 if bond_type == 'double' else 2
                ),
                showlegend=False,
                hoverinfo='skip'
            ))
        
        # 设置布局
        fig.update_layout(
            title=title or f"3D结构: {structure.amino_acid_id}",
            scene=dict(
                xaxis_title='X (Å)',
                yaxis_title='Y (Å)',
                zaxis_title='Z (Å)',
                aspectmode='cube'
            ),
            width=800,
            height=600
        )
        
        return fig
    
    def create_2d_projection(self, structures: List[Structure3D], 
                           method: str = 'pca') -> Optional[Any]:
        """创建2D投影图"""
        if not PLOTLY_AVAILABLE or not SKLEARN_AVAILABLE:
            print("警告: 缺少必要库，无法创建2D投影")
            return None
        
        # 提取特征
        features = []
        labels = []
        
        analyzer = Structure3DAnalyzer()
        
        for structure in structures:
            # 提取几何特征
            geometric_features = analyzer.extract_geometric_features(structure)
            
            feature_vector = [
                structure.molecular_volume,
                structure.surface_area,
                len(structure.atoms),
                len(structure.bonds),
                len(geometric_features.ring_systems),
                len(geometric_features.chiral_centers),
                np.mean(geometric_features.bond_lengths) if geometric_features.bond_lengths else 0,
                np.std(geometric_features.bond_lengths) if geometric_features.bond_lengths else 0,
                np.mean(geometric_features.bond_angles) if geometric_features.bond_angles else 0,
                np.std(geometric_features.bond_angles) if geometric_features.bond_angles else 0
            ]
            
            features.append(feature_vector)
            labels.append(structure.amino_acid_id)
        
        features = np.array(features)
        
        # 降维
        if method == 'pca':
            reducer = PCA(n_components=2)
        elif method == 'tsne':
            reducer = TSNE(n_components=2, random_state=42)
        else:
            raise ValueError(f"不支持的降维方法: {method}")
        
        coords_2d = reducer.fit_transform(features)
        
        # 创建散点图
        fig = px.scatter(
            x=coords_2d[:, 0],
            y=coords_2d[:, 1],
            text=labels,
            title=f"氨基酸化学空间 ({method.upper()})",
            labels={'x': f'{method.upper()}_1', 'y': f'{method.upper()}_2'}
        )
        
        fig.update_traces(textposition="top center")
        fig.update_layout(width=800, height=600)
        
        return fig
    
    def create_property_distribution(self, records: List[AminoAcidRecord], 
                                   property_name: str) -> Optional[Any]:
        """创建性质分布图"""
        if not PLOTLY_AVAILABLE:
            print("警告: Plotly未安装，无法创建分布图")
            return None
        
        # 提取性质值
        values = []
        labels = []
        
        for record in records:
            if property_name == 'molecular_weight':
                values.append(record.molecular_weight)
            elif property_name == 'atom_count':
                values.append(sum(record.atom_composition.values()))
            elif property_name == 'carbon_count':
                values.append(record.atom_composition.get('C', 0))
            elif property_name == 'nitrogen_count':
                values.append(record.atom_composition.get('N', 0))
            elif property_name == 'oxygen_count':
                values.append(record.atom_composition.get('O', 0))
            else:
                continue
            
            labels.append(record.id)
        
        # 创建直方图
        fig = px.histogram(
            x=values,
            nbins=20,
            title=f"{property_name} 分布",
            labels={'x': property_name, 'y': '频次'}
        )
        
        fig.update_layout(width=800, height=400)
        
        return fig

class DrugLikenessAnalyzer:
    """药物相似性分析器"""
    
    def __init__(self):
        # Lipinski五规则阈值
        self.lipinski_thresholds = {
            'molecular_weight': 500,
            'logp': 5,
            'hbd': 5,
            'hba': 10
        }
    
    def calculate_drug_likeness(self, record: AminoAcidRecord, 
                              structure: Optional[Structure3D] = None) -> DrugLikenessScore:
        """计算药物相似性评分"""
        
        # 分子量
        mw = record.molecular_weight
        
        # 估算LogP（简化计算）
        logp = self._estimate_logp(record)
        
        # 氢键供体和受体
        hbd = self._count_hydrogen_bond_donors(record)
        hba = self._count_hydrogen_bond_acceptors(record)
        
        # 拓扑极性表面积（简化估算）
        tpsa = self._estimate_tpsa(record)
        
        # 可旋转键数量
        rotatable_bonds = self._count_rotatable_bonds(record, structure)
        
        # Lipinski违规计数
        violations = 0
        if mw > self.lipinski_thresholds['molecular_weight']:
            violations += 1
        if logp > self.lipinski_thresholds['logp']:
            violations += 1
        if hbd > self.lipinski_thresholds['hbd']:
            violations += 1
        if hba > self.lipinski_thresholds['hba']:
            violations += 1
        
        # 总体评分
        overall_score = max(0, 1.0 - violations * 0.2)
        
        # 药物相似性判断
        drug_like = violations <= 1  # 允许一个违规
        
        return DrugLikenessScore(
            lipinski_violations=violations,
            molecular_weight=mw,
            logp=logp,
            hbd=hbd,
            hba=hba,
            tpsa=tpsa,
            rotatable_bonds=rotatable_bonds,
            overall_score=overall_score,
            drug_like=drug_like
        )
    
    def _estimate_logp(self, record: AminoAcidRecord) -> float:
        """估算LogP值"""
        # 简化的LogP估算
        logp = 0.0
        
        # 基于原子组成
        c_count = record.atom_composition.get('C', 0)
        n_count = record.atom_composition.get('N', 0)
        o_count = record.atom_composition.get('O', 0)
        s_count = record.atom_composition.get('S', 0)
        
        # 疏水性贡献
        logp += c_count * 0.5
        logp += s_count * 0.3
        
        # 亲水性贡献
        logp -= n_count * 0.7
        logp -= o_count * 0.8
        
        # 基于关键特征调整
        if 'carboxyl_group' in record.key_features:
            logp -= 1.5
        if 'amino_group' in record.key_features:
            logp -= 1.0
        if 'hydroxyl_group' in record.key_features:
            logp -= 0.5
        if 'aromatic_ring' in record.key_features:
            logp += 1.0
        
        return round(logp, 2)
    
    def _count_hydrogen_bond_donors(self, record: AminoAcidRecord) -> int:
        """计算氢键供体数量"""
        hbd = 0
        
        # 基于关键特征
        if 'amino_group' in record.key_features:
            hbd += 2  # -NH2
        if 'hydroxyl_group' in record.key_features:
            hbd += 1  # -OH
        if 'carboxyl_group' in record.key_features:
            hbd += 1  # -COOH
        
        return hbd
    
    def _count_hydrogen_bond_acceptors(self, record: AminoAcidRecord) -> int:
        """计算氢键受体数量"""
        hba = 0
        
        # 基于原子组成
        n_count = record.atom_composition.get('N', 0)
        o_count = record.atom_composition.get('O', 0)
        
        hba += n_count + o_count
        
        return hba
    
    def _estimate_tpsa(self, record: AminoAcidRecord) -> float:
        """估算拓扑极性表面积"""
        tpsa = 0.0
        
        # 基于原子组成的简化估算
        n_count = record.atom_composition.get('N', 0)
        o_count = record.atom_composition.get('O', 0)
        
        tpsa += n_count * 23.79  # 氮原子贡献
        tpsa += o_count * 23.06  # 氧原子贡献
        
        return round(tpsa, 2)
    
    def _count_rotatable_bonds(self, record: AminoAcidRecord, 
                             structure: Optional[Structure3D] = None) -> int:
        """计算可旋转键数量"""
        if not structure:
            # 基于分子量的简化估算
            return max(0, int(record.molecular_weight / 50) - 2)
        
        # 基于3D结构的精确计算
        rotatable = 0
        
        for bond in structure.bonds:
            atom1_idx, atom2_idx, bond_type = bond
            atom1 = structure.atoms[atom1_idx]
            atom2 = structure.atoms[atom2_idx]
            
            # 单键且不在环中
            if bond_type == 'single':
                # 简化判断：非氢原子的单键
                if atom1['element'] != 'H' and atom2['element'] != 'H':
                    rotatable += 1
        
        return rotatable
    
    def batch_analyze_drug_likeness(self, records: List[AminoAcidRecord],
                                  structures: Optional[Dict[str, Structure3D]] = None) -> List[DrugLikenessScore]:
        """批量分析药物相似性"""
        results = []
        
        for record in records:
            structure = structures.get(record.id) if structures else None
            score = self.calculate_drug_likeness(record, structure)
            results.append(score)
        
        return results

class ChemicalSpaceAnalyzer:
    """化学空间分析器"""
    
    def __init__(self):
        self.feature_extractors = {
            'basic': self._extract_basic_features,
            'geometric': self._extract_geometric_features,
            'topological': self._extract_topological_features
        }
    
    def analyze_chemical_space(self, records: List[AminoAcidRecord],
                             structures: Optional[Dict[str, Structure3D]] = None,
                             feature_type: str = 'basic') -> Dict[str, Any]:
        """分析化学空间"""
        if not SKLEARN_AVAILABLE:
            raise ImportError("scikit-learn未安装，无法进行化学空间分析")
        
        # 提取特征
        features = []
        labels = []
        
        for record in records:
            structure = structures.get(record.id) if structures else None
            feature_vector = self.feature_extractors[feature_type](record, structure)
            features.append(feature_vector)
            labels.append(record.id)
        
        features = np.array(features)
        
        # PCA分析
        pca = PCA()
        pca_features = pca.fit_transform(features)
        
        # 聚类分析
        optimal_clusters = self._find_optimal_clusters(features)
        kmeans = KMeans(n_clusters=optimal_clusters, random_state=42)
        cluster_labels = kmeans.fit_predict(features)
        
        # DBSCAN聚类
        dbscan = DBSCAN(eps=0.5, min_samples=2)
        dbscan_labels = dbscan.fit_predict(features)
        
        return {
            'features': features,
            'labels': labels,
            'pca': {
                'components': pca_features,
                'explained_variance_ratio': pca.explained_variance_ratio_,
                'cumulative_variance': np.cumsum(pca.explained_variance_ratio_)
            },
            'clustering': {
                'kmeans': {
                    'labels': cluster_labels,
                    'n_clusters': optimal_clusters,
                    'silhouette_score': silhouette_score(features, cluster_labels)
                },
                'dbscan': {
                    'labels': dbscan_labels,
                    'n_clusters': len(set(dbscan_labels)) - (1 if -1 in dbscan_labels else 0),
                    'n_noise': list(dbscan_labels).count(-1)
                }
            }
        }
    
    def _extract_basic_features(self, record: AminoAcidRecord, 
                              structure: Optional[Structure3D] = None) -> List[float]:
        """提取基础特征"""
        features = [
            record.molecular_weight,
            len(record.atom_composition),
            sum(record.atom_composition.values()),
            record.atom_composition.get('C', 0),
            record.atom_composition.get('N', 0),
            record.atom_composition.get('O', 0),
            record.atom_composition.get('S', 0),
            len(record.key_features)
        ]
        
        # 特征编码
        feature_encoding = [0] * 8  # 8种关键特征
        feature_names = ['aromatic_ring', 'carboxyl_group', 'amino_group', 'hydroxyl_group',
                        'methoxy_group', 'double_bond', 'triple_bond', 'sulfur_group']
        
        for i, feature_name in enumerate(feature_names):
            if feature_name in record.key_features:
                feature_encoding[i] = 1
        
        features.extend(feature_encoding)
        
        return features
    
    def _extract_geometric_features(self, record: AminoAcidRecord,
                                  structure: Optional[Structure3D] = None) -> List[float]:
        """提取几何特征"""
        basic_features = self._extract_basic_features(record, structure)
        
        if structure:
            analyzer = Structure3DAnalyzer()
            geometric_features = analyzer.extract_geometric_features(structure)
            
            geometric_stats = [
                structure.molecular_volume,
                structure.surface_area,
                len(structure.bonds),
                len(geometric_features.ring_systems),
                len(geometric_features.chiral_centers),
                np.mean(geometric_features.bond_lengths) if geometric_features.bond_lengths else 0,
                np.std(geometric_features.bond_lengths) if geometric_features.bond_lengths else 0,
                np.mean(geometric_features.bond_angles) if geometric_features.bond_angles else 0,
                np.std(geometric_features.bond_angles) if geometric_features.bond_angles else 0
            ]
            
            basic_features.extend(geometric_stats)
        else:
            # 如果没有3D结构，用零填充
            basic_features.extend([0] * 9)
        
        return basic_features
    
    def _extract_topological_features(self, record: AminoAcidRecord,
                                    structure: Optional[Structure3D] = None) -> List[float]:
        """提取拓扑特征"""
        features = self._extract_basic_features(record, structure)
        
        if structure:
            # 拓扑描述符
            adjacency_matrix = self._build_adjacency_matrix(structure)
            
            # 图论特征
            features.extend([
                np.sum(adjacency_matrix),  # 总边数
                np.max(np.sum(adjacency_matrix, axis=1)),  # 最大度数
                np.mean(np.sum(adjacency_matrix, axis=1)),  # 平均度数
                self._calculate_wiener_index(adjacency_matrix)  # Wiener指数
            ])
        else:
            features.extend([0] * 4)
        
        return features
    
    def _build_adjacency_matrix(self, structure: Structure3D) -> np.ndarray:
        """构建邻接矩阵"""
        n_atoms = len(structure.atoms)
        adjacency = np.zeros((n_atoms, n_atoms))
        
        for bond in structure.bonds:
            atom1_idx, atom2_idx, _ = bond
            adjacency[atom1_idx, atom2_idx] = 1
            adjacency[atom2_idx, atom1_idx] = 1
        
        return adjacency
    
    def _calculate_wiener_index(self, adjacency_matrix: np.ndarray) -> float:
        """计算Wiener指数"""
        # 简化的Wiener指数计算
        n = adjacency_matrix.shape[0]
        distances = np.full((n, n), np.inf)
        
        # 初始化距离矩阵
        for i in range(n):
            distances[i, i] = 0
            for j in range(n):
                if adjacency_matrix[i, j] == 1:
                    distances[i, j] = 1
        
        # Floyd-Warshall算法
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    distances[i, j] = min(distances[i, j], 
                                        distances[i, k] + distances[k, j])
        
        # 计算Wiener指数
        wiener_index = 0
        for i in range(n):
            for j in range(i + 1, n):
                if distances[i, j] != np.inf:
                    wiener_index += distances[i, j]
        
        return wiener_index
    
    def _find_optimal_clusters(self, features: np.ndarray, max_clusters: int = 10) -> int:
        """寻找最优聚类数"""
        if len(features) < 2:
            return 1
        
        max_clusters = min(max_clusters, len(features) - 1)
        silhouette_scores = []
        
        for n_clusters in range(2, max_clusters + 1):
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            cluster_labels = kmeans.fit_predict(features)
            score = silhouette_score(features, cluster_labels)
            silhouette_scores.append(score)
        
        # 返回轮廓系数最高的聚类数
        optimal_idx = np.argmax(silhouette_scores)
        return optimal_idx + 2
