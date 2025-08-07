"""
重构的分子分析器
专注于分子结构分析，提供清晰的接口
"""

from typing import Dict, List, Optional, Tuple, Any

# 重用现有的分子结构分析器，但提供更清晰的接口
from ..molecular_structure_analyzer import MolecularStructureAnalyzer as _BaseAnalyzer
from ...core.models import MolecularAnalysis


class MolecularAnalyzer:
    """
    重构的分子分析器
    
    专注于分子结构分析：
    1. SMILES验证和解析
    2. 芳香性检测
    3. 环系统分析
    4. 官能团识别
    5. 立体化学分析
    """
    
    def __init__(self):
        """初始化分子分析器"""
        self._base_analyzer = _BaseAnalyzer()
    
    def analyze(self, smiles: str) -> MolecularAnalysis:
        """
        分析分子结构
        
        Args:
            smiles: SMILES字符串
            
        Returns:
            分子分析结果
        """
        # 使用基础分析器
        return self._base_analyzer.analyze(smiles)
    
    def is_valid_smiles(self, smiles: str) -> Tuple[bool, str]:
        """
        验证SMILES有效性
        
        Returns:
            (是否有效, 错误信息)
        """
        try:
            analysis = self.analyze(smiles)
            if analysis.is_valid:
                return True, "Valid SMILES"
            else:
                return False, "Invalid SMILES structure"
        except Exception as e:
            return False, f"SMILES parsing error: {str(e)}"
    
    def detect_aromatic_systems(self, smiles: str) -> Dict[str, any]:
        """
        检测芳香体系
        
        Returns:
            芳香体系信息
        """
        analysis = self.analyze(smiles)
        
        is_aromatic, confidence = self._base_analyzer.is_aromatic(smiles)
        
        return {
            'is_aromatic': is_aromatic,
            'confidence': confidence,
            'aromatic_atoms': analysis.aromatic_atoms,
            'aromatic_systems_count': len([ring for ring in analysis.ring_systems 
                                         if any(atom in analysis.aromatic_atoms 
                                               for atom in ring)]) if analysis.ring_systems else 0
        }
    
    def detect_ring_systems(self, smiles: str) -> Dict[str, any]:
        """
        检测环系统
        
        Returns:
            环系统信息
        """
        analysis = self.analyze(smiles)
        has_rings, ring_count, confidence = self._base_analyzer.has_rings(smiles)
        
        return {
            'has_rings': has_rings,
            'ring_count': ring_count,
            'confidence': confidence,
            'ring_systems': analysis.ring_systems,
            'ring_analysis': self._analyze_ring_types(analysis.ring_systems)
        }
    
    def _analyze_ring_types(self, ring_systems: List[List[int]]) -> Dict[str, any]:
        """分析环的类型"""
        ring_types = {}
        for ring in ring_systems:
            size = len(ring)
            ring_type = f"{size}_membered"
            ring_types[ring_type] = ring_types.get(ring_type, 0) + 1
        
        return ring_types
    
    def identify_functional_groups(self, smiles: str) -> Dict[str, List[int]]:
        """
        识别官能团
        
        Returns:
            官能团字典 {官能团名称: [原子索引列表]}
        """
        analysis = self.analyze(smiles)
        return analysis.functional_groups
    
    def analyze_stereochemistry(self, smiles: str) -> Dict[str, any]:
        """
        分析立体化学
        
        Returns:
            立体化学信息
        """
        return self._base_analyzer.detect_stereochemistry(smiles)
    
    def calculate_molecular_descriptors(self, smiles: str) -> Dict[str, float]:
        """
        计算分子描述符
        
        Returns:
            分子描述符字典
        """
        analysis = self.analyze(smiles)
        return analysis.molecular_descriptors
    
    def compare_structures(self, smiles1: str, smiles2: str) -> Dict[str, any]:
        """
        比较两个分子结构
        
        Returns:
            结构比较结果
        """
        analysis1 = self.analyze(smiles1)
        analysis2 = self.analyze(smiles2)
        
        if not (analysis1.is_valid and analysis2.is_valid):
            return {
                'similarity': 0.0,
                'error': 'One or both SMILES are invalid'
            }
        
        # 简单的结构相似性比较
        similarity_factors = []
        
        # 比较芳香原子数量
        aromatic_similarity = 1.0 - abs(len(analysis1.aromatic_atoms) - len(analysis2.aromatic_atoms)) / max(
            len(analysis1.aromatic_atoms) + len(analysis2.aromatic_atoms), 1)
        similarity_factors.append(aromatic_similarity)
        
        # 比较环数量
        ring_similarity = 1.0 - abs(len(analysis1.ring_systems) - len(analysis2.ring_systems)) / max(
            len(analysis1.ring_systems) + len(analysis2.ring_systems), 1)
        similarity_factors.append(ring_similarity)
        
        # 比较手性中心数量
        chiral_similarity = 1.0 - abs(len(analysis1.chiral_centers) - len(analysis2.chiral_centers)) / max(
            len(analysis1.chiral_centers) + len(analysis2.chiral_centers), 1)
        similarity_factors.append(chiral_similarity)
        
        # 计算综合相似性
        overall_similarity = sum(similarity_factors) / len(similarity_factors)
        
        return {
            'similarity': overall_similarity,
            'detailed_comparison': {
                'aromatic_similarity': aromatic_similarity,
                'ring_similarity': ring_similarity,
                'chiral_similarity': chiral_similarity
            },
            'structure1_features': {
                'aromatic_atoms': len(analysis1.aromatic_atoms),
                'rings': len(analysis1.ring_systems),
                'chiral_centers': len(analysis1.chiral_centers)
            },
            'structure2_features': {
                'aromatic_atoms': len(analysis2.aromatic_atoms), 
                'rings': len(analysis2.ring_systems),
                'chiral_centers': len(analysis2.chiral_centers)
            }
        }
    
    def get_analysis_summary(self, smiles: str) -> Dict[str, any]:
        """
        获取完整的分析摘要
        
        Returns:
            完整分析摘要
        """
        analysis = self.analyze(smiles)
        
        if not analysis.is_valid:
            return {
                'valid': False,
                'error': 'Invalid SMILES structure'
            }
        
        # 芳香性分析
        aromatic_info = self.detect_aromatic_systems(smiles)
        
        # 环系统分析
        ring_info = self.detect_ring_systems(smiles)
        
        # 立体化学分析
        stereo_info = self.analyze_stereochemistry(smiles)
        
        return {
            'valid': True,
            'smiles': smiles,
            'basic_features': {
                'is_aromatic': aromatic_info['is_aromatic'],
                'has_rings': ring_info['has_rings'],
                'aromatic_atoms_count': len(analysis.aromatic_atoms),
                'ring_count': ring_info['ring_count'],
                'chiral_centers_count': len(analysis.chiral_centers),
                'functional_groups_count': len(analysis.functional_groups)
            },
            'detailed_analysis': {
                'aromatic_analysis': aromatic_info,
                'ring_analysis': ring_info,
                'stereochemistry': stereo_info,
                'functional_groups': analysis.functional_groups,
                'backbone_analysis': analysis.backbone_analysis
            },
            'molecular_descriptors': analysis.molecular_descriptors
        }
