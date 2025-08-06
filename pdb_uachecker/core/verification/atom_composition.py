"""
原子组成验证器
基于原子组成进行氨基酸验证，采用智能氢原子处理策略
"""

from typing import Dict, Any

from .base import BaseVerifier
from ..models import ResidueInfo, AminoAcidInfo, VerificationMethod
from ...utils.chemistry import ChemistryUtils


class AtomCompositionVerifier(BaseVerifier):
    """原子组成验证器"""
    
    def __init__(self, threshold: float = 0.9):
        super().__init__(threshold)
        self.chemistry_utils = ChemistryUtils()
    
    def get_method(self) -> VerificationMethod:
        return VerificationMethod.ATOM_COMPOSITION
    
    def calculate_score(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> float:
        """
        计算原子组成验证分数
        
        验证逻辑：
        1. 检测氢原子情况
        2. 选择合适的比较策略
        3. 计算Jaccard相似性
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            验证分数 (0.0-1.0)
        """
        residue_atoms = residue.atom_composition or {}
        candidate_atoms = amino_acid.atom_composition or {}
        
        if not residue_atoms or not candidate_atoms:
            return 0.0
        
        # 智能氢原子处理策略
        comp1, comp2 = self._select_comparison_strategy(residue_atoms, candidate_atoms)
        
        # 计算Jaccard相似性
        return self.chemistry_utils.calculate_jaccard_similarity(comp1, comp2)
    
    def _select_comparison_strategy(self, residue_atoms: Dict[str, int], 
                                  candidate_atoms: Dict[str, int]) -> tuple:
        """
        选择比较策略
        
        Args:
            residue_atoms: 残基原子组成
            candidate_atoms: 候选氨基酸原子组成
        
        Returns:
            (比较用的残基组成, 比较用的候选组成)
        """
        # 检测氢原子情况
        residue_has_h = 'H' in residue_atoms
        candidate_has_h = 'H' in candidate_atoms
        
        if residue_has_h == candidate_has_h:
            # 氢原子情况一致，直接比较
            return residue_atoms, candidate_atoms
        else:
            # 氢原子情况不一致，忽略氢原子比较
            comp1 = {k: v for k, v in residue_atoms.items() if k != 'H'}
            comp2 = {k: v for k, v in candidate_atoms.items() if k != 'H'}
            return comp1, comp2
    
    def _get_verification_details(self, residue: ResidueInfo, amino_acid: AminoAcidInfo, score: float) -> Dict[str, Any]:
        """获取原子组成验证详情"""
        residue_atoms = residue.atom_composition or {}
        candidate_atoms = amino_acid.atom_composition or {}
        
        details = super()._get_verification_details(residue, amino_acid, score)
        
        # 分析氢原子情况
        residue_has_h = 'H' in residue_atoms
        candidate_has_h = 'H' in candidate_atoms
        hydrogen_strategy = self._get_hydrogen_strategy(residue_has_h, candidate_has_h)
        
        # 获取比较用的组成
        comp1, comp2 = self._select_comparison_strategy(residue_atoms, candidate_atoms)
        
        details.update({
            'residue_composition': residue_atoms,
            'candidate_composition': candidate_atoms,
            'comparison_composition_residue': comp1,
            'comparison_composition_candidate': comp2,
            'hydrogen_strategy': hydrogen_strategy,
            'residue_has_hydrogen': residue_has_h,
            'candidate_has_hydrogen': candidate_has_h,
            'jaccard_similarity': score
        })
        
        return details
    
    def _get_hydrogen_strategy(self, residue_has_h: bool, candidate_has_h: bool) -> str:
        """获取氢原子处理策略"""
        if residue_has_h and candidate_has_h:
            return "include_hydrogen"  # 包含氢原子
        elif not residue_has_h and not candidate_has_h:
            return "exclude_hydrogen"  # 排除氢原子
        else:
            return "ignore_hydrogen_mismatch"  # 忽略氢原子不匹配
    
    def calculate_weighted_similarity(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> float:
        """
        计算加权相似性（考虑元素重要性）
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            加权相似性分数 (0.0-1.0)
        """
        residue_atoms = residue.atom_composition or {}
        candidate_atoms = amino_acid.atom_composition or {}
        
        if not residue_atoms or not candidate_atoms:
            return 0.0
        
        # 选择比较策略
        comp1, comp2 = self._select_comparison_strategy(residue_atoms, candidate_atoms)
        
        # 计算加权Jaccard相似性
        return self.chemistry_utils.calculate_weighted_jaccard_similarity(comp1, comp2)
    
    def analyze_composition_difference(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """
        分析原子组成差异
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            差异分析结果
        """
        residue_atoms = residue.atom_composition or {}
        candidate_atoms = amino_acid.atom_composition or {}
        
        # 选择比较策略
        comp1, comp2 = self._select_comparison_strategy(residue_atoms, candidate_atoms)
        
        # 分析差异
        all_elements = set(comp1.keys()) | set(comp2.keys())
        differences = {}
        missing_in_residue = {}
        missing_in_candidate = {}
        extra_in_residue = {}
        extra_in_candidate = {}
        
        for element in all_elements:
            count1 = comp1.get(element, 0)
            count2 = comp2.get(element, 0)
            
            if count1 != count2:
                differences[element] = {
                    'residue': count1,
                    'candidate': count2,
                    'difference': count1 - count2
                }
            
            if count1 == 0 and count2 > 0:
                missing_in_residue[element] = count2
            elif count1 > 0 and count2 == 0:
                missing_in_candidate[element] = count1
            elif count1 > count2:
                extra_in_residue[element] = count1 - count2
            elif count2 > count1:
                extra_in_candidate[element] = count2 - count1
        
        return {
            'differences': differences,
            'missing_in_residue': missing_in_residue,
            'missing_in_candidate': missing_in_candidate,
            'extra_in_residue': extra_in_residue,
            'extra_in_candidate': extra_in_candidate,
            'total_differences': len(differences)
        }