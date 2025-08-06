"""
分子式验证器
基于分子式进行氨基酸验证，采用智能fallback策略
"""

from typing import Dict, Any

from .base import BaseVerifier
from ..models import ResidueInfo, AminoAcidInfo, VerificationMethod
from ...utils.chemistry import ChemistryUtils


class MolecularFormulaVerifier(BaseVerifier):
    """分子式验证器"""
    
    def __init__(self, threshold: float = 1.0):
        super().__init__(threshold)
        self.chemistry_utils = ChemistryUtils()
    
    def get_method(self) -> VerificationMethod:
        return VerificationMethod.MOLECULAR_FORMULA
    
    def calculate_score(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> float:
        """
        计算分子式验证分数
        
        验证逻辑：
        1. 优先完整分子式比较（包含氢原子）
        2. Fallback到重原子分子式比较（不含氢原子）
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            验证分数 (0.0, 0.9, 1.0)
        """
        residue_formula = residue.molecular_formula or ""
        candidate_formula = amino_acid.molecular_formula or ""
        
        if not residue_formula or not candidate_formula:
            return 0.0
        
        # 策略1: 完整分子式比较（包含氢原子）
        if residue_formula == candidate_formula:
            return 1.0  # 完美匹配
        
        # 策略2: 重原子分子式比较（fallback）
        residue_no_h = self.chemistry_utils.remove_hydrogen_from_formula(residue_formula)
        candidate_no_h = self.chemistry_utils.remove_hydrogen_from_formula(candidate_formula)
        
        if residue_no_h and candidate_no_h and residue_no_h == candidate_no_h:
            return 0.9  # 稍微降低分数，因为不是完整匹配
        
        return 0.0  # 不匹配
    
    def _get_verification_details(self, residue: ResidueInfo, amino_acid: AminoAcidInfo, score: float) -> Dict[str, Any]:
        """获取分子式验证详情"""
        residue_formula = residue.molecular_formula or ""
        candidate_formula = amino_acid.molecular_formula or ""
        
        details = super()._get_verification_details(residue, amino_acid, score)
        details.update({
            'residue_formula': residue_formula,
            'candidate_formula': candidate_formula,
            'match_type': self._get_match_type(residue_formula, candidate_formula, score)
        })
        
        return details
    
    def _get_match_type(self, residue_formula: str, candidate_formula: str, score: float) -> str:
        """获取匹配类型"""
        if score == 1.0:
            return "complete_match"  # 完整匹配
        elif score == 0.9:
            return "heavy_atom_match"  # 重原子匹配
        else:
            return "no_match"  # 不匹配
    
    def verify_with_hydrogen_handling(self, residue: ResidueInfo, amino_acid: AminoAcidInfo, 
                                    prefer_complete: bool = True) -> Dict[str, Any]:
        """
        带氢原子处理策略的验证
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
            prefer_complete: 是否优先完整匹配
        
        Returns:
            详细验证结果
        """
        residue_formula = residue.molecular_formula or ""
        candidate_formula = amino_acid.molecular_formula or ""
        
        result = {
            'residue_formula': residue_formula,
            'candidate_formula': candidate_formula,
            'complete_match': False,
            'heavy_atom_match': False,
            'final_score': 0.0,
            'match_strategy': 'none'
        }
        
        if not residue_formula or not candidate_formula:
            return result
        
        # 完整匹配检查
        if residue_formula == candidate_formula:
            result.update({
                'complete_match': True,
                'final_score': 1.0,
                'match_strategy': 'complete_formula'
            })
            return result
        
        # 重原子匹配检查
        residue_no_h = self.chemistry_utils.remove_hydrogen_from_formula(residue_formula)
        candidate_no_h = self.chemistry_utils.remove_hydrogen_from_formula(candidate_formula)
        
        if residue_no_h and candidate_no_h and residue_no_h == candidate_no_h:
            result.update({
                'heavy_atom_match': True,
                'final_score': 0.9,
                'match_strategy': 'heavy_atom_formula',
                'residue_heavy_formula': residue_no_h,
                'candidate_heavy_formula': candidate_no_h
            })
        
        return result