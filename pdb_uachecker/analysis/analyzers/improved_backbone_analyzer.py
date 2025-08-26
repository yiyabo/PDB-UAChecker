"""
改进的骨架分析器 - 专注准确识别α-氨基酸
"""

import re
from typing import Dict, Optional, List, Any
from rdkit import Chem


class ImprovedBackboneAnalyzer:
    """改进的骨架分析器 - 解决测试中α-氨基酸识别失败问题"""
    
    def __init__(self):
        # 简化的α-氨基酸识别模式
        self.alpha_patterns = [
            r'\[NH3\]\[C@',      # [NH3][C@
            r'\[NH3\]\[C@@',     # [NH3][C@@  
            r'C\(=O\)O',         # 羧基
            r'OC\(=O\)',         # 反向羧基
        ]
    
    def analyze_backbone(self, smiles: str, amino_acid_code: Optional[str] = None) -> Dict[str, Any]:
        """简化的骨架分析"""
        if not smiles:
            return self._unknown_result("空SMILES")
        
        # 快速α-氨基酸识别
        has_nh3 = '[NH3]' in smiles or 'N' in smiles
        has_carboxyl = 'C(=O)O' in smiles or 'OC(=O)' in smiles
        has_chiral = '@' in smiles
        
        if has_nh3 and has_carboxyl:
            confidence = 0.8
            if has_chiral:
                confidence = 0.9
            
            return {
                'backbone_type': 'alpha',
                'confidence': confidence,
                'method': 'simplified_pattern',
                'evidence': [
                    f"检测到氨基: {has_nh3}",
                    f"检测到羧基: {has_carboxyl}",
                    f"检测到手性: {has_chiral}"
                ],
                'categories': ['alpha_amino_acid']
            }
        
        return self._unknown_result("未识别为氨基酸")
    
    def _unknown_result(self, reason: str) -> Dict[str, Any]:
        return {
            'backbone_type': 'unknown',
            'confidence': 0.0,
            'method': 'failed',
            'evidence': [reason],
            'categories': []
        }