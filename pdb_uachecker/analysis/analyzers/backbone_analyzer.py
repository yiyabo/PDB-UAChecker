"""
骨架分析器
专门负责氨基酸骨架类型的判断（α、β、γ等）
"""

import re
from typing import Dict, List, Optional, Tuple, Any

from ...core.models import AminoAcidInfo


class BackboneAnalyzer:
    """
    氨基酸骨架分析器
    
    判断氨基酸的骨架类型：
    - α-氨基酸（标准）
    - β-氨基酸  
    - γ-氨基酸
    - δ-氨基酸等
    """
    
    def __init__(self):
        """初始化骨架分析器"""
        # 骨架模式定义
        self.backbone_patterns = {
            'alpha': {
                'smiles_patterns': [
                    r'N\[C@@?H?\]\([^)]*\)C\(=O\)O',  # 标准α模式
                    r'N.*\[C@@?H?\].*C\(=O\)O',       # 宽松α模式
                    r'N\[.*?\].*C\(=O\)O'              # 最宽松α模式
                ],
                'name_indicators': ['alpha', 'α'],
                'confidence': 0.95
            },
            'beta': {
                'smiles_patterns': [
                    r'N\[.*?\]C\[.*?\]C\(=O\)O',      # β-氨基酸模式
                    r'NCC\[.*?\]C\(=O\)O',            # 简化β模式
                    r'N.*CC.*C\(=O\)O'                # 宽松β模式
                ],
                'name_indicators': ['beta', 'β', 'b-', '3-amino'],
                'confidence': 0.90
            },
            'gamma': {
                'smiles_patterns': [
                    r'N\[.*?\]CC\[.*?\]C\(=O\)O',     # γ-氨基酸模式
                    r'NCCC\[.*?\]C\(=O\)O',           # 简化γ模式
                    r'N.*CCC.*C\(=O\)O'               # 宽松γ模式
                ],
                'name_indicators': ['gamma', 'γ', 'g-', '4-amino'],
                'confidence': 0.90
            },
            'delta': {
                'smiles_patterns': [
                    r'N\[.*?\]CCC\[.*?\]C\(=O\)O',    # δ-氨基酸模式
                    r'NCCCC\[.*?\]C\(=O\)O',          # 简化δ模式
                    r'N.*CCCC.*C\(=O\)O'              # 宽松δ模式
                ],
                'name_indicators': ['delta', 'δ', 'd-', '5-amino'],
                'confidence': 0.85
            }
        }
        
        # 特殊情况处理
        self.special_cases = {
            'GLY': 'alpha',  # 甘氨酸特殊处理
            'PRO': 'alpha',  # 脯氨酸特殊处理
        }
    
    def analyze(self, smiles: str, amino_acid_code: Optional[str] = None) -> Dict[str, Any]:
        """
        分析氨基酸骨架类型
        
        Args:
            smiles: SMILES字符串
            amino_acid_code: 氨基酸代码（可选）
            
        Returns:
            骨架分析结果
        """
        if not smiles:
            return self._create_unknown_result("No SMILES provided")
        
        # 检查特殊情况
        if amino_acid_code and amino_acid_code in self.special_cases:
            backbone_type = self.special_cases[amino_acid_code]
            return self._create_result(backbone_type, 1.0, "special_case", 
                                     [f"特殊情况: {amino_acid_code}"])
        
        # 按优先级检查各种骨架类型
        results = []
        
        for backbone_type, patterns in self.backbone_patterns.items():
            result = self._match_backbone_patterns(smiles, backbone_type, patterns)
            if result:
                results.append(result)
        
        # 选择最佳匹配
        if results:
            best_result = max(results, key=lambda x: x['confidence'])
            return best_result
        
        # 如果都不匹配，尝试智能推断
        inferred_result = self._intelligent_backbone_inference(smiles)
        if inferred_result:
            return inferred_result
        
        return self._create_unknown_result("No backbone pattern matched")
    
    def _match_backbone_patterns(self, smiles: str, backbone_type: str, patterns: Dict) -> Optional[Dict[str, Any]]:
        """匹配骨架模式"""
        evidence = []
        max_confidence = 0
        
        # SMILES模式匹配
        for pattern in patterns['smiles_patterns']:
            if re.search(pattern, smiles, re.IGNORECASE):
                confidence = patterns['confidence']
                evidence.append(f"SMILES模式匹配: {pattern}")
                max_confidence = max(max_confidence, confidence)
        
        if max_confidence > 0:
            return self._create_result(backbone_type, max_confidence, "pattern_matching", evidence)
        
        return None
    
    def analyze_by_name(self, name: str) -> Dict[str, Any]:
        """
        基于名称分析骨架类型
        
        Args:
            name: 氨基酸名称
            
        Returns:
            骨架分析结果
        """
        name_lower = name.lower()
        
        for backbone_type, patterns in self.backbone_patterns.items():
            for indicator in patterns['name_indicators']:
                if indicator in name_lower:
                    return self._create_result(
                        backbone_type, 0.85, "name_analysis", 
                        [f"名称包含指示符: {indicator}"]
                    )
        
        # 检查数字指示符
        if re.search(r'\b3[-\s]?amino\b', name_lower):
            return self._create_result('beta', 0.80, "name_analysis", ["名称含'3-amino'"])
        elif re.search(r'\b4[-\s]?amino\b', name_lower):
            return self._create_result('gamma', 0.80, "name_analysis", ["名称含'4-amino'"])
        elif re.search(r'\b5[-\s]?amino\b', name_lower):
            return self._create_result('delta', 0.75, "name_analysis", ["名称含'5-amino'"])
        
        return self._create_unknown_result("No backbone indicators in name")
    
    def _intelligent_backbone_inference(self, smiles: str) -> Optional[Dict[str, Any]]:
        """智能骨架推断"""
        evidence = []
        
        # 计算氨基和羧基之间的距离（简化方法）
        amino_pos = self._find_amino_group_position(smiles)
        carboxyl_pos = self._find_carboxyl_group_position(smiles)
        
        if amino_pos is not None and carboxyl_pos is not None:
            # 简单的距离估算
            carbon_count = self._estimate_carbon_chain_length(smiles)
            
            if carbon_count == 1:
                backbone_type = 'alpha'
                confidence = 0.75
                evidence.append("推断: 氨基和羧基间1个碳")
            elif carbon_count == 2:
                backbone_type = 'beta'
                confidence = 0.70
                evidence.append("推断: 氨基和羧基间2个碳")
            elif carbon_count == 3:
                backbone_type = 'gamma'
                confidence = 0.70
                evidence.append("推断: 氨基和羧基间3个碳")
            elif carbon_count == 4:
                backbone_type = 'delta'
                confidence = 0.65
                evidence.append("推断: 氨基和羧基间4个碳")
            else:
                return None
            
            return self._create_result(backbone_type, confidence, "intelligent_inference", evidence)
        
        return None
    
    def _find_amino_group_position(self, smiles: str) -> Optional[int]:
        """查找氨基的位置"""
        # 简化实现：查找N的位置
        for i, char in enumerate(smiles):
            if char == 'N':
                return i
        return None
    
    def _find_carboxyl_group_position(self, smiles: str) -> Optional[int]:
        """查找羧基的位置"""
        # 简化实现：查找C(=O)O模式
        match = re.search(r'C\(=O\)O', smiles)
        if match:
            return match.start()
        return None
    
    def _estimate_carbon_chain_length(self, smiles: str) -> int:
        """估算碳链长度"""
        # 简化方法：计算C的数量，排除羧基碳
        carbon_count = smiles.count('C')
        if 'C(=O)O' in smiles:
            carbon_count -= 1  # 排除羧基碳
        
        # 进一步简化：假设主链碳数
        if carbon_count <= 2:
            return 1  # alpha
        elif carbon_count <= 4:
            return 2  # beta  
        elif carbon_count <= 6:
            return 3  # gamma
        else:
            return 4  # delta
    
    def _create_result(self, backbone_type: str, confidence: float, method: str, evidence: List[str]) -> Dict[str, Any]:
        """创建分析结果"""
        return {
            'backbone_type': backbone_type,
            'confidence': confidence,
            'method': method,
            'evidence': evidence,
            'categories': [f"{backbone_type}_amino_acid"],
            'details': {
                'backbone_classification': backbone_type,
                'analysis_method': method
            }
        }
    
    def _create_unknown_result(self, reason: str) -> Dict[str, Any]:
        """创建未知结果"""
        return {
            'backbone_type': 'unknown',
            'confidence': 0.0,
            'method': 'none',
            'evidence': [reason],
            'categories': ['unknown_backbone'],
            'details': {
                'backbone_classification': 'unknown',
                'reason': reason
            }
        }
    
    def batch_analyze(self, amino_acids: List[Tuple[str, str]]) -> Dict[str, Dict[str, Any]]:
        """
        批量分析骨架类型
        
        Args:
            amino_acids: [(smiles, code), ...]
            
        Returns:
            {code: analysis_result, ...}
        """
        results = {}
        for smiles, code in amino_acids:
            results[code] = self.analyze(smiles, code)
        return results
    
    def get_supported_backbone_types(self) -> List[str]:
        """获取支持的骨架类型"""
        return list(self.backbone_patterns.keys())
    
    def validate_backbone_assignment(self, smiles: str, claimed_type: str) -> Dict[str, Any]:
        """
        验证骨架分配的准确性
        
        Args:
            smiles: SMILES字符串
            claimed_type: 声称的骨架类型
            
        Returns:
            验证结果
        """
        analysis = self.analyze(smiles)
        detected_type = analysis['backbone_type']
        
        is_consistent = (detected_type == claimed_type or 
                        detected_type == 'unknown')  # unknown时不算错误
        
        return {
            'is_consistent': is_consistent,
            'detected_type': detected_type,
            'claimed_type': claimed_type,
            'confidence': analysis['confidence'],
            'evidence': analysis['evidence'],
            'recommendation': 'accept' if is_consistent else 'review_required'
        }
