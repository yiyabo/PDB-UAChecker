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
                    # 高精度α-氨基酸模式
                    r'N\[C@@?H?\]\([^)]*\)C\(=O\)O',                    # 标准α模式: N[C@@H](...)C(=O)O
                    r'\[C@@?H?\]\([^)]*C\(=O\)O[^)]*\)\[NH3?\+?\]',     # 紧密连接α模式: [C@@H](C(=O)O)[NH3]
                    r'\[NH3?\+?\]\[C@@?H?\].*C\(=O\)O',                 # 质子化氨基α模式: [NH3][C@@H]...C(=O)O
                    r'C\[C@@?H?\].*C\(=O\)O.*\[NH3?\+?\]',              # 侧链-手性碳-羧基-氨基模式
                    r'.*\[C@@?H?\]\([^)]*C\(=O\)O[^)]*\).*',           # 手性碳包含羧基的通用模式
                    
                    # 支持各种离子形式和手性标记
                    r'\[NH3?\+?\].*\[C@@?H?\].*C\(=O\).*O[-]?',        # 离子形式支持
                    r'N.*\[C@@?H?\].*C\(=O\)O',                        # 基础α模式（宽松）
                    r'N.*C\(.*\)C\(=O\)O',                             # 无手性标记的α模式（如丙氨酸简化形式）
                    
                    # 反向连接和特殊排列
                    r'C\(=O\)O.*\[C@@?H?\].*N',                        # 反向排列
                    r'.*C\[C@@?H?\].*C\(=O\)O',                        # 常见取代α-氨基酸模式（如苯丙氨酸类）
                ],
                'name_indicators': ['alpha', 'α'],
                'confidence': 0.95
            },
            'beta': {
                'smiles_patterns': [
                    r'NCC\[C@@?H?\].*C\(=O\)O',        # 严格β-氨基酸模式
                    r'\[NH3?\+?\]CC\[.*?\]C\(=O\)O',   # 质子化氨基β模式
                    r'NCC\[.*?\]C\(=O\)O',            # 简化β模式
                    # 移除过宽松的模式 r'N.*CC.*C\(=O\)O' 避免误匹配
                ],
                'name_indicators': ['beta', 'β', 'b-', '3-amino'],
                'confidence': 0.90
            },
            'gamma': {
                'smiles_patterns': [
                    r'NCCCC\(=O\)O',                  # 严格γ模式: NCCCC(=O)O  
                    r'N\[.*?\]CCC\[.*?\]C\(=O\)O',   # 带取代基的γ模式
                    # 移除过宽松的模式 r'N.*CCC.*C\(=O\)O' 避免匹配芳香环
                ],
                'name_indicators': ['gamma', 'γ', 'g-', '4-amino'],
                'confidence': 0.90
            },
            'delta': {
                'smiles_patterns': [
                    r'NCCCCC\(=O\)O',                 # 严格δ模式: NCCCCC(=O)O
                    r'N\[.*?\]CCCC\[.*?\]C\(=O\)O',  # 带取代基的δ模式  
                    # 移除过宽松的模式 r'N.*CCCC.*C\(=O\)O' 避免匹配芳香环
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
        
        # 按优先级检查各种骨架类型（β、γ、δ优先于α，避免α的宽松模式误匹配）
        results = []
        
        # 优先级顺序：先检查更具体的β、γ、δ模式，最后检查α
        priority_order = ['beta', 'gamma', 'delta', 'alpha']
        
        for backbone_type in priority_order:
            if backbone_type in self.backbone_patterns:
                patterns = self.backbone_patterns[backbone_type]
                result = self._match_backbone_patterns(smiles, backbone_type, patterns)
                if result:
                    results.append(result)
        
        # 选择最佳匹配（由于已按优先级排序，直接选第一个）
        if results:
            # 由于results是按优先级顺序(beta, gamma, delta, alpha)构建的
            # 第一个结果就是最高优先级的匹配
            return results[0]
        
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
        """智能骨架推断（改进版）"""
        evidence = []
        
        # 首先尝试基于化学结构的智能识别
        structure_analysis = self._analyze_chemical_structure(smiles)
        if structure_analysis:
            return structure_analysis
        
        # 计算氨基和羧基之间的距离（改进方法）
        amino_pos = self._find_amino_group_position(smiles)
        carboxyl_pos = self._find_carboxyl_group_position(smiles)
        
        if amino_pos is not None and carboxyl_pos is not None:
            # 使用改进的距离分析
            distance_analysis = self._improved_distance_analysis(smiles)
            
            if distance_analysis['backbone_type']:
                return self._create_result(
                    distance_analysis['backbone_type'], 
                    distance_analysis['confidence'],
                    "improved_intelligent_inference", 
                    distance_analysis['evidence']
                )
        
        # 最后的退化：基于SMILES特征进行保守估计
        fallback_analysis = self._fallback_structure_analysis(smiles)
        if fallback_analysis:
            return fallback_analysis
        
        return None
    
    def _analyze_chemical_structure(self, smiles: str) -> Optional[Dict[str, Any]]:
        """基于化学结构的智能分析"""
        evidence = []
        
        # 识别典型的氨基酸骨架特征
        if re.search(r'\[C@@?H?\].*C\(=O\)O', smiles):
            # 发现手性碳，很可能是α-氨基酸
            if self._has_amino_on_chiral_carbon(smiles):
                evidence.append("识别到手性碳直接连接氨基和羧基的α-氨基酸结构")
                return self._create_result('alpha', 0.85, "structure_analysis", evidence)
        
        # 识别简单的氨基酸模式
        if re.search(r'NC.*C\(=O\)O', smiles):
            # N-C-COOH 基本结构，很可能是α-氨基酸
            if not self._has_extended_chain(smiles):
                evidence.append("识别到N-C-COOH基本α-氨基酸结构")
                return self._create_result('alpha', 0.80, "structure_analysis", evidence)
        
        return None
    
    def _has_amino_on_chiral_carbon(self, smiles: str) -> bool:
        """检查氨基是否直接连在手性碳上"""
        # 检查氨基和手性碳的连接模式
        patterns = [
            r'N.*\[C@@?H?\]',               # N...[C@@H]
            r'\[C@@?H?\].*N',               # [C@@H]...N  
            r'\[NH3?\+?\].*\[C@@?H?\]',     # [NH3]...[C@@H]
            r'\[C@@?H?\].*\[NH3?\+?\]'      # [C@@H]...[NH3]
        ]
        
        for pattern in patterns:
            if re.search(pattern, smiles):
                return True
        return False
    
    def _has_extended_chain(self, smiles: str) -> bool:
        """检查是否有延长的碳链（β、γ等）"""
        # 检查明显的β、γ延长链特征
        extended_patterns = [
            r'NCC+.*C\(=O\)O',     # NCC...C(=O)O (β或更长)
            r'NCCC+.*C\(=O\)O',    # NCCC...C(=O)O (γ或更长)
        ]
        
        for pattern in extended_patterns:
            if re.search(pattern, smiles):
                return True
        return False
    
    def _improved_distance_analysis(self, smiles: str) -> Dict[str, Any]:
        """改进的距离分析"""
        # 基于精确距离分析的结果
        precise_distance = self._analyze_amino_carboxyl_distance(smiles)
        
        if precise_distance:
            distance_map = {
                1: {'type': 'alpha', 'confidence': 0.85},
                2: {'type': 'beta', 'confidence': 0.80}, 
                3: {'type': 'gamma', 'confidence': 0.75},
                4: {'type': 'delta', 'confidence': 0.70}
            }
            
            if precise_distance in distance_map:
                mapping = distance_map[precise_distance]
                return {
                    'backbone_type': mapping['type'],
                    'confidence': mapping['confidence'],
                    'evidence': [f"精确距离分析：氨基到羧基{precise_distance}个碳原子"]
                }
        
        # 如果精确分析失败，使用改进的估算方法
        estimated_distance = self._improved_carbon_chain_estimation(smiles)
        if estimated_distance:
            return {
                'backbone_type': estimated_distance['type'],
                'confidence': estimated_distance['confidence'] * 0.9,  # 估算的置信度略低
                'evidence': estimated_distance['evidence']
            }
        
        return {'backbone_type': None, 'confidence': 0, 'evidence': ["无法确定骨架类型"]}
    
    def _improved_carbon_chain_estimation(self, smiles: str) -> Optional[Dict[str, Any]]:
        """改进的碳链估算方法"""
        evidence = []
        
        # 保守的α-氨基酸判断：优先考虑α
        if self._looks_like_alpha_amino_acid(smiles):
            evidence.append("结构特征表明这是α-氨基酸")
            return {
                'type': 'alpha',
                'confidence': 0.75,
                'evidence': evidence
            }
        
        # 明显的β-氨基酸特征
        if self._looks_like_beta_amino_acid(smiles):
            evidence.append("结构特征表明这是β-氨基酸") 
            return {
                'type': 'beta',
                'confidence': 0.70,
                'evidence': evidence
            }
        
        return None
    
    def _looks_like_alpha_amino_acid(self, smiles: str) -> bool:
        """判断是否看起来像α-氨基酸"""
        # 大多数氨基酸都是α-氨基酸，如果没有明显的β、γ特征，倾向于α
        alpha_indicators = [
            r'\[C@@?H?\]',                    # 有手性碳
            r'NC\([^)]*\)C\(=O\)O',          # N-C(...)-COOH结构
            r'N.*C.*C\(=O\)O',               # 基础N-C-COOH结构
        ]
        
        # 排除明显的非α特征
        non_alpha_indicators = [
            r'NCC\[C@@?H?\]',                # β-氨基酸特征
            r'NCCCC\(=O\)O',                 # γ-氨基酸特征
            r'NCCCCC\(=O\)O',                # δ-氨基酸特征
        ]
        
        # 首先检查是否有非α特征
        for pattern in non_alpha_indicators:
            if re.search(pattern, smiles):
                return False
        
        # 然后检查是否有α特征
        for pattern in alpha_indicators:
            if re.search(pattern, smiles):
                return True
        
        return False
    
    def _looks_like_beta_amino_acid(self, smiles: str) -> bool:
        """判断是否看起来像β-氨基酸"""
        beta_indicators = [
            r'NCC\[C@@?H?\]',               # NCC[C@@H] - 明显β特征
            r'\[NH3?\+?\]CC\[',             # [NH3]CC[ - 质子化β
            r'NCCC\(=O\)O',                 # NCCC(=O)O - 简单β
        ]
        
        for pattern in beta_indicators:
            if re.search(pattern, smiles):
                return True
        return False
    
    def _fallback_structure_analysis(self, smiles: str) -> Optional[Dict[str, Any]]:
        """最后的退化结构分析"""
        # 如果所有方法都失败，基于最基本的启发式规则
        if 'N' in smiles and 'C(=O)O' in smiles:
            # 至少有氨基和羧基，默认假设为α-氨基酸（最常见）
            return self._create_result('alpha', 0.60, "fallback_heuristic", 
                                     ["启发式判断：含有氨基和羧基，默认为α-氨基酸"])
        
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
        """估算主链碳链长度（改进版本）"""
        # 首先尝试精确分析氨基和羧基的连接
        amino_carboxyl_distance = self._analyze_amino_carboxyl_distance(smiles)
        if amino_carboxyl_distance is not None:
            return amino_carboxyl_distance
        
        # 降级到简化分析，但改进计数逻辑
        # 排除明显的非主链碳（甲氧基、苯环等）
        cleaned_smiles = self._remove_side_chain_carbons(smiles)
        
        carbon_count = cleaned_smiles.count('C')
        if 'C(=O)O' in cleaned_smiles:
            carbon_count -= 1  # 排除羧基碳
        
        # 更保守的分类，偏向α
        if carbon_count <= 1:
            return 1  # alpha
        elif carbon_count <= 2:
            return 2  # beta  
        elif carbon_count <= 3:
            return 3  # gamma
        else:
            return 4  # delta
    
    def _analyze_amino_carboxyl_distance(self, smiles: str) -> Optional[int]:
        """分析氨基和羧基之间的精确距离"""
        # 增强的α-氨基酸识别模式
        alpha_patterns = [
            # 直接连接模式 - 氨基直接连在α-碳上
            r'N\[C@@?H?\].*C\(=O\)O',                            # N[C@@H]...C(=O)O
            r'\[NH3?\+?\]\[C@@?H?\].*C\(=O\)O',                  # [NH3][C@@H]...C(=O)O 
            r'\[C@@?H?\]\([^)]*C\(=O\)O[^)]*\)\[NH3?\+?\]',     # [C@@H](C(=O)O)[NH3]
            r'C\(=O\)O\[NH3?\+?\]',                              # C(=O)O[NH3] (反向连接)
            
            # 侧链取代的α-氨基酸（典型的氨基酸结构）
            r'.*\[C@@?H?\]\([^)]*C\(=O\)O[^)]*\).*',            # 手性碳包含羧基的模式
            r'.*C\[C@@?H?\].*C\(=O\)O',                         # 取代基-手性碳-羧基（如苯丙氨酸）
            r'N.*C\(.*R.*\)C\(=O\)O',                           # 通用α-氨基酸模式 N-C(R)-COOH
            
            # 无手性标记的α-氨基酸
            r'NC\([^)]*\)C\(=O\)O',                             # NC(...)C(=O)O
            r'N.*C.*C\(=O\)O',                                  # N...C...C(=O)O (最宽松α模式，但排除明显的β、γ模式)
        ]
        
        for pattern in alpha_patterns:
            if re.search(pattern, smiles, re.IGNORECASE):
                # 进一步验证不是β或γ模式
                if not self._is_beta_or_gamma_pattern(smiles):
                    return 1  # α位
        
        # 检查β-氨基酸模式（氨基-C-C-羧基）
        beta_patterns = [
            r'NCC\[C@@?H?\].*C\(=O\)O',         # NCC[C@@H]...C(=O)O
            r'\[NH3?\+?\]CC\[.*?\]C\(=O\)O',    # [NH3]CC[...]C(=O)O
            r'NCCC\(=O\)O',                     # 简单β模式: NCCC(=O)O
            r'NCC.*C\(=O\)O',                   # 广义β模式（但要排除α模式）
        ]
        
        for pattern in beta_patterns:
            if re.search(pattern, smiles, re.IGNORECASE):
                return 2  # β位
        
        # 检查γ-氨基酸模式
        gamma_patterns = [
            r'NCCCC\(=O\)O',                    # NCCCC(=O)O
            r'N.*CCC.*C\(=O\)O'                 # N...CCC...C(=O)O（但要避免匹配芳香环）
        ]
        
        for pattern in gamma_patterns:
            if re.search(pattern, smiles, re.IGNORECASE):
                # 确保不是芳香环中的ccc
                if not re.search(r'c{3,}', smiles):  # 避免匹配芳香环
                    return 3  # γ位
        
        return None  # 无法确定
    
    def _is_beta_or_gamma_pattern(self, smiles: str) -> bool:
        """检查是否是明显的β或γ模式，避免α模式误匹配"""
        # 明显的β模式特征
        beta_indicators = [
            r'NCC\[C@@?H?\]',      # NCC[C@@H] - β氨基酸特征
            r'\[NH3\]CC\[',        # [NH3]CC[ - 质子化β氨基酸
            r'NCCC\(=O\)O'         # NCCC(=O)O - 简单β氨基酸
        ]
        
        for pattern in beta_indicators:
            if re.search(pattern, smiles):
                return True
        
        # 明显的γ模式特征
        gamma_indicators = [
            r'NCCCC\(=O\)O',       # NCCCC(=O)O - 简单γ氨基酸
            r'N\[.*?\]CCC\[.*?\]C\(=O\)O'  # 复杂γ氨基酸
        ]
        
        for pattern in gamma_indicators:
            if re.search(pattern, smiles):
                return True
        
        return False
    
    def _remove_side_chain_carbons(self, smiles: str) -> str:
        """移除明显的侧链碳原子"""
        cleaned = smiles
        
        # 移除甲氧基等小取代基中的碳
        cleaned = re.sub(r'COc', 'Oc', cleaned)  # COc -> Oc
        cleaned = re.sub(r'OC(?![(@])', 'O', cleaned)  # OC -> O (不在立体标记中)
        
        # 注意：不移除苯环碳，因为它们在SMILES中是小写c
        # 苯环碳不会影响主链计数
        
        return cleaned
    
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
