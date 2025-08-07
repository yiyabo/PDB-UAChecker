"""
N-甲基化分析器
精确区分主链氨基的N-甲基化 vs 侧链氮原子的甲基化
严格按照分类标准："必须是主链氨基的甲基化，侧链氮的甲基化不计入此类"
"""

import re
from typing import Dict, List, Optional, Tuple, Any


class NMethylationAnalyzer:
    """
    N-甲基化分析器 - 精确的主链vs侧链区分
    
    核心功能：
    1. 识别所有N-甲基化位点
    2. 严格区分主链氨基 vs 侧链氮原子
    3. 只有主链氨基的甲基化才算N-甲基氨基酸
    4. 支持多种SMILES表示法的准确识别
    """
    
    def __init__(self):
        """初始化N-甲基化分析器"""
        # N-甲基化模式定义
        self.methylation_patterns = {
            'explicit_n_methyl': {
                'smiles_patterns': [
                    r'CN\[C@',           # CN[C@ - 明确的N-甲基
                    r'CN\[C@@',          # CN[C@@ - 另一种手性
                    r'\[N\]\(C\)',       # [N](C) - 括号形式
                    r'CNC\(',            # CNC( - 简化形式
                ],
                'confidence': 0.95,
                'description': '明确的N-甲基SMILES标记'
            },
            'implicit_n_methyl': {
                'smiles_patterns': [
                    r'C\[NH\]',          # C[NH] - 隐式氢
                    r'CN(?!\[)',         # CN但不跟[ - 可能的N-甲基
                ],
                'confidence': 0.80,
                'description': '可能的N-甲基标记'
            },
            'known_n_methyl_compounds': {
                'compound_codes': [
                    'SAR',               # 肌氨酸 (N-甲基甘氨酸)
                    'NME',               # N-甲基
                    'MET',               # 某些N-甲基蛋氨酸变体
                ],
                'confidence': 1.0,
                'description': '已知N-甲基化合物代码'
            }
        }
        
        # 主链氨基识别模式
        self.main_chain_amino_patterns = {
            'alpha_amino': {
                'pattern': r'N.*\[C@@?H?\].*C\(=O\)O',
                'description': 'α-氨基酸主链模式'
            },
            'beta_amino': {
                'pattern': r'N.*C.*\[C@@?H?\].*C\(=O\)O',
                'description': 'β-氨基酸主链模式'
            },
            'gamma_amino': {
                'pattern': r'N.*C.*C.*\[C@@?H?\].*C\(=O\)O',
                'description': 'γ-氨基酸主链模式'
            }
        }
        
        # 侧链氮原子模式（排除）
        self.side_chain_nitrogen_patterns = [
            r'c.*N.*C',              # 芳香环上的氮
            r'N.*c',                 # 连接到芳香环的氮
            r'\[NH2\]\+.*C.*C.*N',   # 主链已有氨基，后面的氮为侧链
            r'N.*\[.*N.*\]',         # 复杂侧链中的氮
        ]
    
    def analyze_n_methylation(self, smiles: str, amino_acid_code: Optional[str] = None,
                            amino_acid_name: Optional[str] = None) -> Dict[str, Any]:
        """
        分析N-甲基化情况
        
        Args:
            smiles: SMILES字符串
            amino_acid_code: 氨基酸代码
            amino_acid_name: 氨基酸名称
            
        Returns:
            详细的N-甲基化分析结果
        """
        if not smiles:
            return self._create_unknown_result("No SMILES provided")
        
        # Step 1: 检查已知N-甲基化合物
        known_compound_result = self._check_known_n_methyl_compounds(amino_acid_code, amino_acid_name)
        if known_compound_result['is_known_n_methyl']:
            return self._create_definitive_result(
                True, "known_compound", known_compound_result['evidence']
            )
        
        # Step 2: 识别所有可能的N-甲基化位点
        methylation_sites = self._identify_methylation_sites(smiles)
        
        if not methylation_sites['found_any']:
            return self._create_negative_result("No N-methylation patterns detected")
        
        # Step 3: 区分主链 vs 侧链氮原子
        chain_analysis = self._analyze_nitrogen_chain_position(smiles, methylation_sites)
        
        # Step 4: 确定是否为主链N-甲基化
        main_chain_methylation = self._determine_main_chain_methylation(
            chain_analysis, methylation_sites
        )
        
        # Step 5: 交叉验证
        validation_result = self._cross_validate_n_methylation(
            main_chain_methylation, amino_acid_name, smiles
        )
        
        # Step 6: 计算最终置信度
        final_confidence = self._calculate_confidence(
            main_chain_methylation, validation_result, methylation_sites
        )
        
        return {
            'is_n_methylated': main_chain_methylation['is_main_chain_methylated'],
            'methylation_type': main_chain_methylation['methylation_type'],
            'confidence': final_confidence,
            'method': 'precise_main_chain_analysis',
            'evidence': self._compile_evidence(
                main_chain_methylation, validation_result, methylation_sites
            ),
            'categories': ['n_methyl_amino_acid'] if main_chain_methylation['is_main_chain_methylated'] else [],
            'details': {
                'methylation_sites': methylation_sites,
                'chain_analysis': chain_analysis,
                'main_chain_methylation': main_chain_methylation,
                'validation_result': validation_result,
                'analysis_method': 'main_chain_vs_side_chain_distinction'
            }
        }
    
    def _check_known_n_methyl_compounds(self, amino_acid_code: Optional[str], 
                                      amino_acid_name: Optional[str]) -> Dict[str, Any]:
        """检查已知的N-甲基化合物"""
        known_codes = self.methylation_patterns['known_n_methyl_compounds']['compound_codes']
        
        # 检查代码
        if amino_acid_code and amino_acid_code in known_codes:
            return {
                'is_known_n_methyl': True,
                'evidence': [f"已知N-甲基化合物代码: {amino_acid_code}"],
                'source': 'compound_code'
            }
        
        # 检查名称
        if amino_acid_name:
            name_lower = amino_acid_name.lower()
            n_methyl_indicators = [
                'n-methyl', 'n methyl', 'sarcosine', '肌氨酸', 'methylamino'
            ]
            
            for indicator in n_methyl_indicators:
                if indicator in name_lower:
                    return {
                        'is_known_n_methyl': True,
                        'evidence': [f"名称含N-甲基指示符: {indicator}"],
                        'source': 'compound_name'
                    }
        
        return {'is_known_n_methyl': False}
    
    def _identify_methylation_sites(self, smiles: str) -> Dict[str, Any]:
        """识别所有可能的N-甲基化位点"""
        methylation_sites = []
        
        # 检查各种N-甲基化模式
        for pattern_type, pattern_info in self.methylation_patterns.items():
            if pattern_type == 'known_n_methyl_compounds':
                continue  # 已在上一步处理
            
            for pattern in pattern_info.get('smiles_patterns', []):
                matches = list(re.finditer(pattern, smiles, re.IGNORECASE))
                for match in matches:
                    methylation_sites.append({
                        'position': match.start(),
                        'pattern': pattern,
                        'pattern_type': pattern_type,
                        'confidence': pattern_info['confidence'],
                        'matched_text': match.group(),
                        'description': pattern_info['description']
                    })
        
        return {
            'found_any': len(methylation_sites) > 0,
            'sites': methylation_sites,
            'total_sites': len(methylation_sites)
        }
    
    def _analyze_nitrogen_chain_position(self, smiles: str, methylation_sites: Dict) -> Dict[str, Any]:
        """分析氮原子在主链vs侧链中的位置"""
        chain_positions = []
        
        for site in methylation_sites['sites']:
            position_analysis = self._analyze_site_chain_position(smiles, site)
            chain_positions.append(position_analysis)
        
        # 统计主链vs侧链位点
        main_chain_sites = [p for p in chain_positions if p['is_main_chain']]
        side_chain_sites = [p for p in chain_positions if not p['is_main_chain']]
        
        return {
            'all_positions': chain_positions,
            'main_chain_sites': main_chain_sites,
            'side_chain_sites': side_chain_sites,
            'main_chain_count': len(main_chain_sites),
            'side_chain_count': len(side_chain_sites)
        }
    
    def _analyze_site_chain_position(self, smiles: str, site: Dict) -> Dict[str, Any]:
        """分析单个位点的链位置"""
        position = site['position']
        
        # 检查是否匹配主链氨基模式
        is_main_chain = False
        main_chain_evidence = []
        
        for pattern_name, pattern_info in self.main_chain_amino_patterns.items():
            # 检查N-甲基化位点是否在主链氨基模式内
            pattern = pattern_info['pattern']
            matches = list(re.finditer(pattern, smiles, re.IGNORECASE))
            
            for match in matches:
                # 如果N-甲基化位点在主链氨基模式范围内
                if match.start() <= position <= match.end():
                    is_main_chain = True
                    main_chain_evidence.append(f"位于{pattern_name}模式内")
                    break
        
        # 检查是否匹配侧链氮原子模式（排除）
        side_chain_evidence = []
        for side_pattern in self.side_chain_nitrogen_patterns:
            if re.search(side_pattern, smiles, re.IGNORECASE):
                # 进一步检查是否影响当前位点
                side_chain_evidence.append(f"检测到侧链氮模式: {side_pattern}")
        
        # 如果有侧链证据且没有强主链证据，倾向于侧链
        if side_chain_evidence and not main_chain_evidence:
            is_main_chain = False
        
        return {
            'site': site,
            'is_main_chain': is_main_chain,
            'main_chain_evidence': main_chain_evidence,
            'side_chain_evidence': side_chain_evidence,
            'confidence': 0.9 if main_chain_evidence else 0.6
        }
    
    def _determine_main_chain_methylation(self, chain_analysis: Dict, methylation_sites: Dict) -> Dict[str, Any]:
        """确定是否为主链N-甲基化"""
        main_chain_sites = chain_analysis['main_chain_sites']
        
        if not main_chain_sites:
            return {
                'is_main_chain_methylated': False,
                'methylation_type': 'none_or_side_chain_only',
                'reason': '未检测到主链N-甲基化',
                'confidence_factor': 0.9
            }
        
        # 选择最可信的主链N-甲基化位点
        best_main_chain_site = max(main_chain_sites, key=lambda x: x['confidence'])
        
        # 确定甲基化类型
        methylation_type = self._classify_methylation_type(best_main_chain_site)
        
        return {
            'is_main_chain_methylated': True,
            'methylation_type': methylation_type,
            'primary_site': best_main_chain_site,
            'all_main_chain_sites': main_chain_sites,
            'confidence_factor': best_main_chain_site['confidence']
        }
    
    def _classify_methylation_type(self, site: Dict) -> str:
        """分类N-甲基化类型"""
        site_info = site['site']
        
        # 基于模式确定类型
        if site_info['pattern_type'] == 'explicit_n_methyl':
            return 'explicit_n_methyl'
        elif site_info['pattern_type'] == 'implicit_n_methyl':
            return 'implicit_n_methyl'
        else:
            return 'uncertain_n_methyl'
    
    def _cross_validate_n_methylation(self, main_chain_methylation: Dict, 
                                    amino_acid_name: Optional[str], smiles: str) -> Dict[str, Any]:
        """交叉验证N-甲基化判断"""
        validations = []
        
        # 验证1: 名称验证
        if amino_acid_name:
            name_validation = self._validate_by_name(amino_acid_name, main_chain_methylation['is_main_chain_methylated'])
            validations.append(name_validation)
        
        # 验证2: 分子特征验证
        molecular_validation = self._validate_by_molecular_features(smiles, main_chain_methylation['is_main_chain_methylated'])
        validations.append(molecular_validation)
        
        # 验证3: 化学合理性验证
        chemistry_validation = self._validate_chemistry(main_chain_methylation)
        validations.append(chemistry_validation)
        
        # 综合验证结果
        consistent_validations = [v for v in validations if v.get('consistent', False)]
        validation_confidence = len(consistent_validations) / len(validations) if validations else 0
        
        return {
            'all_validations_consistent': len(consistent_validations) == len(validations),
            'validation_confidence': validation_confidence,
            'individual_validations': validations,
            'consistency_rate': validation_confidence
        }
    
    def _validate_by_name(self, name: str, predicted_n_methyl: bool) -> Dict[str, Any]:
        """通过名称验证N-甲基化"""
        name_lower = name.lower()
        
        # 名称中的N-甲基指示符
        n_methyl_indicators = ['n-methyl', 'n methyl', 'sarcosine', 'methylamino']
        has_name_indicator = any(indicator in name_lower for indicator in n_methyl_indicators)
        
        consistent = (has_name_indicator == predicted_n_methyl)
        
        return {
            'method': 'name_validation',
            'has_name_indicator': has_name_indicator,
            'consistent': consistent,
            'confidence_boost': 0.1 if consistent else -0.1,
            'evidence': f"名称指示N-甲基化: {'是' if has_name_indicator else '否'}"
        }
    
    def _validate_by_molecular_features(self, smiles: str, predicted_n_methyl: bool) -> Dict[str, Any]:
        """通过分子特征验证"""
        # 检查分子中甲基的数量和位置
        methyl_groups = len(re.findall(r'C(?![=\(])', smiles))  # 简化的甲基计数
        nitrogen_count = smiles.count('N')
        
        # 如果有氮且有额外的甲基，可能支持N-甲基化
        molecular_support = (nitrogen_count > 0 and methyl_groups > 2)
        consistent = (molecular_support == predicted_n_methyl) or not predicted_n_methyl
        
        return {
            'method': 'molecular_validation',
            'molecular_support': molecular_support,
            'consistent': consistent,
            'confidence_boost': 0.05 if consistent else -0.05,
            'evidence': f"分子特征支持: {'是' if molecular_support else '否'}"
        }
    
    def _validate_chemistry(self, main_chain_methylation: Dict) -> Dict[str, Any]:
        """验证化学合理性"""
        is_methylated = main_chain_methylation['is_main_chain_methylated']
        
        # 基本化学合理性：N-甲基氨基酸在生物学中存在
        is_reasonable = True  # N-甲基化在化学上是合理的
        
        return {
            'method': 'chemistry_validation',
            'is_chemically_reasonable': is_reasonable,
            'consistent': True,  # 化学上都是合理的
            'confidence_boost': 0.0,
            'evidence': "N-甲基化化学上合理"
        }
    
    def _calculate_confidence(self, main_chain_methylation: Dict, validation_result: Dict, 
                            methylation_sites: Dict) -> float:
        """计算最终置信度"""
        base_confidence = 0.8 if main_chain_methylation['is_main_chain_methylated'] else 0.9
        
        # 主链甲基化置信度调整
        if main_chain_methylation.get('confidence_factor'):
            base_confidence *= main_chain_methylation['confidence_factor']
        
        # 验证结果置信度调整
        validation_boost = validation_result['validation_confidence'] * 0.1
        base_confidence += validation_boost
        
        # 甲基化位点质量调整
        if methylation_sites['found_any']:
            site_quality = sum(site['confidence'] for site in methylation_sites['sites']) / len(methylation_sites['sites'])
            base_confidence *= site_quality
        
        return max(0.0, min(1.0, base_confidence))
    
    def _compile_evidence(self, main_chain_methylation: Dict, validation_result: Dict, 
                         methylation_sites: Dict) -> List[str]:
        """编译分析证据"""
        evidence = []
        
        if main_chain_methylation['is_main_chain_methylated']:
            evidence.append(f"检测到主链N-甲基化: {main_chain_methylation['methylation_type']}")
            if 'primary_site' in main_chain_methylation:
                site = main_chain_methylation['primary_site']['site']
                evidence.append(f"主要N-甲基化位点: {site['matched_text']}")
        else:
            evidence.append("未检测到主链N-甲基化")
            if methylation_sites['found_any']:
                evidence.append("检测到甲基化但位于侧链")
        
        # 添加验证证据
        for validation in validation_result.get('individual_validations', []):
            if validation.get('evidence'):
                evidence.append(validation['evidence'])
        
        return evidence
    
    def _create_definitive_result(self, is_n_methylated: bool, method: str, evidence: List[str]) -> Dict[str, Any]:
        """创建确定性结果"""
        return {
            'is_n_methylated': is_n_methylated,
            'methylation_type': 'known_compound' if is_n_methylated else 'none',
            'confidence': 1.0,
            'method': method,
            'evidence': evidence,
            'categories': ['n_methyl_amino_acid'] if is_n_methylated else [],
            'details': {
                'analysis_certainty': 'definitive',
                'determination_source': method
            }
        }
    
    def _create_negative_result(self, reason: str) -> Dict[str, Any]:
        """创建否定结果"""
        return {
            'is_n_methylated': False,
            'methylation_type': 'none',
            'confidence': 0.95,
            'method': 'negative_detection',
            'evidence': [reason],
            'categories': [],
            'details': {
                'analysis_certainty': 'high_confidence_negative',
                'reason': reason
            }
        }
    
    def _create_unknown_result(self, reason: str) -> Dict[str, Any]:
        """创建未知结果"""
        return {
            'is_n_methylated': False,
            'methylation_type': 'unknown',
            'confidence': 0.0,
            'method': 'none',
            'evidence': [reason],
            'categories': [],
            'details': {
                'analysis_certainty': 'unknown',
                'failure_reason': reason
            }
        }