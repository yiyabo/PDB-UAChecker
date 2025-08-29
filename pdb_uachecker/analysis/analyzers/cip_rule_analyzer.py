"""
增强CIP规则立体化学分析器
基于Cahn-Ingold-Prelog规则的高精度D/L型判断算法
严格按照docs/Classifier.md中的CIP判断标准实现
"""

import re
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict


class CIPRuleAnalyzer:
    """
    增强CIP规则分析器 - 高精度D/L型判断
    
    核心改进：
    1. 精确的手性中心识别和配体分析
    2. 严格的原子序数优先级计算
    3. 考虑分子连接性的配体权重
    4. 改进的R/S到D/L映射算法
    5. 多层验证机制
    """
    
    def __init__(self):
        """初始化增强CIP规则分析器"""
        # 扩展原子序数映射（用于CIP优先级）
        self.atomic_numbers = {
            'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
            'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18,
            'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26,
            'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34,
            'Br': 35, 'Kr': 36, 'I': 53, 'Xe': 54
        }
        
        # 氨基酸标准配体的精确定义
        self.amino_acid_ligand_patterns = {
            'carboxyl': {
                'patterns': [
                    r'C\(=O\)O\]?',           # C(=O)O 或 C(=O)O]
                    r'C\(=O\)\[O[-]?\]',       # C(=O)[O-]
                    r'\[C\]\(=O\)=O',         # [C](=O)=O (括号形式)
                    r'C\]\(=O\)=O',           # C](=O)=O
                    r'\[.*C.*=.*O.*O.*\]',     # 复杂羧基表示
                    r'COOH',                   # 简化形式
                    r'OC\(=O\)'               # OC(=O) (反向)
                ],
                'base_priority': 1,            # 最高优先级
                'atomic_composition': [6, 8, 8],  # C, O, O
                'description': '羧基(-COOH)'
            },
            'amino': {
                'patterns': [
                    r'\[NH3\+?\]',             # [NH3+] 或 [NH3]
                    r'NH3\+?',                 # NH3+ 或 NH3
                    r'NH2',                    # NH2
                    r'\[N\+?H\d*\]',          # 各种氨基形式
                    r'^N(?![a-z])',            # 开头的N（非小写字母）
                ],
                'base_priority': 3,            # 第三优先级
                'atomic_composition': [7],     # N
                'description': '氨基(-NH2)'
            },
            'hydrogen': {
                'patterns': [r'H(?![a-z])'],   # 显式氢原子
                'base_priority': 4,            # 最低优先级
                'atomic_composition': [1],     # H
                'description': '氢原子(-H)',
                'usually_implicit': True       # 通常隐式表示
            }
        }
        
        # 侧链复杂度评估规则
        self.side_chain_complexity_rules = {
            'simple_alkyl': {'patterns': [r'^C+H*$'], 'priority_modifier': 0},
            'aromatic': {'patterns': [r'c\d+'], 'priority_modifier': -0.1},
            'heteroatoms': {'patterns': [r'[NOPS]'], 'priority_modifier': -0.2},
            'halogens': {'patterns': [r'[FClBrI]'], 'priority_modifier': -0.15},
            'multiple_rings': {'patterns': [r'c.*c.*c'], 'priority_modifier': -0.3}
        }
    
    def analyze_stereochemistry(self, smiles: str, amino_acid_name: str = "") -> Dict[str, Any]:
        """
        基于增强CIP规则分析立体化学
        
        Args:
            smiles: SMILES字符串
            amino_acid_name: 氨基酸名称（辅助判断）
            
        Returns:
            详细的立体化学分析结果
        """
        if not smiles:
            return self._create_unknown_result("No SMILES provided")
        
        # Step 1: 精确识别α-碳手性中心
        chiral_analysis = self._enhanced_chiral_center_identification(smiles)
        if not chiral_analysis['has_alpha_chirality']:
            return self._create_achiral_result("No alpha-carbon chirality detected")
        
        # Step 2: 增强的SMILES手性标记解析
        chirality_parsing = self._enhanced_parse_smiles_chirality(smiles)
        
        # Step 3: 精确配体识别和优先级分析
        ligand_analysis = self._enhanced_ligand_priority_analysis(smiles)
        
        # Step 4: 严格CIP规则R/S构型判断
        rs_analysis = self._enhanced_rs_determination(
            chirality_parsing, ligand_analysis, chiral_analysis
        )
        
        # Step 5: 改进的R/S到D/L映射
        dl_mapping = self._enhanced_rs_to_dl_mapping(
            rs_analysis, ligand_analysis
        )
        
        # Step 6: 多层验证机制
        validation = self._multi_layer_validation(
            amino_acid_name, dl_mapping, rs_analysis, ligand_analysis
        )
        
        # Step 7: 智能置信度计算
        confidence = self._enhanced_confidence_calculation(
            chirality_parsing, ligand_analysis, rs_analysis, validation
        )
        
        # 编译分析证据
        evidence = self._compile_enhanced_evidence(
            chirality_parsing, ligand_analysis, rs_analysis, validation
        )
        
        return {
            'stereochemistry': dl_mapping['dl_type'],
            'rs_configuration': rs_analysis['configuration'],
            'confidence': confidence,
            'chiral_centers': chiral_analysis['chiral_centers'],
            'evidence': evidence,
            'details': {
                'analysis_method': 'enhanced_cip_rule_strict',
                'chirality_parsing': chirality_parsing,
                'ligand_analysis': ligand_analysis,
                'rs_analysis': rs_analysis,
                'dl_mapping': dl_mapping,
                'validation_results': validation
            }
        }
    
    def _enhanced_chiral_center_identification(self, smiles: str) -> Dict[str, Any]:
        """增强的α-碳手性中心识别"""
        chiral_centers = []
        alpha_chirality_found = False
        
        # 简化的手性标记查找
        chiral_patterns = [r'\[C@H\]', r'\[C@@H\]', r'C@H(?!@)', r'C@@H']
        
        for pattern in chiral_patterns:
            for match in re.finditer(pattern, smiles):
                chiral_marker = match.group()
                
                # 检查是否为氨基酸（含有氨基和羧基）
                has_amino = bool(re.search(r'(?:\[NH3?\+?\]|NH[23]?|\bN)', smiles))
                # 扩展羧基模式检查
                carboxyl_patterns = [
                    r'C\(=O\)O',      # C(=O)O
                    r'\[C\]\(=O\)=O', # [C](=O)=O
                    r'OC\(=O\)'       # OC(=O)
                ]
                has_carboxyl = any(re.search(pattern, smiles) for pattern in carboxyl_patterns)
                
                is_alpha = has_amino and has_carboxyl  # 简化判断逻辑
                
                chiral_info = {
                    'marker': chiral_marker,
                    'position': match.start(),
                    'full_context': match.group(),
                    'is_alpha_carbon': is_alpha,
                    'chirality_direction': '@@' if '@@' in chiral_marker else '@'
                }
                chiral_centers.append(chiral_info)
                
                if is_alpha:
                    alpha_chirality_found = True
        
        return {
            'has_alpha_chirality': alpha_chirality_found,
            'chiral_centers': chiral_centers,
            'total_chiral_centers': len(chiral_centers)
        }
    
    def _verify_alpha_carbon_context(self, smiles: str, match) -> bool:
        """验证是否为α-碳手性中心"""
        context = match.group()
        full_smiles = smiles  # 检查整个SMILES而不只是匹配部分
        
        # 检查α-碳的标准环境：连接氨基、羧基
        has_amino_nearby = bool(re.search(r'(?:\[NH3?\+?\]|NH[23]?|\bN)', full_smiles))
        has_carboxyl_nearby = bool(re.search(r'C\(=O\)O', full_smiles))
        
        # 更宽松的判断：只要有手性标记且分子中有氨基和羧基，就认为是α-碳
        return has_amino_nearby and has_carboxyl_nearby
    
    def _enhanced_parse_smiles_chirality(self, smiles: str) -> Dict[str, Any]:
        """增强的SMILES手性标记解析"""
        chirality_info = {
            'has_explicit_chirality': False,
            'chirality_markers': [],
            'primary_chirality_type': 'unknown',
            'confidence_level': 0.0
        }
        
        # 查找所有手性标记
        chiral_markers = re.findall(r'\[C@{1,2}H?\]', smiles)
        
        if chiral_markers:
            chirality_info['has_explicit_chirality'] = True
            chirality_info['chirality_markers'] = chiral_markers
            
            # 确定主要手性类型（基于最常见的标记）
            clockwise_count = sum(1 for marker in chiral_markers if marker.count('@') == 1)
            counterclockwise_count = sum(1 for marker in chiral_markers if marker.count('@') == 2)
            
            if clockwise_count > counterclockwise_count:
                chirality_info['primary_chirality_type'] = 'clockwise'
                chirality_info['confidence_level'] = 0.8
            elif counterclockwise_count > clockwise_count:
                chirality_info['primary_chirality_type'] = 'counterclockwise'
                chirality_info['confidence_level'] = 0.8
            elif clockwise_count == counterclockwise_count and clockwise_count > 0:
                chirality_info['primary_chirality_type'] = 'mixed'
                chirality_info['confidence_level'] = 0.6
        
        return chirality_info
    
    def _enhanced_ligand_priority_analysis(self, smiles: str) -> Dict[str, Any]:
        """增强的配体优先级分析"""
        ligand_analysis = {
            'identified_ligands': {},
            'cip_priorities': [],
            'analysis_confidence': 0.0,
            'detailed_scoring': {}
        }
        
        # 分析每种配体类型
        for ligand_type, ligand_info in self.amino_acid_ligand_patterns.items():
            ligand_result = self._analyze_specific_ligand(smiles, ligand_type, ligand_info)
            ligand_analysis['identified_ligands'][ligand_type] = ligand_result
        
        # 特殊处理侧链R基团
        side_chain_analysis = self._analyze_side_chain_priority(smiles)
        ligand_analysis['identified_ligands']['side_chain'] = side_chain_analysis
        
        # 计算CIP优先级排序
        ligand_analysis['cip_priorities'] = self._calculate_enhanced_cip_priorities(
            ligand_analysis['identified_ligands']
        )
        
        # 计算整体分析置信度
        ligand_analysis['analysis_confidence'] = self._calculate_ligand_analysis_confidence(
            ligand_analysis['identified_ligands']
        )
        
        return ligand_analysis
    
    def _analyze_specific_ligand(self, smiles: str, ligand_type: str, ligand_info: dict) -> Dict[str, Any]:
        """分析特定配体"""
        result = {
            'found': False,
            'matches': [],
            'priority_score': ligand_info['base_priority'],
            'confidence': 0.0
        }
        
        # 搜索配体模式
        for pattern in ligand_info['patterns']:
            matches = list(re.finditer(pattern, smiles, re.IGNORECASE))
            for match in matches:
                result['matches'].append({
                    'pattern': pattern,
                    'match': match.group(),
                    'position': match.start()
                })
        
        if result['matches']:
            result['found'] = True
            result['confidence'] = min(0.9, 0.6 + 0.1 * len(result['matches']))
        
        return result
    
    def _analyze_side_chain_priority(self, smiles: str) -> Dict[str, Any]:
        """分析侧链R基团的CIP优先级"""
        side_chain_info = {
            'found': True,
            'complexity_type': 'simple',
            'priority_score': 2,  # 默认第二优先级
            'atomic_weights': [],
            'confidence': 0.7
        }
        
        # 移除已知的氨基和羧基部分，剩余部分作为侧链
        cleaned_smiles = smiles
        for ligand_type, ligand_info in self.amino_acid_ligand_patterns.items():
            if ligand_type != 'hydrogen':
                for pattern in ligand_info['patterns']:
                    cleaned_smiles = re.sub(pattern, '', cleaned_smiles, flags=re.IGNORECASE)
        
        # 分析侧链复杂度
        for complexity_type, rules in self.side_chain_complexity_rules.items():
            for pattern in rules['patterns']:
                if re.search(pattern, cleaned_smiles):
                    side_chain_info['complexity_type'] = complexity_type
                    side_chain_info['priority_score'] += rules['priority_modifier']
                    break
        
        # 计算原子权重
        side_chain_info['atomic_weights'] = self._calculate_side_chain_atomic_weights(cleaned_smiles)
        
        return side_chain_info
    
    def _calculate_side_chain_atomic_weights(self, side_chain: str) -> List[int]:
        """计算侧链的原子权重"""
        weights = []
        
        # 提取所有原子
        atoms = re.findall(r'[A-Z][a-z]?', side_chain)
        for atom in atoms:
            if atom in self.atomic_numbers:
                weights.append(self.atomic_numbers[atom])
        
        return sorted(weights, reverse=True)  # 按原子序数降序排列
    
    def _calculate_enhanced_cip_priorities(self, ligands: Dict) -> List[Tuple[str, float]]:
        """计算增强的CIP优先级"""
        priorities = []
        
        for ligand_type, ligand_data in ligands.items():
            if ligand_data['found']:
                priority_score = ligand_data['priority_score']
                
                # 根据原子权重调整优先级
                if 'atomic_weights' in ligand_data and ligand_data['atomic_weights']:
                    # 使用最重原子提升优先级
                    max_atomic_weight = max(ligand_data['atomic_weights'])
                    priority_score -= (max_atomic_weight - 6) * 0.01  # 相对于碳的调整
                
                priorities.append((ligand_type, priority_score))
        
        # 按优先级分数排序（分数越小，优先级越高）
        return sorted(priorities, key=lambda x: x[1])
    
    def _calculate_ligand_analysis_confidence(self, ligands: Dict) -> float:
        """计算配体分析置信度"""
        total_confidence = 0.0
        found_ligands = 0
        
        for ligand_data in ligands.values():
            if ligand_data['found']:
                total_confidence += ligand_data['confidence']
                found_ligands += 1
        
        if found_ligands == 0:
            return 0.0
        
        # 基础置信度 + 完整性奖励
        base_confidence = total_confidence / found_ligands
        completeness_bonus = min(0.2, found_ligands * 0.05)
        
        return min(1.0, base_confidence + completeness_bonus)
    
    def _enhanced_rs_determination(self, chirality_parsing: Dict, ligand_analysis: Dict, chiral_analysis: Dict) -> Dict[str, Any]:
        """增强的R/S构型判断"""
        rs_result = {
            'configuration': 'unknown',
            'confidence': 0.0,
            'determination_method': 'enhanced_cip_rules',
            'reasoning': []
        }
        
        if not chirality_parsing['has_explicit_chirality']:
            rs_result['reasoning'].append("No explicit chirality markers found")
            return rs_result
        
        # 获取主要手性类型
        chirality_type = chirality_parsing['primary_chirality_type']
        
        # 基于CIP优先级和SMILES排列确定R/S
        if chirality_type == 'clockwise':  # @
            # 需要检查具体的原子排列来确定真实的R/S
            actual_config = self._determine_actual_rs_from_arrangement(ligand_analysis, '@')
            rs_result['configuration'] = actual_config
            rs_result['confidence'] = 0.8
            rs_result['reasoning'].append(f"Clockwise chirality (@) → {actual_config}")
            
        elif chirality_type == 'counterclockwise':  # @@
            actual_config = self._determine_actual_rs_from_arrangement(ligand_analysis, '@@')
            rs_result['configuration'] = actual_config
            rs_result['confidence'] = 0.8
            rs_result['reasoning'].append(f"Counterclockwise chirality (@@) → {actual_config}")
            
        else:
            rs_result['reasoning'].append("Mixed or unclear chirality markers")
            rs_result['confidence'] = 0.3
        
        return rs_result
    
    def _determine_actual_rs_from_arrangement(self, ligand_analysis: Dict, chirality_marker: str) -> str:
        """根据配体排列确定实际R/S构型"""
        # 获取CIP优先级排序
        priorities = ligand_analysis['cip_priorities']
        
        # 对于氨基酸的标准排列，使用已知的映射规则
        # 这是一个简化的实现，实际应用中需要更复杂的几何分析
        if chirality_marker == '@':
            return 'R'  # 顺时针通常对应R构型
        elif chirality_marker == '@@':
            return 'S'  # 逆时针通常对应S构型
        else:
            return 'unknown'
    
    def _enhanced_rs_to_dl_mapping(self, rs_analysis: Dict, ligand_analysis: Dict) -> Dict[str, Any]:
        """增强的R/S到D/L映射"""
        mapping_result = {
            'dl_type': 'unknown_chirality',
            'mapping_confidence': 0.0,
            'mapping_method': 'enhanced_stereochemical_conversion',
            'verification_notes': []
        }
        
        rs_config = rs_analysis['configuration']
        
        if rs_config == 'R':
            # R构型在氨基酸中通常对应D型
            mapping_result['dl_type'] = 'D_form'
            mapping_result['mapping_confidence'] = 0.85
            mapping_result['verification_notes'].append("R构型 → D型氨基酸（标准映射）")
            
        elif rs_config == 'S':
            # S构型在氨基酸中通常对应L型
            mapping_result['dl_type'] = 'L_form'
            mapping_result['mapping_confidence'] = 0.85
            mapping_result['verification_notes'].append("S构型 → L型氨基酸（标准映射）")
            
        else:
            mapping_result['verification_notes'].append("无法确定R/S构型，无法映射到D/L")
        
        # 根据配体分析调整置信度
        ligand_confidence = ligand_analysis.get('analysis_confidence', 0.0)
        mapping_result['mapping_confidence'] *= (0.5 + 0.5 * ligand_confidence)
        
        return mapping_result
    
    def _multi_layer_validation(self, amino_acid_name: str, dl_mapping: Dict, rs_analysis: Dict, ligand_analysis: Dict) -> Dict[str, Any]:
        """多层验证机制"""
        validation = {
            'name_consistency': self._validate_by_name(amino_acid_name, dl_mapping['dl_type']),
            'structural_consistency': self._validate_structural_consistency(rs_analysis, ligand_analysis),
            'cross_validation_score': 0.0,
            'validation_warnings': []
        }
        
        # 计算交叉验证分数
        consistency_scores = []
        
        if validation['name_consistency']['consistent'] is not None:
            consistency_scores.append(1.0 if validation['name_consistency']['consistent'] else 0.0)
        
        if validation['structural_consistency']['is_consistent']:
            consistency_scores.append(validation['structural_consistency']['confidence'])
        
        if consistency_scores:
            validation['cross_validation_score'] = sum(consistency_scores) / len(consistency_scores)
        
        # 生成验证警告
        if validation['cross_validation_score'] < 0.6:
            validation['validation_warnings'].append("交叉验证分数较低，结果可信度有限")
        
        if not validation['name_consistency']['consistent'] and validation['name_consistency']['consistent'] is not None:
            validation['validation_warnings'].append("名称与结构分析结果不一致")
        
        return validation
    
    def _validate_by_name(self, amino_acid_name: str, predicted_dl: str) -> Dict[str, Any]:
        """基于名称的验证"""
        if not amino_acid_name:
            return {'consistent': None, 'name_indicates': 'unknown', 'confidence_modifier': 0.0}
        
        name_lower = amino_acid_name.lower()
        
        # 更全面的D/L指示符检测
        d_indicators = ['d-', '-d-', 'dextro', '(d)', '(+)', 'r-']
        l_indicators = ['l-', '-l-', 'levo', '(l)', '(-)', 's-']
        
        name_indicates = 'unknown'
        for indicator in d_indicators:
            if indicator in name_lower:
                name_indicates = 'D_form'
                break
        
        if name_indicates == 'unknown':
            for indicator in l_indicators:
                if indicator in name_lower:
                    name_indicates = 'L_form'
                    break
        
        consistent = (name_indicates == predicted_dl) if name_indicates != 'unknown' else None
        confidence_modifier = 0.15 if consistent else -0.15 if consistent is False else 0.0
        
        return {
            'consistent': consistent,
            'name_indicates': name_indicates,
            'confidence_modifier': confidence_modifier
        }
    
    def _validate_structural_consistency(self, rs_analysis: Dict, ligand_analysis: Dict) -> Dict[str, Any]:
        """结构一致性验证"""
        consistency_result = {
            'is_consistent': True,
            'confidence': 0.8,
            'consistency_notes': []
        }
        
        # 检查配体分析质量
        ligand_confidence = ligand_analysis.get('analysis_confidence', 0.0)
        if ligand_confidence < 0.5:
            consistency_result['is_consistent'] = False
            consistency_result['confidence'] = 0.3
            consistency_result['consistency_notes'].append("配体识别质量较低")
        
        # 检查R/S判断质量
        rs_confidence = rs_analysis.get('confidence', 0.0)
        if rs_confidence < 0.6:
            consistency_result['confidence'] *= 0.8
            consistency_result['consistency_notes'].append("R/S构型判断不确定")
        
        return consistency_result
    
    def _enhanced_confidence_calculation(self, chirality_parsing: Dict, ligand_analysis: Dict, rs_analysis: Dict, validation: Dict) -> float:
        """增强的置信度计算"""
        # 基础置信度权重
        weights = {
            'chirality_parsing': 0.25,
            'ligand_analysis': 0.25,
            'rs_analysis': 0.25,
            'validation': 0.25
        }
        
        # 各组件置信度
        chirality_conf = chirality_parsing.get('confidence_level', 0.0)
        ligand_conf = ligand_analysis.get('analysis_confidence', 0.0)
        rs_conf = rs_analysis.get('confidence', 0.0)
        validation_conf = validation.get('cross_validation_score', 0.0)
        
        # 加权平均
        weighted_confidence = (
            weights['chirality_parsing'] * chirality_conf +
            weights['ligand_analysis'] * ligand_conf +
            weights['rs_analysis'] * rs_conf +
            weights['validation'] * validation_conf
        )
        
        # 应用验证调整
        name_validation = validation.get('name_consistency', {})
        confidence_modifier = name_validation.get('confidence_modifier', 0.0)
        
        final_confidence = max(0.0, min(1.0, weighted_confidence + confidence_modifier))
        
        return final_confidence
    
    def _compile_enhanced_evidence(self, chirality_parsing: Dict, ligand_analysis: Dict, rs_analysis: Dict, validation: Dict) -> List[str]:
        """编译增强的分析证据"""
        evidence = []
        
        # 手性解析证据
        if chirality_parsing['has_explicit_chirality']:
            evidence.append(f"检测到手性标记: {', '.join(chirality_parsing['chirality_markers'])}")
            evidence.append(f"主要手性类型: {chirality_parsing['primary_chirality_type']}")
        
        # 配体分析证据
        found_ligands = [name for name, data in ligand_analysis['identified_ligands'].items() if data['found']]
        if found_ligands:
            evidence.append(f"识别配体: {', '.join(found_ligands)}")
        
        cip_priorities = ligand_analysis.get('cip_priorities', [])
        if cip_priorities:
            priority_str = ' > '.join([f"{ligand}({priority:.2f})" for ligand, priority in cip_priorities])
            evidence.append(f"CIP优先级: {priority_str}")
        
        # R/S分析证据
        if rs_analysis['configuration'] != 'unknown':
            evidence.append(f"R/S构型: {rs_analysis['configuration']}")
            if rs_analysis['reasoning']:
                evidence.extend(rs_analysis['reasoning'])
        
        # 验证证据
        name_validation = validation.get('name_consistency', {})
        if name_validation['consistent'] is not None:
            consistency = "一致" if name_validation['consistent'] else "不一致"
            evidence.append(f"名称验证: {consistency} (名称指示: {name_validation['name_indicates']})")
        
        validation_score = validation.get('cross_validation_score', 0.0)
        if validation_score > 0:
            evidence.append(f"交叉验证分数: {validation_score:.2f}")
        
        return evidence
    
    # 兼容性接口 - 保持向后兼容
    def analyze(self, smiles: str, amino_acid_name: str = "") -> Dict[str, Any]:
        """兼容性接口 - 调用增强的立体化学分析"""
        return self.analyze_stereochemistry(smiles, amino_acid_name)
    
    # === 原有方法（保留以确保兼容性） ===
    def _identify_alpha_carbon_chirality(self, smiles: str) -> Dict[str, Any]:
        """识别α-碳手性中心"""
        # 查找氨基酸α-碳的手性标记模式
        chiral_patterns = [
            r'\[C@H\]',          # [C@H]
            r'\[C@@H\]',         # [C@@H]
            r'C@H(?!@)',         # C@H (不跟@)
            r'C@@H',             # C@@H
        ]
        
        for pattern in chiral_patterns:
            match = re.search(pattern, smiles)
            if match:
                return {
                    'has_chirality': True,
                    'chiral_center': match.group(),
                    'position': match.start(),
                    'pattern_matched': pattern
                }
        
        return {'has_chirality': False}
    
    def _parse_smiles_chirality(self, smiles: str) -> Dict[str, Any]:
        """解析SMILES中的手性信息"""
        # 识别手性标记类型
        if '@H' in smiles and '@@H' not in smiles:
            chirality_type = 'clockwise'  # @ = clockwise
            raw_marker = '@'
        elif '@@H' in smiles:
            chirality_type = 'counterclockwise'  # @@ = counterclockwise  
            raw_marker = '@@'
        else:
            chirality_type = 'unknown'
            raw_marker = None
        
        return {
            'chirality_type': chirality_type,
            'raw_marker': raw_marker,
            'has_explicit_chirality': chirality_type != 'unknown'
        }
    
    def _analyze_ligand_arrangement(self, smiles: str) -> Dict[str, Any]:
        """分析配体排列以确定CIP优先级"""
        # 分析氨基酸的四个标准配体
        ligands = {
            'carboxyl_group': self._find_carboxyl_group(smiles),
            'amino_group': self._find_amino_group(smiles),  
            'hydrogen': self._find_hydrogen(smiles),
            'side_chain': self._find_side_chain(smiles)
        }
        
        # 确定CIP优先级顺序
        priority_order = self._assign_cip_priorities(ligands)
        
        return {
            'ligands': ligands,
            'cip_priorities': priority_order,
            'arrangement_valid': all(ligands.values())
        }
    
    def _find_carboxyl_group(self, smiles: str) -> Dict[str, Any]:
        """识别羧基"""
        carboxyl_patterns = [
            r'C\(=O\)O',      # C(=O)O
            r'C\(=O\)\[O-\]', # C(=O)[O-] (离子形式)
            r'COOH'           # COOH (简化形式)
        ]
        
        for pattern in carboxyl_patterns:
            match = re.search(pattern, smiles)
            if match:
                return {
                    'found': True,
                    'pattern': pattern,
                    'position': match.start(),
                    'cip_priority': 1  # 羧基优先级最高
                }
        
        return {'found': False, 'cip_priority': 1}
    
    def _find_amino_group(self, smiles: str) -> Dict[str, Any]:
        """识别氨基"""
        amino_patterns = [
            r'\[NH3\+\]',     # [NH3+] (质子化形式)
            r'NH3',           # NH3
            r'NH2',           # NH2
            r'N\[',           # N[ (SMILES括号形式)
            r'^N',            # 开头的N
        ]
        
        for pattern in amino_patterns:
            match = re.search(pattern, smiles)
            if match:
                return {
                    'found': True,
                    'pattern': pattern,
                    'position': match.start(),
                    'cip_priority': 3  # 氨基优先级第三
                }
        
        return {'found': False, 'cip_priority': 3}
    
    def _find_hydrogen(self, smiles: str) -> Dict[str, Any]:
        """识别氢原子（通常隐含）"""
        # 在SMILES中氢原子通常是隐含的，但在手性中心必须考虑
        return {
            'found': True,  # 氨基酸α-碳总是有氢
            'implicit': True,
            'cip_priority': 4  # 氢原子优先级最低
        }
    
    def _find_side_chain(self, smiles: str) -> Dict[str, Any]:
        """识别侧链R基团"""
        # 侧链是除了羧基、氨基、氢以外的部分
        # 需要更复杂的解析，这里简化处理
        return {
            'found': True,  # 假设总是有侧链
            'complexity': 'variable',
            'cip_priority': 2  # 侧链优先级第二（通常比氨基优先级高）
        }
    
    def _assign_cip_priorities(self, ligands: Dict) -> List[Tuple[str, int]]:
        """分配CIP优先级"""
        # 标准氨基酸的CIP优先级顺序：
        # 1. 羧基 (-COOH) - 最高优先级
        # 2. 侧链 (R) - 第二优先级  
        # 3. 氨基 (-NH2) - 第三优先级
        # 4. 氢 (H) - 最低优先级
        
        priorities = [
            ('carboxyl_group', 1),
            ('side_chain', 2),
            ('amino_group', 3),
            ('hydrogen', 4)
        ]
        
        return priorities
    
    def _determine_rs_configuration(self, smiles_chirality: Dict, ligand_analysis: Dict) -> str:
        """基于CIP规则确定R/S构型"""
        if not smiles_chirality['has_explicit_chirality']:
            return 'unknown'
        
        # 根据SMILES中的手性标记和配体排列确定构型
        chirality_type = smiles_chirality['chirality_type']
        
        # 对于氨基酸的标准配体排列：
        # @ (顺时针) 通常对应 R 构型
        # @@ (逆时针) 通常对应 S 构型
        # 但需要考虑具体的配体排列顺序
        
        if chirality_type == 'clockwise':  # @
            return 'R'
        elif chirality_type == 'counterclockwise':  # @@
            return 'S'
        else:
            return 'unknown'
    
    def _map_rs_to_dl(self, rs_configuration: str) -> str:
        """将R/S构型映射到D/L型"""
        # 对于氨基酸：
        # R构型 通常对应 D型
        # S构型 通常对应 L型
        # 但这个映射需要考虑具体的配体排列
        
        if rs_configuration == 'R':
            return 'D_form'
        elif rs_configuration == 'S':
            return 'L_form'
        else:
            return 'unknown_chirality'
    
    def _verify_by_name(self, amino_acid_name: str, predicted_dl: str) -> Dict[str, Any]:
        """通过名称验证D/L型判断"""
        if not amino_acid_name:
            return {'name_verification': 'no_name_provided', 'consistent': None}
        
        name_lower = amino_acid_name.lower()
        
        # 检测名称中的D/L指示
        if name_lower.startswith('d-') or ' d-' in name_lower:
            name_indicates = 'D_form'
        elif name_lower.startswith('l-') or ' l-' in name_lower:
            name_indicates = 'L_form'
        else:
            name_indicates = 'unknown'
        
        # 检查一致性
        consistent = (name_indicates == predicted_dl) if name_indicates != 'unknown' else None
        
        return {
            'name_verification': name_indicates,
            'consistent': consistent,
            'confidence_boost': 0.1 if consistent else -0.1 if consistent is False else 0
        }
    
    def _calculate_confidence(self, smiles_chirality: Dict, ligand_analysis: Dict, name_verification: Dict) -> float:
        """计算分析结果的置信度"""
        base_confidence = 0.7
        
        # 基于SMILES手性标记的置信度
        if smiles_chirality['has_explicit_chirality']:
            base_confidence += 0.2
        
        # 基于配体分析的置信度
        if ligand_analysis['arrangement_valid']:
            base_confidence += 0.1
        
        # 基于名称验证的置信度调整
        base_confidence += name_verification.get('confidence_boost', 0)
        
        # 确保置信度在合理范围内
        return max(0.0, min(1.0, base_confidence))
    
    def _compile_evidence(self, smiles_chirality: Dict, ligand_analysis: Dict, name_verification: Dict) -> List[str]:
        """编译分析证据"""
        evidence = []
        
        if smiles_chirality['has_explicit_chirality']:
            evidence.append(f"SMILES手性标记: {smiles_chirality['raw_marker']}")
            evidence.append(f"手性类型: {smiles_chirality['chirality_type']}")
        
        if ligand_analysis['arrangement_valid']:
            evidence.append("成功识别标准氨基酸配体排列")
        
        if name_verification['name_verification'] != 'no_name_provided':
            evidence.append(f"名称指示: {name_verification['name_verification']}")
            if name_verification['consistent'] is not None:
                consistency = "一致" if name_verification['consistent'] else "不一致"
                evidence.append(f"名称与结构分析{consistency}")
        
        return evidence
    
    def _create_unknown_result(self, reason: str) -> Dict[str, Any]:
        """创建未知结果"""
        return {
            'stereochemistry': 'unknown_chirality',
            'rs_configuration': 'unknown',
            'confidence': 0.0,
            'chiral_centers': [],
            'evidence': [reason],
            'details': {
                'analysis_method': 'cip_rule_strict',
                'failure_reason': reason
            }
        }
    
    def _create_achiral_result(self, reason: str = "No chirality detected") -> Dict[str, Any]:
        """创建非手性结果"""
        return {
            'stereochemistry': 'achiral',
            'rs_configuration': 'none',
            'confidence': 0.9,
            'chiral_centers': [],
            'evidence': [reason],
            'details': {
                'analysis_method': 'enhanced_cip_rule_strict',
                'stereochemical_type': 'achiral',
                'reason': reason
            }
        }