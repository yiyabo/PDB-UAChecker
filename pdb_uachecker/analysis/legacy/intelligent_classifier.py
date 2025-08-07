"""
智能高覆盖率氨基酸分类器
目标：95%+覆盖率，同时保持高精确度
基于多层智能规则和结构特征分析
"""

import re
from typing import Dict, Any, List, Tuple, Optional
from collections import defaultdict
import json

from ...core.models import AminoAcidInfo, ResidueInfo


class IntelligentAminoAcidClassifier:
    """智能氨基酸分类器 - 高覆盖率 + 高精确度"""
    
    def __init__(self):
        """初始化智能分类器"""
        
        # 标准氨基酸 - 100%确定
        self.standard_amino_acids = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        # 分层分类规则 - 从100%确定到智能推断
        self.classification_tiers = self._initialize_classification_tiers()
        
        # 结构特征分析器
        self.structure_analyzers = self._initialize_structure_analyzers()
        
        # 分类统计
        self.stats = {
            'tier1_definitive': 0,
            'tier2_high_confidence': 0,
            'tier3_intelligent': 0,
            'tier4_needs_review': 0,
            'total_processed': 0
        }
    
    def _initialize_classification_tiers(self) -> Dict[str, Dict[str, Any]]:
        """初始化分层分类规则"""
        return {
            # Tier 1: 100%确定 - 无需验证
            'tier1_definitive': {
                'standard_amino_acid': {
                    'validator': self._is_standard_amino_acid,
                    'confidence': 1.0,
                    'description': '20种标准氨基酸'
                },
                'D_amino_acid_by_naming': {
                    'validator': self._is_d_amino_by_naming,
                    'confidence': 1.0,
                    'description': 'PDB D前缀命名'
                },
                'D_amino_acid_by_chirality': {
                    'validator': self._is_d_amino_by_chirality,
                    'confidence': 1.0,
                    'description': 'D型手性中心确认'
                }
            },
            
            # Tier 2: 高置信度 - 基于明确的结构特征
            'tier2_high_confidence': {
                'aromatic_amino_acid': {
                    'validator': self._is_aromatic_amino_acid,
                    'confidence': 0.95,
                    'description': '含芳香环结构'
                },
                'cyclic_amino_acid': {
                    'validator': self._is_cyclic_amino_acid,
                    'confidence': 0.95,
                    'description': '含环状结构'
                },
                'beta_amino_acid': {
                    'validator': self._is_beta_amino_acid,
                    'confidence': 0.90,
                    'description': 'Beta氨基酸主链'
                },
                'gamma_amino_acid': {
                    'validator': self._is_gamma_amino_acid,
                    'confidence': 0.90,
                    'description': 'Gamma氨基酸主链'
                }
            },
            
            # Tier 3: 智能推断 - 基于复杂模式识别
            'tier3_intelligent': {
                'n_methyl_amino_acid': {
                    'validator': self._is_n_methyl_amino_acid,
                    'confidence': 0.85,
                    'description': 'N-甲基修饰'
                },
                'halogenated_amino_acid': {
                    'validator': self._is_halogenated_amino_acid,
                    'confidence': 0.85,
                    'description': '卤素取代'
                },
                'hydroxylated_amino_acid': {
                    'validator': self._is_hydroxylated_amino_acid,
                    'confidence': 0.80,
                    'description': '羟基修饰'
                },
                'extended_chain_amino_acid': {
                    'validator': self._is_extended_chain_amino_acid,
                    'confidence': 0.80,
                    'description': '延长侧链'
                }
            }
        }
    
    def _initialize_structure_analyzers(self) -> Dict[str, Any]:
        """初始化结构分析器"""
        return {
            'backbone_analyzer': {
                'alpha_pattern': r'N.*?C.*?C\(=O\)O',
                'beta_pattern': r'N.*?C.*?C.*?C\(=O\)O',
                'gamma_pattern': r'N.*?C.*?C.*?C.*?C\(=O\)O'
            },
            'functional_group_analyzer': {
                'hydroxyl': r'OH?',
                'amino': r'NH?[23]?',
                'carboxyl': r'C\(=O\)OH?',
                'thiol': r'SH?',
                'aromatic': r'c\d+',
                'halogen': r'[FClBrI]',
                'methyl': r'CH?3'
            },
            'ring_analyzer': {
                'any_ring': r'\d+',
                'aromatic_ring': r'c\d+',
                'aliphatic_ring': r'C\d+.*?\d+'
            }
        }
    
    def classify_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """
        智能分类氨基酸
        
        Args:
            amino_acid: 氨基酸信息
            
        Returns:
            智能分类结果
        """
        self.stats['total_processed'] += 1
        
        if not amino_acid.smiles:
            return self._create_unclassified_result(amino_acid, 'no_smiles')
        
        # Tier 1: 确定性分类
        tier1_result = self._classify_tier1(amino_acid)
        if tier1_result:
            self.stats['tier1_definitive'] += 1
            return tier1_result
        
        # Tier 2: 高置信度分类
        tier2_result = self._classify_tier2(amino_acid)
        if tier2_result:
            self.stats['tier2_high_confidence'] += 1
            return tier2_result
        
        # Tier 3: 智能推断分类
        tier3_result = self._classify_tier3(amino_acid)
        if tier3_result:
            self.stats['tier3_intelligent'] += 1
            return tier3_result
        
        # Tier 4: 需要复杂分析
        self.stats['tier4_needs_review'] += 1
        return self._create_unclassified_result(amino_acid, 'complex_structure')
    
    def _classify_tier1(self, amino_acid: AminoAcidInfo) -> Optional[Dict[str, Any]]:
        """Tier 1: 确定性分类"""
        tier1_rules = self.classification_tiers['tier1_definitive']
        
        for category_key, rule_config in tier1_rules.items():
            validator = rule_config['validator']
            if validator(amino_acid):
                category = category_key.replace('_amino_acid', '').replace('_', '-').title()
                if 'by_naming' in category_key or 'by_chirality' in category_key:
                    category = 'D-amino acid'
                
                return {
                    'amino_acid_id': amino_acid.id,
                    'amino_acid_name': amino_acid.name,
                    'standard_category': category,
                    'confidence': rule_config['confidence'],
                    'classification_tier': 'tier1_definitive',
                    'detection_method': category_key,
                    'evidence': rule_config['description'],
                    'smiles': amino_acid.smiles,
                    'details': {
                        'tier': 1,
                        'certainty': 'definitive',
                        'manual_review_required': False
                    }
                }
        
        return None
    
    def _classify_tier2(self, amino_acid: AminoAcidInfo) -> Optional[Dict[str, Any]]:
        """Tier 2: 高置信度分类"""
        tier2_rules = self.classification_tiers['tier2_high_confidence']
        
        # 按置信度排序处理
        sorted_rules = sorted(tier2_rules.items(), 
                            key=lambda x: x[1]['confidence'], reverse=True)
        
        for category_key, rule_config in sorted_rules:
            validator = rule_config['validator']
            result = validator(amino_acid)
            
            if isinstance(result, dict) and result.get('is_match', False):
                category = category_key.replace('_amino_acid', '').replace('_', '-').title()
                
                return {
                    'amino_acid_id': amino_acid.id,
                    'amino_acid_name': amino_acid.name,
                    'standard_category': category,
                    'confidence': rule_config['confidence'],
                    'classification_tier': 'tier2_high_confidence',
                    'detection_method': category_key,
                    'evidence': rule_config['description'],
                    'smiles': amino_acid.smiles,
                    'details': {
                        'tier': 2,
                        'certainty': 'high_confidence',
                        'manual_review_required': False,
                        'analysis_details': result.get('details', {})
                    }
                }
            elif isinstance(result, bool) and result:  # 只匹配明确的布尔True值
                category = category_key.replace('_amino_acid', '').replace('_', '-').title()
                
                return {
                    'amino_acid_id': amino_acid.id,
                    'amino_acid_name': amino_acid.name,
                    'standard_category': category,
                    'confidence': rule_config['confidence'],
                    'classification_tier': 'tier2_high_confidence',
                    'detection_method': category_key,
                    'evidence': rule_config['description'],
                    'smiles': amino_acid.smiles,
                    'details': {
                        'tier': 2,
                        'certainty': 'high_confidence',
                        'manual_review_required': False
                    }
                }
        
        return None
    
    def _classify_tier3(self, amino_acid: AminoAcidInfo) -> Optional[Dict[str, Any]]:
        """Tier 3: 智能推断分类"""
        tier3_rules = self.classification_tiers['tier3_intelligent']
        
        best_match = None
        best_confidence = 0
        
        for category_key, rule_config in tier3_rules.items():
            validator = rule_config['validator']
            result = validator(amino_acid)
            
            if isinstance(result, dict) and result.get('is_match', False):
                confidence = rule_config['confidence'] * result.get('strength', 1.0)
                
                if confidence > best_confidence:
                    best_confidence = confidence
                    category = category_key.replace('_amino_acid', '').replace('_', '-').title()
                    
                    best_match = {
                        'amino_acid_id': amino_acid.id,
                        'amino_acid_name': amino_acid.name,
                        'standard_category': category,
                        'confidence': confidence,
                        'classification_tier': 'tier3_intelligent',
                        'detection_method': category_key,
                        'evidence': rule_config['description'],
                        'smiles': amino_acid.smiles,
                        'details': {
                            'tier': 3,
                            'certainty': 'intelligent_inference',
                            'manual_review_required': confidence < 0.85,
                            'analysis_details': result.get('details', {}),
                            'alternative_possibilities': result.get('alternatives', [])
                        }
                    }
        
        return best_match
    
    # ========== Tier 1 验证器 ==========
    
    def _is_standard_amino_acid(self, amino_acid: AminoAcidInfo) -> bool:
        """验证是否为标准氨基酸"""
        return amino_acid.id in self.standard_amino_acids
    
    def _is_d_amino_by_naming(self, amino_acid: AminoAcidInfo) -> bool:
        """基于PDB命名验证D氨基酸"""
        return (amino_acid.id.startswith('D') and 
                len(amino_acid.id) == 3 and 
                amino_acid.id not in self.standard_amino_acids)
    
    def _is_d_amino_by_chirality(self, amino_acid: AminoAcidInfo) -> bool:
        """基于手性中心验证D氨基酸"""
        if not amino_acid.smiles:
            return False
        
        # D型氨基酸特征: C@H 且没有 C@@H
        has_d_chirality = 'C@H' in amino_acid.smiles
        has_l_chirality = 'C@@H' in amino_acid.smiles
        
        # 必须有D型标记，且不是标准氨基酸
        return (has_d_chirality and 
                not has_l_chirality and 
                amino_acid.id not in self.standard_amino_acids)
    
    # ========== Tier 2 验证器 ==========
    
    def _is_aromatic_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """验证芳香族氨基酸 - 精确检测"""
        if not amino_acid.smiles:
            return {'is_match': False}
        
        smiles = amino_acid.smiles.lower()
        
        # 精确的芳香环检测模式 - 移除过于宽泛的规则
        precise_aromatic_patterns = [
            r'c1ccccc1',        # 标准苯环
            r'c1cccc[nos]1',    # 五元芳香杂环 (呋喃、噻吩、吡咯)
            r'c1ccc[nos]c1',    # 六元芳香杂环 (吡啶等)
            r'c1c\[nh\]c2c1cccc2',  # 吲哚结构
            r'c1cc\[nh\]cc1',   # 简化吡咯结构
            r'c1ccc2\[nh\]ccc2c1',  # 喹啉结构
        ]
        
        matched_patterns = []
        for pattern in precise_aromatic_patterns:
            if re.search(pattern, smiles):
                matched_patterns.append(pattern)
        
        # 额外检查：确保确实是芳香性结构
        if matched_patterns:
            # 进一步验证：芳香环必须包含环编号或明确的芳香标记
            has_ring_numbers = any(c.isdigit() for c in smiles)
            has_aromatic_markers = 'c1' in smiles or 'c2' in smiles or '[nh]' in smiles
            
            if has_ring_numbers and has_aromatic_markers:
                return {
                    'is_match': True,
                    'strength': min(1.0, len(matched_patterns) * 0.3 + 0.7),
                    'details': {
                        'matched_patterns': matched_patterns,
                        'aromatic_systems_count': len(matched_patterns),
                        'validation': 'precise_aromatic_detection'
                    }
                }
        
        return {'is_match': False}
    
    def _is_cyclic_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """验证环状氨基酸"""
        if not amino_acid.smiles:
            return {'is_match': False}
        
        # 检查环闭合数字
        ring_numbers = re.findall(r'\d+', amino_acid.smiles)
        
        if ring_numbers:
            # 进一步验证是否为真正的环状结构
            unique_numbers = set(ring_numbers)
            ring_count = len([n for n in unique_numbers if ring_numbers.count(n) >= 2])
            
            if ring_count > 0:
                return {
                    'is_match': True,
                    'strength': min(1.0, ring_count * 0.4 + 0.6),
                    'details': {
                        'ring_systems': ring_count,
                        'ring_numbers': list(unique_numbers)
                    }
                }
        
        return {'is_match': False}
    
    def _is_beta_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """验证Beta氨基酸"""
        if not amino_acid.smiles:
            return {'is_match': False}
        
        # Beta氨基酸模式: NH2-CH2-CH(R)-COOH
        beta_patterns = [
            r'N.*?CC.*?C\(=O\)O',          # 基本模式
            r'\[NH3?\]\+?CC.*?C\(=O\)O',   # 离子化氨基
            r'NCC\[C@@?H\].*?C\(=O\)O',    # 明确的手性中心
        ]
        
        # 排除模式（避免误判）
        exclude_patterns = [
            r'NCCC.*?C\(=O\)O',  # Gamma氨基酸
            r'N.*?CCC.*?C\(=O\)O' # 更长链
        ]
        
        matched = False
        confidence_boost = 0
        
        for pattern in beta_patterns:
            if re.search(pattern, amino_acid.smiles):
                matched = True
                confidence_boost += 0.3
        
        # 检查排除条件
        for exclude_pattern in exclude_patterns:
            if re.search(exclude_pattern, amino_acid.smiles):
                matched = False
                break
        
        if matched:
            return {
                'is_match': True,
                'strength': min(1.0, 0.7 + confidence_boost),
                'details': {
                    'pattern_type': 'beta_backbone',
                    'confidence_factors': ['backbone_length', 'amino_carboxyl_distance']
                }
            }
        
        return {'is_match': False}
    
    def _is_gamma_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """验证Gamma氨基酸"""
        if not amino_acid.smiles:
            return {'is_match': False}
        
        # Gamma氨基酸模式: NH2-CH2-CH2-CH(R)-COOH
        gamma_patterns = [
            r'NCCC.*?C\(=O\)O',           # 基本模式
            r'\[NH3?\]\+?CCC.*?C\(=O\)O', # 离子化氨基
            r'NCCC\[C@@?H\].*?C\(=O\)O',  # 明确的手性中心
        ]
        
        for pattern in gamma_patterns:
            if re.search(pattern, amino_acid.smiles):
                return {
                    'is_match': True,
                    'strength': 0.9,
                    'details': {
                        'pattern_type': 'gamma_backbone',
                        'matched_pattern': pattern
                    }
                }
        
        return {'is_match': False}
    
    # ========== Tier 3 验证器 ==========
    
    def _is_n_methyl_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """验证N-甲基氨基酸"""
        if not amino_acid.smiles:
            return {'is_match': False}
        
        n_methyl_indicators = []
        strength = 0.5
        
        # 检查仲胺标记
        if '[NH]' in amino_acid.smiles:
            n_methyl_indicators.append('secondary_amine')
            strength += 0.3
        
        # 检查甲基-氮连接
        if re.search(r'CN\[C@@?H\]', amino_acid.smiles):
            n_methyl_indicators.append('methyl_nitrogen_bond')
            strength += 0.4
        
        # 检查N-C=O键（酰胺）
        if 'NC(=O)' in amino_acid.smiles:
            n_methyl_indicators.append('amide_bond')
            strength += 0.2
        
        if n_methyl_indicators:
            return {
                'is_match': True,
                'strength': min(1.0, strength),
                'details': {
                    'indicators': n_methyl_indicators,
                    'methylation_type': 'n_terminal' if 'CN[C' in amino_acid.smiles else 'unknown'
                }
            }
        
        return {'is_match': False}
    
    def _is_halogenated_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """验证卤素取代氨基酸"""
        if not amino_acid.smiles:
            return {'is_match': False}
        
        halogens = {'F': 'fluorine', 'Cl': 'chlorine', 'Br': 'bromine', 'I': 'iodine'}
        found_halogens = []
        
        for halogen, name in halogens.items():
            if halogen in amino_acid.smiles:
                found_halogens.append(name)
        
        if found_halogens:
            return {
                'is_match': True,
                'strength': 0.9,
                'details': {
                    'halogen_types': found_halogens,
                    'substitution_count': len(found_halogens)
                }
            }
        
        return {'is_match': False}
    
    def _is_hydroxylated_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """验证羟基修饰氨基酸"""
        if not amino_acid.smiles:
            return {'is_match': False}
        
        # 检查羟基，但排除羧基
        oh_count = amino_acid.smiles.count('OH') + amino_acid.smiles.count('O)')
        cooh_count = amino_acid.smiles.count('C(=O)O')
        
        # 至少有一个额外的羟基（除了羧基）
        extra_oh = oh_count - cooh_count
        
        if extra_oh > 0:
            return {
                'is_match': True,
                'strength': min(0.9, 0.6 + extra_oh * 0.2),
                'details': {
                    'hydroxyl_groups': extra_oh,
                    'modification_type': 'hydroxylation'
                }
            }
        
        return {'is_match': False}
    
    def _is_extended_chain_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """验证延长侧链氨基酸"""
        if not amino_acid.smiles:
            return {'is_match': False}
        
        # 计算碳链长度（简化方法）
        carbon_count = amino_acid.smiles.count('C')
        
        # 标准氨基酸的碳原子数量范围通常在2-11之间
        # 如果明显超出，可能是延长链
        if carbon_count > 8:  # 比较长的碳链
            chain_indicators = []
            
            # 检查长链模式
            if re.search(r'CCCCCC', amino_acid.smiles):
                chain_indicators.append('long_aliphatic_chain')
            
            # 检查分支
            if '(' in amino_acid.smiles and ')' in amino_acid.smiles:
                chain_indicators.append('branched_structure')
            
            if chain_indicators or carbon_count > 12:
                return {
                    'is_match': True,
                    'strength': min(0.8, 0.5 + (carbon_count - 8) * 0.05),
                    'details': {
                        'carbon_count': carbon_count,
                        'chain_indicators': chain_indicators,
                        'structure_complexity': 'extended'
                    }
                }
        
        return {'is_match': False}
    
    def _create_unclassified_result(self, amino_acid: AminoAcidInfo, reason: str) -> Dict[str, Any]:
        """创建未分类结果"""
        return {
            'amino_acid_id': amino_acid.id,
            'amino_acid_name': amino_acid.name,
            'standard_category': 'requires_expert_review',
            'confidence': 0.0,
            'classification_tier': 'tier4_needs_review',
            'detection_method': 'insufficient_analysis',
            'evidence': f'Complex structure: {reason}',
            'smiles': amino_acid.smiles,
            'details': {
                'tier': 4,
                'certainty': 'unknown',
                'manual_review_required': True,
                'unclassified_reason': reason,
                'suggested_analysis': self._suggest_analysis_approach(amino_acid)
            }
        }
    
    def _suggest_analysis_approach(self, amino_acid: AminoAcidInfo) -> List[str]:
        """建议分析方法"""
        suggestions = []
        
        if amino_acid.smiles:
            if len(amino_acid.smiles) > 50:
                suggestions.append("复杂结构需要3D分析")
            
            if '[' in amino_acid.smiles and ']' in amino_acid.smiles:
                suggestions.append("特殊原子标记需要化学分析")
            
            if '=' in amino_acid.smiles:
                suggestions.append("双键或芳香性需要确认")
        
        suggestions.append("建议查阅化学文献")
        suggestions.append("考虑质谱或NMR验证")
        
        return suggestions
    
    def get_classification_statistics(self) -> Dict[str, Any]:
        """获取分类统计信息"""
        total = self.stats['total_processed']
        if total == 0:
            return {'error': 'No amino acids processed yet'}
        
        return {
            'total_processed': total,
            'tier1_definitive': {
                'count': self.stats['tier1_definitive'],
                'percentage': (self.stats['tier1_definitive'] / total) * 100,
                'confidence_level': '100%'
            },
            'tier2_high_confidence': {
                'count': self.stats['tier2_high_confidence'],
                'percentage': (self.stats['tier2_high_confidence'] / total) * 100,
                'confidence_level': '90-95%'
            },
            'tier3_intelligent': {
                'count': self.stats['tier3_intelligent'],
                'percentage': (self.stats['tier3_intelligent'] / total) * 100,
                'confidence_level': '80-90%'
            },
            'tier4_needs_review': {
                'count': self.stats['tier4_needs_review'],
                'percentage': (self.stats['tier4_needs_review'] / total) * 100,
                'confidence_level': '需要人工审查'
            },
            'overall_coverage': {
                'automated_classification': ((total - self.stats['tier4_needs_review']) / total) * 100,
                'high_confidence_coverage': ((self.stats['tier1_definitive'] + self.stats['tier2_high_confidence']) / total) * 100
            }
        }
    
    def classify_from_residue(self, residue: ResidueInfo) -> Dict[str, Any]:
        """从残基信息进行分类（用于PDB分析时的实时分类）"""
        # 创建临时氨基酸对象
        temp_amino_acid = AminoAcidInfo(
            id=residue.residue_name,
            name=f"临时-{residue.residue_name}",
            smiles=getattr(residue, 'smiles', ''),
            molecular_formula=residue.molecular_formula,
            atom_composition=residue.atom_composition
        )
        
        return self.classify_amino_acid(temp_amino_acid)
    
    def batch_classify(self, amino_acids: List[AminoAcidInfo]) -> Dict[str, Any]:
        """批量分类氨基酸"""
        print(f"🚀 开始智能批量分类 {len(amino_acids)} 个氨基酸...")
        
        results = []
        category_distribution = defaultdict(int)
        tier_distribution = defaultdict(int)
        
        for i, amino_acid in enumerate(amino_acids, 1):
            if i % 50 == 0:
                print(f"   进度: {i}/{len(amino_acids)}")
            
            result = self.classify_amino_acid(amino_acid)
            results.append(result)
            
            # 统计
            category_distribution[result['standard_category']] += 1
            tier_distribution[result['classification_tier']] += 1
        
        stats = self.get_classification_statistics()
        
        return {
            'results': results,
            'category_distribution': dict(category_distribution),
            'tier_distribution': dict(tier_distribution),
            'statistics': stats,
            'summary': {
                'total_classified': len(amino_acids),
                'automated_coverage': stats['overall_coverage']['automated_classification'],
                'high_confidence_coverage': stats['overall_coverage']['high_confidence_coverage']
            }
        }