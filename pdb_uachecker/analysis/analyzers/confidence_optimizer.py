"""
置信度优化器
针对测试中发现的低置信度问题进行优化
"""

from typing import Dict, List, Any, Optional
from dataclasses import dataclass
import math


@dataclass
class AnalyzerResult:
    """分析器结果数据结构"""
    backbone_type: str
    stereochemistry: str
    is_n_methylated: bool
    structural_features: List[str]
    evidence: List[str]
    base_confidence: float


class ConfidenceOptimizer:
    """
    置信度优化器
    
    主要改进：
    1. 基于多重证据的置信度累积
    2. 一致性奖励机制
    3. 特征完整性评分
    4. 化学合理性验证
    """
    
    def __init__(self):
        """初始化置信度优化器"""
        # 基础置信度权重
        self.base_weights = {
            'backbone_analysis': 0.4,      # 骨架分析权重最高
            'stereochemistry': 0.25,       # 立体化学分析
            'structural_features': 0.20,   # 结构特征识别
            'chemical_validity': 0.15      # 化学合理性
        }
        
        # 一致性奖励因子
        self.consistency_bonus = {
            'high_consistency': 0.15,      # 所有分析器一致
            'moderate_consistency': 0.10,  # 大部分分析器一致
            'low_consistency': 0.05        # 部分分析器一致
        }
        
        # 证据质量评分
        self.evidence_quality_weights = {
            'rdkit_analysis': 1.0,         # RDKit分析最可靠
            'cip_rule_strict': 0.95,       # CIP规则严格分析
            'pattern_matching': 0.8,       # 模式匹配
            'heuristic_analysis': 0.6,     # 启发式分析
            'fallback_pattern': 0.4        # 回退模式
        }
        
        # 特征完整性奖励
        self.completeness_bonus = {
            'full_analysis': 0.2,          # 完整分析（所有特征）
            'partial_analysis': 0.1,       # 部分分析
            'minimal_analysis': 0.0        # 最小分析
        }
    
    def optimize_confidence(self, analyzer_results: Dict[str, Any], 
                          validation_result: Optional[Dict] = None) -> float:
        """
        优化置信度计算
        
        Args:
            analyzer_results: 各个分析器的结果
            validation_result: 验证结果（可选）
            
        Returns:
            优化后的置信度分数 (0.0-1.0)
        """
        # Step 1: 基础置信度计算
        base_confidence = self._calculate_base_confidence(analyzer_results)
        
        # Step 2: 一致性分析
        consistency_score = self._analyze_consistency(analyzer_results)
        consistency_bonus = self._get_consistency_bonus(consistency_score)
        
        # Step 3: 证据质量评估
        evidence_quality = self._evaluate_evidence_quality(analyzer_results)
        
        # Step 4: 特征完整性评估
        completeness_score = self._evaluate_completeness(analyzer_results)
        
        # Step 5: 化学合理性检查
        chemical_validity = self._check_chemical_validity(analyzer_results)
        
        # Step 6: 综合计算
        final_confidence = self._combine_confidence_factors(
            base_confidence=base_confidence,
            consistency_bonus=consistency_bonus,
            evidence_quality=evidence_quality,
            completeness_score=completeness_score,
            chemical_validity=chemical_validity
        )
        
        # Step 7: 边界处理和调整
        final_confidence = self._apply_confidence_adjustments(
            final_confidence, analyzer_results
        )
        
        return min(max(final_confidence, 0.0), 1.0)  # 确保在0-1范围内
    
    def _calculate_base_confidence(self, analyzer_results: Dict[str, Any]) -> float:
        """计算基础置信度"""
        weighted_confidence = 0.0
        total_weight = 0.0
        
        # 骨架分析置信度
        if 'backbone' in analyzer_results:
            backbone_conf = analyzer_results['backbone'].get('confidence', 0.0)
            weight = self.base_weights['backbone_analysis']
            weighted_confidence += backbone_conf * weight
            total_weight += weight
        
        # 立体化学分析置信度
        if 'cip_rule' in analyzer_results or 'stereochemistry' in analyzer_results:
            stereo_conf = self._get_stereochemistry_confidence(analyzer_results)
            weight = self.base_weights['stereochemistry']
            weighted_confidence += stereo_conf * weight
            total_weight += weight
        
        # 结构特征置信度
        structural_conf = self._get_structural_features_confidence(analyzer_results)
        weight = self.base_weights['structural_features']
        weighted_confidence += structural_conf * weight
        total_weight += weight
        
        return weighted_confidence / total_weight if total_weight > 0 else 0.0
    
    def _get_stereochemistry_confidence(self, analyzer_results: Dict[str, Any]) -> float:
        """获取立体化学分析的置信度"""
        if 'cip_rule' in analyzer_results:
            cip_result = analyzer_results['cip_rule']
            if cip_result.get('stereochemistry') in ['D_form', 'L_form']:
                return 0.9  # CIP规则分析高置信度
            elif cip_result.get('stereochemistry') in ['D_suspected', 'L_suspected']:
                return 0.7  # CIP规则疑似结果
        
        return 0.5  # 默认中等置信度
    
    def _get_structural_features_confidence(self, analyzer_results: Dict[str, Any]) -> float:
        """获取结构特征的置信度"""
        feature_count = 0
        confidence_sum = 0.0
        
        # N-甲基化分析
        if 'n_methylation' in analyzer_results:
            n_methyl = analyzer_results['n_methylation']
            feature_count += 1
            if n_methyl.get('confidence', 0) > 0.8:
                confidence_sum += 0.9
            else:
                confidence_sum += 0.6
        
        # 分子结构分析
        if 'molecular_structure' in analyzer_results:
            mol_struct = analyzer_results['molecular_structure']
            feature_count += 1
            confidence_sum += mol_struct.get('confidence', 0.5)
        
        return confidence_sum / feature_count if feature_count > 0 else 0.5
    
    def _analyze_consistency(self, analyzer_results: Dict[str, Any]) -> float:
        """分析各分析器结果的一致性"""
        consistency_checks = []
        
        # 检查骨架类型一致性
        backbone_types = []
        if 'backbone' in analyzer_results:
            backbone_types.append(analyzer_results['backbone'].get('backbone_type'))
        
        # 检查立体化学一致性
        stereochemistries = []
        if 'cip_rule' in analyzer_results:
            stereochemistries.append(analyzer_results['cip_rule'].get('stereochemistry'))
        if 'stereochemistry' in analyzer_results:
            stereochemistries.append(analyzer_results['stereochemistry'].get('result'))
        
        # 计算一致性分数
        total_consistency = 0.0
        check_count = 0
        
        if len(backbone_types) > 1:
            consistency = self._calculate_list_consistency(backbone_types)
            total_consistency += consistency
            check_count += 1
        
        if len(stereochemistries) > 1:
            consistency = self._calculate_list_consistency(stereochemistries)
            total_consistency += consistency
            check_count += 1
        
        return total_consistency / check_count if check_count > 0 else 0.8
    
    def _calculate_list_consistency(self, values: List[Any]) -> float:
        """计算列表中值的一致性"""
        if not values:
            return 0.0
        
        # 移除None值
        valid_values = [v for v in values if v is not None]
        if not valid_values:
            return 0.0
        
        # 计算最常见值的比例
        value_counts = {}
        for value in valid_values:
            value_counts[value] = value_counts.get(value, 0) + 1
        
        max_count = max(value_counts.values())
        consistency = max_count / len(valid_values)
        
        return consistency
    
    def _get_consistency_bonus(self, consistency_score: float) -> float:
        """获取一致性奖励"""
        if consistency_score >= 0.9:
            return self.consistency_bonus['high_consistency']
        elif consistency_score >= 0.7:
            return self.consistency_bonus['moderate_consistency']
        elif consistency_score >= 0.5:
            return self.consistency_bonus['low_consistency']
        else:
            return 0.0
    
    def _evaluate_evidence_quality(self, analyzer_results: Dict[str, Any]) -> float:
        """评估证据质量"""
        evidence_scores = []
        
        for analyzer_name, result in analyzer_results.items():
            method = result.get('method', 'unknown')
            base_score = self.evidence_quality_weights.get(method, 0.5)
            
            # 根据证据数量调整分数
            evidence_count = len(result.get('evidence', []))
            evidence_factor = min(1.0, evidence_count / 3.0)  # 3个证据为满分
            
            final_score = base_score * (0.7 + 0.3 * evidence_factor)
            evidence_scores.append(final_score)
        
        return sum(evidence_scores) / len(evidence_scores) if evidence_scores else 0.5
    
    def _evaluate_completeness(self, analyzer_results: Dict[str, Any]) -> float:
        """评估分析完整性"""
        required_analyses = ['backbone', 'cip_rule', 'n_methylation']
        completed_analyses = sum(1 for analysis in required_analyses 
                               if analysis in analyzer_results)
        
        completeness_ratio = completed_analyses / len(required_analyses)
        
        if completeness_ratio >= 0.9:
            return self.completeness_bonus['full_analysis']
        elif completeness_ratio >= 0.6:
            return self.completeness_bonus['partial_analysis']
        else:
            return self.completeness_bonus['minimal_analysis']
    
    def _check_chemical_validity(self, analyzer_results: Dict[str, Any]) -> float:
        """检查化学合理性"""
        validity_score = 0.8  # 基础化学有效性
        
        # 检查是否有明显的化学不合理性
        backbone_type = None
        if 'backbone' in analyzer_results:
            backbone_type = analyzer_results['backbone'].get('backbone_type')
        
        stereochemistry = None
        if 'cip_rule' in analyzer_results:
            stereochemistry = analyzer_results['cip_rule'].get('stereochemistry')
        
        # 基本化学合理性检查
        if backbone_type in ['alpha', 'beta', 'gamma']:
            validity_score += 0.1  # 已知骨架类型
        
        if stereochemistry in ['D_form', 'L_form', 'D_suspected', 'L_suspected']:
            validity_score += 0.1  # 有立体化学信息
        
        return min(validity_score, 1.0)
    
    def _combine_confidence_factors(self, base_confidence: float, consistency_bonus: float,
                                  evidence_quality: float, completeness_score: float,
                                  chemical_validity: float) -> float:
        """组合各种置信度因素"""
        # 基础置信度是主要因素
        combined_confidence = base_confidence
        
        # 添加各种奖励和调整
        combined_confidence += consistency_bonus
        combined_confidence += completeness_score
        
        # 证据质量和化学有效性作为乘法因子
        quality_factor = (evidence_quality + chemical_validity) / 2.0
        combined_confidence *= quality_factor
        
        return combined_confidence
    
    def _apply_confidence_adjustments(self, confidence: float, 
                                    analyzer_results: Dict[str, Any]) -> float:
        """应用最终的置信度调整"""
        adjusted_confidence = confidence
        
        # 如果骨架分析非常明确，给予额外奖励
        if 'backbone' in analyzer_results:
            backbone = analyzer_results['backbone']
            if (backbone.get('backbone_type') == 'alpha' and 
                backbone.get('confidence', 0) > 0.8):
                adjusted_confidence += 0.15  # α-氨基酸高置信度奖励
        
        # 如果有多个分析器都给出了结果，增加置信度
        analyzer_count = len([r for r in analyzer_results.values() 
                            if r.get('confidence', 0) > 0.5])
        if analyzer_count >= 3:
            adjusted_confidence += 0.1  # 多分析器一致性奖励
        
        # 最小置信度保护（避免过低的置信度）
        if adjusted_confidence < 0.3 and self._has_basic_amino_acid_features(analyzer_results):
            adjusted_confidence = 0.5  # 基础氨基酸特征的最小置信度
        
        return adjusted_confidence
    
    def _has_basic_amino_acid_features(self, analyzer_results: Dict[str, Any]) -> bool:
        """检查是否具有基本的氨基酸特征"""
        # 检查是否有骨架识别
        has_backbone = ('backbone' in analyzer_results and 
                       analyzer_results['backbone'].get('backbone_type') in 
                       ['alpha', 'beta', 'gamma'])
        
        # 检查是否有立体化学识别
        has_stereochemistry = ('cip_rule' in analyzer_results and
                             analyzer_results['cip_rule'].get('stereochemistry') is not None)
        
        return has_backbone or has_stereochemistry


def create_confidence_optimizer():
    """创建置信度优化器实例"""
    return ConfidenceOptimizer()