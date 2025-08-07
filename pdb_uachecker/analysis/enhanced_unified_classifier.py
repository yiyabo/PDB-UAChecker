"""
增强统一分类器
整合所有新算法的高精度分类系统
基于严格的分类标准和多层验证
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from enum import Enum

# 导入核心组件
from .classification_validator import ClassificationValidator, ValidationLevel
from .analyzers.cip_rule_analyzer import CIPRuleAnalyzer
from .analyzers.robust_chirality_analyzer import RobustChiralityAnalyzer
from .analyzers.enhanced_backbone_analyzer import EnhancedBackboneAnalyzer
from .analyzers.n_methylation_analyzer import NMethylationAnalyzer

# 导入核心模型
from ..core.models import AminoAcidInfo


class ClassificationTier(Enum):
    """分类层级"""
    DEFINITIVE = "definitive"           # 100%确定 (>0.95)
    HIGH_CONFIDENCE = "high_confidence"  # 高置信度 (>0.85)
    MODERATE = "moderate"               # 中等置信度 (>0.70)
    LOW_CONFIDENCE = "low_confidence"   # 低置信度 (>0.50)
    UNCERTAIN = "uncertain"             # 不确定 (<=0.50)


@dataclass
class EnhancedClassificationResult:
    """增强分类结果"""
    amino_acid_id: str
    amino_acid_name: str
    categories: List[str]
    confidence: float
    classification_tier: ClassificationTier
    validation_passed: bool
    
    # 详细信息
    stereochemistry: Optional[str] = None
    backbone_type: Optional[str] = None
    is_n_methylated: bool = False
    structural_features: List[str] = None
    
    # 分析详情
    classification_method: str = "enhanced_multi_analyzer"
    evidence: List[str] = None
    inconsistencies: List[str] = None
    recommendations: List[str] = None
    
    # 原始分析结果
    raw_analyzer_results: Dict[str, Any] = None
    validation_details: Dict[str, Any] = None


class EnhancedUnifiedClassifier:
    """
    增强统一分类器
    
    核心改进：
    1. 基于严格CIP规则的D/L型判断
    2. 精确的主链氨基骨架识别
    3. 准确的N-甲基化判断（区分主链vs侧链）
    4. 多层验证确保分类一致性
    5. 智能冲突解决机制
    """
    
    def __init__(self, validation_level: ValidationLevel = ValidationLevel.MODERATE):
        """初始化增强分类器"""
        self.validation_level = validation_level
        
        # 核心分析器
        self.cip_analyzer = CIPRuleAnalyzer()
        self.robust_chirality_analyzer = RobustChiralityAnalyzer()  # 高级手性分析器
        self.backbone_analyzer = EnhancedBackboneAnalyzer()
        self.n_methylation_analyzer = NMethylationAnalyzer()
        
        # 验证器
        self.validator = ClassificationValidator(validation_level)
        
        # 性能统计
        self.stats = {
            'total_classified': 0,
            'definitive_results': 0,
            'high_confidence_results': 0,
            'moderate_confidence_results': 0,
            'low_confidence_results': 0,
            'uncertain_results': 0,
            'validation_passed': 0,
            'validation_failed': 0
        }
    
    def classify_amino_acid(self, amino_acid: AminoAcidInfo) -> EnhancedClassificationResult:
        """
        分类氨基酸（增强版）
        
        Args:
            amino_acid: 氨基酸信息
            
        Returns:
            增强分类结果
        """
        self.stats['total_classified'] += 1
        
        # Step 1: 运行分类验证器（包含所有分析器）
        validation_result = self.validator.validate_classification(
            amino_acid.smiles, 
            amino_acid.id, 
            amino_acid.name
        )
        
        # Step 2: 确定分类层级
        classification_tier = self._determine_classification_tier(validation_result.confidence_score)
        
        # Step 3: 提取核心分类信息
        core_classifications = self._extract_core_classifications(validation_result)
        
        # Step 4: 生成最终分类结果
        result = EnhancedClassificationResult(
            amino_acid_id=amino_acid.id,
            amino_acid_name=amino_acid.name,
            categories=validation_result.final_categories,
            confidence=validation_result.confidence_score,
            classification_tier=classification_tier,
            validation_passed=validation_result.is_consistent,
            
            # 详细信息
            stereochemistry=core_classifications.get('stereochemistry'),
            backbone_type=core_classifications.get('backbone_type'),
            is_n_methylated=core_classifications.get('is_n_methylated', False),
            structural_features=core_classifications.get('structural_features', []),
            
            # 分析详情
            evidence=self._compile_classification_evidence(validation_result),
            inconsistencies=validation_result.inconsistencies,
            recommendations=validation_result.recommendations,
            
            # 原始结果
            raw_analyzer_results=validation_result.validation_details.get('analyzer_results', {}),
            validation_details=validation_result.validation_details
        )
        
        # Step 5: 更新统计信息
        self._update_statistics(result)
        
        return result
    
    def _determine_classification_tier(self, confidence: float) -> ClassificationTier:
        """确定分类层级"""
        if confidence > 0.95:
            return ClassificationTier.DEFINITIVE
        elif confidence > 0.85:
            return ClassificationTier.HIGH_CONFIDENCE
        elif confidence > 0.70:
            return ClassificationTier.MODERATE
        elif confidence > 0.50:
            return ClassificationTier.LOW_CONFIDENCE
        else:
            return ClassificationTier.UNCERTAIN
    
    def _extract_core_classifications(self, validation_result) -> Dict[str, Any]:
        """提取核心分类信息"""
        core_info = {}
        
        # 从验证结果中提取信息
        validation_details = validation_result.validation_details
        analyzer_results = validation_details.get('analyzer_results', {})
        conflict_resolution = validation_details.get('conflict_resolution', {})
        resolved_categories = conflict_resolution.get('resolved_categories', {})
        
        # 提取立体化学
        if 'stereochemistry' in resolved_categories:
            core_info['stereochemistry'] = resolved_categories['stereochemistry']
        elif analyzer_results.get('cip_rule', {}).get('stereochemistry'):
            core_info['stereochemistry'] = analyzer_results['cip_rule']['stereochemistry']
        
        # 提取骨架类型
        if 'backbone_type' in resolved_categories:
            core_info['backbone_type'] = resolved_categories['backbone_type']
        elif analyzer_results.get('backbone', {}).get('backbone_type'):
            core_info['backbone_type'] = analyzer_results['backbone']['backbone_type']
        
        # 提取N-甲基化
        n_methyl_result = analyzer_results.get('n_methylation', {})
        core_info['is_n_methylated'] = n_methyl_result.get('is_n_methylated', False)
        
        # 提取结构特征
        if 'structural_features' in resolved_categories:
            core_info['structural_features'] = resolved_categories['structural_features']
        else:
            core_info['structural_features'] = []
        
        return core_info
    
    def _compile_classification_evidence(self, validation_result) -> List[str]:
        """编译分类证据"""
        evidence = []
        
        # 从各个分析器提取证据
        analyzer_results = validation_result.validation_details.get('analyzer_results', {})
        
        for analyzer_name, result in analyzer_results.items():
            if 'evidence' in result and result['evidence']:
                analyzer_display_name = {
                    'cip_rule': 'CIP规则分析',
                    'backbone': '骨架分析',
                    'n_methylation': 'N-甲基化分析',
                    'stereochemistry': '立体化学分析',
                    'molecular_structure': '分子结构分析'
                }.get(analyzer_name, analyzer_name)
                
                evidence.append(f"【{analyzer_display_name}】")
                evidence.extend(result['evidence'])
        
        # 添加验证结果证据
        if validation_result.is_consistent:
            evidence.append("【验证结果】所有分析器结果一致")
        else:
            evidence.append("【验证结果】存在不一致性，已通过冲突解决")
        
        return evidence
    
    def _update_statistics(self, result: EnhancedClassificationResult):
        """更新统计信息"""
        # 更新层级统计
        tier_stats_map = {
            ClassificationTier.DEFINITIVE: 'definitive_results',
            ClassificationTier.HIGH_CONFIDENCE: 'high_confidence_results',
            ClassificationTier.MODERATE: 'moderate_confidence_results',
            ClassificationTier.LOW_CONFIDENCE: 'low_confidence_results',
            ClassificationTier.UNCERTAIN: 'uncertain_results'
        }
        
        stat_key = tier_stats_map.get(result.classification_tier)
        if stat_key:
            self.stats[stat_key] += 1
        
        # 更新验证统计
        if result.validation_passed:
            self.stats['validation_passed'] += 1
        else:
            self.stats['validation_failed'] += 1
    
    def batch_classify(self, amino_acids: List[AminoAcidInfo]) -> List[EnhancedClassificationResult]:
        """批量分类"""
        results = []
        for amino_acid in amino_acids:
            result = self.classify_amino_acid(amino_acid)
            results.append(result)
        return results
    
    def get_classification_statistics(self) -> Dict[str, Any]:
        """获取分类统计信息"""
        total = self.stats['total_classified']
        if total == 0:
            return self.stats
        
        return {
            **self.stats,
            'accuracy_metrics': {
                'definitive_rate': self.stats['definitive_results'] / total,
                'high_confidence_rate': self.stats['high_confidence_results'] / total,
                'reliable_rate': (self.stats['definitive_results'] + self.stats['high_confidence_results']) / total,
                'validation_success_rate': self.stats['validation_passed'] / total
            },
            'confidence_distribution': {
                'definitive': f"{self.stats['definitive_results'] / total:.1%}",
                'high_confidence': f"{self.stats['high_confidence_results'] / total:.1%}",
                'moderate': f"{self.stats['moderate_confidence_results'] / total:.1%}",
                'low_confidence': f"{self.stats['low_confidence_results'] / total:.1%}",
                'uncertain': f"{self.stats['uncertain_results'] / total:.1%}"
            }
        }
    
    def generate_classification_report(self, results: List[EnhancedClassificationResult]) -> str:
        """生成分类报告"""
        if not results:
            return "无分类结果"
        
        report_lines = [
            "# 氨基酸分类报告",
            f"总计分析: {len(results)}个氨基酸",
            "",
            "## 分类统计",
        ]
        
        # 统计各个层级
        tier_counts = {}
        for result in results:
            tier = result.classification_tier.value
            tier_counts[tier] = tier_counts.get(tier, 0) + 1
        
        for tier, count in tier_counts.items():
            percentage = count / len(results) * 100
            report_lines.append(f"- {tier}: {count}个 ({percentage:.1f}%)")
        
        # 统计验证情况
        validation_passed = sum(1 for r in results if r.validation_passed)
        validation_rate = validation_passed / len(results) * 100
        report_lines.extend([
            "",
            f"验证通过率: {validation_passed}/{len(results)} ({validation_rate:.1f}%)",
            "",
            "## 分类详情"
        ])
        
        # 详细结果
        for result in results:
            confidence_str = f"{result.confidence:.1%}"
            categories_str = ", ".join(result.categories) if result.categories else "无特殊分类"
            
            report_lines.append(f"- {result.amino_acid_id}: {categories_str} (置信度: {confidence_str})")
            
            if result.inconsistencies:
                report_lines.append(f"  ⚠️ 不一致性: {'; '.join(result.inconsistencies[:2])}")
            
            if result.recommendations:
                report_lines.append(f"  💡 建议: {result.recommendations[0]}")
        
        return "\n".join(report_lines)
    
    def export_detailed_results(self, results: List[EnhancedClassificationResult], 
                               output_format: str = "json") -> str:
        """导出详细结果"""
        if output_format.lower() == "json":
            import json
            
            export_data = []
            for result in results:
                export_item = {
                    'amino_acid_id': result.amino_acid_id,
                    'amino_acid_name': result.amino_acid_name,
                    'final_categories': result.categories,
                    'confidence': result.confidence,
                    'classification_tier': result.classification_tier.value,
                    'validation_passed': result.validation_passed,
                    'detailed_classification': {
                        'stereochemistry': result.stereochemistry,
                        'backbone_type': result.backbone_type,
                        'is_n_methylated': result.is_n_methylated,
                        'structural_features': result.structural_features
                    },
                    'analysis_info': {
                        'evidence': result.evidence,
                        'inconsistencies': result.inconsistencies,
                        'recommendations': result.recommendations
                    }
                }
                export_data.append(export_item)
            
            return json.dumps(export_data, indent=2, ensure_ascii=False)
        
        elif output_format.lower() == "csv":
            import csv
            from io import StringIO
            
            output = StringIO()
            writer = csv.writer(output)
            
            # 写入标题行
            headers = [
                'amino_acid_id', 'amino_acid_name', 'categories', 'confidence',
                'classification_tier', 'validation_passed', 'stereochemistry',
                'backbone_type', 'is_n_methylated', 'structural_features',
                'inconsistencies', 'recommendations'
            ]
            writer.writerow(headers)
            
            # 写入数据行
            for result in results:
                row = [
                    result.amino_acid_id,
                    result.amino_acid_name,
                    '; '.join(result.categories),
                    f"{result.confidence:.3f}",
                    result.classification_tier.value,
                    result.validation_passed,
                    result.stereochemistry or '',
                    result.backbone_type or '',
                    result.is_n_methylated,
                    '; '.join(result.structural_features) if result.structural_features else '',
                    '; '.join(result.inconsistencies) if result.inconsistencies else '',
                    '; '.join(result.recommendations) if result.recommendations else ''
                ]
                writer.writerow(row)
            
            return output.getvalue()
        
        else:
            raise ValueError(f"不支持的导出格式: {output_format}")


def create_enhanced_classifier(validation_level: str = "moderate") -> EnhancedUnifiedClassifier:
    """创建增强分类器的工厂函数"""
    level_map = {
        'strict': ValidationLevel.STRICT,
        'moderate': ValidationLevel.MODERATE,
        'lenient': ValidationLevel.LENIENT
    }
    
    level = level_map.get(validation_level.lower(), ValidationLevel.MODERATE)
    return EnhancedUnifiedClassifier(level)