"""
分类验证层
确保氨基酸分类结果的一致性和准确性
集成所有分析器的结果，进行交叉验证和一致性检查
"""

from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass
from enum import Enum

# 导入新的分析器
from .analyzers.cip_rule_analyzer import CIPRuleAnalyzer
from .analyzers.robust_chirality_analyzer import RobustChiralityAnalyzer
from .analyzers.enhanced_backbone_analyzer import EnhancedBackboneAnalyzer
from .analyzers.n_methylation_analyzer import NMethylationAnalyzer

# 导入现有分析器
from .analyzers.stereochemistry_analyzer import StereochemistryAnalyzer
from .molecular_structure_analyzer import MolecularStructureAnalyzer


class ValidationLevel(Enum):
    """验证等级"""
    STRICT = "strict"           # 严格验证，所有分析器必须一致
    MODERATE = "moderate"       # 中等验证，多数分析器一致
    LENIENT = "lenient"        # 宽松验证，允许部分不一致


@dataclass
class ValidationResult:
    """验证结果"""
    is_consistent: bool
    confidence_score: float
    inconsistencies: List[str]
    recommendations: List[str]
    final_categories: List[str]
    validation_details: Dict[str, Any]


class ClassificationValidator:
    """
    分类验证层
    
    核心功能：
    1. 整合所有分析器结果
    2. 检测分类不一致性
    3. 解决分类冲突
    4. 提供最终可信分类结果
    5. 生成验证报告
    """
    
    def __init__(self, validation_level: ValidationLevel = ValidationLevel.MODERATE):
        """初始化分类验证器"""
        self.validation_level = validation_level
        
        # 初始化所有分析器
        self.analyzers = {
            'cip_rule': CIPRuleAnalyzer(),
            'robust_chirality': RobustChiralityAnalyzer(),  # 高级手性分析器
            'backbone': EnhancedBackboneAnalyzer(),
            'n_methylation': NMethylationAnalyzer(),
            'stereochemistry': StereochemistryAnalyzer(),
            'molecular_structure': MolecularStructureAnalyzer()
        }
        
        # 分类标准权重（根据可信度调整）
        self.analyzer_weights = {
            'cip_rule': 0.25,              # CIP规则权重最高
            'backbone': 0.25,              # 骨架分析权重最高
            'n_methylation': 0.20,         # N-甲基化分析
            'stereochemistry': 0.15,       # 立体化学分析（作为CIP的补充）
            'molecular_structure': 0.15    # 分子结构分析
        }
        
        # 分类一致性规则
        self.consistency_rules = {
            'stereochemistry_consistency': {
                'description': 'CIP规则和立体化学分析应该一致',
                'weight': 0.3
            },
            'backbone_consistency': {
                'description': '骨架分析和分子结构分析应该一致',
                'weight': 0.3
            },
            'methylation_consistency': {
                'description': 'N-甲基化分析应该与其他结构特征一致',
                'weight': 0.2
            },
            'overall_coherence': {
                'description': '所有分类应该在化学上合理',
                'weight': 0.2
            }
        }
    
    def validate_classification(self, smiles: str, amino_acid_code: Optional[str] = None,
                              amino_acid_name: Optional[str] = None) -> ValidationResult:
        """
        验证氨基酸分类
        
        Args:
            smiles: SMILES字符串
            amino_acid_code: 氨基酸代码
            amino_acid_name: 氨基酸名称
            
        Returns:
            完整的验证结果
        """
        if not smiles:
            return self._create_failed_validation("No SMILES provided")
        
        # Step 1: 运行所有分析器
        analyzer_results = self._run_all_analyzers(smiles, amino_acid_code, amino_acid_name)
        
        # Step 2: 提取分类信息
        classification_summary = self._extract_classifications(analyzer_results)
        
        # Step 3: 检测不一致性
        inconsistency_analysis = self._detect_inconsistencies(classification_summary)
        
        # Step 4: 解决冲突
        conflict_resolution = self._resolve_conflicts(classification_summary, inconsistency_analysis)
        
        # Step 5: 计算最终置信度
        final_confidence = self._calculate_final_confidence(
            analyzer_results, inconsistency_analysis, conflict_resolution
        )
        
        # Step 6: 生成最终分类
        final_categories = self._generate_final_categories(conflict_resolution)
        
        # Step 7: 生成建议
        recommendations = self._generate_recommendations(
            inconsistency_analysis, conflict_resolution, final_confidence
        )
        
        return ValidationResult(
            is_consistent=(len(inconsistency_analysis['major_inconsistencies']) == 0),
            confidence_score=final_confidence,
            inconsistencies=inconsistency_analysis['all_inconsistencies'],
            recommendations=recommendations,
            final_categories=final_categories,
            validation_details={
                'analyzer_results': analyzer_results,
                'classification_summary': classification_summary,
                'inconsistency_analysis': inconsistency_analysis,
                'conflict_resolution': conflict_resolution,
                'validation_level': self.validation_level.value
            }
        )
    
    def _run_all_analyzers(self, smiles: str, amino_acid_code: Optional[str], 
                          amino_acid_name: Optional[str]) -> Dict[str, Any]:
        """运行所有分析器"""
        results = {}
        
        # 智能手性分析：先尝试基础CIP，如果不可靠则使用强化分析器
        results['cip_rule'] = self._analyze_chirality_intelligent(smiles, amino_acid_name or "", amino_acid_code or "")
        
        try:
            # 骨架分析
            results['backbone'] = self.analyzers['backbone'].analyze_backbone_precise(
                smiles, amino_acid_code, amino_acid_name
            )
        except Exception as e:
            results['backbone'] = {'error': str(e), 'backbone_type': 'unknown'}
        
        try:
            # N-甲基化分析
            results['n_methylation'] = self.analyzers['n_methylation'].analyze_n_methylation(
                smiles, amino_acid_code, amino_acid_name
            )
        except Exception as e:
            results['n_methylation'] = {'error': str(e), 'is_n_methylated': False}
        
        try:
            # 立体化学分析（作为CIP的补充）
            results['stereochemistry'] = self.analyzers['stereochemistry'].analyze(smiles)
        except Exception as e:
            results['stereochemistry'] = {'error': str(e), 'stereochemistry': 'unknown'}
        
        try:
            # 分子结构分析
            analysis_result = self.analyzers['molecular_structure'].analyze(smiles)
            # 将MolecularAnalysis对象转换为字典
            results['molecular_structure'] = {
                'is_valid': analysis_result.is_valid,
                'aromatic_atoms': analysis_result.aromatic_atoms,
                'ring_systems': analysis_result.ring_systems,
                'chiral_centers': analysis_result.chiral_centers,
                'functional_groups': analysis_result.functional_groups,
                'backbone_analysis': analysis_result.backbone_analysis,
                'molecular_descriptors': analysis_result.molecular_descriptors
            }
        except Exception as e:
            results['molecular_structure'] = {'error': str(e), 'is_valid': False}
        
        return results
    
    def _extract_classifications(self, analyzer_results: Dict[str, Any]) -> Dict[str, Any]:
        """提取和标准化分类信息"""
        classifications = {
            'stereochemistry': [],
            'backbone_type': [],
            'n_methylation': [],
            'structural_features': [],
            'all_categories': []
        }
        
        # 提取立体化学分类
        cip_result = analyzer_results.get('cip_rule', {})
        stereo_result = analyzer_results.get('stereochemistry', {})
        
        if cip_result.get('stereochemistry') and cip_result['stereochemistry'] != 'unknown_chirality':
            classifications['stereochemistry'].append({
                'source': 'cip_rule',
                'value': cip_result['stereochemistry'],
                'confidence': cip_result.get('confidence', 0)
            })
        
        if stereo_result.get('stereochemistry') and stereo_result['stereochemistry'] != 'unknown':
            classifications['stereochemistry'].append({
                'source': 'stereochemistry',
                'value': stereo_result['stereochemistry'],
                'confidence': stereo_result.get('confidence', 0)
            })
        
        # 提取骨架类型
        backbone_result = analyzer_results.get('backbone', {})
        if backbone_result.get('backbone_type') and backbone_result['backbone_type'] != 'unknown':
            classifications['backbone_type'].append({
                'source': 'backbone',
                'value': backbone_result['backbone_type'],
                'confidence': backbone_result.get('confidence', 0)
            })
        
        # 提取N-甲基化
        n_methyl_result = analyzer_results.get('n_methylation', {})
        if n_methyl_result.get('is_n_methylated'):
            classifications['n_methylation'].append({
                'source': 'n_methylation',
                'value': 'n_methyl',
                'confidence': n_methyl_result.get('confidence', 0)
            })
        
        # 提取结构特征
        molecular_result = analyzer_results.get('molecular_structure', {})
        if molecular_result.get('is_valid'):
            # 添加检测到的结构特征
            if molecular_result.get('aromatic_atoms'):
                classifications['structural_features'].append({
                    'source': 'molecular_structure',
                    'value': 'aromatic',
                    'confidence': 0.9
                })
            if molecular_result.get('ring_systems'):
                classifications['structural_features'].append({
                    'source': 'molecular_structure', 
                    'value': 'cyclic',
                    'confidence': 0.9
                })
        
        # 收集所有类别
        for category_type, items in classifications.items():
            if category_type != 'all_categories':
                for item in items:
                    classifications['all_categories'].append(item)
        
        return classifications
    
    def _detect_inconsistencies(self, classification_summary: Dict[str, Any]) -> Dict[str, Any]:
        """检测分类不一致性"""
        inconsistencies = {
            'major_inconsistencies': [],
            'minor_inconsistencies': [],
            'all_inconsistencies': [],
            'consistency_score': 1.0
        }
        
        # 检查立体化学一致性
        stereo_classes = classification_summary['stereochemistry']
        if len(stereo_classes) > 1:
            values = [c['value'] for c in stereo_classes]
            if len(set(values)) > 1:
                inconsistency = f"立体化学不一致: {values}"
                inconsistencies['major_inconsistencies'].append(inconsistency)
                inconsistencies['consistency_score'] -= 0.3
        
        # 检查骨架类型一致性
        backbone_classes = classification_summary['backbone_type']
        if len(backbone_classes) > 1:
            values = [c['value'] for c in backbone_classes]
            if len(set(values)) > 1:
                inconsistency = f"骨架类型不一致: {values}"
                inconsistencies['major_inconsistencies'].append(inconsistency)
                inconsistencies['consistency_score'] -= 0.3
        
        # 检查化学合理性
        chemistry_check = self._check_chemical_reasonableness(classification_summary)
        if not chemistry_check['is_reasonable']:
            inconsistencies['minor_inconsistencies'].extend(chemistry_check['issues'])
            inconsistencies['consistency_score'] -= 0.1
        
        # 合并所有不一致性
        inconsistencies['all_inconsistencies'] = (
            inconsistencies['major_inconsistencies'] + 
            inconsistencies['minor_inconsistencies']
        )
        
        return inconsistencies
    
    def _check_chemical_reasonableness(self, classification_summary: Dict[str, Any]) -> Dict[str, Any]:
        """检查化学合理性"""
        issues = []
        
        # 检查分类组合的合理性
        all_values = [item['value'] for item in classification_summary['all_categories']]
        
        # 例如：D型和L型同时存在是不合理的
        if 'D_form' in all_values and 'L_form' in all_values:
            issues.append("同时检测到D型和L型，不合理")
        
        # 检查骨架类型的合理性
        backbone_types = [item['value'] for item in classification_summary['backbone_type']]
        unusual_backbones = [b for b in backbone_types if b not in ['alpha', 'beta', 'gamma']]
        if unusual_backbones:
            issues.append(f"不常见骨架类型: {unusual_backbones}")
        
        return {
            'is_reasonable': len(issues) == 0,
            'issues': issues
        }
    
    def _resolve_conflicts(self, classification_summary: Dict[str, Any], 
                          inconsistency_analysis: Dict[str, Any]) -> Dict[str, Any]:
        """解决分类冲突"""
        resolved_categories = {}
        resolution_methods = []
        
        # 解决立体化学冲突
        stereo_classes = classification_summary['stereochemistry']
        if len(stereo_classes) > 1:
            # 选择置信度最高的
            best_stereo = max(stereo_classes, key=lambda x: x['confidence'])
            resolved_categories['stereochemistry'] = best_stereo['value']
            resolution_methods.append(f"立体化学冲突解决：选择最高置信度 {best_stereo['value']}")
        elif len(stereo_classes) == 1:
            resolved_categories['stereochemistry'] = stereo_classes[0]['value']
        
        # 解决骨架类型冲突
        backbone_classes = classification_summary['backbone_type']
        if len(backbone_classes) > 1:
            best_backbone = max(backbone_classes, key=lambda x: x['confidence'])
            resolved_categories['backbone_type'] = best_backbone['value']
            resolution_methods.append(f"骨架类型冲突解决：选择最高置信度 {best_backbone['value']}")
        elif len(backbone_classes) == 1:
            resolved_categories['backbone_type'] = backbone_classes[0]['value']
        
        # N-甲基化不会有冲突（只有一个分析器）
        n_methyl_classes = classification_summary['n_methylation']
        if n_methyl_classes:
            resolved_categories['n_methylation'] = 'n_methyl'
        
        # 结构特征
        structural_features = list(set(item['value'] for item in classification_summary['structural_features']))
        if structural_features:
            resolved_categories['structural_features'] = structural_features
        
        return {
            'resolved_categories': resolved_categories,
            'resolution_methods': resolution_methods,
            'conflict_count': len(resolution_methods)
        }
    
    def _calculate_final_confidence(self, analyzer_results: Dict[str, Any], 
                                   inconsistency_analysis: Dict[str, Any],
                                   conflict_resolution: Dict[str, Any]) -> float:
        """计算最终置信度"""
        # 基础置信度：各分析器置信度加权平均
        weighted_confidences = []
        
        for analyzer_name, weight in self.analyzer_weights.items():
            result = analyzer_results.get(analyzer_name, {})
            confidence = result.get('confidence', 0)
            if 'error' not in result and confidence > 0:
                weighted_confidences.append(confidence * weight)
        
        base_confidence = sum(weighted_confidences) if weighted_confidences else 0.5
        
        # 一致性调整
        consistency_penalty = (1 - inconsistency_analysis['consistency_score']) * 0.2
        final_confidence = base_confidence - consistency_penalty
        
        # 冲突解决调整
        conflict_penalty = conflict_resolution['conflict_count'] * 0.05
        final_confidence -= conflict_penalty
        
        return max(0.0, min(1.0, final_confidence))
    
    def _generate_final_categories(self, conflict_resolution: Dict[str, Any]) -> List[str]:
        """生成最终分类类别"""
        categories = []
        
        resolved = conflict_resolution['resolved_categories']
        
        # 添加骨架类型
        if 'backbone_type' in resolved:
            categories.append(f"{resolved['backbone_type']}_amino_acid")
        
        # 添加立体化学
        if 'stereochemistry' in resolved:
            stereo = resolved['stereochemistry']
            if stereo in ['D_form', 'L_form']:
                categories.append(stereo.lower())
        
        # 添加N-甲基化
        if 'n_methylation' in resolved:
            categories.append('n_methyl_amino_acid')
        
        # 添加结构特征
        if 'structural_features' in resolved:
            categories.extend([f"{feature}_amino_acid" for feature in resolved['structural_features']])
        
        # 去重并排序
        return sorted(list(set(categories)))
    
    def _generate_recommendations(self, inconsistency_analysis: Dict[str, Any],
                                conflict_resolution: Dict[str, Any], 
                                final_confidence: float) -> List[str]:
        """生成建议"""
        recommendations = []
        
        # 基于置信度的建议
        if final_confidence > 0.9:
            recommendations.append("分类结果高度可信，建议直接采用")
        elif final_confidence > 0.7:
            recommendations.append("分类结果较为可信，建议采用但需注意不确定性")
        elif final_confidence > 0.5:
            recommendations.append("分类结果不确定性较高，建议人工复核")
        else:
            recommendations.append("分类结果不可信，强烈建议人工分析")
        
        # 基于不一致性的建议
        if inconsistency_analysis['major_inconsistencies']:
            recommendations.append("检测到主要不一致性，建议检查原始数据和分析过程")
        
        if conflict_resolution['conflict_count'] > 0:
            recommendations.append("存在分类冲突，建议查看冲突解决过程和原因")
        
        # 基于验证等级的建议
        if self.validation_level == ValidationLevel.STRICT and not inconsistency_analysis['all_inconsistencies']:
            recommendations.append("通过严格验证，分类结果高度可靠")
        
        return recommendations
    
    def _create_failed_validation(self, reason: str) -> ValidationResult:
        """创建失败的验证结果"""
        return ValidationResult(
            is_consistent=False,
            confidence_score=0.0,
            inconsistencies=[reason],
            recommendations=["验证失败，无法进行分类"],
            final_categories=[],
            validation_details={'failure_reason': reason}
        )
    
    def _analyze_chirality_intelligent(self, smiles: str, amino_acid_name: str, amino_acid_code: str) -> Dict[str, Any]:
        """
        智能手性分析：层级选择机制
        
        1. 首先使用基础CIP分析器（快速、稳定）
        2. 如果结果不可靠（置信度低或出错），则升级到强化分析器
        3. 强化分析器提供多方法验证和一致性检测
        """
        try:
            # Step 1: 尝试基础CIP规则分析
            basic_result = self.analyzers['cip_rule'].analyze_stereochemistry(smiles, amino_acid_name)
            
            # Step 2: 评估基础分析结果的可靠性
            needs_robust_analysis = self._should_use_robust_chirality_analysis(basic_result, smiles)
            
            if not needs_robust_analysis:
                # 基础分析结果可靠，直接使用
                basic_result['analysis_method'] = 'basic_cip'
                return basic_result
            
            # Step 3: 升级到强化分析器
            try:
                robust_result = self.analyzers['robust_chirality'].analyze_chirality_robust(
                    smiles, amino_acid_name, amino_acid_code
                )
                
                # 添加升级原因到证据中
                if 'evidence' not in robust_result:
                    robust_result['evidence'] = []
                robust_result['evidence'].append("因基础CIP分析不可靠，升级到强化分析器")
                
                return robust_result
                
            except Exception as robust_error:
                # 强化分析失败，回退到基础分析并标记
                basic_result['fallback_reason'] = f"强化分析失败: {str(robust_error)}"
                basic_result['analysis_method'] = 'basic_cip_fallback'
                return basic_result
                
        except Exception as basic_error:
            # 基础分析完全失败，尝试强化分析器
            try:
                robust_result = self.analyzers['robust_chirality'].analyze_chirality_robust(
                    smiles, amino_acid_name, amino_acid_code
                )
                robust_result['analysis_method'] = 'robust_emergency'
                if 'evidence' not in robust_result:
                    robust_result['evidence'] = []
                robust_result['evidence'].append(f"基础CIP分析失败({str(basic_error)})，直接使用强化分析器")
                return robust_result
                
            except Exception as robust_error:
                # 两个分析器都失败
                return {
                    'error': f"所有手性分析方法失败 - 基础: {str(basic_error)}, 强化: {str(robust_error)}", 
                    'stereochemistry': 'unknown_chirality',
                    'confidence': 0.0,
                    'analysis_method': 'failed'
                }
    
    def _should_use_robust_chirality_analysis(self, basic_result: Dict[str, Any], smiles: str) -> bool:
        """
        判断是否需要使用强化手性分析器
        
        升级条件：
        1. 基础分析出错
        2. 置信度太低 (<0.6)
        3. 结果不确定 (unknown_chirality)
        4. 分子复杂度高 (多个手性中心)
        5. 包含争议性结构特征
        """
        # 1. 基础分析出错
        if basic_result.get('error') or not basic_result.get('stereochemistry'):
            return True
        
        # 2. 置信度太低
        confidence = basic_result.get('confidence', 0)
        if confidence < 0.6:
            return True
        
        # 3. 结果不确定
        if basic_result.get('stereochemistry') == 'unknown_chirality':
            return True
        
        # 4. 分子复杂度评估
        complexity_score = self._assess_chirality_complexity(smiles)
        if complexity_score > 5:  # 高复杂度阈值
            return True
        
        # 5. 检测争议性特征（多个@标记、复杂分支等）
        if '@' in smiles and smiles.count('@') > 1:  # 多个手性中心
            return True
        
        # 其他情况使用基础分析器
        return False
    
    def _assess_chirality_complexity(self, smiles: str) -> int:
        """
        评估手性复杂度分数
        
        考虑因素：
        - 手性中心数量 (每个+2分)
        - 环结构 (每个+1分)  
        - 分支结构 (每个+1分)
        - 芳香原子 (每个+0.5分)
        """
        score = 0
        
        # 手性中心
        score += smiles.count('@') * 2
        
        # 环结构（数字1-9表示环闭合）
        for digit in '123456789':
            score += smiles.count(digit)
        
        # 分支结构
        score += smiles.count('(')
        
        # 芳香原子（小写字母）
        aromatic_count = sum(1 for c in smiles if c.islower() and c.isalpha())
        score += int(aromatic_count * 0.5)
        
        return score