"""
验证引擎
整合四重验证算法，提供统一的验证接口
"""

import time
from typing import List, Dict, Any, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

from .molecular_formula import MolecularFormulaVerifier
from .atom_composition import AtomCompositionVerifier
from .fingerprint_similarity import FingerprintSimilarityVerifier
from .structure_3d import Structure3DVerifier
from ..models import (
    ResidueInfo, AminoAcidInfo, VerificationResult, VerificationScore,
    VerificationMethod, VerificationError
)
from ...utils.config import Config


class VerificationEngine:
    """四重验证引擎"""
    
    def __init__(self, config: Optional[Config] = None):
        if config is None:
            from ...utils.config import default_config
            config = default_config
        
        self.config = config
        self.thresholds = config.thresholds
        
        # 初始化验证器
        self.verifiers = {
            VerificationMethod.MOLECULAR_FORMULA: MolecularFormulaVerifier(
                self.thresholds.molecular_formula
            ),
            VerificationMethod.ATOM_COMPOSITION: AtomCompositionVerifier(
                self.thresholds.atom_composition
            ),
            VerificationMethod.FINGERPRINT_SIMILARITY: FingerprintSimilarityVerifier(
                self.thresholds.fingerprint_similarity
            ),
            VerificationMethod.STRUCTURE_3D: Structure3DVerifier(
                self.thresholds.structure_3d
            )
        }
        
        self.enable_parallel = config.performance.enable_parallel
        self.max_workers = config.performance.max_workers
    
    def verify_amino_acid(self, residue: ResidueInfo, amino_acid: AminoAcidInfo,
                         enabled_methods: Optional[List[VerificationMethod]] = None) -> VerificationResult:
        """
        对单个氨基酸执行四重验证
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
            enabled_methods: 启用的验证方法列表，None表示全部启用
        
        Returns:
            验证结果
        """
        if enabled_methods is None:
            enabled_methods = list(self.verifiers.keys())
        
        start_time = time.time()
        
        try:
            # 执行验证
            if self.enable_parallel and len(enabled_methods) > 1:
                scores = self._verify_parallel(residue, amino_acid, enabled_methods)
            else:
                scores = self._verify_sequential(residue, amino_acid, enabled_methods)
            
            # 计算综合置信度
            overall_confidence = self._calculate_overall_confidence(scores)
            
            # 统计通过的验证数量
            passed_count = sum(1 for score in scores if score.passed)
            
            verification_time = time.time() - start_time
            
            return VerificationResult(
                amino_acid_id=amino_acid.id,
                amino_acid_name=amino_acid.name,
                scores=scores,
                overall_confidence=overall_confidence,
                passed_verifications=passed_count,
                total_verifications=len(scores)
            )
        
        except Exception as e:
            raise VerificationError(f"验证失败: {e}")
    
    def _verify_parallel(self, residue: ResidueInfo, amino_acid: AminoAcidInfo,
                        methods: List[VerificationMethod]) -> List[VerificationScore]:
        """并行执行验证"""
        scores = []
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # 提交验证任务
            future_to_method = {
                executor.submit(self.verifiers[method].verify, residue, amino_acid): method
                for method in methods if method in self.verifiers
            }
            
            # 收集结果
            for future in as_completed(future_to_method):
                try:
                    score = future.result(timeout=self.config.performance.timeout_seconds)
                    scores.append(score)
                except Exception as e:
                    method = future_to_method[future]
                    # 创建失败的验证分数
                    failed_score = VerificationScore(
                        method=method,
                        score=0.0,
                        passed=False,
                        details={'error': str(e)}
                    )
                    scores.append(failed_score)
        
        # 按方法顺序排序
        method_order = {method: i for i, method in enumerate(methods)}
        scores.sort(key=lambda s: method_order.get(s.method, 999))
        
        return scores
    
    def _verify_sequential(self, residue: ResidueInfo, amino_acid: AminoAcidInfo,
                          methods: List[VerificationMethod]) -> List[VerificationScore]:
        """顺序执行验证"""
        scores = []
        
        for method in methods:
            if method in self.verifiers:
                try:
                    score = self.verifiers[method].verify(residue, amino_acid)
                    scores.append(score)
                except Exception as e:
                    # 创建失败的验证分数
                    failed_score = VerificationScore(
                        method=method,
                        score=0.0,
                        passed=False,
                        details={'error': str(e)}
                    )
                    scores.append(failed_score)
        
        return scores
    
    def _calculate_overall_confidence(self, scores: List[VerificationScore]) -> float:
        """
        计算综合置信度
        
        采用分层验证策略：核心验证优先，辅助验证加分
        
        Args:
            scores: 验证分数列表
        
        Returns:
            综合置信度 (0.0-1.0)
        """
        if not scores:
            return 0.0
        
        # 分层验证策略
        core_methods = {VerificationMethod.MOLECULAR_FORMULA, VerificationMethod.ATOM_COMPOSITION}
        auxiliary_methods = {VerificationMethod.FINGERPRINT_SIMILARITY, VerificationMethod.STRUCTURE_3D}
        
        core_scores = [s for s in scores if s.method in core_methods]
        auxiliary_scores = [s for s in scores if s.method in auxiliary_methods]
        
        # 核心验证评估
        core_confidence = 0.0
        if core_scores:
            core_passed = [s for s in core_scores if s.passed]
            if len(core_passed) >= 1:  # 至少一个核心验证通过
                # 计算核心验证的加权平均
                core_weights = {
                    VerificationMethod.MOLECULAR_FORMULA: 0.6,
                    VerificationMethod.ATOM_COMPOSITION: 0.4
                }
                
                total_weight = 0.0
                weighted_sum = 0.0
                
                for score in core_scores:
                    weight = core_weights.get(score.method, 0.5)
                    weighted_sum += score.score * weight
                    total_weight += weight
                
                core_confidence = weighted_sum / total_weight if total_weight > 0 else 0.0
        
        # 辅助验证加分
        auxiliary_bonus = 0.0
        if auxiliary_scores:
            auxiliary_passed = [s for s in auxiliary_scores if s.passed]
            if auxiliary_passed:
                # 每个通过的辅助验证提供额外加分
                bonus_per_method = 0.1
                auxiliary_bonus = min(len(auxiliary_passed) * bonus_per_method, 0.2)
        
        # 最终置信度
        final_confidence = min(core_confidence + auxiliary_bonus, 1.0)
        
        # 如果没有任何验证通过，返回0
        if not any(score.passed for score in scores):
            return 0.0
        
        return final_confidence
    
    def batch_verify(self, residue: ResidueInfo, amino_acids: List[AminoAcidInfo],
                    enabled_methods: Optional[List[VerificationMethod]] = None) -> List[VerificationResult]:
        """
        批量验证多个氨基酸
        
        Args:
            residue: 残基信息
            amino_acids: 氨基酸信息列表
            enabled_methods: 启用的验证方法列表
        
        Returns:
            验证结果列表
        """
        results = []
        
        for amino_acid in amino_acids:
            try:
                result = self.verify_amino_acid(residue, amino_acid, enabled_methods)
                results.append(result)
            except Exception as e:
                print(f"⚠️ 验证氨基酸 {amino_acid.id} 失败: {e}")
        
        # 按置信度排序
        results.sort(key=lambda r: r.overall_confidence, reverse=True)
        
        return results
    
    def find_best_matches(self, residue: ResidueInfo, amino_acids: List[AminoAcidInfo],
                         top_k: int = 5, min_confidence: float = 0.0,
                         enabled_methods: Optional[List[VerificationMethod]] = None) -> List[VerificationResult]:
        """
        找到最佳匹配的氨基酸
        
        Args:
            residue: 残基信息
            amino_acids: 氨基酸信息列表
            top_k: 返回前k个结果
            min_confidence: 最小置信度阈值
            enabled_methods: 启用的验证方法列表
        
        Returns:
            最佳匹配结果列表
        """
        # 批量验证
        all_results = self.batch_verify(residue, amino_acids, enabled_methods)
        
        # 过滤低置信度结果
        filtered_results = [r for r in all_results if r.overall_confidence >= min_confidence]
        
        # 返回前k个结果
        return filtered_results[:top_k]
    
    def get_verification_statistics(self, results: List[VerificationResult]) -> Dict[str, Any]:
        """
        获取验证统计信息
        
        Args:
            results: 验证结果列表
        
        Returns:
            统计信息
        """
        if not results:
            return {}
        
        stats = {
            'total_verifications': len(results),
            'successful_matches': sum(1 for r in results if r.is_match),
            'average_confidence': sum(r.overall_confidence for r in results) / len(results),
            'method_success_rates': {},
            'confidence_distribution': {
                'high': 0,    # >= 0.8
                'medium': 0,  # 0.6-0.8
                'low': 0,     # 0.4-0.6
                'very_low': 0 # < 0.4
            }
        }
        
        # 计算各方法成功率
        method_counts = {}
        method_successes = {}
        
        for result in results:
            for score in result.scores:
                method = score.method.value
                method_counts[method] = method_counts.get(method, 0) + 1
                if score.passed:
                    method_successes[method] = method_successes.get(method, 0) + 1
        
        for method, count in method_counts.items():
            success_count = method_successes.get(method, 0)
            stats['method_success_rates'][method] = success_count / count if count > 0 else 0.0
        
        # 置信度分布
        for result in results:
            confidence = result.overall_confidence
            if confidence >= 0.8:
                stats['confidence_distribution']['high'] += 1
            elif confidence >= 0.6:
                stats['confidence_distribution']['medium'] += 1
            elif confidence >= 0.4:
                stats['confidence_distribution']['low'] += 1
            else:
                stats['confidence_distribution']['very_low'] += 1
        
        return stats
    
    def update_thresholds(self, **thresholds):
        """
        更新验证阈值
        
        Args:
            **thresholds: 阈值参数
        """
        for method_name, threshold in thresholds.items():
            if hasattr(VerificationMethod, method_name.upper()):
                method = getattr(VerificationMethod, method_name.upper())
                if method in self.verifiers:
                    self.verifiers[method].set_threshold(threshold)
                    print(f"✅ 更新 {method_name} 阈值为 {threshold}")
    
    def get_current_thresholds(self) -> Dict[str, float]:
        """获取当前阈值设置"""
        return {
            method.value: verifier.get_threshold()
            for method, verifier in self.verifiers.items()
        }
    
    def validate_configuration(self) -> Dict[str, Any]:
        """
        验证配置有效性
        
        Returns:
            验证结果
        """
        validation = {
            'valid': True,
            'errors': [],
            'warnings': [],
            'verifiers_status': {}
        }
        
        # 检查各验证器
        for method, verifier in self.verifiers.items():
            method_name = method.value
            try:
                threshold = verifier.get_threshold()
                if not 0.0 <= threshold <= 1.0:
                    validation['errors'].append(f"{method_name} 阈值超出范围: {threshold}")
                    validation['valid'] = False
                
                validation['verifiers_status'][method_name] = {
                    'available': True,
                    'threshold': threshold
                }
            except Exception as e:
                validation['errors'].append(f"{method_name} 验证器错误: {e}")
                validation['valid'] = False
                validation['verifiers_status'][method_name] = {
                    'available': False,
                    'error': str(e)
                }
        
        # 检查RDKit可用性（影响指纹相似性验证）
        if not self.config.is_rdkit_available():
            validation['warnings'].append("RDKit不可用，指纹相似性验证功能受限")
        
        return validation