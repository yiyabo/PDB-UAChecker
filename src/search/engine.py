#!/usr/bin/env python3
"""
主搜索引擎模块
提供模块化可扩展的搜索引擎功能
"""

import time
from typing import Dict, List, Any, Optional
from dataclasses import dataclass

from ..core.models import AminoAcidRecord, SearchResult
from ..core.database import AminoAcidDatabase
from ..core.exceptions import SearchError, InvalidQueryError
from .indexing import IndexManager
from .strategies import (
    ResidueNameMatcher,
    BasicFingerprintMatcher,
    MolecularFormulaSearcher,
    MolecularWeightSearcher,
    FeatureSearcher
)

# 导入同分异构体识别器
try:
    from isomer_identifier import IsomerIdentifier
    ISOMER_IDENTIFIER_AVAILABLE = True
except ImportError:
    ISOMER_IDENTIFIER_AVAILABLE = False

# 导入原子级分析器
try:
    from atomic_analyzer import AtomicAnalyzer
    ATOMIC_ANALYZER_AVAILABLE = True
except ImportError:
    ATOMIC_ANALYZER_AVAILABLE = False


class ScalableSearchEngine:
    """模块化可扩展搜索引擎"""

    def __init__(self, database_path: str = "amino_acids.db"):
        # 初始化核心组件
        self.database = AminoAcidDatabase(database_path)
        self.index_manager = IndexManager(self.database)
        
        # 初始化搜索策略组件
        self.residue_matcher = ResidueNameMatcher(self.index_manager)
        self.fingerprint_matcher = BasicFingerprintMatcher(self.database)
        self.formula_searcher = MolecularFormulaSearcher(self.index_manager, self.database)
        self.weight_searcher = MolecularWeightSearcher(self.index_manager, self.database)
        self.feature_searcher = FeatureSearcher(self.index_manager)
        
        # 初始化同分异构体识别器
        if ISOMER_IDENTIFIER_AVAILABLE:
            self.isomer_identifier = IsomerIdentifier()
            print("✓ 同分异构体识别器已启用")
        else:
            self.isomer_identifier = None
            print("⚠ 同分异构体识别器未启用，将使用传统方法")

        # 初始化原子级分析器
        if ATOMIC_ANALYZER_AVAILABLE:
            self.atomic_analyzer = AtomicAnalyzer()
            print("✓ 原子级分析器已启用")
        else:
            self.atomic_analyzer = None
            print("⚠ 原子级分析器未启用")

        # 搜索策略配置
        self.search_strategies = {
            'residue_name': self._search_by_residue_name,
            'molecular_formula': self._search_by_molecular_formula,
            'atom_composition': self._search_by_atom_composition,
            'fingerprint_similarity': self._search_by_fingerprint,
            'molecular_weight': self._search_by_molecular_weight,
            'features': self._search_by_features
        }
        
        # 如果启用了同分异构体识别，添加新的搜索策略
        if ISOMER_IDENTIFIER_AVAILABLE:
            self.search_strategies.update({
                'ecfp_similarity': self._search_by_ecfp_similarity,
                'structural_similarity': self._search_by_structural_similarity,
                'isomer_aware': self._search_isomer_aware
            })

        # 置信度权重
        self.confidence_weights = {
            'residue_name': 1.0,
            'molecular_formula': 0.95,
            'atom_composition': 0.85,
            'fingerprint_similarity': 0.8,
            'molecular_weight': 0.7,
            'features': 0.75
        }
        
        # 同分异构体识别相关的置信度权重
        if ISOMER_IDENTIFIER_AVAILABLE:
            self.confidence_weights.update({
                'ecfp_similarity': 0.90,
                'structural_similarity': 0.88,
                'isomer_aware': 0.92
            })

        print(f"搜索引擎初始化完成，支持 {self.database.get_amino_acid_count()} 种氨基酸")

    def search(self, query_data: Dict[str, Any],
               methods: List[str] = None,
               max_results: int = 10) -> List[SearchResult]:
        """执行多策略搜索"""
        if methods is None:
            methods = ['residue_name', 'molecular_formula', 'atom_composition']

        start_time = time.time()
        all_results = []

        # 执行各种搜索策略
        for method in methods:
            if method in self.search_strategies:
                try:
                    results = self.search_strategies[method](query_data)
                    all_results.extend(results)
                except Exception as e:
                    print(f"搜索方法 {method} 执行失败: {e}")

        # 合并和排序结果
        merged_results = self._merge_and_rank_results(all_results)

        search_time = time.time() - start_time
        print(f"搜索完成，耗时 {search_time:.3f}s，找到 {len(merged_results)} 个结果")

        return merged_results[:max_results]

    def _search_by_residue_name(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """残基名搜索策略"""
        results = []
        residue_name = query_data.get('residue_name', '')

        if not residue_name:
            return results

        try:
            # 精确匹配
            amino_id = self.residue_matcher.exact_match(residue_name)
            if amino_id:
                record = self.database.get_amino_acid(amino_id)
                if record:
                    results.append(SearchResult(
                        amino_acid_id=amino_id,
                        match_method='residue_name',
                        confidence_score=1.0,
                        amino_acid_record=record,
                        additional_info={'match_type': 'exact'}
                    ))

            # 模糊匹配
            fuzzy_matches = self.residue_matcher.fuzzy_match(residue_name, max_distance=1)
            for match_id, similarity in fuzzy_matches[:3]:
                if match_id != amino_id:
                    record = self.database.get_amino_acid(match_id)
                    if record:
                        results.append(SearchResult(
                            amino_acid_id=match_id,
                            match_method='residue_name',
                            confidence_score=similarity * 0.8,
                            amino_acid_record=record,
                            additional_info={'match_type': 'fuzzy', 'similarity': similarity}
                        ))
        except Exception as e:
            raise SearchError(f"残基名搜索失败: {e}")

        return results

    def _search_by_molecular_formula(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """分子式搜索策略"""
        results = []
        molecular_formula = query_data.get('molecular_formula', '')

        if not molecular_formula:
            return results

        try:
            amino_ids = self.formula_searcher.search_by_formula(molecular_formula)

            for amino_id in amino_ids:
                record = self.database.get_amino_acid(amino_id)
                if record:
                    results.append(SearchResult(
                        amino_acid_id=amino_id,
                        match_method='molecular_formula',
                        confidence_score=0.95,
                        amino_acid_record=record,
                        additional_info={'formula': molecular_formula}
                    ))
        except Exception as e:
            raise SearchError(f"分子式搜索失败: {e}")

        return results

    def _search_by_atom_composition(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """原子组成搜索策略"""
        results = []
        atom_composition = query_data.get('atom_composition', {})

        if not atom_composition:
            return results

        try:
            composition_results = self.formula_searcher.search_by_composition(atom_composition)

            for amino_id, similarity in composition_results:
                record = self.database.get_amino_acid(amino_id)
                if record:
                    # 增强的同分异构体识别
                    enhanced_similarity = self._enhance_with_isomer_detection(
                        query_data, record, similarity
                    )
                    
                    results.append(SearchResult(
                        amino_acid_id=amino_id,
                        match_method='atom_composition',
                        confidence_score=enhanced_similarity * 0.85,
                        amino_acid_record=record,
                        additional_info={
                            'composition_similarity': similarity,
                            'enhanced_similarity': enhanced_similarity
                        }
                    ))
        except Exception as e:
            raise SearchError(f"原子组成搜索失败: {e}")

        return results

    def _search_by_fingerprint(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """指纹相似性搜索策略"""
        results = []

        try:
            # 需要从查询数据构建临时记录
            if 'smiles' in query_data or 'atom_composition' in query_data:
                query_record = self._create_query_record(query_data)
                similar_acids = self.fingerprint_matcher.find_similar_amino_acids(
                    query_record, threshold=0.7
                )

                for amino_id, similarity in similar_acids:
                    record = self.database.get_amino_acid(amino_id)
                    if record:
                        results.append(SearchResult(
                            amino_acid_id=amino_id,
                            match_method='fingerprint_similarity',
                            confidence_score=similarity * 0.8,
                            amino_acid_record=record,
                            additional_info={'fingerprint_similarity': similarity}
                        ))
        except Exception as e:
            raise SearchError(f"指纹相似性搜索失败: {e}")

        return results

    def _search_by_molecular_weight(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """分子量搜索策略"""
        results = []
        molecular_weight = query_data.get('molecular_weight', 0)
        tolerance = query_data.get('weight_tolerance', 5.0)

        if molecular_weight <= 0:
            return results

        try:
            weight_results = self.weight_searcher.search_by_weight_range(molecular_weight, tolerance)

            for amino_id, score in weight_results:
                record = self.database.get_amino_acid(amino_id)
                if record:
                    results.append(SearchResult(
                        amino_acid_id=amino_id,
                        match_method='molecular_weight',
                        confidence_score=score * 0.7,
                        amino_acid_record=record,
                        additional_info={
                            'target_weight': molecular_weight,
                            'actual_weight': record.molecular_weight,
                            'weight_diff': abs(record.molecular_weight - molecular_weight)
                        }
                    ))
        except Exception as e:
            raise SearchError(f"分子量搜索失败: {e}")

        return results

    def _search_by_features(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """特征搜索策略"""
        results = []
        features = query_data.get('features', [])
        match_all = query_data.get('match_all_features', False)

        if not features:
            return results

        try:
            amino_ids = self.feature_searcher.search_by_features(features, match_all)

            for amino_id in amino_ids:
                record = self.database.get_amino_acid(amino_id)
                if record:
                    # 计算特征匹配度
                    feature_score = self._calculate_feature_score(features, record.key_features, match_all)
                    
                    results.append(SearchResult(
                        amino_acid_id=amino_id,
                        match_method='features',
                        confidence_score=feature_score * 0.75,
                        amino_acid_record=record,
                        additional_info={
                            'query_features': features,
                            'matched_features': list(set(features) & set(record.key_features)),
                            'feature_score': feature_score
                        }
                    ))
        except Exception as e:
            raise SearchError(f"特征搜索失败: {e}")

        return results

    def _search_by_ecfp_similarity(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """ECFP相似性搜索策略（需要同分异构体识别器）"""
        results = []
        
        if not ISOMER_IDENTIFIER_AVAILABLE or not self.isomer_identifier:
            return results

        query_smiles = query_data.get('smiles', '')
        if not query_smiles:
            return results

        try:
            # 使用ECFP进行相似性搜索
            all_records = self.database.get_all_amino_acids()
            
            for record in all_records:
                if record.smiles:
                    similarity = self.isomer_identifier.calculate_ecfp_similarity(
                        query_smiles, record.smiles
                    )
                    
                    if similarity >= 0.7:
                        results.append(SearchResult(
                            amino_acid_id=record.id,
                            match_method='ecfp_similarity',
                            confidence_score=similarity * 0.9,
                            amino_acid_record=record,
                            additional_info={
                                'ecfp_similarity': similarity,
                                'query_smiles': query_smiles
                            }
                        ))
        except Exception as e:
            raise SearchError(f"ECFP相似性搜索失败: {e}")

        return results

    def _search_by_structural_similarity(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """结构相似性搜索策略（需要同分异构体识别器）"""
        results = []
        
        if not ISOMER_IDENTIFIER_AVAILABLE or not self.isomer_identifier:
            return results

        query_smiles = query_data.get('smiles', '')
        if not query_smiles:
            return results

        try:
            # 使用结构指纹进行相似性搜索
            all_records = self.database.get_all_amino_acids()
            
            for record in all_records:
                if record.smiles:
                    similarity = self.isomer_identifier.calculate_structural_similarity(
                        query_smiles, record.smiles
                    )
                    
                    if similarity >= 0.6:
                        results.append(SearchResult(
                            amino_acid_id=record.id,
                            match_method='structural_similarity',
                            confidence_score=similarity * 0.88,
                            amino_acid_record=record,
                            additional_info={
                                'structural_similarity': similarity,
                                'query_smiles': query_smiles
                            }
                        ))
        except Exception as e:
            raise SearchError(f"结构相似性搜索失败: {e}")

        return results

    def _search_isomer_aware(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """同分异构体感知搜索策略（需要同分异构体识别器）"""
        results = []
        
        if not ISOMER_IDENTIFIER_AVAILABLE or not self.isomer_identifier:
            return results

        query_smiles = query_data.get('smiles', '')
        molecular_formula = query_data.get('molecular_formula', '')
        
        if not query_smiles or not molecular_formula:
            return results

        try:
            # 首先找到所有相同分子式的化合物
            formula_matches = self.formula_searcher.search_by_formula(molecular_formula)
            
            for amino_id in formula_matches:
                record = self.database.get_amino_acid(amino_id)
                if record and record.smiles:
                    # 进行同分异构体分析
                    isomer_result = self.isomer_identifier.identify_isomer_relationship(
                        query_smiles, record.smiles
                    )
                    
                    if isomer_result['are_isomers']:
                        confidence = 0.92
                        if isomer_result['isomer_type'] == 'identical':
                            confidence = 1.0
                        elif isomer_result['isomer_type'] == 'stereoisomer':
                            confidence = 0.95
                        elif isomer_result['isomer_type'] == 'structural':
                            confidence = 0.88
                        
                        results.append(SearchResult(
                            amino_acid_id=amino_id,
                            match_method='isomer_aware',
                            confidence_score=confidence,
                            amino_acid_record=record,
                            additional_info={
                                'isomer_type': isomer_result['isomer_type'],
                                'structural_similarity': isomer_result.get('structural_similarity', 0),
                                'query_smiles': query_smiles
                            }
                        ))
        except Exception as e:
            raise SearchError(f"同分异构体感知搜索失败: {e}")

        return results

    def _merge_and_rank_results(self, all_results: List[SearchResult]) -> List[SearchResult]:
        """合并和排序搜索结果"""
        # 按氨基酸ID分组
        result_groups = {}
        for result in all_results:
            amino_id = result.amino_acid_id
            if amino_id not in result_groups:
                result_groups[amino_id] = []
            result_groups[amino_id].append(result)

        # 合并每组结果
        merged_results = []
        for amino_id, group_results in result_groups.items():
            if len(group_results) == 1:
                merged_results.append(group_results[0])
            else:
                # 合并多个结果
                merged_result = self._merge_result_group(group_results)
                merged_results.append(merged_result)

        # 按置信度排序
        merged_results.sort(key=lambda x: x.confidence_score, reverse=True)
        return merged_results

    def _merge_result_group(self, group_results: List[SearchResult]) -> SearchResult:
        """合并同一氨基酸的多个搜索结果"""
        # 选择置信度最高的结果作为基础
        base_result = max(group_results, key=lambda x: x.confidence_score)
        
        # 合并匹配方法
        methods = [result.match_method for result in group_results]
        combined_method = "+".join(sorted(set(methods)))
        
        # 计算加权平均置信度
        weighted_confidence = 0.0
        total_weight = 0.0
        
        for result in group_results:
            weight = self.confidence_weights.get(result.match_method, 0.5)
            weighted_confidence += result.confidence_score * weight
            total_weight += weight
        
        final_confidence = weighted_confidence / total_weight if total_weight > 0 else base_result.confidence_score
        
        # 合并附加信息
        combined_info = {}
        for result in group_results:
            if result.additional_info:
                combined_info.update(result.additional_info)
        combined_info['matched_methods'] = methods
        
        return SearchResult(
            amino_acid_id=base_result.amino_acid_id,
            match_method=combined_method,
            confidence_score=min(final_confidence, 1.0),
            amino_acid_record=base_result.amino_acid_record,
            additional_info=combined_info
        )

    def _create_query_record(self, query_data: Dict[str, Any]) -> AminoAcidRecord:
        """从查询数据创建临时记录"""
        return AminoAcidRecord(
            id="query_temp",
            name=query_data.get('name', 'Query'),
            molecular_formula=query_data.get('molecular_formula', ''),
            molecular_weight=query_data.get('molecular_weight', 0.0),
            smiles=query_data.get('smiles', ''),
            atom_composition=query_data.get('atom_composition', {}),
            key_features=query_data.get('features', [])
        )

    def _enhance_with_isomer_detection(self, query_data: Dict[str, Any], 
                                     record: AminoAcidRecord, 
                                     base_similarity: float) -> float:
        """使用同分异构体检测增强相似度"""
        if not ISOMER_IDENTIFIER_AVAILABLE or not self.isomer_identifier:
            return base_similarity

        query_smiles = query_data.get('smiles', '')
        if not query_smiles or not record.smiles:
            return base_similarity

        try:
            # 计算结构相似性
            structural_similarity = self.isomer_identifier.calculate_structural_similarity(
                query_smiles, record.smiles
            )
            
            # 加权平均
            enhanced_similarity = (base_similarity * 0.6 + structural_similarity * 0.4)
            return min(enhanced_similarity, 1.0)
        except Exception:
            return base_similarity

    def _calculate_feature_score(self, query_features: List[str], 
                               record_features: List[str], 
                               match_all: bool) -> float:
        """计算特征匹配分数"""
        if not query_features:
            return 0.0

        matched_features = set(query_features) & set(record_features)
        
        if match_all:
            return 1.0 if len(matched_features) == len(query_features) else 0.0
        else:
            return len(matched_features) / len(query_features)

    def get_search_statistics(self) -> Dict[str, Any]:
        """获取搜索引擎统计信息"""
        stats = {
            'database_size': self.database.get_amino_acid_count(),
            'index_stats': self.index_manager.get_index_statistics(),
            'residue_matcher_stats': self.residue_matcher.get_match_statistics(),
            'fingerprint_matcher_stats': self.fingerprint_matcher.get_fingerprint_statistics(),
            'available_strategies': list(self.search_strategies.keys()),
            'isomer_identifier_enabled': ISOMER_IDENTIFIER_AVAILABLE,
            'atomic_analyzer_enabled': ATOMIC_ANALYZER_AVAILABLE
        }
        
        return stats

    def clear_caches(self):
        """清空所有缓存"""
        self.residue_matcher.clear_cache()
        self.fingerprint_matcher.clear_cache()
        self.database.clear_cache()
        
        print("所有缓存已清空")

    def rebuild_indices(self):
        """重建所有索引"""
        self.index_manager.rebuild_indices()
        print("所有索引已重建")

    def add_amino_acid(self, record: AminoAcidRecord):
        """添加新的氨基酸记录"""
        try:
            self.database.add_amino_acid(record)
            self.index_manager.add_record_to_indices(record)
            print(f"成功添加氨基酸: {record.id}")
        except Exception as e:
            raise SearchError(f"添加氨基酸失败: {e}")

    def remove_amino_acid(self, amino_acid_id: str):
        """移除氨基酸记录"""
        try:
            self.database.remove_amino_acid(amino_acid_id)
            self.index_manager.remove_record_from_indices(amino_acid_id)
            print(f"成功移除氨基酸: {amino_acid_id}")
        except Exception as e:
            raise SearchError(f"移除氨基酸失败: {e}")
