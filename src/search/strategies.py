#!/usr/bin/env python3
"""
搜索策略模块
包含各种搜索方法和匹配器
"""

import hashlib
from typing import Dict, List, Set, Optional, Tuple, Any
from collections import defaultdict

from ..core.models import AminoAcidRecord, SearchResult
from ..core.database import AminoAcidDatabase
from ..core.exceptions import SearchError, InvalidQueryError
from .indexing import IndexManager


class ResidueNameMatcher:
    """残基名匹配器，支持精确和模糊匹配"""

    def __init__(self, index_manager: IndexManager):
        self.index_manager = index_manager
        self.exact_match_cache = {}
        self.fuzzy_match_cache = {}

    def exact_match(self, residue_name: str) -> Optional[str]:
        """精确匹配残基名 O(1)"""
        if not residue_name:
            return None
            
        if residue_name in self.exact_match_cache:
            return self.exact_match_cache[residue_name]

        result = self.index_manager.find_by_residue_name(residue_name)
        self.exact_match_cache[residue_name] = result
        return result

    def fuzzy_match(self, residue_name: str, max_distance: int = 1) -> List[Tuple[str, float]]:
        """模糊匹配残基名"""
        if not residue_name:
            return []
            
        cache_key = f"{residue_name}_{max_distance}"
        if cache_key in self.fuzzy_match_cache:
            return self.fuzzy_match_cache[cache_key]

        candidates = []

        try:
            # 获取所有残基名
            all_residue_names = list(self.index_manager.residue_name_index.keys())

            for candidate in all_residue_names:
                distance = self._calculate_edit_distance(residue_name, candidate)
                if distance <= max_distance:
                    similarity = 1.0 - (distance / max(len(residue_name), len(candidate)))
                    candidates.append((candidate, similarity))

            # 按相似度排序
            candidates.sort(key=lambda x: x[1], reverse=True)

            self.fuzzy_match_cache[cache_key] = candidates
            return candidates
        except Exception as e:
            raise SearchError(f"模糊匹配失败: {e}")

    def _calculate_edit_distance(self, s1: str, s2: str) -> int:
        """计算编辑距离（Levenshtein距离）"""
        if len(s1) < len(s2):
            s1, s2 = s2, s1

        if len(s2) == 0:
            return len(s1)

        previous_row = list(range(len(s2) + 1))
        for i, c1 in enumerate(s1):
            current_row = [i + 1]
            for j, c2 in enumerate(s2):
                insertions = previous_row[j + 1] + 1
                deletions = current_row[j] + 1
                substitutions = previous_row[j] + (c1 != c2)
                current_row.append(min(insertions, deletions, substitutions))
            previous_row = current_row

        return previous_row[-1]

    def batch_match(self, residue_names: List[str]) -> Dict[str, Optional[str]]:
        """批量匹配残基名"""
        results = {}
        for name in residue_names:
            results[name] = self.exact_match(name)
        return results

    def get_match_statistics(self) -> Dict[str, int]:
        """获取匹配统计信息"""
        return {
            'exact_cache_size': len(self.exact_match_cache),
            'fuzzy_cache_size': len(self.fuzzy_match_cache),
            'total_residue_names': len(self.index_manager.residue_name_index)
        }

    def clear_cache(self):
        """清空缓存"""
        self.exact_match_cache.clear()
        self.fuzzy_match_cache.clear()


class BasicFingerprintMatcher:
    """基础分子指纹匹配器"""

    def __init__(self, database: AminoAcidDatabase):
        self.database = database
        self.fingerprint_cache = {}
        self.similarity_cache = {}

    def generate_basic_fingerprint(self, record: AminoAcidRecord) -> List[int]:
        """生成基础分子指纹"""
        if record.id in self.fingerprint_cache:
            return self.fingerprint_cache[record.id]

        try:
            # 基础指纹：基于原子组成和关键特征
            fingerprint = [0] * 64  # 64位指纹

            # 原子组成特征 (0-31位)
            atom_features = {
                'C': 0, 'H': 1, 'O': 2, 'N': 3, 'S': 4, 'P': 5,
                'F': 6, 'Cl': 7, 'Br': 8, 'I': 9
            }

            for atom, count in record.atom_composition.items():
                if atom in atom_features:
                    bit_pos = atom_features[atom]
                    # 设置原子存在位和数量位
                    for i in range(min(count, 8)):  # 最多8个同类原子
                        if bit_pos * 3 + i < 32:
                            fingerprint[bit_pos * 3 + i] = 1

            # 关键特征 (32-63位)
            feature_bits = {
                'aromatic_ring': 32, 'carboxyl_group': 33, 'amino_group': 34,
                'hydroxyl_group': 35, 'methoxy_group': 36, 'double_bond': 37,
                'triple_bond': 38, 'sulfur_containing': 39, 'phosphate': 40,
                'halogen': 41, 'nitro': 42, 'alpha_amino': 43, 'beta_amino': 44
            }

            for feature in record.key_features:
                if feature in feature_bits:
                    fingerprint[feature_bits[feature]] = 1

            self.fingerprint_cache[record.id] = fingerprint
            return fingerprint
        except Exception as e:
            raise SearchError(f"生成分子指纹失败: {e}")

    def calculate_tanimoto_similarity(self, fp1: List[int], fp2: List[int]) -> float:
        """计算Tanimoto相似性"""
        if len(fp1) != len(fp2):
            return 0.0

        try:
            intersection = sum(a & b for a, b in zip(fp1, fp2))
            union = sum(a | b for a, b in zip(fp1, fp2))
            return intersection / union if union > 0 else 0.0
        except Exception as e:
            raise SearchError(f"计算Tanimoto相似性失败: {e}")

    def find_similar_amino_acids(self, query_record: AminoAcidRecord,
                                threshold: float = 0.7) -> List[Tuple[str, float]]:
        """查找相似的氨基酸"""
        query_fp = self.generate_basic_fingerprint(query_record)
        similarities = []

        try:
            all_records = self.database.get_all_amino_acids()

            for record in all_records:
                if record.id == query_record.id:
                    continue

                record_fp = self.generate_basic_fingerprint(record)
                similarity = self.calculate_tanimoto_similarity(query_fp, record_fp)

                if similarity >= threshold:
                    similarities.append((record.id, similarity))

            # 按相似度排序
            similarities.sort(key=lambda x: x[1], reverse=True)
            return similarities
        except Exception as e:
            raise SearchError(f"查找相似氨基酸失败: {e}")

    def get_fingerprint_statistics(self) -> Dict[str, Any]:
        """获取指纹统计信息"""
        if not self.fingerprint_cache:
            return {'cache_size': 0, 'average_bits_set': 0}

        total_bits = sum(sum(fp) for fp in self.fingerprint_cache.values())
        total_fingerprints = len(self.fingerprint_cache)

        return {
            'cache_size': total_fingerprints,
            'average_bits_set': total_bits / total_fingerprints if total_fingerprints > 0 else 0,
            'fingerprint_length': 64
        }

    def clear_cache(self):
        """清空缓存"""
        self.fingerprint_cache.clear()
        self.similarity_cache.clear()


class MolecularFormulaSearcher:
    """分子式搜索器"""

    def __init__(self, index_manager: IndexManager, database: AminoAcidDatabase):
        self.index_manager = index_manager
        self.database = database

    def search_by_formula(self, formula: str) -> List[str]:
        """按分子式搜索"""
        if not formula:
            raise InvalidQueryError("分子式不能为空")

        try:
            return self.index_manager.find_by_molecular_formula(formula)
        except Exception as e:
            raise SearchError(f"分子式搜索失败: {e}")

    def search_by_composition(self, composition: Dict[str, int]) -> List[Tuple[str, float]]:
        """按原子组成搜索"""
        if not composition:
            raise InvalidQueryError("原子组成不能为空")

        try:
            results = []
            all_records = self.database.get_all_amino_acids()

            for record in all_records:
                similarity = self._calculate_composition_similarity(composition, record.atom_composition)
                if similarity > 0.5:  # 相似度阈值
                    results.append((record.id, similarity))

            # 按相似度排序
            results.sort(key=lambda x: x[1], reverse=True)
            return results
        except Exception as e:
            raise SearchError(f"原子组成搜索失败: {e}")

    def _calculate_composition_similarity(self, comp1: Dict[str, int], comp2: Dict[str, int]) -> float:
        """计算原子组成相似度"""
        if not comp1 or not comp2:
            return 0.0

        all_atoms = set(comp1.keys()) | set(comp2.keys())
        total_diff = 0
        total_atoms = 0

        for atom in all_atoms:
            count1 = comp1.get(atom, 0)
            count2 = comp2.get(atom, 0)
            total_diff += abs(count1 - count2)
            total_atoms += max(count1, count2)

        if total_atoms == 0:
            return 0.0

        return 1.0 - (total_diff / total_atoms)


class MolecularWeightSearcher:
    """分子量搜索器"""

    def __init__(self, index_manager: IndexManager, database: AminoAcidDatabase):
        self.index_manager = index_manager
        self.database = database

    def search_by_weight_range(self, target_weight: float, tolerance: float = 5.0) -> List[Tuple[str, float]]:
        """按分子量范围搜索"""
        if target_weight <= 0:
            raise InvalidQueryError("分子量必须大于0")

        try:
            min_weight = target_weight - tolerance
            max_weight = target_weight + tolerance

            amino_acid_ids = self.index_manager.find_by_molecular_weight_range(min_weight, max_weight)
            results = []

            for amino_id in amino_acid_ids:
                record = self.database.get_amino_acid(amino_id)
                if record:
                    weight_diff = abs(record.molecular_weight - target_weight)
                    # 计算权重分数（越接近目标分子量分数越高）
                    score = 1.0 - (weight_diff / tolerance) if tolerance > 0 else 1.0
                    results.append((amino_id, score))

            # 按分数排序
            results.sort(key=lambda x: x[1], reverse=True)
            return results
        except Exception as e:
            raise SearchError(f"分子量搜索失败: {e}")


class FeatureSearcher:
    """特征搜索器"""

    def __init__(self, index_manager: IndexManager):
        self.index_manager = index_manager

    def search_by_features(self, features: List[str], match_all: bool = False) -> Set[str]:
        """按特征搜索"""
        if not features:
            raise InvalidQueryError("特征列表不能为空")

        try:
            return self.index_manager.find_by_features(features, match_all)
        except Exception as e:
            raise SearchError(f"特征搜索失败: {e}")

    def get_available_features(self) -> List[str]:
        """获取可用的特征列表"""
        return list(self.index_manager.feature_index.keys())

    def get_feature_statistics(self) -> Dict[str, int]:
        """获取特征统计信息"""
        stats = {}
        for feature, amino_ids in self.index_manager.feature_index.items():
            stats[feature] = len(amino_ids)
        return stats