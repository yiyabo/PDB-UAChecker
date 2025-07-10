#!/usr/bin/env python3
"""
索引管理模块
提供高效的数据索引和查询功能
"""

import time
import bisect
from typing import Dict, List, Set, Optional
from collections import defaultdict

from ..core.models import AminoAcidRecord
from ..core.database import AminoAcidDatabase
from ..core.exceptions import SearchError


class IndexManager:
    """高效索引管理系统"""

    def __init__(self, database: AminoAcidDatabase):
        self.database = database
        self.residue_name_index = {}  # 残基名哈希索引 O(1)
        self.molecular_formula_index = {}  # 分子式索引
        self.molecular_weight_index = []  # 分子量B树索引 O(log n)
        self.feature_index = defaultdict(set)  # 特征倒排索引

        self._build_indices()

    def _build_indices(self):
        """构建所有索引"""
        print("正在构建索引...")
        start_time = time.time()

        try:
            amino_acids = self.database.get_all_amino_acids()

            for record in amino_acids:
                self._add_to_indices(record)

            # 对分子量索引排序以支持二分查找
            self.molecular_weight_index.sort(key=lambda x: x[0])

            build_time = time.time() - start_time
            print(f"索引构建完成，耗时 {build_time:.3f}s，包含 {len(amino_acids)} 条记录")
            
        except Exception as e:
            raise SearchError(f"索引构建失败: {e}")

    def _add_to_indices(self, record: AminoAcidRecord):
        """将记录添加到所有索引"""
        # 残基名索引 (哈希表)
        self.residue_name_index[record.id] = record.id

        # 分子式索引
        if record.molecular_formula not in self.molecular_formula_index:
            self.molecular_formula_index[record.molecular_formula] = []
        self.molecular_formula_index[record.molecular_formula].append(record.id)

        # 分子量索引 (用于范围查询)
        self.molecular_weight_index.append((record.molecular_weight, record.id))

        # 特征索引 (倒排索引)
        for feature in record.key_features:
            self.feature_index[feature].add(record.id)

    def find_by_residue_name(self, residue_name: str) -> Optional[str]:
        """通过残基名查找 O(1)"""
        return self.residue_name_index.get(residue_name)

    def find_by_molecular_formula(self, formula: str) -> List[str]:
        """通过分子式查找 O(1)"""
        return self.molecular_formula_index.get(formula, [])

    def find_by_molecular_weight_range(self, min_weight: float,
                                     max_weight: float) -> List[str]:
        """通过分子量范围查找 O(log n)"""
        try:
            # 使用二分查找找到范围
            left_idx = bisect.bisect_left(self.molecular_weight_index, (min_weight, ''))
            right_idx = bisect.bisect_right(self.molecular_weight_index, (max_weight, 'zzz'))

            return [item[1] for item in self.molecular_weight_index[left_idx:right_idx]]
        except Exception as e:
            raise SearchError(f"分子量范围查询失败: {e}")

    def find_by_features(self, features: List[str],
                        match_all: bool = False) -> Set[str]:
        """通过特征查找"""
        if not features:
            return set()

        try:
            if match_all:
                # 交集：必须包含所有特征
                result = self.feature_index[features[0]].copy()
                for feature in features[1:]:
                    result &= self.feature_index[feature]
                return result
            else:
                # 并集：包含任一特征
                result = set()
                for feature in features:
                    result |= self.feature_index[feature]
                return result
        except Exception as e:
            raise SearchError(f"特征查询失败: {e}")

    def rebuild_indices(self):
        """重建所有索引"""
        self.residue_name_index.clear()
        self.molecular_formula_index.clear()
        self.molecular_weight_index.clear()
        self.feature_index.clear()
        
        self._build_indices()

    def add_record_to_indices(self, record: AminoAcidRecord):
        """添加新记录到索引"""
        self._add_to_indices(record)
        # 重新排序分子量索引
        self.molecular_weight_index.sort(key=lambda x: x[0])

    def remove_record_from_indices(self, amino_acid_id: str):
        """从索引中移除记录"""
        try:
            # 从残基名索引移除
            if amino_acid_id in self.residue_name_index:
                del self.residue_name_index[amino_acid_id]

            # 从分子式索引移除
            for formula, ids in self.molecular_formula_index.items():
                if amino_acid_id in ids:
                    ids.remove(amino_acid_id)
                    if not ids:  # 如果列表为空，删除该分子式项
                        del self.molecular_formula_index[formula]
                    break

            # 从分子量索引移除
            self.molecular_weight_index = [
                (weight, id_) for weight, id_ in self.molecular_weight_index 
                if id_ != amino_acid_id
            ]

            # 从特征索引移除
            for feature, ids in self.feature_index.items():
                ids.discard(amino_acid_id)

        except Exception as e:
            raise SearchError(f"从索引移除记录失败: {e}")

    def get_index_statistics(self) -> Dict[str, int]:
        """获取索引统计信息"""
        return {
            'residue_name_count': len(self.residue_name_index),
            'molecular_formula_count': len(self.molecular_formula_index),
            'molecular_weight_count': len(self.molecular_weight_index),
            'feature_count': len(self.feature_index),
            'total_features': sum(len(ids) for ids in self.feature_index.values())
        }