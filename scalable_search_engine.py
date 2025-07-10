#!/usr/bin/env python3
"""
可扩展非天然氨基酸PDB搜索引擎 - 第一阶段核心架构
支持从5种氨基酸扩展到300-400种的模块化搜索系统

核心组件：
- AminoAcidDatabase: 可扩展数据库管理
- IndexManager: 高效索引系统
- ResidueNameMatcher: 残基名匹配器
- BasicFingerprintMatcher: 基础指纹匹配
- ScalableSearchEngine: 模块化搜索引擎
- 新增：同分异构体识别支持
"""

import os
import json
import sqlite3
import hashlib
import time
from typing import Dict, List, Set, Optional, Tuple, Any, Union
from collections import defaultdict, Counter
from dataclasses import dataclass, asdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import bisect
import pickle

# 导入原子级分析器
try:
    from atomic_analyzer import AtomicAnalyzer, CarbonAnalysis
    ATOMIC_ANALYZER_AVAILABLE = True
except ImportError:
    ATOMIC_ANALYZER_AVAILABLE = False
    print("警告: 原子分析器不可用")

# 导入同分异构体识别器
try:
    from isomer_identifier import IsomerIdentifier, ECFPGenerator, StructuralFingerprint
    ISOMER_IDENTIFIER_AVAILABLE = True
except ImportError:
    ISOMER_IDENTIFIER_AVAILABLE = False
    print("警告: 同分异构体识别器不可用，将使用传统方法")

@dataclass
class AminoAcidRecord:
    """氨基酸记录数据类"""
    id: str
    name: str
    molecular_formula: str
    molecular_weight: float
    smiles: str
    atom_composition: Dict[str, int]
    key_features: List[str]
    structure_data: Optional[Dict] = None
    fingerprints: Optional[Dict] = None
    metadata: Optional[Dict] = None
    
    def to_dict(self) -> Dict:
        """转换为字典格式"""
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'AminoAcidRecord':
        """从字典创建记录"""
        return cls(**data)

class AminoAcidDatabase:
    """可扩展的氨基酸数据库管理器"""
    
    def __init__(self, database_path: str = "amino_acids.db", 
                 cache_size: int = 1000):
        self.database_path = database_path
        self.cache_size = cache_size
        self._cache = {}
        self._cache_order = []
        
        # 原子量表
        self.atomic_weights = {
            'C': 12.01, 'H': 1.008, 'O': 16.00, 'N': 14.01, 
            'S': 32.07, 'P': 30.97, 'F': 19.00, 'Cl': 35.45,
            'Br': 79.90, 'I': 126.90
        }
        
        self._initialize_database()
        self._load_default_amino_acids()
    
    def _initialize_database(self):
        """初始化数据库结构"""
        with sqlite3.connect(self.database_path) as conn:
            conn.execute('''
                CREATE TABLE IF NOT EXISTS amino_acids (
                    id TEXT PRIMARY KEY,
                    name TEXT NOT NULL,
                    molecular_formula TEXT,
                    molecular_weight REAL,
                    smiles TEXT,
                    atom_composition TEXT,
                    key_features TEXT,
                    structure_data TEXT,
                    fingerprints TEXT,
                    metadata TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            
            # 创建索引
            conn.execute('CREATE INDEX IF NOT EXISTS idx_molecular_formula ON amino_acids(molecular_formula)')
            conn.execute('CREATE INDEX IF NOT EXISTS idx_molecular_weight ON amino_acids(molecular_weight)')
            conn.execute('CREATE INDEX IF NOT EXISTS idx_name ON amino_acids(name)')
            
            conn.commit()
    
    def _load_default_amino_acids(self):
        """加载默认的5种氨基酸数据"""
        default_amino_acids = {
            "0A1": {
                "name": "4-甲氧基苯丙氨酸",
                "smiles": "COc1ccc(cc1)C[C@@H](N)C(=O)O",
                "key_features": ["aromatic_ring", "carboxyl_group", "amino_group", "methoxy_group"]
            },
            "0AF": {
                "name": "5-羟基色氨酸", 
                "smiles": "N[C@@H](Cc1c[nH]c2c1cccc2O)C(=O)O",
                "key_features": ["aromatic_ring", "carboxyl_group", "amino_group", "hydroxyl_group"]
            },
            "0BN": {
                "name": "4-胍基苯丙氨酸",
                "smiles": "N[C@@H](Cc1ccc(cc1)C(=N)N)C(=O)O",
                "key_features": ["aromatic_ring", "carboxyl_group", "amino_group"]
            },
            "2AG": {
                "name": "烯丙基甘氨酸",
                "smiles": "N[C@@H](CC=C)C(=O)O",
                "key_features": ["carboxyl_group", "amino_group", "double_bond"]
            },
            "2AS": {
                "name": "天冬氨酸衍生物",
                "smiles": "N[C@@H]([C@H](C)C(=O)O)C(=O)O",
                "key_features": ["carboxyl_group", "amino_group"]
            }
        }
        
        # 检查是否已经加载过数据
        if self.get_amino_acid_count() == 0:
            print("正在加载默认氨基酸数据...")
            for amino_id, data in default_amino_acids.items():
                record = self._create_amino_acid_record(amino_id, data)
                self.add_amino_acid(record)
            print(f"已加载 {len(default_amino_acids)} 种默认氨基酸")
    
    def _create_amino_acid_record(self, amino_id: str, data: Dict) -> AminoAcidRecord:
        """创建氨基酸记录"""
        atom_composition = self._parse_smiles_composition(data["smiles"])
        molecular_formula = self._get_molecular_formula(atom_composition)
        molecular_weight = self._calculate_molecular_weight(atom_composition)
        
        return AminoAcidRecord(
            id=amino_id,
            name=data["name"],
            molecular_formula=molecular_formula,
            molecular_weight=molecular_weight,
            smiles=data["smiles"],
            atom_composition=atom_composition,
            key_features=data["key_features"],
            metadata={"source": "default", "version": "1.0"}
        )
    
    def _parse_smiles_composition(self, smiles: str) -> Dict[str, int]:
        """从SMILES解析原子组成"""
        import re
        composition = defaultdict(int)
        
        # 移除立体化学和电荷符号
        clean_smiles = re.sub(r'[@\[\]()=+\-#]', '', smiles)
        
        # 查找原子符号和数量
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, clean_smiles)
        
        for atom, count in matches:
            count = int(count) if count else 1
            composition[atom] += count
        
        return dict(composition)
    
    def _get_molecular_formula(self, composition: Dict[str, int]) -> str:
        """生成分子式"""
        formula_parts = []
        
        # 按照标准顺序排列：C, H, 其他按字母顺序
        if 'C' in composition:
            count = composition['C']
            formula_parts.append(f"C{count}" if count > 1 else "C")
        
        if 'H' in composition:
            count = composition['H']
            formula_parts.append(f"H{count}" if count > 1 else "H")
        
        for atom in sorted(composition.keys()):
            if atom not in ['C', 'H']:
                count = composition[atom]
                formula_parts.append(f"{atom}{count}" if count > 1 else atom)
        
        return ''.join(formula_parts)
    
    def _calculate_molecular_weight(self, composition: Dict[str, int]) -> float:
        """计算分子量"""
        weight = 0.0
        for atom, count in composition.items():
            weight += self.atomic_weights.get(atom, 0) * count
        return round(weight, 4)
    
    def add_amino_acid(self, record: AminoAcidRecord) -> bool:
        """添加氨基酸记录"""
        try:
            with sqlite3.connect(self.database_path) as conn:
                conn.execute('''
                    INSERT OR REPLACE INTO amino_acids 
                    (id, name, molecular_formula, molecular_weight, smiles, 
                     atom_composition, key_features, structure_data, fingerprints, metadata)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    record.id,
                    record.name,
                    record.molecular_formula,
                    record.molecular_weight,
                    record.smiles,
                    json.dumps(record.atom_composition),
                    json.dumps(record.key_features),
                    json.dumps(record.structure_data) if record.structure_data else None,
                    json.dumps(record.fingerprints) if record.fingerprints else None,
                    json.dumps(record.metadata) if record.metadata else None
                ))
                conn.commit()
            
            # 更新缓存
            self._update_cache(record.id, record)
            return True
            
        except Exception as e:
            print(f"添加氨基酸记录失败: {e}")
            return False
    
    def get_amino_acid(self, amino_id: str) -> Optional[AminoAcidRecord]:
        """获取氨基酸记录"""
        # 先检查缓存
        if amino_id in self._cache:
            self._update_cache_order(amino_id)
            return self._cache[amino_id]
        
        # 从数据库查询
        with sqlite3.connect(self.database_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(
                'SELECT * FROM amino_acids WHERE id = ?', (amino_id,)
            )
            row = cursor.fetchone()
            
            if row:
                record = self._row_to_record(row)
                self._update_cache(amino_id, record)
                return record
        
        return None
    
    def _row_to_record(self, row: sqlite3.Row) -> AminoAcidRecord:
        """将数据库行转换为记录对象"""
        return AminoAcidRecord(
            id=row['id'],
            name=row['name'],
            molecular_formula=row['molecular_formula'],
            molecular_weight=row['molecular_weight'],
            smiles=row['smiles'],
            atom_composition=json.loads(row['atom_composition']) if row['atom_composition'] else {},
            key_features=json.loads(row['key_features']) if row['key_features'] else [],
            structure_data=json.loads(row['structure_data']) if row['structure_data'] else None,
            fingerprints=json.loads(row['fingerprints']) if row['fingerprints'] else None,
            metadata=json.loads(row['metadata']) if row['metadata'] else None
        )
    
    def _update_cache(self, amino_id: str, record: AminoAcidRecord):
        """更新缓存"""
        if amino_id in self._cache:
            self._cache_order.remove(amino_id)
        elif len(self._cache) >= self.cache_size:
            # 移除最旧的记录
            oldest_id = self._cache_order.pop(0)
            del self._cache[oldest_id]
        
        self._cache[amino_id] = record
        self._cache_order.append(amino_id)
    
    def _update_cache_order(self, amino_id: str):
        """更新缓存访问顺序"""
        if amino_id in self._cache_order:
            self._cache_order.remove(amino_id)
            self._cache_order.append(amino_id)
    
    def get_all_amino_acids(self) -> List[AminoAcidRecord]:
        """获取所有氨基酸记录"""
        records = []
        with sqlite3.connect(self.database_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute('SELECT * FROM amino_acids ORDER BY id')
            
            for row in cursor:
                records.append(self._row_to_record(row))
        
        return records
    
    def get_amino_acid_count(self) -> int:
        """获取氨基酸数量"""
        with sqlite3.connect(self.database_path) as conn:
            cursor = conn.execute('SELECT COUNT(*) FROM amino_acids')
            return cursor.fetchone()[0]
    
    def search_by_molecular_formula(self, formula: str) -> List[AminoAcidRecord]:
        """按分子式搜索"""
        records = []
        with sqlite3.connect(self.database_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(
                'SELECT * FROM amino_acids WHERE molecular_formula = ?', (formula,)
            )
            
            for row in cursor:
                records.append(self._row_to_record(row))
        
        return records
    
    def search_by_molecular_weight_range(self, min_weight: float, 
                                       max_weight: float) -> List[AminoAcidRecord]:
        """按分子量范围搜索"""
        records = []
        with sqlite3.connect(self.database_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(
                'SELECT * FROM amino_acids WHERE molecular_weight BETWEEN ? AND ?',
                (min_weight, max_weight)
            )
            
            for row in cursor:
                records.append(self._row_to_record(row))
        
        return records
    
    def get_statistics(self) -> Dict[str, Any]:
        """获取数据库统计信息"""
        with sqlite3.connect(self.database_path) as conn:
            cursor = conn.execute('''
                SELECT 
                    COUNT(*) as total_count,
                    AVG(molecular_weight) as avg_molecular_weight,
                    MIN(molecular_weight) as min_molecular_weight,
                    MAX(molecular_weight) as max_molecular_weight
                FROM amino_acids
            ''')
            stats = cursor.fetchone()
            
            return {
                'total_amino_acids': stats[0],
                'average_molecular_weight': round(stats[1], 2) if stats[1] else 0,
                'molecular_weight_range': (stats[2], stats[3]),
                'cache_size': len(self._cache),
                'cache_hit_ratio': self._calculate_cache_hit_ratio()
            }
    
    def _calculate_cache_hit_ratio(self) -> float:
        """计算缓存命中率（简化实现）"""
        return len(self._cache) / max(self.cache_size, 1)

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

        amino_acids = self.database.get_all_amino_acids()

        for record in amino_acids:
            self._add_to_indices(record)

        # 对分子量索引排序以支持二分查找
        self.molecular_weight_index.sort(key=lambda x: x[0])

        build_time = time.time() - start_time
        print(f"索引构建完成，耗时 {build_time:.3f}s，包含 {len(amino_acids)} 条记录")

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
        # 使用二分查找找到范围
        left_idx = bisect.bisect_left(self.molecular_weight_index, (min_weight, ''))
        right_idx = bisect.bisect_right(self.molecular_weight_index, (max_weight, 'zzz'))

        return [item[1] for item in self.molecular_weight_index[left_idx:right_idx]]

    def find_by_features(self, features: List[str],
                        match_all: bool = False) -> Set[str]:
        """通过特征查找"""
        if not features:
            return set()

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

    def get_index_statistics(self) -> Dict[str, Any]:
        """获取索引统计信息"""
        return {
            'residue_name_index_size': len(self.residue_name_index),
            'molecular_formula_index_size': len(self.molecular_formula_index),
            'molecular_weight_index_size': len(self.molecular_weight_index),
            'feature_index_size': len(self.feature_index),
            'total_features': sum(len(ids) for ids in self.feature_index.values())
        }

    def rebuild_indices(self):
        """重建所有索引"""
        self.residue_name_index.clear()
        self.molecular_formula_index.clear()
        self.molecular_weight_index.clear()
        self.feature_index.clear()

        self._build_indices()

class ResidueNameMatcher:
    """重构的残基名匹配器，支持大规模数据集"""

    def __init__(self, index_manager: IndexManager):
        self.index_manager = index_manager
        self.exact_match_cache = {}
        self.fuzzy_match_cache = {}

    def exact_match(self, residue_name: str) -> Optional[str]:
        """精确匹配残基名 O(1)"""
        if residue_name in self.exact_match_cache:
            return self.exact_match_cache[residue_name]

        result = self.index_manager.find_by_residue_name(residue_name)
        self.exact_match_cache[residue_name] = result
        return result

    def fuzzy_match(self, residue_name: str, max_distance: int = 1) -> List[Tuple[str, float]]:
        """模糊匹配残基名"""
        cache_key = f"{residue_name}_{max_distance}"
        if cache_key in self.fuzzy_match_cache:
            return self.fuzzy_match_cache[cache_key]

        candidates = []

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
                # 使用原子数量的二进制表示
                for i in range(min(count, 8)):  # 最多8位
                    if bit_pos * 3 + i < 32:
                        fingerprint[bit_pos * 3 + i] = 1

        # 关键特征 (32-63位)
        feature_bits = {
            'aromatic_ring': 32,
            'carboxyl_group': 33,
            'amino_group': 34,
            'hydroxyl_group': 35,
            'methoxy_group': 36,
            'double_bond': 37,
            'triple_bond': 38,
            'sulfur_group': 39
        }

        for feature in record.key_features:
            if feature in feature_bits:
                fingerprint[feature_bits[feature]] = 1

        # 分子量范围特征 (40-47位)
        weight_ranges = [
            (0, 100), (100, 150), (150, 200), (200, 250),
            (250, 300), (300, 400), (400, 500), (500, float('inf'))
        ]

        for i, (min_w, max_w) in enumerate(weight_ranges):
            if min_w <= record.molecular_weight < max_w:
                fingerprint[40 + i] = 1
                break

        self.fingerprint_cache[record.id] = fingerprint
        return fingerprint

    def calculate_tanimoto_similarity(self, fp1: List[int], fp2: List[int]) -> float:
        """计算Tanimoto相似性系数"""
        if len(fp1) != len(fp2):
            raise ValueError("指纹长度不匹配")

        intersection = sum(a & b for a, b in zip(fp1, fp2))
        union = sum(a | b for a, b in zip(fp1, fp2))

        return intersection / union if union > 0 else 0.0

    def find_similar_amino_acids(self, query_record: AminoAcidRecord,
                                threshold: float = 0.7) -> List[Tuple[str, float]]:
        """查找相似的氨基酸"""
        query_fp = self.generate_basic_fingerprint(query_record)
        similar_acids = []

        all_records = self.database.get_all_amino_acids()

        for record in all_records:
            if record.id == query_record.id:
                continue

            record_fp = self.generate_basic_fingerprint(record)
            similarity = self.calculate_tanimoto_similarity(query_fp, record_fp)

            if similarity >= threshold:
                similar_acids.append((record.id, similarity))

        # 按相似度排序
        similar_acids.sort(key=lambda x: x[1], reverse=True)
        return similar_acids

    def batch_similarity_search(self, query_records: List[AminoAcidRecord],
                               threshold: float = 0.7) -> Dict[str, List[Tuple[str, float]]]:
        """批量相似性搜索"""
        results = {}

        for query_record in query_records:
            results[query_record.id] = self.find_similar_amino_acids(
                query_record, threshold
            )

        return results

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

@dataclass
class SearchResult:
    """搜索结果数据类"""
    amino_acid_id: str
    match_method: str
    confidence_score: float
    amino_acid_record: AminoAcidRecord
    additional_info: Dict[str, Any] = None

    def to_dict(self) -> Dict:
        """转换为字典格式"""
        result = asdict(self)
        result['amino_acid_record'] = self.amino_acid_record.to_dict()
        return result

class ScalableSearchEngine:
    """模块化可扩展搜索引擎"""

    def __init__(self, database_path: str = "amino_acids.db"):
        # 初始化核心组件
        self.database = AminoAcidDatabase(database_path)
        self.index_manager = IndexManager(self.database)
        self.residue_matcher = ResidueNameMatcher(self.index_manager)
        self.fingerprint_matcher = BasicFingerprintMatcher(self.database)
        
        # 初始化同分异构体识别器
        if ISOMER_IDENTIFIER_AVAILABLE:
            self.isomer_identifier = IsomerIdentifier()
            print("✓ 同分异构体识别器已启用")
        else:
            self.isomer_identifier = None
            print("⚠ 同分异构体识别器未启用，将使用传统方法")

        # 搜索策略配置
        self.search_strategies = {
            'residue_name': self._search_by_residue_name,
            'molecular_formula': self._search_by_molecular_formula,
            'atom_composition': self._search_by_atom_composition,
            'fingerprint_similarity': self._search_by_fingerprint,
            'molecular_weight': self._search_by_molecular_weight
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
            'molecular_weight': 0.7
        }
        
        # 同分异构体识别相关的置信度权重
        if ISOMER_IDENTIFIER_AVAILABLE:
            self.confidence_weights.update({
                'ecfp_similarity': 0.90,        # ECFP相似性
                'structural_similarity': 0.88,  # 结构相似性
                'isomer_aware': 0.92            # 同分异构体感知搜索
            })

        # 初始化原子级分析器
        if ATOMIC_ANALYZER_AVAILABLE:
            self.atomic_analyzer = AtomicAnalyzer()
            print(f"搜索引擎初始化完成，支持 {self.database.get_amino_acid_count()} 种氨基酸")
            print("原子级分析器已启用")
        else:
            self.atomic_analyzer = None
            print(f"搜索引擎初始化完成，支持 {self.database.get_amino_acid_count()} 种氨基酸")
            print("原子级分析器不可用")

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
        for match_id, similarity in fuzzy_matches[:3]:  # 最多3个模糊匹配
            if match_id != amino_id:  # 避免重复
                record = self.database.get_amino_acid(match_id)
                if record:
                    results.append(SearchResult(
                        amino_acid_id=match_id,
                        match_method='residue_name',
                        confidence_score=similarity * 0.8,  # 模糊匹配降低置信度
                        amino_acid_record=record,
                        additional_info={'match_type': 'fuzzy', 'similarity': similarity}
                    ))

        return results

    def _search_by_molecular_formula(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """分子式搜索策略"""
        results = []
        molecular_formula = query_data.get('molecular_formula', '')

        if not molecular_formula:
            return results

        amino_ids = self.index_manager.find_by_molecular_formula(molecular_formula)

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

        return results

    def _search_by_atom_composition(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """原子组成搜索策略"""
        results = []
        atom_composition = query_data.get('atom_composition', {})

        if not atom_composition:
            return results

        all_records = self.database.get_all_amino_acids()

        for record in all_records:
            # 传递SMILES信息以支持同分异构体识别
            query_smiles = query_data.get('smiles', '')
            target_smiles = record.smiles if record.smiles else ''
            
            similarity = self._calculate_composition_similarity(
                atom_composition, record.atom_composition,
                query_smiles, target_smiles
            )

            if similarity >= 0.8:  # 相似度阈值
                results.append(SearchResult(
                    amino_acid_id=record.id,
                    match_method='atom_composition',
                    confidence_score=similarity * 0.85,
                    amino_acid_record=record,
                    additional_info={
                        'composition_similarity': similarity,
                        'structural_enhanced': bool(query_smiles and target_smiles and ISOMER_IDENTIFIER_AVAILABLE)
                    }
                ))

        return results

    def _search_by_fingerprint(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """指纹相似性搜索策略"""
        results = []

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

        return results

    def _search_by_molecular_weight(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """分子量搜索策略"""
        results = []
        molecular_weight = query_data.get('molecular_weight', 0)
        tolerance = query_data.get('weight_tolerance', 5.0)  # 默认±5Da

        if molecular_weight <= 0:
            return results

        min_weight = molecular_weight - tolerance
        max_weight = molecular_weight + tolerance

        amino_ids = self.index_manager.find_by_molecular_weight_range(min_weight, max_weight)

        for amino_id in amino_ids:
            record = self.database.get_amino_acid(amino_id)
            if record:
                weight_diff = abs(record.molecular_weight - molecular_weight)
                confidence = max(0.0, 1.0 - (weight_diff / tolerance)) * 0.7

                results.append(SearchResult(
                    amino_acid_id=amino_id,
                    match_method='molecular_weight',
                    confidence_score=confidence,
                    amino_acid_record=record,
                    additional_info={
                        'weight_difference': weight_diff,
                        'query_weight': molecular_weight
                    }
                ))

        return results

    def _calculate_composition_similarity(self, comp1: Dict[str, int],
                                        comp2: Dict[str, int], 
                                        smiles1: str = None,
                                        smiles2: str = None) -> float:
        """计算原子组成相似度（增强版：支持同分异构体识别）"""
        if not comp1 or not comp2:
            return 0.0

        # 计算基础的原子组成相似度
        basic_similarity = self._calculate_basic_composition_similarity(comp1, comp2)
        
        # 如果有SMILES信息且启用了同分异构体识别，则使用增强算法
        if (ISOMER_IDENTIFIER_AVAILABLE and smiles1 and smiles2 and 
            hasattr(self, 'isomer_identifier')):
            try:
                return self.isomer_identifier.enhanced_composition_similarity(
                    smiles1, smiles2, basic_similarity
                )
            except Exception as e:
                print(f"同分异构体识别失败，使用基础算法: {e}")
                return basic_similarity
        
        return basic_similarity
    
    def _calculate_basic_composition_similarity(self, comp1: Dict[str, int],
                                             comp2: Dict[str, int]) -> float:
        """计算基础原子组成相似度"""
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

    def _create_query_record(self, query_data: Dict[str, Any]) -> AminoAcidRecord:
        """从查询数据创建临时记录"""
        return AminoAcidRecord(
            id='QUERY',
            name='Query Record',
            molecular_formula=query_data.get('molecular_formula', ''),
            molecular_weight=query_data.get('molecular_weight', 0.0),
            smiles=query_data.get('smiles', ''),
            atom_composition=query_data.get('atom_composition', {}),
            key_features=query_data.get('key_features', [])
        )

    def _merge_and_rank_results(self, results: List[SearchResult]) -> List[SearchResult]:
        """合并和排序搜索结果"""
        # 按氨基酸ID分组
        grouped_results = defaultdict(list)
        for result in results:
            grouped_results[result.amino_acid_id].append(result)

        # 合并同一氨基酸的多个匹配结果
        merged_results = []
        for amino_id, group_results in grouped_results.items():
            if len(group_results) == 1:
                merged_results.append(group_results[0])
            else:
                # 选择置信度最高的结果，并合并信息
                best_result = max(group_results, key=lambda x: x.confidence_score)

                # 合并匹配方法信息
                all_methods = [r.match_method for r in group_results]
                all_info = {}
                for r in group_results:
                    if r.additional_info:
                        all_info.update(r.additional_info)

                best_result.additional_info = all_info
                best_result.additional_info['all_match_methods'] = all_methods

                merged_results.append(best_result)

        # 按置信度排序
        merged_results.sort(key=lambda x: x.confidence_score, reverse=True)
        return merged_results

    def search_pdb_residues(self, pdb_residues: List[Dict[str, Any]]) -> List[SearchResult]:
        """搜索PDB残基（保持向后兼容）"""
        all_results = []

        for residue_data in pdb_residues:
            query_data = {
                'residue_name': residue_data.get('residue_name', ''),
                'atom_composition': residue_data.get('composition', {}),
                'molecular_weight': self._calculate_molecular_weight_from_composition(
                    residue_data.get('composition', {})
                )
            }

            results = self.search(query_data, methods=['residue_name', 'atom_composition'])
            all_results.extend(results)

        return all_results

    def _calculate_molecular_weight_from_composition(self, composition: Dict[str, int]) -> float:
        """从原子组成计算分子量"""
        weight = 0.0
        for atom, count in composition.items():
            weight += self.database.atomic_weights.get(atom, 0) * count
        return round(weight, 4)
    
    # 新增：同分异构体识别相关的搜索方法
    def _search_by_ecfp_similarity(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """基于ECFP指纹的相似性搜索"""
        results = []
        query_smiles = query_data.get('smiles', '')
        
        if not query_smiles or not self.isomer_identifier:
            return results
        
        # 生成查询分子的ECFP指纹
        query_fp = self.isomer_identifier.ecfp_generator.generate_structural_fingerprint(query_smiles)
        
        # 与数据库中的所有氨基酸进行比较
        all_amino_acids = self.database.get_all_amino_acids()
        
        for record in all_amino_acids:
            if not record.smiles:
                continue
                
            try:
                # 生成目标分子的ECFP指纹
                target_fp = self.isomer_identifier.ecfp_generator.generate_structural_fingerprint(record.smiles)
                
                # 计算相似性
                similarity = self.isomer_identifier.ecfp_generator.calculate_structural_similarity(query_fp, target_fp)
                
                if similarity > 0.7:  # 相似性阈值
                    confidence = self.confidence_weights.get('ecfp_similarity', 0.9) * similarity
                    
                    results.append(SearchResult(
                        amino_acid_id=record.id,
                        match_method='ecfp_similarity',
                        confidence_score=confidence,
                        amino_acid_record=record,
                        additional_info={
                            'ecfp_similarity': similarity,
                            'query_smiles': query_smiles,
                            'target_smiles': record.smiles
                        }
                    ))
                    
            except Exception as e:
                print(f"ECFP比较失败 ({record.id}): {e}")
                continue
        
        return results
    
    def _search_by_structural_similarity(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """基于结构相似性的搜索"""
        results = []
        query_smiles = query_data.get('smiles', '')
        
        if not query_smiles or not self.isomer_identifier:
            return results
        
        all_amino_acids = self.database.get_all_amino_acids()
        
        for record in all_amino_acids:
            if not record.smiles:
                continue
                
            try:
                # 分析同分异构体关系
                isomer_result = self.isomer_identifier.analyze_isomers(query_smiles, record.smiles)
                
                # 如果是同分异构体或结构相似，加入结果
                if (isomer_result.is_isomer or 
                    isomer_result.similarity_score > 0.8 or 
                    isomer_result.isomer_type in ['structural', 'stereoisomer']):
                    
                    confidence = self.confidence_weights.get('structural_similarity', 0.88) * isomer_result.similarity_score
                    
                    results.append(SearchResult(
                        amino_acid_id=record.id,
                        match_method='structural_similarity',
                        confidence_score=confidence,
                        amino_acid_record=record,
                        additional_info={
                            'structural_similarity': isomer_result.similarity_score,
                            'isomer_type': isomer_result.isomer_type,
                            'is_isomer': isomer_result.is_isomer,
                            'confidence': isomer_result.confidence,
                            'structural_differences': isomer_result.structural_differences,
                            'query_smiles': query_smiles,
                            'target_smiles': record.smiles
                        }
                    ))
                    
            except Exception as e:
                print(f"结构相似性分析失败 ({record.id}): {e}")
                continue
        
        return results
    
    def _search_isomer_aware(self, query_data: Dict[str, Any]) -> List[SearchResult]:
        """同分异构体感知搜索（综合方法）"""
        results = []
        query_smiles = query_data.get('smiles', '')
        
        if not query_smiles or not self.isomer_identifier:
            return results
        
        # 首先尝试精确匹配
        exact_results = self._search_by_residue_name(query_data)
        if exact_results:
            return exact_results
        
        # 然后尝试ECFP相似性搜索
        ecfp_results = self._search_by_ecfp_similarity(query_data)
        
        # 最后尝试结构相似性搜索
        struct_results = self._search_by_structural_similarity(query_data)
        
        # 合并结果，去重并重新计算置信度
        all_results = ecfp_results + struct_results
        merged_results = self._merge_and_rank_results(all_results)
        
        # 为合并后的结果重新计算置信度
        for result in merged_results:
            # 提升同分异构体感知搜索的置信度
            result.confidence_score = min(1.0, result.confidence_score * 1.1)
            result.match_method = 'isomer_aware'
        
        return merged_results
    
    def detect_isomers_in_database(self) -> Dict[str, List[str]]:
        """检测数据库中的同分异构体组"""
        if not self.isomer_identifier:
            return {}
        
        all_amino_acids = self.database.get_all_amino_acids()
        smiles_list = [(record.id, record.smiles) for record in all_amino_acids if record.smiles]
        
        isomer_groups = defaultdict(list)
        
        for i, (id1, smiles1) in enumerate(smiles_list):
            for j, (id2, smiles2) in enumerate(smiles_list[i+1:], i+1):
                try:
                    result = self.isomer_identifier.analyze_isomers(smiles1, smiles2)
                    if result.is_isomer:
                        # 使用分子式作为分组键
                        record1 = self.database.get_amino_acid(id1)
                        if record1:
                            group_key = record1.molecular_formula
                            isomer_groups[group_key].extend([id1, id2])
                except Exception as e:
                    print(f"检测同分异构体失败 ({id1}, {id2}): {e}")
                    continue
        
        # 去重
        for group_key in isomer_groups:
            isomer_groups[group_key] = list(set(isomer_groups[group_key]))
        
        return dict(isomer_groups)

    def add_amino_acid(self, amino_acid_data: Dict[str, Any]) -> bool:
        """添加新的氨基酸（支持扩展）"""
        try:
            # 创建记录
            record = AminoAcidRecord(
                id=amino_acid_data['id'],
                name=amino_acid_data['name'],
                molecular_formula=amino_acid_data.get('molecular_formula', ''),
                molecular_weight=amino_acid_data.get('molecular_weight', 0.0),
                smiles=amino_acid_data.get('smiles', ''),
                atom_composition=amino_acid_data.get('atom_composition', {}),
                key_features=amino_acid_data.get('key_features', []),
                structure_data=amino_acid_data.get('structure_data'),
                fingerprints=amino_acid_data.get('fingerprints'),
                metadata=amino_acid_data.get('metadata', {})
            )

            # 添加到数据库
            success = self.database.add_amino_acid(record)

            if success:
                # 重建索引
                self.index_manager.rebuild_indices()
                print(f"成功添加氨基酸: {record.id} - {record.name}")

            return success

        except Exception as e:
            print(f"添加氨基酸失败: {e}")
            return False

    def get_system_statistics(self) -> Dict[str, Any]:
        """获取系统统计信息"""
        return {
            'database_stats': self.database.get_statistics(),
            'index_stats': self.index_manager.get_index_statistics(),
            'residue_matcher_stats': self.residue_matcher.get_match_statistics(),
            'fingerprint_stats': self.fingerprint_matcher.get_fingerprint_statistics(),
            'supported_search_methods': list(self.search_strategies.keys())
        }

    def benchmark_search_performance(self, num_queries: int = 100) -> Dict[str, float]:
        """性能基准测试"""
        import random

        all_records = self.database.get_all_amino_acids()
        if not all_records:
            return {'error': 'No amino acids in database'}

        # 生成随机查询
        test_queries = []
        for _ in range(num_queries):
            record = random.choice(all_records)
            test_queries.append({
                'residue_name': record.id,
                'molecular_formula': record.molecular_formula,
                'atom_composition': record.atom_composition
            })

        # 测试各种搜索方法
        performance_results = {}

        for method in self.search_strategies.keys():
            start_time = time.time()

            for query in test_queries:
                self.search(query, methods=[method], max_results=5)

            total_time = time.time() - start_time
            avg_time = (total_time / num_queries) * 1000  # 转换为毫秒
            performance_results[method] = round(avg_time, 3)

        return performance_results

    def analyze_atomic_structure(self, amino_acid_id: str) -> Optional[Dict[str, Any]]:
        """分析氨基酸的原子级结构"""
        if not self.atomic_analyzer:
            return {'error': '原子分析器不可用'}

        # 获取氨基酸记录
        record = self.database.get_amino_acid(amino_acid_id)
        if not record:
            return {'error': f'未找到氨基酸: {amino_acid_id}'}

        try:
            # 使用SMILES进行原子级分析
            analysis = self.atomic_analyzer.analyze_amino_acid(
                record.smiles, record.name
            )

            return {
                'amino_acid_id': amino_acid_id,
                'name': record.name,
                'analysis': analysis,
                'success': True
            }

        except Exception as e:
            return {
                'amino_acid_id': amino_acid_id,
                'error': str(e),
                'success': False
            }

    def get_carbon_connectivity_report(self, amino_acid_id: str) -> str:
        """获取碳原子连接性报告"""
        if not self.atomic_analyzer:
            return "错误: 原子分析器不可用"

        record = self.database.get_amino_acid(amino_acid_id)
        if not record:
            return f"错误: 未找到氨基酸 {amino_acid_id}"

        try:
            report = self.atomic_analyzer.get_carbon_connectivity_report(record.smiles)
            return f"氨基酸: {record.name} ({amino_acid_id})\n{report}"
        except Exception as e:
            return f"错误: 分析失败 - {e}"

    def search_by_carbon_pattern(self, carbon_pattern: Dict[str, Any]) -> List[SearchResult]:
        """根据碳原子模式搜索氨基酸"""
        if not self.atomic_analyzer:
            return []

        results = []
        all_records = self.database.get_all_amino_acids()

        for record in all_records:
            try:
                analysis = self.atomic_analyzer.analyze_molecule(record.smiles)
                carbon_analyses = analysis.get('carbon_analysis', [])

                # 检查是否匹配碳原子模式
                if self._matches_carbon_pattern(carbon_analyses, carbon_pattern):
                    confidence = self._calculate_pattern_confidence(
                        carbon_analyses, carbon_pattern
                    )

                    result = SearchResult(
                        amino_acid_id=record.id,
                        amino_acid_record=record,
                        match_method='carbon_pattern',
                        confidence_score=confidence,
                        additional_info={
                            'carbon_analysis': carbon_analyses,
                            'pattern_match': carbon_pattern
                        }
                    )
                    results.append(result)

            except Exception as e:
                continue

        # 按置信度排序
        results.sort(key=lambda x: x.confidence_score, reverse=True)
        return results

    def _matches_carbon_pattern(self, carbon_analyses: List, pattern: Dict[str, Any]) -> bool:
        """检查碳原子分析是否匹配指定模式"""
        # 检查碳原子数量
        if 'carbon_count' in pattern:
            if len(carbon_analyses) != pattern['carbon_count']:
                return False

        # 检查特定碳原子的连接
        if 'carbon_connections' in pattern:
            for carbon_num, expected_atoms in pattern['carbon_connections'].items():
                carbon_idx = int(carbon_num.replace('C', '')) - 1
                if carbon_idx < len(carbon_analyses):
                    actual_atoms = carbon_analyses[carbon_idx].attached_atoms
                    if not all(atom in actual_atoms for atom in expected_atoms):
                        return False

        # 检查功能基团
        if 'functional_groups' in pattern:
            all_groups = []
            for carbon in carbon_analyses:
                all_groups.extend(carbon.functional_groups)

            for required_group in pattern['functional_groups']:
                if required_group not in all_groups:
                    return False

        return True

    def _calculate_pattern_confidence(self, carbon_analyses: List, pattern: Dict[str, Any]) -> float:
        """计算模式匹配的置信度"""
        confidence = 0.0
        total_checks = 0

        # 碳原子数量匹配
        if 'carbon_count' in pattern:
            total_checks += 1
            if len(carbon_analyses) == pattern['carbon_count']:
                confidence += 1.0

        # 连接性匹配
        if 'carbon_connections' in pattern:
            for carbon_num, expected_atoms in pattern['carbon_connections'].items():
                total_checks += 1
                carbon_idx = int(carbon_num.replace('C', '')) - 1
                if carbon_idx < len(carbon_analyses):
                    actual_atoms = carbon_analyses[carbon_idx].attached_atoms
                    match_ratio = len(set(expected_atoms) & set(actual_atoms)) / len(expected_atoms)
                    confidence += match_ratio

        # 功能基团匹配
        if 'functional_groups' in pattern:
            all_groups = []
            for carbon in carbon_analyses:
                all_groups.extend(carbon.functional_groups)

            for required_group in pattern['functional_groups']:
                total_checks += 1
                if required_group in all_groups:
                    confidence += 1.0

        return confidence / total_checks if total_checks > 0 else 0.0

# 向后兼容接口
class CompatibilityWrapper:
    """向后兼容包装器，保持与原有接口的兼容性"""

    def __init__(self, database_path: str = "amino_acids.db"):
        self.engine = ScalableSearchEngine(database_path)

    def search_pdb_files(self, pdb_files: List[str],
                        methods: List[str] = None) -> List[Dict[str, Any]]:
        """兼容原有的PDB文件搜索接口"""
        if methods is None:
            methods = ['residue_name', 'molecular_formula', 'atom_composition']

        all_results = []

        for pdb_file in pdb_files:
            residues = self._extract_residues_from_pdb(pdb_file)
            pdb_id = os.path.basename(pdb_file).replace('.pdb', '')

            for residue_key, residue_data in residues.items():
                residue_name = residue_key.split('_')[0]
                chain_id = residue_key.split('_')[1] if len(residue_key.split('_')) > 1 else ''
                residue_number = residue_key.split('_')[2] if len(residue_key.split('_')) > 2 else ''

                # 跳过标准氨基酸和核酸
                if self._is_standard_amino_acid(residue_name) or self._is_nucleic_acid(residue_name, residue_data):
                    continue

                query_data = {
                    'residue_name': residue_name,
                    'atom_composition': dict(residue_data['composition']),
                    'molecular_weight': self.engine._calculate_molecular_weight_from_composition(
                        residue_data['composition']
                    )
                }

                search_results = self.engine.search(query_data, methods=methods, max_results=1)

                if search_results:
                    result = search_results[0]
                    # 转换为原有格式
                    all_results.append({
                        'pdb_id': pdb_id,
                        'amino_acid_id': result.amino_acid_id,
                        'residue_name': residue_name,
                        'chain_id': chain_id,
                        'residue_number': residue_number,
                        'match_method': result.match_method,
                        'confidence_score': result.confidence_score,
                        'atom_composition': dict(residue_data['composition']),
                        'additional_info': result.additional_info or {}
                    })

        return all_results

    def _extract_residues_from_pdb(self, pdb_file: str) -> Dict[str, Dict]:
        """从PDB文件提取残基信息（简化版）"""
        residues = defaultdict(lambda: {
            'atoms': [],
            'coordinates': [],
            'composition': defaultdict(int),
            'atom_details': []
        })

        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    residue_name = line[17:20].strip()
                    chain_id = line[21:22].strip()
                    residue_seq = line[22:26].strip()
                    element = line[76:78].strip()

                    if not element:
                        element = self._infer_element_from_atom_name(line[12:16].strip())

                    key = f"{residue_name}_{chain_id}_{residue_seq}"

                    if element:
                        residues[key]['composition'][element] += 1

        return dict(residues)

    def _infer_element_from_atom_name(self, atom_name: str) -> str:
        """从原子名推断元素符号"""
        atom_name = atom_name.strip()
        if atom_name.startswith('C'):
            return 'C'
        elif atom_name.startswith('N'):
            return 'N'
        elif atom_name.startswith('O'):
            return 'O'
        elif atom_name.startswith('S'):
            return 'S'
        elif atom_name.startswith('P'):
            return 'P'
        elif atom_name.startswith('H'):
            return 'H'
        else:
            return atom_name[0] if atom_name else ''

    def _is_standard_amino_acid(self, residue_name: str) -> bool:
        """检查是否为标准氨基酸"""
        standard_aa = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        return residue_name in standard_aa

    def _is_nucleic_acid(self, residue_name: str, residue_data: Dict) -> bool:
        """检查是否为核酸残基"""
        nucleic_acids = {
            'A', 'T', 'G', 'C', 'U', 'DA', 'DT', 'DG', 'DC',
            'RA', 'RT', 'RG', 'RC', 'RU'
        }

        if residue_name in nucleic_acids:
            return True

        # 检查是否含有磷原子
        composition = residue_data.get('composition', {})
        if 'P' in composition:
            return True

        return False

def main():
    """主程序 - 演示新架构的使用"""
    print("=" * 80)
    print("可扩展非天然氨基酸PDB搜索引擎 - 第一阶段")
    print("=" * 80)

    # 初始化搜索引擎
    engine = ScalableSearchEngine()

    # 显示系统统计信息
    stats = engine.get_system_statistics()
    print(f"\n系统统计信息:")
    print(f"  数据库中的氨基酸数量: {stats['database_stats']['total_amino_acids']}")
    print(f"  索引大小: {stats['index_stats']['residue_name_index_size']}")
    print(f"  支持的搜索方法: {', '.join(stats['supported_search_methods'])}")

    # 示例搜索
    print(f"\n示例搜索:")

    # 1. 残基名搜索
    print(f"\n1. 残基名搜索 '0A1':")
    results = engine.search({'residue_name': '0A1'})
    for result in results:
        print(f"   找到: {result.amino_acid_id} - {result.amino_acid_record.name}")
        print(f"   置信度: {result.confidence_score:.3f}")

    # 2. 分子式搜索
    print(f"\n2. 分子式搜索 'C10H13NO3':")
    results = engine.search({'molecular_formula': 'C10H13NO3'})
    for result in results:
        print(f"   找到: {result.amino_acid_id} - {result.amino_acid_record.name}")
        print(f"   置信度: {result.confidence_score:.3f}")

    # 3. 原子组成搜索
    print(f"\n3. 原子组成搜索:")
    atom_comp = {'C': 6, 'H': 11, 'N': 1, 'O': 2}
    results = engine.search({'atom_composition': atom_comp})
    for result in results:
        print(f"   找到: {result.amino_acid_id} - {result.amino_acid_record.name}")
        print(f"   置信度: {result.confidence_score:.3f}")

    # 性能测试
    print(f"\n性能基准测试 (100次查询):")
    perf_results = engine.benchmark_search_performance(100)
    for method, avg_time in perf_results.items():
        print(f"  {method}: {avg_time:.3f}ms 平均")

    # 向后兼容测试
    print(f"\n向后兼容性测试:")
    wrapper = CompatibilityWrapper()

    # 查找test目录中的PDB文件
    test_dir = "test"
    if os.path.exists(test_dir):
        pdb_files = [os.path.join(test_dir, f) for f in os.listdir(test_dir) if f.endswith('.pdb')]
        if pdb_files:
            print(f"  搜索 {len(pdb_files)} 个PDB文件...")
            compat_results = wrapper.search_pdb_files(pdb_files)
            print(f"  找到 {len(compat_results)} 个匹配")

            for result in compat_results:
                print(f"    PDB: {result['pdb_id']} | 氨基酸: {result['amino_acid_id']} | "
                      f"置信度: {result['confidence_score']:.3f}")

    print(f"\n" + "=" * 80)
    print("第一阶段核心架构演示完成！")
    print("=" * 80)

if __name__ == "__main__":
    main()
