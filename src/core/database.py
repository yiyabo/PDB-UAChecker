#!/usr/bin/env python3
"""
数据库管理模块
负责氨基酸数据的存储、查询和管理
"""

import os
import json
import sqlite3
import time
from typing import Dict, List, Optional, Any
from collections import defaultdict

from .models import AminoAcidRecord
from .exceptions import DatabaseError, DatabaseConnectionError, DatabaseCorruptionError


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
        try:
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
        except sqlite3.Error as e:
            raise DatabaseConnectionError(f"数据库初始化失败: {e}")
    
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
            
        except sqlite3.Error as e:
            raise DatabaseError(f"添加氨基酸记录失败: {e}")
    
    def get_amino_acid(self, amino_id: str) -> Optional[AminoAcidRecord]:
        """获取氨基酸记录"""
        # 先检查缓存
        if amino_id in self._cache:
            self._update_cache_order(amino_id)
            return self._cache[amino_id]
        
        # 从数据库查询
        try:
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
        except sqlite3.Error as e:
            raise DatabaseError(f"查询氨基酸记录失败: {e}")
        
        return None
    
    def get_all_amino_acids(self) -> List[AminoAcidRecord]:
        """获取所有氨基酸记录"""
        records = []
        try:
            with sqlite3.connect(self.database_path) as conn:
                conn.row_factory = sqlite3.Row
                cursor = conn.execute('SELECT * FROM amino_acids ORDER BY id')
                
                for row in cursor:
                    records.append(self._row_to_record(row))
        except sqlite3.Error as e:
            raise DatabaseError(f"查询所有氨基酸记录失败: {e}")
        
        return records
    
    def get_amino_acid_count(self) -> int:
        """获取氨基酸数量"""
        try:
            with sqlite3.connect(self.database_path) as conn:
                cursor = conn.execute('SELECT COUNT(*) FROM amino_acids')
                return cursor.fetchone()[0]
        except sqlite3.Error as e:
            raise DatabaseError(f"查询氨基酸数量失败: {e}")
    
    def search_by_formula(self, formula: str) -> List[AminoAcidRecord]:
        """按分子式搜索"""
        records = []
        try:
            with sqlite3.connect(self.database_path) as conn:
                conn.row_factory = sqlite3.Row
                cursor = conn.execute(
                    'SELECT * FROM amino_acids WHERE molecular_formula = ?', (formula,)
                )
                
                for row in cursor:
                    records.append(self._row_to_record(row))
        except sqlite3.Error as e:
            raise DatabaseError(f"按分子式搜索失败: {e}")
        
        return records
    
    def search_by_weight_range(self, min_weight: float, max_weight: float) -> List[AminoAcidRecord]:
        """按分子量范围搜索"""
        records = []
        try:
            with sqlite3.connect(self.database_path) as conn:
                conn.row_factory = sqlite3.Row
                cursor = conn.execute(
                    'SELECT * FROM amino_acids WHERE molecular_weight BETWEEN ? AND ?',
                    (min_weight, max_weight)
                )
                
                for row in cursor:
                    records.append(self._row_to_record(row))
        except sqlite3.Error as e:
            raise DatabaseError(f"按分子量范围搜索失败: {e}")
        
        return records
    
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
    
    def clear_cache(self):
        """清空缓存"""
        self._cache.clear()
        self._cache_order.clear()
    
    def get_cache_stats(self) -> Dict[str, Any]:
        """获取缓存统计信息"""
        return {
            'cache_size': len(self._cache),
            'max_cache_size': self.cache_size,
            'cache_usage': len(self._cache) / self.cache_size if self.cache_size > 0 else 0
        }