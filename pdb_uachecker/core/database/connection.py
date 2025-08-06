"""
数据库连接管理
"""

import sqlite3
import json
import os
from pathlib import Path
from typing import List, Dict, Optional, Iterator
from contextlib import contextmanager

from ..models import AminoAcidInfo, DatabaseError
from ...utils.config import Config


class DatabaseManager:
    """数据库管理器"""
    
    def __init__(self, config: Optional[Config] = None):
        if config is None:
            from ...utils.config import default_config
            config = default_config
        
        self.config = config
        self.db_path = config.get_database_path()
        self._ensure_database_exists()
    
    def _ensure_database_exists(self):
        """确保数据库文件存在"""
        if not os.path.exists(self.db_path):
            if self.config.database.auto_create:
                self._create_database()
            else:
                raise DatabaseError(f"数据库文件不存在: {self.db_path}")
    
    def _create_database(self):
        """创建数据库"""
        # 确保目录存在
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            # 创建氨基酸表
            cursor.execute('''
                CREATE TABLE IF NOT EXISTS amino_acids (
                    id TEXT PRIMARY KEY,
                    name TEXT NOT NULL,
                    molecular_formula TEXT,
                    molecular_weight REAL,
                    smiles TEXT,
                    atom_composition TEXT,
                    key_features TEXT,
                    fingerprints TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            ''')
            
            # 创建索引
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_molecular_formula ON amino_acids(molecular_formula)')
            cursor.execute('CREATE INDEX IF NOT EXISTS idx_name ON amino_acids(name)')
            
            conn.commit()
            print(f"✅ 数据库已创建: {self.db_path}")
    
    @contextmanager
    def get_connection(self) -> Iterator[sqlite3.Connection]:
        """获取数据库连接（上下文管理器）"""
        conn = None
        try:
            conn = sqlite3.connect(self.db_path)
            conn.row_factory = sqlite3.Row  # 使结果可以按列名访问
            yield conn
        except Exception as e:
            if conn:
                conn.rollback()
            raise DatabaseError(f"数据库操作失败: {e}")
        finally:
            if conn:
                conn.close()
    
    def get_amino_acid_by_id(self, amino_acid_id: str) -> Optional[AminoAcidInfo]:
        """根据ID获取氨基酸信息"""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT id, name, molecular_formula, molecular_weight, smiles,
                       atom_composition, key_features, fingerprints
                FROM amino_acids WHERE id = ?
            ''', (amino_acid_id,))
            
            row = cursor.fetchone()
            if row:
                return self._row_to_amino_acid_info(row)
            return None
    
    def get_amino_acids_by_formula(self, molecular_formula: str) -> List[AminoAcidInfo]:
        """根据分子式获取氨基酸列表"""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT id, name, molecular_formula, molecular_weight, smiles,
                       atom_composition, key_features, fingerprints
                FROM amino_acids WHERE molecular_formula = ?
            ''', (molecular_formula,))
            
            return [self._row_to_amino_acid_info(row) for row in cursor.fetchall()]
    
    def get_all_amino_acids(self) -> List[AminoAcidInfo]:
        """获取所有氨基酸信息"""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute('''
                SELECT id, name, molecular_formula, molecular_weight, smiles,
                       atom_composition, key_features, fingerprints
                FROM amino_acids ORDER BY id
            ''')
            
            return [self._row_to_amino_acid_info(row) for row in cursor.fetchall()]
    
    def search_amino_acids(self, **criteria) -> List[AminoAcidInfo]:
        """根据条件搜索氨基酸"""
        conditions = []
        params = []
        
        for key, value in criteria.items():
            if key in ['id', 'name', 'molecular_formula', 'smiles']:
                conditions.append(f"{key} = ?")
                params.append(value)
            elif key == 'molecular_weight_range':
                min_weight, max_weight = value
                conditions.append("molecular_weight BETWEEN ? AND ?")
                params.extend([min_weight, max_weight])
        
        if not conditions:
            return self.get_all_amino_acids()
        
        query = f'''
            SELECT id, name, molecular_formula, molecular_weight, smiles,
                   atom_composition, key_features, fingerprints
            FROM amino_acids WHERE {' AND '.join(conditions)}
            ORDER BY id
        '''
        
        with self.get_connection() as conn:
            cursor = conn.cursor()
            cursor.execute(query, params)
            return [self._row_to_amino_acid_info(row) for row in cursor.fetchall()]
    
    def add_amino_acid(self, amino_acid: AminoAcidInfo) -> bool:
        """添加氨基酸信息"""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute('''
                    INSERT OR REPLACE INTO amino_acids 
                    (id, name, molecular_formula, molecular_weight, smiles,
                     atom_composition, key_features, fingerprints)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    amino_acid.id,
                    amino_acid.name,
                    amino_acid.molecular_formula,
                    amino_acid.molecular_weight,
                    amino_acid.smiles,
                    json.dumps(amino_acid.atom_composition),
                    json.dumps(amino_acid.key_features),
                    json.dumps(amino_acid.fingerprints)
                ))
                conn.commit()
                return True
        except Exception as e:
            print(f"⚠️ 添加氨基酸失败 {amino_acid.id}: {e}")
            return False
    
    def update_amino_acid(self, amino_acid: AminoAcidInfo) -> bool:
        """更新氨基酸信息"""
        return self.add_amino_acid(amino_acid)  # INSERT OR REPLACE
    
    def delete_amino_acid(self, amino_acid_id: str) -> bool:
        """删除氨基酸信息"""
        try:
            with self.get_connection() as conn:
                cursor = conn.cursor()
                cursor.execute('DELETE FROM amino_acids WHERE id = ?', (amino_acid_id,))
                conn.commit()
                return cursor.rowcount > 0
        except Exception as e:
            print(f"⚠️ 删除氨基酸失败 {amino_acid_id}: {e}")
            return False
    
    def get_database_stats(self) -> Dict[str, int]:
        """获取数据库统计信息"""
        with self.get_connection() as conn:
            cursor = conn.cursor()
            
            # 总氨基酸数量
            cursor.execute('SELECT COUNT(*) FROM amino_acids')
            total_count = cursor.fetchone()[0]
            
            # 有SMILES的数量
            cursor.execute('SELECT COUNT(*) FROM amino_acids WHERE smiles IS NOT NULL AND smiles != ""')
            smiles_count = cursor.fetchone()[0]
            
            # 有指纹的数量
            cursor.execute('SELECT COUNT(*) FROM amino_acids WHERE fingerprints IS NOT NULL AND fingerprints != "{}"')
            fingerprint_count = cursor.fetchone()[0]
            
            return {
                'total_amino_acids': total_count,
                'with_smiles': smiles_count,
                'with_fingerprints': fingerprint_count
            }
    
    def _row_to_amino_acid_info(self, row: sqlite3.Row) -> AminoAcidInfo:
        """将数据库行转换为AminoAcidInfo对象"""
        return AminoAcidInfo(
            id=row['id'],
            name=row['name'],
            molecular_formula=row['molecular_formula'] or '',
            molecular_weight=row['molecular_weight'] or 0.0,
            smiles=row['smiles'] or '',
            atom_composition=json.loads(row['atom_composition']) if row['atom_composition'] else {},
            key_features=json.loads(row['key_features']) if row['key_features'] else [],
            fingerprints=json.loads(row['fingerprints']) if row['fingerprints'] else {}
        )
    
    def migrate_from_legacy_database(self, legacy_db_path: str) -> bool:
        """从旧版数据库迁移数据"""
        if not os.path.exists(legacy_db_path):
            print(f"⚠️ 旧版数据库不存在: {legacy_db_path}")
            return False
        
        try:
            # 连接旧版数据库
            legacy_conn = sqlite3.connect(legacy_db_path)
            legacy_conn.row_factory = sqlite3.Row
            legacy_cursor = legacy_conn.cursor()
            
            # 查询旧版数据
            legacy_cursor.execute('''
                SELECT id, name, molecular_formula, molecular_weight, smiles,
                       atom_composition, key_features, fingerprints
                FROM amino_acids
            ''')
            
            migrated_count = 0
            for row in legacy_cursor.fetchall():
                amino_acid = self._row_to_amino_acid_info(row)
                if self.add_amino_acid(amino_acid):
                    migrated_count += 1
            
            legacy_conn.close()
            print(f"✅ 成功迁移 {migrated_count} 个氨基酸数据")
            return True
            
        except Exception as e:
            print(f"❌ 数据迁移失败: {e}")
            return False