"""
简化的数据库管理器
专注于非天然氨基酸数据的高效查询
"""

import sqlite3
import json
from typing import List, Dict, Optional
from pathlib import Path
import logging

from .models_simple import AminoAcidInfo


class NNADatabase:
    """
    非天然氨基酸数据库
    简化版本，专注于核心功能
    """
    
    def __init__(self, db_path: Optional[str] = None):
        """
        初始化数据库
        
        Args:
            db_path: 数据库路径，None则使用默认路径
        """
        if db_path is None:
            db_path = Path(__file__).parent.parent.parent / "data" / "non_natural_amino_acids.db"
        
        self.db_path = Path(db_path)
        self.connection = None
        
        # 确保数据库存在
        self._ensure_database_exists()
    
    def _ensure_database_exists(self):
        """确保数据库文件存在"""
        if not self.db_path.exists():
            self.db_path.parent.mkdir(parents=True, exist_ok=True)
            self._create_database()
    
    def _create_database(self):
        """创建数据库表结构"""
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                
                # 创建氨基酸表
                cursor.execute('''
                CREATE TABLE IF NOT EXISTS amino_acids (
                    id TEXT PRIMARY KEY,
                    name TEXT NOT NULL,
                    molecular_formula TEXT NOT NULL,
                    molecular_weight REAL NOT NULL,
                    smiles TEXT NOT NULL,
                    atom_composition TEXT NOT NULL,
                    fingerprints TEXT,
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
                ''')
                
                # 创建索引
                cursor.execute('CREATE INDEX IF NOT EXISTS idx_molecular_formula ON amino_acids(molecular_formula)')
                cursor.execute('CREATE INDEX IF NOT EXISTS idx_molecular_weight ON amino_acids(molecular_weight)')
                cursor.execute('CREATE INDEX IF NOT EXISTS idx_smiles ON amino_acids(smiles)')
                
                conn.commit()
                logging.info(f"数据库创建完成: {self.db_path}")
                
        except Exception as e:
            logging.error(f"数据库创建失败: {e}")
            raise
    
    def connect(self) -> sqlite3.Connection:
        """获取数据库连接"""
        if self.connection is None:
            self.connection = sqlite3.connect(self.db_path)
            self.connection.row_factory = sqlite3.Row  # 启用字典式访问
        return self.connection
    
    def close(self):
        """关闭数据库连接"""
        if self.connection:
            self.connection.close()
            self.connection = None
    
    def add_amino_acid(self, amino_acid: AminoAcidInfo) -> bool:
        """
        添加氨基酸到数据库
        
        Args:
            amino_acid: 氨基酸信息
            
        Returns:
            是否添加成功
        """
        try:
            conn = self.connect()
            cursor = conn.cursor()
            
            cursor.execute('''
            INSERT OR REPLACE INTO amino_acids 
            (id, name, molecular_formula, molecular_weight, smiles, atom_composition, fingerprints)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ''', (
                amino_acid.id,
                amino_acid.name,
                amino_acid.molecular_formula,
                amino_acid.molecular_weight,
                amino_acid.smiles,
                json.dumps(amino_acid.atom_composition),
                json.dumps({})  # 指纹稍后计算
            ))
            
            conn.commit()
            return True
            
        except Exception as e:
            logging.error(f"添加氨基酸失败 {amino_acid.id}: {e}")
            return False
    
    def get_amino_acid_by_id(self, amino_acid_id: str) -> Optional[AminoAcidInfo]:
        """根据ID获取氨基酸"""
        try:
            conn = self.connect()
            cursor = conn.cursor()
            
            cursor.execute('SELECT * FROM amino_acids WHERE id = ?', (amino_acid_id,))
            row = cursor.fetchone()
            
            if row:
                return self._row_to_amino_acid(row)
            return None
            
        except Exception as e:
            logging.error(f"获取氨基酸失败 {amino_acid_id}: {e}")
            return None
    
    def get_all_amino_acids(self) -> List[AminoAcidInfo]:
        """获取所有氨基酸"""
        try:
            conn = self.connect()
            cursor = conn.cursor()
            
            cursor.execute('SELECT * FROM amino_acids ORDER BY name')
            rows = cursor.fetchall()
            
            return [self._row_to_amino_acid(row) for row in rows]
            
        except Exception as e:
            logging.error(f"获取所有氨基酸失败: {e}")
            return []
    
    def search_by_molecular_formula(self, formula: str) -> List[AminoAcidInfo]:
        """根据分子式搜索"""
        try:
            conn = self.connect()
            cursor = conn.cursor()
            
            cursor.execute('SELECT * FROM amino_acids WHERE molecular_formula = ?', (formula,))
            rows = cursor.fetchall()
            
            return [self._row_to_amino_acid(row) for row in rows]
            
        except Exception as e:
            logging.error(f"分子式搜索失败 {formula}: {e}")
            return []
    
    def search_by_molecular_weight(self, weight: float, tolerance: float = 0.1) -> List[AminoAcidInfo]:
        """根据分子量搜索"""
        try:
            conn = self.connect()
            cursor = conn.cursor()
            
            min_weight = weight - tolerance
            max_weight = weight + tolerance
            
            cursor.execute('''
            SELECT * FROM amino_acids 
            WHERE molecular_weight BETWEEN ? AND ?
            ORDER BY ABS(molecular_weight - ?)
            ''', (min_weight, max_weight, weight))
            
            rows = cursor.fetchall()
            return [self._row_to_amino_acid(row) for row in rows]
            
        except Exception as e:
            logging.error(f"分子量搜索失败 {weight}: {e}")
            return []
    
    def search_by_atom_composition(self, composition: Dict[str, int], exact: bool = True) -> List[AminoAcidInfo]:
        """根据原子组成搜索"""
        try:
            all_amino_acids = self.get_all_amino_acids()
            matches = []
            
            for aa in all_amino_acids:
                if exact:
                    # 精确匹配
                    if aa.atom_composition == composition:
                        matches.append(aa)
                else:
                    # 包含匹配（候选分子包含所有指定元素）
                    if all(aa.atom_composition.get(elem, 0) >= count 
                           for elem, count in composition.items()):
                        matches.append(aa)
            
            return matches
            
        except Exception as e:
            logging.error(f"原子组成搜索失败: {e}")
            return []
    
    def get_candidates_for_residue(self, molecular_formula: str, 
                                  molecular_weight: Optional[float] = None,
                                  weight_tolerance: float = 1.0) -> List[AminoAcidInfo]:
        """
        为残基获取候选氨基酸
        
        Args:
            molecular_formula: 分子式
            molecular_weight: 分子量（可选）
            weight_tolerance: 分子量容差
            
        Returns:
            候选氨基酸列表
        """
        candidates = []
        
        # 首先按分子式精确匹配
        exact_matches = self.search_by_molecular_formula(molecular_formula)
        candidates.extend(exact_matches)
        
        # 如果提供了分子量，按分子量搜索作为补充
        if molecular_weight is not None and not exact_matches:
            weight_matches = self.search_by_molecular_weight(molecular_weight, weight_tolerance)
            candidates.extend(weight_matches)
        
        # 去重
        seen_ids = set()
        unique_candidates = []
        for candidate in candidates:
            if candidate.id not in seen_ids:
                unique_candidates.append(candidate)
                seen_ids.add(candidate.id)
        
        return unique_candidates
    
    def _row_to_amino_acid(self, row: sqlite3.Row) -> AminoAcidInfo:
        """将数据库行转换为氨基酸对象"""
        return AminoAcidInfo(
            id=row['id'],
            name=row['name'],
            molecular_formula=row['molecular_formula'],
            molecular_weight=row['molecular_weight'],
            smiles=row['smiles'],
            atom_composition=json.loads(row['atom_composition'])
        )
    
    def get_statistics(self) -> Dict[str, any]:
        """获取数据库统计信息"""
        try:
            conn = self.connect()
            cursor = conn.cursor()
            
            # 总数量
            cursor.execute('SELECT COUNT(*) as total FROM amino_acids')
            total = cursor.fetchone()['total']
            
            # 按分子量分布
            cursor.execute('''
            SELECT 
                CASE 
                    WHEN molecular_weight < 200 THEN 'light (<200)'
                    WHEN molecular_weight < 400 THEN 'medium (200-400)'
                    ELSE 'heavy (>400)'
                END as weight_range,
                COUNT(*) as count
            FROM amino_acids
            GROUP BY weight_range
            ''')
            weight_distribution = {row['weight_range']: row['count'] for row in cursor.fetchall()}
            
            return {
                'total_amino_acids': total,
                'weight_distribution': weight_distribution,
                'database_path': str(self.db_path)
            }
            
        except Exception as e:
            logging.error(f"获取统计信息失败: {e}")
            return {}
    
    def batch_add_amino_acids(self, amino_acids: List[AminoAcidInfo]) -> int:
        """批量添加氨基酸"""
        success_count = 0
        
        try:
            conn = self.connect()
            cursor = conn.cursor()
            
            for aa in amino_acids:
                try:
                    cursor.execute('''
                    INSERT OR REPLACE INTO amino_acids 
                    (id, name, molecular_formula, molecular_weight, smiles, atom_composition, fingerprints)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    ''', (
                        aa.id,
                        aa.name,
                        aa.molecular_formula,
                        aa.molecular_weight,
                        aa.smiles,
                        json.dumps(aa.atom_composition),
                        json.dumps({})
                    ))
                    success_count += 1
                    
                except Exception as e:
                    logging.warning(f"添加氨基酸 {aa.id} 失败: {e}")
            
            conn.commit()
            logging.info(f"批量添加完成: {success_count}/{len(amino_acids)} 成功")
            
        except Exception as e:
            logging.error(f"批量添加失败: {e}")
        
        return success_count
    
    def clear_database(self):
        """清空数据库（谨慎使用）"""
        try:
            conn = self.connect()
            cursor = conn.cursor()
            cursor.execute('DELETE FROM amino_acids')
            conn.commit()
            logging.info("数据库已清空")
            
        except Exception as e:
            logging.error(f"清空数据库失败: {e}")


# 数据库单例
_database_instance = None

def get_database(db_path: Optional[str] = None) -> NNADatabase:
    """获取数据库单例"""
    global _database_instance
    if _database_instance is None:
        _database_instance = NNADatabase(db_path)
    return _database_instance