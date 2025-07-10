#!/usr/bin/env python3
"""
核心数据模型
定义氨基酸记录、搜索结果等核心数据结构
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass, asdict

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

@dataclass
class SearchQuery:
    """搜索查询数据类"""
    residue_name: Optional[str] = None
    molecular_formula: Optional[str] = None
    atom_composition: Optional[Dict[str, int]] = None
    molecular_weight: Optional[float] = None
    weight_tolerance: Optional[float] = 5.0
    smiles: Optional[str] = None
    methods: List[str] = None
    max_results: int = 10
    
    def __post_init__(self):
        if self.methods is None:
            self.methods = ['residue_name', 'molecular_formula', 'atom_composition']

@dataclass
class PerformanceConfig:
    """性能配置数据类"""
    cache_size: int = 1000
    enable_lsh: bool = True
    lsh_num_tables: int = 10
    lsh_hash_size: int = 16
    parallel_threads: int = 4
    batch_size: int = 50
    memory_limit_mb: int = 200
    similarity_threshold: float = 0.7