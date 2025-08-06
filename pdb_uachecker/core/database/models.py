"""
数据库模型
定义数据库相关的数据结构和操作
"""

from typing import Dict, List, Optional
from dataclasses import dataclass

from ..models import AminoAcidInfo


@dataclass
class AminoAcidDB:
    """氨基酸数据库模型"""
    amino_acids: Dict[str, AminoAcidInfo]
    fingerprints: Dict[str, Dict[str, str]]
    
    def __init__(self):
        self.amino_acids = {}
        self.fingerprints = {}
    
    def add_amino_acid(self, amino_acid: AminoAcidInfo):
        """添加氨基酸"""
        self.amino_acids[amino_acid.id] = amino_acid
        if amino_acid.fingerprints:
            self.fingerprints[amino_acid.id] = amino_acid.fingerprints
    
    def get_amino_acid(self, amino_acid_id: str) -> Optional[AminoAcidInfo]:
        """获取氨基酸"""
        return self.amino_acids.get(amino_acid_id)
    
    def get_amino_acids_by_formula(self, formula: str) -> List[AminoAcidInfo]:
        """根据分子式获取氨基酸列表"""
        return [
            aa for aa in self.amino_acids.values()
            if aa.molecular_formula == formula
        ]
    
    def get_all_amino_acids(self) -> List[AminoAcidInfo]:
        """获取所有氨基酸"""
        return list(self.amino_acids.values())
    
    def get_fingerprint(self, amino_acid_id: str, fp_type: str = 'ecfp2') -> Optional[str]:
        """获取指纹"""
        if amino_acid_id in self.fingerprints:
            return self.fingerprints[amino_acid_id].get(fp_type)
        return None
    
    @property
    def size(self) -> int:
        """数据库大小"""
        return len(self.amino_acids)