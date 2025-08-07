"""
知识库模块
包含氨基酸化学知识和标准定义
"""

from .chemical_database import ChemicalDatabase
from .amino_acid_registry import StandardAminoAcids

__all__ = [
    'ChemicalDatabase',
    'StandardAminoAcids'
]
