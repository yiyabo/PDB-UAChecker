"""
PDB-UAChecker: 非天然氨基酸识别系统
Advanced Non-Natural Amino Acid Recognition System

基于四重验证策略的高精度氨基酸识别系统
"""

__version__ = "2.0.0"
__author__ = "Xinxiang wang"

from .analysis.analyzer import PDBAnalyzer
from .core.verification.engine import VerificationEngine
from .utils.config import Config

__all__ = ["PDBAnalyzer", "VerificationEngine", "Config"]