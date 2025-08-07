"""
分析模块 - 氨基酸分类和识别
重构后的统一架构
"""

# 🎯 新的统一分类系统 - 推荐使用
from .unified_classifier import UnifiedClassifier, ClassificationTier

# 🧪 知识库和分析器组件
from .knowledge import ChemicalDatabase, StandardAminoAcids
from .analyzers import BackboneAnalyzer, StereochemistryAnalyzer
from .molecular_structure_analyzer import MolecularStructureAnalyzer

# 📊 保留主要分析器
from .analyzer import PDBAnalyzer

# ⚠️ Legacy导入 - 现在位于legacy/目录，不推荐直接使用
# 如需使用历史分类器，请明确导入：
# from .legacy.classifier import AminoAcidClassifier
# from .legacy.expert_classifier import ExpertAminoAcidClassifier  

__all__ = [
    # 新架构 - 推荐使用
    'UnifiedClassifier',
    'ClassificationTier',
    
    # 组件
    'ChemicalDatabase', 
    'StandardAminoAcids',
    'MolecularStructureAnalyzer',
    'BackboneAnalyzer', 
    'StereochemistryAnalyzer',
    
    # 分析器
    'PDBAnalyzer'
]