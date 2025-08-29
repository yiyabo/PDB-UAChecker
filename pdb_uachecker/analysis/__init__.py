"""
分析模块 - 氨基酸分类和识别
增强统一架构
"""

# 🎯 增强统一分类系统 - 推荐使用
try:
    from .enhanced_unified_classifier import EnhancedUnifiedClassifier, ClassificationTier
    from .classification_validator import ClassificationValidator, ValidationLevel
except ImportError:
    pass  # 允许部分导入失败

# 🧪 知识库和分析器组件
try:
    from .molecular_structure_analyzer import MolecularStructureAnalyzer
    from .chemical_knowledge_base import ChemicalDatabase, StandardAminoAcids
except ImportError:
    pass

# 📊 保留主要分析器
try:
    from .analyzer import PDBAnalyzer
except ImportError:
    pass

__all__ = [
    # 新架构 - 推荐使用
    'EnhancedUnifiedClassifier',
    'ClassificationValidator', 
    'ValidationLevel',
    'ClassificationTier',
    
    # 组件
    'MolecularStructureAnalyzer',
    'ChemicalDatabase', 
    'StandardAminoAcids',
    
    # 分析器
    'PDBAnalyzer'
]