#!/usr/bin/env python3
"""
自定义异常类
定义系统中使用的异常类型
"""

class PDBSearchEngineError(Exception):
    """PDB搜索引擎基础异常"""
    pass

class DatabaseError(PDBSearchEngineError):
    """数据库相关异常"""
    pass

class DatabaseConnectionError(DatabaseError):
    """数据库连接异常"""
    pass

class DatabaseCorruptionError(DatabaseError):
    """数据库损坏异常"""
    pass

class SearchError(PDBSearchEngineError):
    """搜索相关异常"""
    pass

class InvalidQueryError(SearchError):
    """无效查询异常"""
    pass

class SearchTimeoutError(SearchError):
    """搜索超时异常"""
    pass

class ChemistryError(PDBSearchEngineError):
    """化学分析相关异常"""
    pass

class InvalidSMILESError(ChemistryError):
    """无效SMILES格式异常"""
    pass

class FingerprintGenerationError(ChemistryError):
    """指纹生成异常"""
    pass

class IsomerAnalysisError(ChemistryError):
    """同分异构体分析异常"""
    pass

class PerformanceError(PDBSearchEngineError):
    """性能相关异常"""
    pass

class MemoryLimitExceededError(PerformanceError):
    """内存限制超出异常"""
    pass

class CacheError(PerformanceError):
    """缓存相关异常"""
    pass

class APIError(PDBSearchEngineError):
    """API相关异常"""
    pass

class ValidationError(APIError):
    """参数验证异常"""
    pass

class RateLimitError(APIError):
    """频率限制异常"""
    pass