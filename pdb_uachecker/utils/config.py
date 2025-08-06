"""
配置管理模块
统一管理系统配置和参数
"""

import os
from pathlib import Path
from typing import Dict, Any, Optional
from dataclasses import dataclass


@dataclass
class VerificationThresholds:
    """验证阈值配置"""
    molecular_formula: float = 1.0      # 分子式验证阈值
    atom_composition: float = 0.9       # 原子组成验证阈值
    fingerprint_similarity: float = 0.7 # 指纹相似性阈值
    structure_3d: float = 0.5           # 3D结构验证阈值


@dataclass
class DatabaseConfig:
    """数据库配置"""
    path: Optional[str] = None
    auto_create: bool = True
    
    def __post_init__(self):
        if self.path is None:
            # 默认数据库路径
            current_dir = Path(__file__).parent.parent
            self.path = str(current_dir / "data" / "amino_acids.db")


@dataclass
class PerformanceConfig:
    """性能配置"""
    enable_parallel: bool = True
    max_workers: int = 4
    cache_size: int = 1000
    timeout_seconds: int = 30


class Config:
    """系统配置管理器"""
    
    def __init__(self, config_file: Optional[str] = None):
        self.thresholds = VerificationThresholds()
        self.database = DatabaseConfig()
        self.performance = PerformanceConfig()
        
        if config_file and os.path.exists(config_file):
            self.load_from_file(config_file)
    
    def load_from_file(self, config_file: str):
        """从配置文件加载配置"""
        import json
        try:
            with open(config_file, 'r', encoding='utf-8') as f:
                config_data = json.load(f)
            
            # 更新阈值配置
            if 'thresholds' in config_data:
                threshold_data = config_data['thresholds']
                self.thresholds = VerificationThresholds(**threshold_data)
            
            # 更新数据库配置
            if 'database' in config_data:
                db_data = config_data['database']
                self.database = DatabaseConfig(**db_data)
            
            # 更新性能配置
            if 'performance' in config_data:
                perf_data = config_data['performance']
                self.performance = PerformanceConfig(**perf_data)
                
        except Exception as e:
            print(f"⚠️ 配置文件加载失败: {e}")
    
    def save_to_file(self, config_file: str):
        """保存配置到文件"""
        import json
        config_data = {
            'thresholds': {
                'molecular_formula': self.thresholds.molecular_formula,
                'atom_composition': self.thresholds.atom_composition,
                'fingerprint_similarity': self.thresholds.fingerprint_similarity,
                'structure_3d': self.thresholds.structure_3d
            },
            'database': {
                'path': self.database.path,
                'auto_create': self.database.auto_create
            },
            'performance': {
                'enable_parallel': self.performance.enable_parallel,
                'max_workers': self.performance.max_workers,
                'cache_size': self.performance.cache_size,
                'timeout_seconds': self.performance.timeout_seconds
            }
        }
        
        try:
            with open(config_file, 'w', encoding='utf-8') as f:
                json.dump(config_data, f, indent=2, ensure_ascii=False)
        except Exception as e:
            print(f"⚠️ 配置文件保存失败: {e}")
    
    def get_database_path(self) -> str:
        """获取数据库路径"""
        return self.database.path
    
    def is_rdkit_available(self) -> bool:
        """检查RDKit是否可用"""
        try:
            import rdkit
            return True
        except ImportError:
            return False
    
    def validate(self) -> bool:
        """验证配置的有效性"""
        # 检查阈值范围
        if not (0.0 <= self.thresholds.molecular_formula <= 1.0):
            return False
        if not (0.0 <= self.thresholds.atom_composition <= 1.0):
            return False
        if not (0.0 <= self.thresholds.fingerprint_similarity <= 1.0):
            return False
        if not (0.0 <= self.thresholds.structure_3d <= 1.0):
            return False
        
        # 检查数据库路径
        if not self.database.path:
            return False
        
        # 检查性能参数
        if self.performance.max_workers <= 0:
            return False
        if self.performance.cache_size <= 0:
            return False
        if self.performance.timeout_seconds <= 0:
            return False
        
        return True


# 默认配置实例
default_config = Config()