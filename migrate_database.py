#!/usr/bin/env python3
"""
数据库迁移脚本
将旧版本的数据库迁移到新的重构版本
"""

import sys
import os
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.core.database import DatabaseManager
from pdb_uachecker.utils.config import Config


def find_legacy_database():
    """查找旧版数据库文件"""
    possible_paths = [
        "core/amino_acids.db",
        "amino_acids.db",
        "data/amino_acids.db",
        "legacy/amino_acids.db"
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    
    return None


def main():
    """主函数"""
    print("🔄 PDB-UAChecker 数据库迁移工具")
    print("=" * 50)
    
    # 查找旧版数据库
    legacy_db_path = find_legacy_database()
    
    if not legacy_db_path:
        print("❌ 未找到旧版数据库文件")
        print("请确保以下路径之一存在数据库文件:")
        print("  - core/amino_acids.db")
        print("  - amino_acids.db")
        print("  - data/amino_acids.db")
        print("  - legacy/amino_acids.db")
        return 1
    
    print(f"📁 找到旧版数据库: {legacy_db_path}")
    
    # 初始化新版数据库管理器
    config = Config()
    db_manager = DatabaseManager(config)
    
    print(f"📁 新版数据库路径: {config.get_database_path()}")
    
    # 执行迁移
    print("🚀 开始数据迁移...")
    
    try:
        success = db_manager.migrate_from_legacy_database(legacy_db_path)
        
        if success:
            print("✅ 数据迁移完成！")
            
            # 显示统计信息
            stats = db_manager.get_database_stats()
            print("\n📊 新数据库统计:")
            for key, value in stats.items():
                print(f"  {key.replace('_', ' ').title()}: {value}")
            
            return 0
        else:
            print("❌ 数据迁移失败")
            return 1
    
    except Exception as e:
        print(f"❌ 迁移过程中出错: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())