#!/usr/bin/env python3
"""
项目清理脚本 v2.0
删除重构后不再需要的文件和目录
"""

import os
import shutil
from pathlib import Path


def safe_remove(path):
    """安全删除文件或目录"""
    try:
        if os.path.isfile(path):
            os.remove(path)
            print(f"✅ 删除文件: {path}")
        elif os.path.isdir(path):
            shutil.rmtree(path)
            print(f"✅ 删除目录: {path}")
        else:
            print(f"⚠️ 路径不存在: {path}")
    except Exception as e:
        print(f"❌ 删除失败 {path}: {e}")


def main():
    """主清理函数"""
    print("🧹 PDB-UAChecker v2.0 项目清理")
    print("=" * 50)
    print("⚠️ 这将删除重构后不再需要的文件和目录")
    
    # 确认操作
    response = input("是否继续? (y/N): ").lower().strip()
    if response != 'y':
        print("❌ 操作已取消")
        return
    
    print("\n🚀 开始清理...")
    
    # 1. 删除旧版本的主要文件
    old_main_files = [
        "analyze_pdb.py",  # 旧版主入口
        "optimized_pdb_analyzer.py",  # 旧版优化分析器
        "amino_acid_classifier_improved.py",  # 已集成到新架构
        "improved_algorithms.py",  # 已重构
        "final_classification_report.py",  # 临时脚本
        "cleanup_project.py",  # 旧版清理脚本
    ]
    
    print("\n📁 删除旧版主要文件:")
    for file in old_main_files:
        if os.path.exists(file):
            safe_remove(file)
    
    # 2. 删除legacy目录（保留一个备份说明）
    print("\n📁 处理legacy目录:")
    if os.path.exists("legacy"):
        # 创建备份说明
        with open("LEGACY_BACKUP_INFO.md", "w", encoding="utf-8") as f:
            f.write("# Legacy代码备份信息\n\n")
            f.write("以下文件已在v2.0重构中被替代，如需查看历史版本请查看git历史:\n\n")
            
            for root, dirs, files in os.walk("legacy"):
                for file in files:
                    if file.endswith('.py'):
                        f.write(f"- {os.path.join(root, file)}\n")
        
        safe_remove("legacy")
        print("✅ 创建了 LEGACY_BACKUP_INFO.md 记录备份信息")
    
    # 3. 删除旧版core目录（新版在pdb_uachecker/core）
    print("\n📁 删除旧版core目录:")
    if os.path.exists("core"):
        safe_remove("core")
    
    # 4. 删除src目录（功能已迁移到pdb_uachecker）
    print("\n📁 删除src目录:")
    if os.path.exists("src"):
        safe_remove("src")
    
    # 5. 删除旧版tests目录（新版在pdb_uachecker/tests）
    print("\n📁 删除旧版tests目录:")
    if os.path.exists("tests"):
        safe_remove("tests")
    
    # 6. 删除临时和生成的文件
    temp_files = [
        "amino_acid_classification_improved.csv",
        "amino_acid_classification_improved.json", 
        "classification_improved_summary.txt",
        "classification_report.txt",
        "amino_acids_for_manual_classification.csv",
        "test_protein.pdb",  # 测试文件
        "test_report.txt",   # 测试报告
        "test_results.json", # 测试结果
        "amino_acids.db",    # 根目录的数据库文件（新版在pdb_uachecker/data/）
    ]
    
    print("\n📁 删除临时和生成文件:")
    for file in temp_files:
        if os.path.exists(file):
            safe_remove(file)
    
    # 7. 删除缓存目录
    cache_dirs = [
        "__pycache__",
        "cache",
        ".pytest_cache",
    ]
    
    print("\n📁 删除缓存目录:")
    for dir_name in cache_dirs:
        if os.path.exists(dir_name):
            safe_remove(dir_name)
    
    # 8. 删除IDE和系统文件
    ide_files = [
        ".DS_Store",
        ".idea",  # PyCharm
    ]
    
    print("\n📁 删除IDE和系统文件:")
    for item in ide_files:
        if os.path.exists(item):
            safe_remove(item)
    
    # 9. 整理文档文件
    print("\n📁 整理文档文件:")
    doc_files_to_move = [
        "ALGORITHM_LOGIC.md",
        "Four-fold_Verification_Algorithm_PPT_Content.md", 
        "指纹相似性计算方法说明.md",
        "四重验证算法PPT内容.md",
        "流程图.png"
    ]
    
    # 确保docs目录存在
    os.makedirs("docs", exist_ok=True)
    
    for file in doc_files_to_move:
        if os.path.exists(file):
            try:
                shutil.move(file, f"docs/{file}")
                print(f"✅ 移动文档: {file} -> docs/{file}")
            except Exception as e:
                print(f"❌ 移动失败 {file}: {e}")
    
    # 10. 创建清理后的项目结构说明
    print("\n📝 创建项目结构说明:")
    with open("PROJECT_STRUCTURE.md", "w", encoding="utf-8") as f:
        f.write("# PDB-UAChecker v2.0 项目结构\n\n")
        f.write("## 核心目录\n")
        f.write("- `pdb_uachecker/` - 主要代码包\n")
        f.write("  - `core/` - 核心算法模块\n")
        f.write("  - `analysis/` - 分析模块\n")
        f.write("  - `utils/` - 工具模块\n")
        f.write("  - `api/` - 接口层\n")
        f.write("  - `tests/` - 测试模块\n\n")
        f.write("## 主要文件\n")
        f.write("- `analyze_pdb_v2.py` - 主入口文件\n")
        f.write("- `migrate_database.py` - 数据库迁移工具\n")
        f.write("- `demo_v2.py` - 功能演示脚本\n")
        f.write("- `README_v2.md` - 项目文档\n")
        f.write("- `requirements.txt` - 依赖列表\n")
        f.write("- `setup.py` - 安装配置\n\n")
        f.write("## 配置和数据\n")
        f.write("- `config_example.json` - 配置文件示例\n")
        f.write("- `pdb_uachecker/data/` - 数据库文件\n\n")
        f.write("## 文档\n")
        f.write("- `docs/` - 技术文档和算法说明\n")
        f.write("- `tools/` - 数据处理工具\n")
    
    print("✅ 创建了 PROJECT_STRUCTURE.md")
    
    # 11. 更新.gitignore
    print("\n📝 更新.gitignore:")
    gitignore_content = """# Python
__pycache__/
*.py[cod]
*$py.class
*.so
.Python
build/
develop-eggs/
dist/
downloads/
eggs/
.eggs/
lib/
lib64/
parts/
sdist/
var/
wheels/
*.egg-info/
.installed.cfg
*.egg

# Virtual Environment
venv/
env/
ENV/

# IDE
.vscode/
.idea/
*.swp
*.swo

# OS
.DS_Store
Thumbs.db

# Project specific
*.pdb
*.db
test_*.txt
test_*.json
*_report.txt
*_results.json
cache/
temp/

# Logs
*.log
logs/

# Data files (except examples)
data/*.csv
data/*.json
!data/example_*
"""
    
    with open(".gitignore", "w", encoding="utf-8") as f:
        f.write(gitignore_content)
    print("✅ 更新了 .gitignore")
    
    print("\n🎉 项目清理完成！")
    print("\n📊 清理后的项目结构:")
    print("```")
    print("PDB-UAChecker/")
    print("├── pdb_uachecker/          # 主要代码包")
    print("├── docs/                   # 文档")
    print("├── tools/                  # 工具脚本")
    print("├── analyze_pdb_v2.py       # 主入口")
    print("├── migrate_database.py     # 数据库迁移")
    print("├── demo_v2.py             # 功能演示")
    print("├── README_v2.md           # 项目文档")
    print("├── requirements.txt       # 依赖")
    print("├── setup.py              # 安装配置")
    print("└── config_example.json   # 配置示例")
    print("```")
    
    print("\n💡 下一步建议:")
    print("1. 检查清理结果是否符合预期")
    print("2. 运行测试确保功能正常: python demo_v2.py")
    print("3. 更新git仓库: git add . && git commit -m 'v2.0 项目重构和清理'")
    print("4. 查看 PROJECT_STRUCTURE.md 了解新的项目结构")


if __name__ == "__main__":
    main()