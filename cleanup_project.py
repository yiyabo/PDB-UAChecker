#!/usr/bin/env python3
"""
项目清理脚本
清理冗余文件，统一项目结构
"""

import os
import shutil
from pathlib import Path

def cleanup_project():
    """清理项目"""
    
    print("🧹 开始项目清理...")
    print("=" * 50)
    
    # 1. 清理冗余的分析器文件
    cleanup_redundant_analyzers()
    
    # 2. 清理重复的数据库文件
    cleanup_duplicate_databases()
    
    # 3. 清理缓存文件
    cleanup_cache_files()
    
    # 4. 处理src目录
    handle_src_directory()
    
    # 5. 移动工具文件
    organize_tool_files()
    
    # 6. 更新项目配置
    update_project_config()
    
    print("\n🎉 项目清理完成！")

def cleanup_redundant_analyzers():
    """清理冗余的分析器文件"""
    
    print("\n📁 清理冗余分析器文件:")
    
    # 保留的核心文件
    keep_files = {
        'analyze_pdb.py',           # 主入口（最新优化版）
        'optimized_pdb_analyzer.py' # 最新分析器
    }
    
    # 需要删除的冗余文件
    redundant_files = [
        'enhanced_pdb_analyzer.py',    # 被optimized版本替代
        'pdb_analyzer_with_3d.py',     # 功能已集成到optimized版本
        'add_fingerprints.py',         # 功能已完成
    ]
    
    for file_name in redundant_files:
        if os.path.exists(file_name):
            print(f"   🗑️  删除冗余文件: {file_name}")
            # 先移动到legacy目录作为备份
            legacy_path = f"legacy/{file_name}"
            if not os.path.exists(legacy_path):
                shutil.move(file_name, legacy_path)
                print(f"      → 已备份到 legacy/{file_name}")
            else:
                os.remove(file_name)
                print(f"      → 已删除（legacy中已有备份）")

def cleanup_duplicate_databases():
    """清理重复的数据库文件"""
    
    print("\n💾 清理重复数据库文件:")
    
    # 检查根目录的数据库文件
    root_db = "amino_acids.db"
    core_db = "core/amino_acids.db"
    
    if os.path.exists(root_db) and os.path.exists(core_db):
        # 比较文件大小，保留较新的
        root_size = os.path.getsize(root_db)
        core_size = os.path.getsize(core_db)
        
        print(f"   📊 根目录数据库: {root_size} bytes")
        print(f"   📊 core目录数据库: {core_size} bytes")
        
        if root_size == core_size:
            print(f"   🗑️  删除根目录重复数据库")
            os.remove(root_db)
        else:
            print(f"   ⚠️  数据库文件大小不同，请手动检查")
    elif os.path.exists(root_db):
        print(f"   ✅ 仅存在根目录数据库")
    elif os.path.exists(core_db):
        print(f"   ✅ 仅存在core目录数据库")

def cleanup_cache_files():
    """清理缓存文件"""
    
    print("\n🗂️  清理缓存文件:")
    
    # 清理__pycache__目录
    pycache_dirs = []
    for root, dirs, files in os.walk('.'):
        for dir_name in dirs:
            if dir_name == '__pycache__':
                pycache_dirs.append(os.path.join(root, dir_name))
    
    for pycache_dir in pycache_dirs:
        print(f"   🗑️  删除缓存目录: {pycache_dir}")
        shutil.rmtree(pycache_dir)
    
    # 清理.pyc文件
    pyc_files = []
    for root, dirs, files in os.walk('.'):
        for file_name in files:
            if file_name.endswith('.pyc'):
                pyc_files.append(os.path.join(root, file_name))
    
    for pyc_file in pyc_files:
        print(f"   🗑️  删除缓存文件: {pyc_file}")
        os.remove(pyc_file)
    
    print(f"   ✅ 清理了 {len(pycache_dirs)} 个缓存目录和 {len(pyc_files)} 个缓存文件")

def handle_src_directory():
    """处理src目录"""
    
    print("\n📦 处理src目录:")
    
    if os.path.exists('src'):
        print("   ❓ 发现src目录，这是一个旧的包结构")
        print("   💡 建议:")
        print("      1. 如果src中有重要功能，需要手动迁移")
        print("      2. 如果src是旧版本，可以删除")
        print("      3. 当前主要功能已在根目录实现")
        
        # 检查src目录是否有活跃使用的代码
        src_files = []
        for root, dirs, files in os.walk('src'):
            for file_name in files:
                if file_name.endswith('.py') and file_name != '__init__.py':
                    src_files.append(os.path.join(root, file_name))
        
        print(f"   📊 src目录包含 {len(src_files)} 个Python文件")
        
        # 暂时不自动删除，让用户决定
        print("   ⚠️  建议手动检查src目录内容后决定是否保留")

def organize_tool_files():
    """整理工具文件"""
    
    print("\n🔧 整理工具文件:")
    
    # 移动工具文件到tools目录
    tool_files = [
        'add_3d_structures.py',
        'reorganize_project.py'
    ]
    
    # 确保tools目录存在
    os.makedirs('tools', exist_ok=True)
    
    for file_name in tool_files:
        if os.path.exists(file_name):
            dest_path = f"tools/{file_name}"
            if not os.path.exists(dest_path):
                print(f"   📁 移动工具文件: {file_name} → tools/")
                shutil.move(file_name, dest_path)
            else:
                print(f"   ✅ 工具文件已存在: tools/{file_name}")

def update_project_config():
    """更新项目配置"""
    
    print("\n⚙️  更新项目配置:")
    
    # 更新setup.py中的入口点
    setup_file = "setup.py"
    if os.path.exists(setup_file):
        print("   📝 需要手动更新setup.py中的入口点:")
        print("      将 'pdb-uachecker=src.api.cli:main' 改为")
        print("      'pdb-uachecker=analyze_pdb:main'")
    
    # 创建.gitignore文件
    gitignore_content = """# Python缓存文件
__pycache__/
*.py[cod]
*$py.class

# 分发/打包
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

# 环境变量
.env
.venv
env/
venv/
ENV/
env.bak/
venv.bak/

# IDE
.vscode/
.idea/
*.swp
*.swo

# 项目特定
cache/
*.cache
*.log
test_results/
"""
    
    with open('.gitignore', 'w') as f:
        f.write(gitignore_content)
    print("   ✅ 创建/更新 .gitignore 文件")

def generate_cleanup_report():
    """生成清理报告"""
    
    print("\n📊 项目清理报告:")
    print("=" * 50)
    
    # 统计当前文件
    python_files = []
    for root, dirs, files in os.walk('.'):
        # 跳过隐藏目录和缓存目录
        dirs[:] = [d for d in dirs if not d.startswith('.') and d != '__pycache__']
        for file_name in files:
            if file_name.endswith('.py'):
                python_files.append(os.path.join(root, file_name))
    
    print(f"📁 当前Python文件数量: {len(python_files)}")
    
    # 按目录分类
    by_directory = {}
    for file_path in python_files:
        dir_name = os.path.dirname(file_path) or '根目录'
        if dir_name not in by_directory:
            by_directory[dir_name] = []
        by_directory[dir_name].append(os.path.basename(file_path))
    
    for dir_name, files in sorted(by_directory.items()):
        print(f"\n📂 {dir_name}:")
        for file_name in sorted(files):
            print(f"   • {file_name}")
    
    print(f"\n🎯 建议的下一步:")
    print("1. 检查legacy目录中的备份文件")
    print("2. 测试主入口: python analyze_pdb.py")
    print("3. 手动检查src目录是否需要保留")
    print("4. 更新setup.py中的入口点")
    print("5. 提交清理后的代码到git")

def main():
    """主函数"""
    
    print("🧹 项目清理工具")
    print("=" * 30)
    
    # 确认操作
    response = input("是否继续清理项目？这将删除一些冗余文件 (y/N): ").lower().strip()
    
    if response == 'y':
        cleanup_project()
        generate_cleanup_report()
    else:
        print("操作已取消")

if __name__ == "__main__":
    main()
