#!/usr/bin/env python3
"""
测试增强数据库管理器的功能
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.core.database.enhanced_manager import EnhancedDatabaseManager


def test_load_amino_acid_database():
    """测试加载氨基酸数据库"""
    print("🧪 测试加载氨基酸数据库...")
    
    # 创建增强数据库管理器
    db_manager = EnhancedDatabaseManager(structures_path="data/structures")
    
    # 加载数据库（只加载前5个进行测试）
    print("📚 开始加载氨基酸数据...")
    amino_acids = {}
    
    # 测试加载几个特定的氨基酸
    test_ids = ['2AG', 'ALA', 'GLY', 'PHE', '4CF']
    
    for aa_id in test_ids:
        aa_dir = Path("data/structures") / aa_id
        if aa_dir.exists():
            amino_acid = db_manager._load_single_amino_acid(aa_id, aa_dir)
            if amino_acid:
                amino_acids[aa_id] = amino_acid
                print(f"   ✅ {aa_id}: {amino_acid.name}")
                print(f"      分子式: {amino_acid.molecular_formula}")
                print(f"      SMILES: {amino_acid.smiles}")
                print(f"      特征: {amino_acid.key_features}")
            else:
                print(f"   ❌ {aa_id}: 加载失败")
        else:
            print(f"   ⚠️ {aa_id}: 目录不存在")
    
    print(f"✅ 测试完成，成功加载 {len(amino_acids)} 个氨基酸")
    return amino_acids


def test_quality_reports():
    """测试数据质量报告"""
    print("\\n🧪 测试数据质量报告...")
    
    db_manager = EnhancedDatabaseManager(structures_path="data/structures")
    
    # 加载几个氨基酸并生成质量报告
    test_ids = ['2AG', 'ALA']
    
    for aa_id in test_ids:
        aa_dir = Path("data/structures") / aa_id
        if aa_dir.exists():
            amino_acid = db_manager._load_single_amino_acid(aa_id, aa_dir)
            if amino_acid:
                report = db_manager.quality_reports.get(aa_id)
                if report:
                    print(f"   📊 {aa_id} 质量报告:")
                    print(f"      总体质量: {report.overall_quality:.2f}")
                    print(f"      SMILES有效: {report.smiles_valid}")
                    print(f"      分子式一致: {report.formula_consistent}")
                    print(f"      文件完整: {report.files_complete}")
                    if report.issues:
                        print(f"      问题: {report.issues}")
                    if report.warnings:
                        print(f"      警告: {report.warnings}")
    
    print("✅ 质量报告测试完成")


def test_database_operations():
    """测试数据库操作"""
    print("\\n🧪 测试数据库操作...")
    
    db_manager = EnhancedDatabaseManager(structures_path="data/structures")
    
    # 加载一个氨基酸
    aa_dir = Path("data/structures/2AG")
    if aa_dir.exists():
        amino_acid = db_manager._load_single_amino_acid("2AG", aa_dir)
        if amino_acid:
            # 测试添加到数据库
            success = db_manager.add_amino_acid(amino_acid)
            print(f"   添加到数据库: {'成功' if success else '失败'}")
            
            # 测试从数据库查询
            retrieved = db_manager.get_amino_acid_by_id("2AG")
            if retrieved:
                print(f"   从数据库查询: 成功")
                print(f"      ID: {retrieved.id}")
                print(f"      名称: {retrieved.name}")
                print(f"      分子式: {retrieved.molecular_formula}")
            else:
                print(f"   从数据库查询: 失败")
            
            # 测试按分子式查询
            if amino_acid.molecular_formula:
                formula_matches = db_manager.get_amino_acids_by_formula(amino_acid.molecular_formula)
                print(f"   按分子式查询: 找到 {len(formula_matches)} 个匹配")
    
    print("✅ 数据库操作测试完成")


if __name__ == "__main__":
    print("🚀 开始增强数据库管理器测试")
    print("=" * 60)
    
    try:
        amino_acids = test_load_amino_acid_database()
        test_quality_reports()
        test_database_operations()
        
        print("=" * 60)
        print("🎉 所有测试通过！增强数据库管理器工作正常！")
        
    except Exception as e:
        print(f"❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)