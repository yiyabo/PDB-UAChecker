#!/usr/bin/env python3
"""
加载完整的氨基酸数据库
从data/structures目录扫描并加载所有氨基酸数据
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.core.database.enhanced_manager import EnhancedDatabaseManager


def main():
    """主函数"""
    print("🚀 开始加载完整的氨基酸数据库")
    print("=" * 60)
    
    try:
        # 创建增强数据库管理器
        db_manager = EnhancedDatabaseManager(structures_path="data/structures")
        
        # 加载所有氨基酸数据
        print("📚 扫描并加载所有氨基酸数据...")
        amino_acids = db_manager.load_amino_acid_database()
        
        # 构建搜索索引
        db_manager.build_search_indices()
        
        # 获取统计信息
        stats = db_manager.get_enhanced_database_stats()
        
        print("\\n📊 数据库统计信息:")
        print(f"   总氨基酸数: {stats['total_amino_acids']}")
        print(f"   有SMILES的: {stats['with_smiles']}")
        print(f"   有指纹的: {stats['with_fingerprints']}")
        print(f"   高质量数据: {stats['high_quality']}")
        print(f"   中等质量数据: {stats['medium_quality']}")
        print(f"   低质量数据: {stats['low_quality']}")
        print(f"   有问题的数据: {stats['with_issues']}")
        print(f"   有警告的数据: {stats['with_warnings']}")
        
        # 显示质量报告摘要
        quality_reports = db_manager.get_quality_report()
        
        print("\\n⚠️ 数据质量问题摘要:")
        issues_count = 0
        warnings_count = 0
        
        for aa_id, report in quality_reports.items():
            if report.issues:
                issues_count += len(report.issues)
                print(f"   ❌ {aa_id}: {', '.join(report.issues)}")
            elif report.warnings:
                warnings_count += len(report.warnings)
                if len(report.warnings) <= 2:  # 只显示主要警告
                    print(f"   ⚠️ {aa_id}: {', '.join(report.warnings[:2])}")
        
        print(f"\\n📈 质量摘要: {issues_count}个问题, {warnings_count}个警告")
        
        # 测试几个查询
        print("\\n🔍 测试数据库查询功能:")
        
        # 测试按ID查询
        test_aa = db_manager.get_amino_acid_by_id("2AG")
        if test_aa:
            print(f"   ✅ 按ID查询: {test_aa.id} - {test_aa.name}")
        
        # 测试按分子式查询
        formula_matches = db_manager.get_amino_acids_by_formula("C5H9NO2")
        print(f"   ✅ 按分子式查询(C5H9NO2): 找到 {len(formula_matches)} 个匹配")
        
        # 测试获取所有氨基酸
        all_aas = db_manager.get_all_amino_acids()
        print(f"   ✅ 获取所有氨基酸: {len(all_aas)} 个")
        
        print("=" * 60)
        print("🎉 数据库加载完成！")
        print(f"✅ 成功加载 {len(amino_acids)} 个氨基酸到数据库")
        
    except Exception as e:
        print(f"❌ 数据库加载失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()