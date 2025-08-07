#!/usr/bin/env python3
"""
测试智能高覆盖率分类器
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.analysis.intelligent_classifier import IntelligentAminoAcidClassifier
from pdb_uachecker.core.database.enhanced_manager import EnhancedDatabaseManager
from pdb_uachecker.utils.config import Config


def test_intelligent_classifier():
    """测试智能分类器"""
    print("🚀 智能高覆盖率分类器测试")
    print("=" * 60)
    
    # 初始化
    config = Config()
    db = EnhancedDatabaseManager(config)
    intelligent_classifier = IntelligentAminoAcidClassifier()
    
    # 获取所有氨基酸
    amino_acids = db.get_all_amino_acids()
    print(f"📊 测试数据: {len(amino_acids)} 个氨基酸")
    
    # 执行批量分类
    batch_result = intelligent_classifier.batch_classify(amino_acids)
    
    print("\n" + "=" * 60)
    print("📈 智能分类结果")
    print("=" * 60)
    
    # 显示统计信息
    stats = batch_result['statistics']
    
    print(f"🎯 分类覆盖率:")
    print(f"  自动化覆盖率: {stats['overall_coverage']['automated_classification']:.1f}%")
    print(f"  高置信度覆盖率: {stats['overall_coverage']['high_confidence_coverage']:.1f}%")
    print(f"  需要人工审查: {stats['tier4_needs_review']['percentage']:.1f}%")
    
    print(f"\n📊 分层分类统计:")
    tier_names = {
        'tier1_definitive': 'Tier 1 - 100%确定',
        'tier2_high_confidence': 'Tier 2 - 高置信度',
        'tier3_intelligent': 'Tier 3 - 智能推断',
        'tier4_needs_review': 'Tier 4 - 需要审查'
    }
    
    for tier_key, tier_name in tier_names.items():
        tier_stats = stats[tier_key]
        print(f"  {tier_name}: {tier_stats['count']} ({tier_stats['percentage']:.1f}%) - {tier_stats['confidence_level']}")
    
    print(f"\n🏷️ 分类分布:")
    for category, count in sorted(batch_result['category_distribution'].items()):
        percentage = (count / len(amino_acids)) * 100
        print(f"  {category}: {count} ({percentage:.1f}%)")
    
    print(f"\n🎯 性能对比:")
    old_coverage = 10.5  # 之前超高精度分类器的覆盖率
    new_coverage = stats['overall_coverage']['automated_classification']
    improvement = new_coverage - old_coverage
    
    print(f"  原覆盖率: {old_coverage}%")
    print(f"  新覆盖率: {new_coverage:.1f}%")
    print(f"  提升幅度: +{improvement:.1f} 个百分点")
    print(f"  人工工作量减少: {(100 - new_coverage) - (100 - old_coverage):.1f}%")
    
    # 显示各层级示例
    print(f"\n📋 分类示例:")
    
    results_by_tier = {}
    for result in batch_result['results']:
        tier = result['classification_tier']
        if tier not in results_by_tier:
            results_by_tier[tier] = []
        results_by_tier[tier].append(result)
    
    for tier_key, tier_name in tier_names.items():
        if tier_key in results_by_tier:
            examples = results_by_tier[tier_key][:3]  # 显示前3个
            print(f"\n  {tier_name}:")
            for example in examples:
                print(f"    {example['amino_acid_id']}: {example['standard_category']} (置信度: {example['confidence']:.2f})")
                print(f"      方法: {example['detection_method']}")
                if len(example['smiles']) > 50:
                    print(f"      SMILES: {example['smiles'][:50]}...")
                else:
                    print(f"      SMILES: {example['smiles']}")
    
    # 显示需要审查的案例
    needs_review = [r for r in batch_result['results'] if r['classification_tier'] == 'tier4_needs_review']
    if needs_review:
        print(f"\n⚠️ 需要专家审查的案例 ({len(needs_review)}个):")
        for i, case in enumerate(needs_review[:5], 1):
            print(f"  {i}. {case['amino_acid_id']}: {case['details']['unclassified_reason']}")
            suggestions = case['details']['suggested_analysis']
            if suggestions:
                print(f"     建议: {suggestions[0]}")
    
    print(f"\n🎉 智能分类器测试完成！")
    print(f"💡 成功实现 {new_coverage:.1f}% 自动化分类覆盖率！")
    
    return batch_result


def random_classification_check():
    """随机抽检分类结果验证准确性"""
    import random
    
    print("\n🔍 随机抽检智能分类结果")
    print("=" * 60)
    
    # 初始化
    config = Config()
    db = EnhancedDatabaseManager(config)
    intelligent_classifier = IntelligentAminoAcidClassifier()
    
    # 获取所有氨基酸
    amino_acids = db.get_all_amino_acids()
    
    # 随机选择12个进行检查
    sample_size = 12
    random_sample = random.sample(amino_acids, sample_size)
    
    print(f"📊 随机抽取 {sample_size} 个氨基酸进行验证:\n")
    
    correct_count = 0
    total_count = 0
    
    for i, amino_acid in enumerate(random_sample, 1):
        print(f"🧪 [{i:2d}] {amino_acid.id} - {amino_acid.name}")
        print(f"     SMILES: {amino_acid.smiles}")
        
        # 运行分类
        classification = intelligent_classifier.classify_amino_acid(amino_acid)
        
        print(f"     🏷️  分类: {classification['standard_category']}")
        print(f"     🎯 置信度: {classification['confidence']:.2f}")
        if 'detection_method' in classification:
            print(f"     🔧 方法: {classification['detection_method']}")
        elif 'method' in classification:
            print(f"     🔧 方法: {classification['method']}")
        if 'classification_tier' in classification:
            print(f"     📈 层级: {classification['classification_tier']}")
        elif 'tier' in classification:
            print(f"     📈 层级: {classification['tier']}")
        
        # 手动验证
        manual_check = ""
        is_correct = False
        
        if "D-amino" in classification['standard_category']:
            if amino_acid.smiles and 'C@H' in amino_acid.smiles and 'C@@H' not in amino_acid.smiles:
                manual_check = "✅ D型手性检测正确"
                is_correct = True
            else:
                manual_check = "❌ 可能误判D型"
        elif "Aromatic" in classification['standard_category']:
            aromatic_patterns = ['c1ccccc1', 'c1cccc', 'c1ccc', 'c1cc', 'c1c[nH]c', 'c1cnc']
            if amino_acid.smiles and any(pattern in amino_acid.smiles for pattern in aromatic_patterns):
                manual_check = "✅ 芳香环检测正确" 
                is_correct = True
            else:
                manual_check = "❓ 需要核实芳香性"
        elif classification['standard_category'] == 'Standard':
            standard_ids = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}
            if amino_acid.id in standard_ids:
                manual_check = "✅ 标准氨基酸正确"
                is_correct = True
            else:
                manual_check = "❓ 非标准氨基酸被分类为标准"
        
        if manual_check:
            print(f"     🔎 验证: {manual_check}")
            if is_correct:
                correct_count += 1
            total_count += 1
        print()
    
    # 统计结果
    if total_count > 0:
        accuracy = (correct_count / total_count) * 100
        print(f"🎯 抽检结果统计:")
        print(f"   验证样本: {total_count} 个")
        print(f"   正确分类: {correct_count} 个") 
        print(f"   准确率: {accuracy:.1f}%")
    
    print("🔍 抽检完成！")


if __name__ == "__main__":
    result = test_intelligent_classifier()
    random_classification_check()