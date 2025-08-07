#!/usr/bin/env python3
"""
增强分类器测试脚本
验证新算法的准确性和性能
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.analysis.enhanced_unified_classifier import EnhancedUnifiedClassifier
from pdb_uachecker.core.models import AminoAcidInfo


def create_test_amino_acids() -> list:
    """创建测试用氨基酸数据"""
    test_cases = [
        # 标准氨基酸测试
        {
            'id': 'ALA',
            'name': 'Alanine', 
            'smiles': 'N[C@@H](C)C(=O)O',
            'expected_categories': ['alpha_amino_acid', 'l_form']
        },
        # D型氨基酸测试
        {
            'id': 'DAL',
            'name': 'D-Alanine',
            'smiles': 'N[C@H](C)C(=O)O', 
            'expected_categories': ['alpha_amino_acid', 'd_form']
        },
        # β-氨基酸测试
        {
            'id': 'BAL',
            'name': 'β-Alanine',
            'smiles': 'NCCC(=O)O',
            'expected_categories': ['beta_amino_acid']
        },
        # N-甲基氨基酸测试
        {
            'id': 'SAR',
            'name': 'Sarcosine',
            'smiles': 'CN[C@@H](C)C(=O)O',
            'expected_categories': ['alpha_amino_acid', 'n_methyl_amino_acid']
        },
        # 芳香族氨基酸测试
        {
            'id': 'PHE',
            'name': 'Phenylalanine',
            'smiles': 'N[C@@H](Cc1ccccc1)C(=O)O',
            'expected_categories': ['alpha_amino_acid', 'aromatic_amino_acid']
        },
        # γ-氨基酸测试
        {
            'id': 'GAB', 
            'name': 'γ-Aminobutyric acid',
            'smiles': 'NCCCC(=O)O',
            'expected_categories': ['gamma_amino_acid']
        },
    ]
    
    amino_acids = []
    for case in test_cases:
        amino_acid = AminoAcidInfo(
            id=case['id'],
            name=case['name'],
            smiles=case['smiles'],
            molecular_formula='',  # 将在分析中计算
            molecular_weight=0.0,  # 将在分析中计算
            atom_composition={}    # 将在分析中计算
        )
        # 添加期望的分类用于验证
        amino_acid.expected_categories = case['expected_categories']
        amino_acids.append(amino_acid)
    
    return amino_acids


def test_enhanced_classifier():
    """测试增强分类器"""
    print("🚀 开始测试增强分类器...")
    print("=" * 60)
    
    # 创建分类器
    classifier = EnhancedUnifiedClassifier()
    
    # 准备测试数据
    test_amino_acids = create_test_amino_acids()
    
    # 执行分类
    print(f"📊 测试数据: {len(test_amino_acids)}个氨基酸")
    print()
    
    results = []
    for i, amino_acid in enumerate(test_amino_acids, 1):
        print(f"🔍 测试 {i}/{len(test_amino_acids)}: {amino_acid.id} ({amino_acid.name})")
        print(f"   SMILES: {amino_acid.smiles}")
        
        try:
            result = classifier.classify_amino_acid(amino_acid)
            results.append(result)
            
            # 显示结果
            print(f"   ✅ 分类成功")
            print(f"   📋 分类结果: {', '.join(result.categories) if result.categories else '无特殊分类'}")
            print(f"   🎯 置信度: {result.confidence:.1%} ({result.classification_tier.value})")
            print(f"   ✔️ 验证通过: {'是' if result.validation_passed else '否'}")
            
            if result.stereochemistry:
                print(f"   🔬 立体化学: {result.stereochemistry}")
            if result.backbone_type:
                print(f"   🦴 骨架类型: {result.backbone_type}")
            if result.is_n_methylated:
                print(f"   🧪 N-甲基化: 是")
            
            # 检查预期结果
            if hasattr(amino_acid, 'expected_categories'):
                expected = set(amino_acid.expected_categories)
                actual = set(result.categories)
                
                if expected.issubset(actual):
                    print(f"   ✅ 预期验证: 通过")
                else:
                    missing = expected - actual
                    print(f"   ⚠️ 预期验证: 部分通过 (缺失: {missing})")
            
            if result.inconsistencies:
                print(f"   ⚠️ 不一致性: {'; '.join(result.inconsistencies[:2])}")
            
            if result.recommendations:
                print(f"   💡 建议: {result.recommendations[0]}")
            
        except Exception as e:
            print(f"   ❌ 分类失败: {e}")
            
        print()
    
    # 生成统计报告
    print("📈 分类统计:")
    print("=" * 60)
    
    stats = classifier.get_classification_statistics()
    print(f"总计分析: {stats['total_classified']}个")
    print(f"验证通过率: {stats['accuracy_metrics']['validation_success_rate']:.1%}")
    print(f"高可信率: {stats['accuracy_metrics']['reliable_rate']:.1%}")
    print()
    
    print("置信度分布:")
    for tier, percentage in stats['confidence_distribution'].items():
        print(f"  {tier}: {percentage}")
    print()
    
    # 生成详细报告
    print("📋 详细报告:")
    print("=" * 60)
    report = classifier.generate_classification_report(results)
    print(report)
    
    return results, stats


def test_individual_analyzers():
    """测试单个分析器"""
    print("\n🔬 测试单个分析器...")
    print("=" * 60)
    
    # 导入单个分析器
    from pdb_uachecker.analysis.analyzers.cip_rule_analyzer import CIPRuleAnalyzer
    from pdb_uachecker.analysis.analyzers.enhanced_backbone_analyzer import EnhancedBackboneAnalyzer
    from pdb_uachecker.analysis.analyzers.n_methylation_analyzer import NMethylationAnalyzer
    
    # 创建分析器
    cip_analyzer = CIPRuleAnalyzer()
    backbone_analyzer = EnhancedBackboneAnalyzer()
    n_methyl_analyzer = NMethylationAnalyzer()
    
    # 测试用例
    test_smiles = {
        'L-Alanine': 'N[C@@H](C)C(=O)O',
        'D-Alanine': 'N[C@H](C)C(=O)O',
        'Sarcosine': 'CN[C@@H](C)C(=O)O',
        'β-Alanine': 'NCCC(=O)O'
    }
    
    for name, smiles in test_smiles.items():
        print(f"\n测试: {name}")
        print(f"SMILES: {smiles}")
        
        # CIP规则分析
        try:
            cip_result = cip_analyzer.analyze_stereochemistry(smiles, name)
            print(f"  CIP规则: {cip_result.get('stereochemistry', 'unknown')} (置信度: {cip_result.get('confidence', 0):.2f})")
        except Exception as e:
            print(f"  CIP规则: 错误 - {e}")
        
        # 骨架分析
        try:
            backbone_result = backbone_analyzer.analyze_backbone_precise(smiles)
            print(f"  骨架类型: {backbone_result.get('backbone_type', 'unknown')} (置信度: {backbone_result.get('confidence', 0):.2f})")
        except Exception as e:
            print(f"  骨架类型: 错误 - {e}")
        
        # N-甲基化分析
        try:
            n_methyl_result = n_methyl_analyzer.analyze_n_methylation(smiles, None, name)
            is_n_methylated = n_methyl_result.get('is_n_methylated', False)
            print(f"  N-甲基化: {'是' if is_n_methylated else '否'} (置信度: {n_methyl_result.get('confidence', 0):.2f})")
        except Exception as e:
            print(f"  N-甲基化: 错误 - {e}")


def main():
    """主函数"""
    print("🌟 增强分类器综合测试")
    print("=" * 80)
    print()
    
    try:
        # 测试增强分类器
        results, stats = test_enhanced_classifier()
        
        # 测试单个分析器
        test_individual_analyzers()
        
        # 导出结果（可选）
        print("\n💾 导出测试结果...")
        classifier = EnhancedUnifiedClassifier()
        
        # JSON导出
        json_output = classifier.export_detailed_results(results, "json")
        with open("classification_test_results.json", "w", encoding="utf-8") as f:
            f.write(json_output)
        print("✅ JSON结果已保存到 classification_test_results.json")
        
        # CSV导出
        csv_output = classifier.export_detailed_results(results, "csv")
        with open("classification_test_results.csv", "w", encoding="utf-8") as f:
            f.write(csv_output)
        print("✅ CSV结果已保存到 classification_test_results.csv")
        
        print("\n🎉 测试完成！")
        
        # 总结
        print("\n📊 测试总结:")
        print(f"  - 总计测试: {len(results)}个氨基酸")
        print(f"  - 成功分类: {len([r for r in results if r.confidence > 0.5])}个")
        print(f"  - 验证通过: {len([r for r in results if r.validation_passed])}个") 
        if results:
            print(f"  - 平均置信度: {sum(r.confidence for r in results) / len(results):.1%}")
        else:
            print(f"  - 平均置信度: N/A (无有效结果)")
        
    except Exception as e:
        print(f"❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()