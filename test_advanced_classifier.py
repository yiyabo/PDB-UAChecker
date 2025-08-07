#!/usr/bin/env python3
"""
测试高级专家分类器
"""

import sys
from pathlib import Path
import time

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.analysis.advanced_expert_classifier import AdvancedExpertAminoAcidClassifier
from pdb_uachecker.analysis.intelligent_classifier import IntelligentAminoAcidClassifier
from pdb_uachecker.core.database.enhanced_manager import EnhancedDatabaseManager
from pdb_uachecker.utils.config import Config


def test_advanced_classifier():
    """测试高级专家分类器"""
    print("🔬 高级专家分类器测试")
    print("=" * 70)
    
    # 初始化
    config = Config()
    db = EnhancedDatabaseManager(config)
    advanced_classifier = AdvancedExpertAminoAcidClassifier()
    
    # 重点测试案例
    critical_test_cases = [
        # D型氨基酸
        '0AF',  # 已知D型（含吲哚环）
        '0BN',  # 已知D型（含苯环）
        
        # 芳香性氨基酸
        'PHE',  # 标准苯丙氨酸
        'TYR',  # 标准酪氨酸
        '004',  # 非标准芳香性
        
        # 之前误判的案例
        'NLE',  # 正亮氨酸
        'ABA',  # α-氨基丁酸
        'CSA',  # 半胱氨酸修饰
        
        # 复杂结构
        'PRO',  # 脯氨酸（环状）
        'HIS',  # 组氨酸（含氮杂环）
    ]
    
    # 从数据库获取测试氨基酸
    all_amino_acids = db.get_all_amino_acids()
    test_amino_acids = []
    
    for case_id in critical_test_cases:
        matches = [aa for aa in all_amino_acids if aa.id == case_id]
        if matches:
            test_amino_acids.append(matches[0])
        else:
            print(f"⚠️ 未找到: {case_id}")
    
    print(f"\\n📊 重点测试 {len(test_amino_acids)} 个关键氨基酸:")
    print("-" * 70)
    
    for i, amino_acid in enumerate(test_amino_acids, 1):
        print(f"\\n🧪 [{i:2d}] {amino_acid.id} - {amino_acid.name}")
        print(f"     SMILES: {amino_acid.smiles}")
        
        # 高级分类
        start_time = time.time()
        result = advanced_classifier.classify_amino_acid(amino_acid)
        analysis_time = time.time() - start_time
        
        print(f"     🏷️  主要分类: {', '.join(result['primary_categories'])}")
        
        if result.get('secondary_features'):
            print(f"     ✨ 次要特征: {', '.join(result['secondary_features'])}")
        
        print(f"     🎯 置信度: {result['confidence']:.2f}")
        print(f"     ⏱️  分析耗时: {analysis_time*1000:.1f}ms")
        
        # 显示详细分析
        if 'chirality_analysis' in result:
            chirality = result['chirality_analysis']
            if chirality.get('has_chirality'):
                d_l_form = chirality.get('d_l_form', 'unknown')
                cip_config = chirality.get('configuration', 'unknown')
                print(f"     🔄 手性分析: {d_l_form}型 ({cip_config}构型)")
        
        if 'structure_analysis' in result:
            structure = result['structure_analysis']
            backbone = structure.get('backbone_analysis', {})
            backbone_type = backbone.get('type', 'unknown')
            backbone_conf = backbone.get('confidence', 0)
            print(f"     🏗️  骨架分析: {backbone_type} (置信度: {backbone_conf:.2f})")
        
        # PDB验证结果
        if 'pdb_validation' in result:
            validation = result['pdb_validation']
            if validation.get('confirmations'):
                print(f"     ✅ PDB验证: {', '.join(validation['confirmations'])}")
            if validation.get('warnings'):
                print(f"     ⚠️  PDB警告: {', '.join(validation['warnings'])}")
        
        # 错误处理
        if 'error_message' in result:
            print(f"     ❌ 错误: {result['error_message']}")


def comprehensive_comparison():
    """全面对比三种分类器"""
    print("\\n\\n🔄 三种分类器全面对比")
    print("=" * 70)
    
    # 初始化所有分类器
    config = Config()
    db = EnhancedDatabaseManager(config)
    
    old_classifier = IntelligentAminoAcidClassifier()
    advanced_classifier = AdvancedExpertAminoAcidClassifier()
    
    # 对比测试案例
    comparison_cases = ['0AF', 'NLE', 'PHE', 'ABA', '004', 'PRO']
    all_amino_acids = db.get_all_amino_acids()
    
    print(f"\\n📋 对比分析 {len(comparison_cases)} 个关键案例:")
    print("-" * 70)
    
    for case_id in comparison_cases:
        matches = [aa for aa in all_amino_acids if aa.id == case_id]
        if not matches:
            continue
        
        aa = matches[0]
        
        print(f"\\n🔍 {aa.id} - {aa.name}")
        print(f"   SMILES: {aa.smiles}")
        
        # 旧分类器
        old_result = old_classifier.classify_amino_acid(aa)
        print(f"   📜 智能分类器: {old_result['standard_category']} (置信度: {old_result['confidence']:.2f})")
        
        # 高级分类器
        advanced_result = advanced_classifier.classify_amino_acid(aa)
        primary = ', '.join(advanced_result['primary_categories'])
        print(f"   🔬 高级专家分类器: {primary} (置信度: {advanced_result['confidence']:.2f})")
        
        # 分析差异
        old_category = old_result['standard_category']
        new_primary = advanced_result['primary_categories'][0] if advanced_result['primary_categories'] else 'None'
        
        if old_category != new_primary:
            print(f"   🔍 分类变化: {old_category} → {new_primary}")
            
            # 分析变化原因
            if 'chirality_analysis' in advanced_result:
                chirality = advanced_result['chirality_analysis']
                if chirality.get('has_chirality') and chirality.get('d_l_form') == 'D':
                    print(f"   💡 原因: CIP分析检测到D型手性")
            
            if 'structure_analysis' in advanced_result:
                structure = advanced_result['structure_analysis']
                backbone = structure.get('backbone_analysis', {})
                if backbone.get('type') in ['beta', 'gamma']:
                    print(f"   💡 原因: 精确骨架分析检测到{backbone['type']}型")


def performance_benchmark():
    """性能基准测试"""
    print("\\n\\n⚡ 性能基准测试")
    print("=" * 70)
    
    config = Config()
    db = EnhancedDatabaseManager(config)
    
    old_classifier = IntelligentAminoAcidClassifier()
    advanced_classifier = AdvancedExpertAminoAcidClassifier()
    
    # 获取测试样本
    all_amino_acids = db.get_all_amino_acids()
    test_sample = all_amino_acids[:50]  # 测试前50个
    
    print(f"\\n🧪 测试样本: {len(test_sample)} 个氨基酸")
    
    # 测试旧分类器
    print("\\n📜 智能分类器性能测试...")
    start_time = time.time()
    old_results = []
    for aa in test_sample:
        result = old_classifier.classify_amino_acid(aa)
        old_results.append(result)
    old_time = time.time() - start_time
    
    # 测试高级分类器
    print("🔬 高级专家分类器性能测试...")
    start_time = time.time()
    advanced_results = []
    for aa in test_sample:
        result = advanced_classifier.classify_amino_acid(aa)
        advanced_results.append(result)
    advanced_time = time.time() - start_time
    
    # 性能对比
    print(f"\\n📊 性能对比结果:")
    print(f"   智能分类器:")
    print(f"     总耗时: {old_time:.2f}s")
    print(f"     平均耗时: {old_time/len(test_sample)*1000:.1f}ms/样本")
    
    print(f"   高级专家分类器:")
    print(f"     总耗时: {advanced_time:.2f}s") 
    print(f"     平均耗时: {advanced_time/len(test_sample)*1000:.1f}ms/样本")
    
    print(f"   性能比较: {advanced_time/old_time:.1f}x 耗时")
    
    # 准确性对比
    print(f"\\n🎯 准确性对比:")
    
    # 统计高置信度分类
    old_high_conf = len([r for r in old_results if r.get('confidence', 0) >= 0.8])
    advanced_high_conf = len([r for r in advanced_results if r.get('confidence', 0) >= 0.8])
    
    print(f"   高置信度分类 (≥0.8):")
    print(f"     智能分类器: {old_high_conf}/{len(test_sample)} ({old_high_conf/len(test_sample)*100:.1f}%)")
    print(f"     高级分类器: {advanced_high_conf}/{len(test_sample)} ({advanced_high_conf/len(test_sample)*100:.1f}%)")


if __name__ == "__main__":
    try:
        test_advanced_classifier()
        comprehensive_comparison() 
        performance_benchmark()
        
        print("\\n\\n🎉 高级专家分类器测试完成！")
        
    except Exception as e:
        print(f"\\n❌ 测试过程中发生错误: {str(e)}")
        import traceback
        traceback.print_exc()