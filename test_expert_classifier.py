#!/usr/bin/env python3
"""
测试专家标准分类器
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.analysis.expert_classifier import ExpertAminoAcidClassifier
from pdb_uachecker.core.database.enhanced_manager import EnhancedDatabaseManager
from pdb_uachecker.utils.config import Config

def test_expert_classifier():
    """测试专家分类器"""
    print("🔬 专家标准氨基酸分类器测试")
    print("=" * 60)
    
    # 初始化
    config = Config()
    db = EnhancedDatabaseManager(config)
    expert_classifier = ExpertAminoAcidClassifier()
    
    # 测试案例
    test_cases = [
        # D型氨基酸测试
        'DAL',  # D-丙氨酸（如果存在）
        '0AF',  # 已知D型
        '0BN',  # 已知D型
        
        # 芳香性氨基酸
        'PHE',  # 标准苯丙氨酸
        'TYR',  # 标准酪氨酸
        '004',  # 芳香性氨基酸
        
        # 可能的Beta/Gamma氨基酸
        'ABA',  # α-氨基丁酸
        'B3A',  # 如果存在Beta-丙氨酸
        
        # N-甲基化
        'SAR',  # 肌氨酸（N-甲基甘氨酸）
        
        # 之前误判的案例
        'NLE',  # 正亮氨酸
        'CSA',  # 半胱氨酸修饰
    ]
    
    # 从数据库获取测试氨基酸
    all_amino_acids = db.get_all_amino_acids()
    test_amino_acids = []
    
    for case_id in test_cases:
        matches = [aa for aa in all_amino_acids if aa.id == case_id]
        if matches:
            test_amino_acids.append(matches[0])
        else:
            print(f"⚠️ 未找到: {case_id}")
    
    print(f"\\n📊 测试 {len(test_amino_acids)} 个氨基酸:")
    print("-" * 60)
    
    for i, amino_acid in enumerate(test_amino_acids, 1):
        print(f"\\n🧪 [{i:2d}] {amino_acid.id} - {amino_acid.name}")
        print(f"     SMILES: {amino_acid.smiles}")
        
        # 专家分类
        result = expert_classifier.classify_amino_acid(amino_acid)
        
        print(f"     🏗️  骨架: {result.get('backbone_type', {}).get('type', 'unknown')}")
        print(f"     🏷️  主要分类: {', '.join(result['primary_categories'])}")
        
        if result.get('secondary_features'):
            print(f"     ✨ 次要特征: {', '.join(result['secondary_features'])}")
        
        print(f"     🎯 置信度: {result['confidence']:.2f}")
        
        # 显示详细分析
        if 'modifications' in result:
            mods = result['modifications']
            
            if mods.get('chirality', {}).get('has_chirality'):
                chirality = mods['chirality']
                print(f"     🔄 手性: {chirality.get('d_l_form', 'unknown')} (置信度: {chirality.get('confidence', 0):.2f})")
            
            if mods.get('n_methylation', {}).get('is_n_methylated'):
                print(f"     🔗 N-甲基化: 是")
            
            if mods.get('aromaticity', {}).get('is_aromatic'):
                print(f"     💍 芳香性: 是")
            
            if mods.get('cyclization', {}).get('is_cyclic'):
                rings = mods['cyclization'].get('ring_count', 0)
                print(f"     🔄 环状结构: {rings}个环")
        
        # PDB验证
        if 'expert_validation' in result:
            validation = result['expert_validation']
            if validation.get('confidence_boost', 0) > 0:
                print(f"     ✅ PDB验证加成: +{validation['confidence_boost']:.2f}")
            if validation.get('warnings'):
                print(f"     ⚠️  验证警告: {', '.join(validation['warnings'])}")

def compare_with_old_classifier():
    """与原分类器对比"""
    print("\\n\\n🔄 新旧分类器对比")
    print("=" * 60)
    
    # 导入原分类器
    from pdb_uachecker.analysis.intelligent_classifier import IntelligentAminoAcidClassifier
    
    config = Config()
    db = EnhancedDatabaseManager(config)
    
    old_classifier = IntelligentAminoAcidClassifier()
    new_classifier = ExpertAminoAcidClassifier()
    
    # 选择几个关键测试案例
    comparison_cases = ['0AF', 'NLE', 'PHE', 'ABA']
    all_amino_acids = db.get_all_amino_acids()
    
    for case_id in comparison_cases:
        matches = [aa for aa in all_amino_acids if aa.id == case_id]
        if not matches:
            continue
        
        aa = matches[0]
        
        print(f"\\n📋 {aa.id} - {aa.name}")
        print(f"   SMILES: {aa.smiles}")
        
        # 旧分类器
        old_result = old_classifier.classify_amino_acid(aa)
        print(f"   📜 旧分类: {old_result['standard_category']} (置信度: {old_result['confidence']:.2f})")
        
        # 新分类器
        new_result = new_classifier.classify_amino_acid(aa)
        print(f"   🔬 新分类: {', '.join(new_result['primary_categories'])} (置信度: {new_result['confidence']:.2f})")
        
        # 对比分析
        new_primary = new_result['primary_categories'][0] if new_result['primary_categories'] else 'None'
        if old_result['standard_category'] != new_primary:
            print(f"   🔍 分类差异: {old_result['standard_category']} → {new_primary}")

if __name__ == "__main__":
    test_expert_classifier()
    compare_with_old_classifier()