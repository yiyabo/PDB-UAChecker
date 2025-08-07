#!/usr/bin/env python3
"""
详细调试分类决策过程
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.analysis.intelligent_classifier import IntelligentAminoAcidClassifier
from pdb_uachecker.core.database.enhanced_manager import EnhancedDatabaseManager
from pdb_uachecker.utils.config import Config

def detailed_debug():
    """详细调试分类过程"""
    config = Config()
    db = EnhancedDatabaseManager(config)
    classifier = IntelligentAminoAcidClassifier()

    # 获取NLE案例
    amino_acids = [aa for aa in db.get_all_amino_acids() if aa.id == 'NLE']
    aa = amino_acids[0]

    print('🔍 详细调试 NLE 分类决策:')
    print(f'SMILES: {aa.smiles}')

    # 逐步检查
    print('\n1️⃣ Tier 1 检查:')
    tier1 = classifier._classify_tier1(aa)
    print(f'   结果: {tier1}')

    print('\n2️⃣ Tier 2 检查:') 
    tier2 = classifier._classify_tier2(aa)
    print(f'   结果: {tier2}')

    print('\n3️⃣ Tier 3 检查:')
    tier3 = classifier._classify_tier3(aa)
    print(f'   结果: {tier3}')
    
    print('\n🎯 完整分类:')
    full_result = classifier.classify_amino_acid(aa)
    print(f'   最终: {full_result}')
    
    # 检查Tier 2中的芳香性验证
    print('\n🔬 芳香性验证详细检查:')
    aromatic_result = classifier._is_aromatic_amino_acid(aa)
    print(f'   芳香性检测: {aromatic_result}')
    
    # 检查Tier 3的所有规则
    print('\n📊 Tier 3 规则逐项检查:')
    tier3_rules = classifier.classification_tiers['tier3_intelligent']
    for rule_name, rule_config in tier3_rules.items():
        validator = rule_config['validator']
        result = validator(aa)
        print(f'   {rule_name}: {result}')

if __name__ == "__main__":
    detailed_debug()