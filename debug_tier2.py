#!/usr/bin/env python3
"""
专门调试Tier 2分类问题
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.analysis.intelligent_classifier import IntelligentAminoAcidClassifier
from pdb_uachecker.core.database.enhanced_manager import EnhancedDatabaseManager
from pdb_uachecker.utils.config import Config

def debug_tier2_logic():
    """专门调试Tier 2分类逻辑"""
    config = Config()
    db = EnhancedDatabaseManager(config)
    classifier = IntelligentAminoAcidClassifier()

    # 获取NLE案例
    amino_acids = [aa for aa in db.get_all_amino_acids() if aa.id == 'NLE']
    aa = amino_acids[0]

    print('🔧 Tier 2 逻辑调试:')
    print(f'SMILES: {aa.smiles}')

    # 获取Tier 2规则
    tier2_rules = classifier.classification_tiers['tier2_high_confidence']
    print(f'\\n📋 Tier 2 规则列表: {list(tier2_rules.keys())}')

    # 按置信度排序（模拟_classify_tier2的逻辑）
    sorted_rules = sorted(tier2_rules.items(), 
                        key=lambda x: x[1]['confidence'], reverse=True)
    
    print('\\n🔍 按置信度顺序检查每个规则:')
    for i, (category_key, rule_config) in enumerate(sorted_rules, 1):
        validator = rule_config['validator']
        result = validator(aa)
        
        print(f'\\n{i}. {category_key} (置信度: {rule_config["confidence"]})')
        print(f'   验证器结果: {result}')
        print(f'   is_match: {result.get("is_match", False) if isinstance(result, dict) else result}')
        
        if isinstance(result, dict) and result.get('is_match', False):
            category = category_key.replace('_amino_acid', '').replace('_', '-').title()
            print(f'   ✅ 匹配！将返回: {category}')
            break
        else:
            print(f'   ❌ 不匹配')

if __name__ == "__main__":
    debug_tier2_logic()