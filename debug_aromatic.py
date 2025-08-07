#!/usr/bin/env python3
"""
调试芳香性分类误判问题
"""

import sys
from pathlib import Path
import re

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.analysis.intelligent_classifier import IntelligentAminoAcidClassifier
from pdb_uachecker.core.database.enhanced_manager import EnhancedDatabaseManager
from pdb_uachecker.utils.config import Config

def debug_aromatic_classification():
    """调试芳香性分类问题"""
    print("🔧 调试芳香性分类误判问题")
    print("=" * 60)
    
    # 初始化
    config = Config()
    db = EnhancedDatabaseManager(config)
    classifier = IntelligentAminoAcidClassifier()
    
    # 获取问题案例
    problem_cases = ['NLE', 'ABA', 'CSA']
    
    print("🧪 分析问题案例:")
    for aa_id in problem_cases:
        amino_acids = [aa for aa in db.get_all_amino_acids() if aa.id == aa_id]
        if amino_acids:
            aa = amino_acids[0]
            print(f"\n📋 {aa.id} - {aa.name}")
            print(f"   SMILES: {aa.smiles}")
            
            # 测试当前芳香性检测规则
            smiles = aa.smiles.lower()
            aromatic_patterns = [
                r'c1ccccc1',    # 苯环
                r'c1cccc[noh]1', # 五元芳香环
                r'c1ccc[noh]c1', # 六元含杂原子芳香环
                r'c\d+c.*c\d+',  # 一般芳香环模式 - 问题规则
            ]
            
            print("   当前规则匹配:")
            for i, pattern in enumerate(aromatic_patterns, 1):
                match = re.search(pattern, smiles)
                if match:
                    print(f"     ✓ 规则{i}: {pattern} -> {match.group()}")
                else:
                    print(f"     ✗ 规则{i}: {pattern}")
            
            # 分析SMILES中的芳香性特征
            print("   SMILES分析:")
            print(f"     小写c数量: {smiles.count('c')}")
            print(f"     环标记: {[c for c in smiles if c.isdigit()]}")
            print(f"     是否包含苯环: {'c1ccccc1' in smiles}")
            print(f"     是否包含杂环: {any(x in smiles for x in ['[nh]', '[oh]', 'n', 'o', 's']) and 'c' in smiles}")
    
    print("\n🔬 真正的芳香性案例检查:")
    aromatic_cases = ['PHE', '004', '0A1', 'TYR']  # 已知芳香性氨基酸
    
    for aa_id in aromatic_cases:
        amino_acids = [aa for aa in db.get_all_amino_acids() if aa.id == aa_id]
        if amino_acids:
            aa = amino_acids[0]
            print(f"\n✅ {aa.id} - {aa.name}")
            print(f"   SMILES: {aa.smiles}")
            smiles = aa.smiles.lower()
            if 'c1ccccc1' in smiles:
                print("   ✓ 含苯环")
            elif any(pattern in smiles for pattern in ['c1cccc', 'c1ccc', 'c1cc']):
                print("   ✓ 含芳香环结构")

def propose_improved_rules():
    """提出改进的芳香性检测规则"""
    print("\n💡 改进的芳香性检测规则:")
    print("=" * 60)
    
    improved_patterns = [
        r'c1ccccc1',        # 苯环 - 精确匹配
        r'c1cccc[noh]1',    # 含杂原子的五元芳香环
        r'c1ccc[noh]c1',    # 含杂原子的六元芳香环  
        r'c1c\[nh\]c[2-9]c1cccc[2-9]',  # 吲哚类结构
        r'c1cc[2-9]c[2-9]c1',  # 简化芳香环 (移除问题规则)
    ]
    
    print("新规则:")
    for i, pattern in enumerate(improved_patterns, 1):
        print(f"  {i}. {pattern}")
    
    print("\n❌ 移除的问题规则:")
    print(f"     r'c\\d+c.*c\\d+' - 过于宽泛，匹配任何含c的结构")

if __name__ == "__main__":
    debug_aromatic_classification()
    propose_improved_rules()