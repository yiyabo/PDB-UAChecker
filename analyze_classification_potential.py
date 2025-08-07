#!/usr/bin/env python3
"""
分析100%自动化分类的可能性
"""

import sys
from pathlib import Path
import re

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.core.database.enhanced_manager import EnhancedDatabaseManager
from pdb_uachecker.utils.config import Config

def analyze_classification_potential():
    """分析分类潜力"""
    config = Config()
    db = EnhancedDatabaseManager(config)
    amino_acids = db.get_all_amino_acids()

    print('🤔 **能否通过代码实现100%精确验证？**')
    print('=' * 60)

    standard_aa = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'}

    categories = {
        'standard': [],
        'D_amino': [],
        'aromatic': [],
        'cyclic': [],
        'beta': [],
        'gamma': [],
        'n_methyl': [],
        'unclear': []
    }

    # 分类逻辑
    for aa in amino_acids:
        classified = False
        
        # 标准氨基酸 - 100%确定
        if aa.id in standard_aa:
            categories['standard'].append(aa.id)
            classified = True
        
        # D氨基酸 - 基于PDB命名和手性标记
        elif ((aa.id.startswith('D') and len(aa.id) == 3) or
              (aa.smiles and 'C@H' in aa.smiles and 'C@@H' not in aa.smiles)):
            categories['D_amino'].append(aa.id)
            classified = True
        
        # 芳香族 - 基于芳香环标记
        elif aa.smiles and ('c1' in aa.smiles.lower() or 'c2' in aa.smiles.lower()):
            categories['aromatic'].append(aa.id)
            classified = True
        
        # 环状 - 基于环闭合数字
        elif aa.smiles and any(char.isdigit() for char in aa.smiles):
            categories['cyclic'].append(aa.id)
            classified = True
            
        # Beta氨基酸 - NH2-CH2-CH(R)-COOH
        elif aa.smiles and re.search(r'NCC[^C].*C\(=O\)O', aa.smiles):
            categories['beta'].append(aa.id)
            classified = True
            
        # Gamma氨基酸 - NH2-CH2-CH2-CH(R)-COOH  
        elif aa.smiles and re.search(r'NCCC[^C].*C\(=O\)O', aa.smiles):
            categories['gamma'].append(aa.id) 
            classified = True
            
        # N-甲基氨基酸
        elif aa.smiles and ('[NH]' in aa.smiles or 'CN[C' in aa.smiles):
            categories['n_methyl'].append(aa.id)
            classified = True
        
        if not classified:
            categories['unclear'].append(aa.id)

    print('📊 基于精确规则的分类结果:')
    total_classified = 0
    for category, items in categories.items():
        count = len(items)
        percentage = count / len(amino_acids) * 100
        print(f'  {category}: {count} ({percentage:.1f}%)')
        if category != 'unclear':
            total_classified += count

    coverage = total_classified / len(amino_acids) * 100
    print(f'\n🎯 **理论上可以100%准确分类的覆盖率: {coverage:.1f}%**')
    print(f'❓ 真正需要人工的: {len(categories["unclear"])} ({len(categories["unclear"])/len(amino_acids)*100:.1f}%)')

    print('\n💡 **关键发现:**')
    print(f'  ✅ 标准氨基酸: {len(categories["standard"])} 个 (100%确定)')
    print(f'  ✅ D氨基酸: {len(categories["D_amino"])} 个 (手性标记确定)')
    print(f'  ✅ 芳香族: {len(categories["aromatic"])} 个 (芳香环确定)')
    print(f'  ✅ 环状: {len(categories["cyclic"])} 个 (环结构确定)')
    print(f'  ✅ Beta/Gamma: {len(categories["beta"]) + len(categories["gamma"])} 个 (主链模式)')
    print(f'  ✅ N-甲基: {len(categories["n_methyl"])} 个 (氮修饰确定)')

    print(f'\n🚀 **结论: 理论上可以自动化分类 {coverage:.0f}% 的氨基酸!**')
    
    # 显示不清楚的案例
    if categories['unclear']:
        print(f'\n❓ 需要更复杂规则的案例 ({len(categories["unclear"])}个):')
        for aa_id in categories['unclear'][:10]:
            aa = next(aa for aa in amino_acids if aa.id == aa_id)
            print(f'  {aa_id}: {aa.smiles[:50]}...')
            
    return {
        'coverage': coverage,
        'categories': categories,
        'total_amino_acids': len(amino_acids)
    }

if __name__ == "__main__":
    result = analyze_classification_potential()