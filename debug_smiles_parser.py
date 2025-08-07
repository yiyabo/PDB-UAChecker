#!/usr/bin/env python3
"""
调试SMILES解析器
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.utils.smiles_parser import AminoAcidStructureAnalyzer

def test_smiles_parser():
    """测试SMILES解析器"""
    print("🔧 SMILES解析器调试")
    print("=" * 50)
    
    analyzer = AminoAcidStructureAnalyzer()
    
    # 测试案例
    test_smiles = [
        ('Simple alpha', 'N[C@@H](CC)C(=O)O'),  # ABA - 简单α氨基酸
        ('D-form', 'N[C@H](C(=O)O)Cc1ccccc1'),  # 简化的D型
        ('Aromatic', 'N[C@@H](Cc1ccccc1)C(=O)O'),  # PHE - 芳香性
        ('Complex', 'N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2O'),  # 0AF - 复杂结构
    ]
    
    for name, smiles in test_smiles:
        print(f"\\n🧪 测试 {name}: {smiles}")
        
        try:
            result = analyzer.analyze_amino_acid_structure(smiles)
            
            print(f"   骨架分析: {result.get('backbone_analysis', {})}")
            print(f"   官能团: {len(result.get('functional_groups', []))} 个")
            print(f"   手性分析: {result.get('chirality_analysis', {})}")
            print(f"   环系统: {len(result.get('ring_systems', []))} 个")
            print(f"   置信度: {result.get('confidence', 0):.2f}")
            
            # 详细显示官能团
            for group in result.get('functional_groups', []):
                print(f"     - {group['type']}: {group}")
                
        except Exception as e:
            print(f"   ❌ 解析错误: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    test_smiles_parser()