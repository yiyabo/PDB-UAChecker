#!/usr/bin/env python3
"""
调试手性分析问题
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.utils.cip_rules import AminoAcidChiralityAnalyzer

def test_chirality_analysis():
    """测试手性分析"""
    print("🔄 手性分析调试")
    print("=" * 50)
    
    analyzer = AminoAcidChiralityAnalyzer()
    
    # 测试案例
    test_cases = [
        ('0AF-D型', 'N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2O'),  # D型
        ('0BN-D型', 'N[C@H](C(=O)O)Cc1ccc(cc1)C(=N)N'),   # D型  
        ('PHE-L型', 'N[C@@H](Cc1ccccc1)C(=O)O'),            # L型标准
        ('ABA-L型', 'N[C@@H](CC)C(=O)O'),                   # L型
    ]
    
    for name, smiles in test_cases:
        print(f"\\n🧪 测试 {name}: {smiles}")
        
        try:
            result = analyzer.analyze_chirality(smiles)
            
            print(f"   有手性: {result.get('has_chirality', False)}")
            print(f"   构型: {result.get('configuration', 'unknown')}")
            print(f"   D/L型: {result.get('d_l_form', 'unknown')}")
            print(f"   置信度: {result.get('confidence', 0):.2f}")
            
            # CIP分析详情
            if 'cip_analysis' in result:
                cip = result['cip_analysis']
                if cip.get('chiral_centers'):
                    for center in cip['chiral_centers']:
                        print(f"   手性中心: 原子{center.get('chiral_center', '?')} - {center.get('configuration', '?')}")
                        
        except Exception as e:
            print(f"   ❌ 分析错误: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    test_chirality_analysis()