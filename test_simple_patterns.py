#!/usr/bin/env python3
"""
测试简单模式匹配
"""

import re

def test_patterns():
    """测试SMILES模式匹配"""
    
    test_cases = [
        ('ABA', 'N[C@@H](CC)C(=O)O'),
        ('D-form', 'N[C@H](C(=O)O)Cc1ccccc1'), 
        ('PHE', 'N[C@@H](Cc1ccccc1)C(=O)O'),
        ('Complex', 'N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2O'),
    ]
    
    # α氨基酸模式
    alpha_patterns = [
        r'N\[C@@?H?\]\([^)]*\)C\(=O\)O',  # N[C@H](R)C(=O)O
        r'N\[C@@?H?\].*C\(=O\)O',        # N[C@H]...C(=O)O
    ]
    
    print("🧪 测试α氨基酸模式匹配:")
    for name, smiles in test_cases:
        print(f"\\n{name}: {smiles}")
        
        for i, pattern in enumerate(alpha_patterns, 1):
            match = re.search(pattern, smiles)
            print(f"  模式{i}: {'✅' if match else '❌'} {pattern}")
            if match:
                print(f"    匹配: {match.group()}")

if __name__ == "__main__":
    test_patterns()