#!/usr/bin/env python3
"""
高级SMILES修正器 - 修复价态错误和其他化学结构问题
"""

import re
import sys
from pathlib import Path

def advanced_smiles_fix(smiles: str) -> str:
    """高级SMILES修正，处理价态错误"""
    
    # 第一轮：基础修正
    fixed = smiles
    
    # 移除质子化状态
    basic_fixes = [
        ('[NH3]', 'N'),
        ('[NH2]', 'N'),
        ('[NH+]', 'N'),
        ('[N+]', 'N'),
        ('[OH2+]', 'O'),
        ('[O-]', 'O'),
        ('[COO-]', 'C(=O)O'),
    ]
    
    for old, new in basic_fixes:
        fixed = fixed.replace(old, new)
    
    # 第二轮：价态错误修正
    valence_fixes = [
        # 碳价态错误
        (r'\[C\]\(=O\)=O', 'C(=O)O'),           # [C](=O)=O -> C(=O)O
        (r'\[C\]\(=N\)=N', 'C(=N)N'),           # [C](=N)=N -> C(=N)N
        (r'\[C\]=\([^)]+\)=O', lambda m: m.group(0).replace('[C]', 'C').replace('=O', 'O')),
        
        # 氮价态错误
        (r'\[NH\]=\[C\]\(=N\)=N', 'NC(=N)N'),   # [NH]=[C](=N)=N -> NC(=N)N
        (r'\[NH\]=\[C\]\(=O\)=N', 'NC(=O)N'),   # [NH]=[C](=O)=N -> NC(=O)N
        
        # 通用价态修正
        (r'\[([CNOS])\]', r'\1'),               # 移除不必要的方括号
        (r'=\[([CNOS])\]=', r'=\1'),            # 移除双键中的方括号
    ]
    
    for pattern, replacement in valence_fixes:
        if callable(replacement):
            fixed = re.sub(pattern, replacement, fixed)
        else:
            fixed = re.sub(pattern, replacement, fixed)
    
    # 第三轮：特殊结构修正
    special_fixes = [
        # 羧酸基团标准化
        ('C(=O)=O', 'C(=O)O'),
        ('C(O)=O', 'C(=O)O'),
        
        # 胍基修正
        ('C(=N)=N', 'C(=N)N'),
        ('C(N)=N', 'C(=N)N'),
        
        # 醛基修正
        ('C=O=O', 'C(=O)O'),
        ('C(=O)=C', 'C(=O)C'),
        
        # 硒原子修正
        ('[SeH]', '[Se]'),
        
        # 氯原子修正
        ('ClC', 'CCl'),
    ]
    
    for old, new in special_fixes:
        fixed = fixed.replace(old, new)
    
    return fixed

def test_problematic_smiles():
    """测试有问题的SMILES修正"""
    
    problematic_cases = [
        "CC[C@H]([C@@H]([C](=O)=O)N)O",           # VAH - 价态错误
        "N[C@H]([C](=O)=O)Cc1c[nH]c2c1cccc2O",    # 0AF - 价态错误
        "N[C@H](C(=O)O)Cc1ccc(cc1)[C](=N)=N",     # 0BN - 胍基错误
        "OC(=O)[C@H](CCC[C](=[C](=O)=O)=O)N",     # 26P - 复杂价态错误
        "O[C@H](C[C@@H](C(=O)O)N)C[NH]=[C](=N)=N", # ARO - 胍基错误
        "N[C@H]([C](=O)=O)C[SeH]",                # CSE - 硒原子错误
    ]
    
    print("🧪 测试有问题的SMILES修正:")
    print("=" * 60)
    
    for i, smiles in enumerate(problematic_cases, 1):
        print(f"\n{i}. 原始: {smiles}")
        fixed = advanced_smiles_fix(smiles)
        print(f"   修正: {fixed}")
        
        # 验证修正结果
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(fixed)
            valid = mol is not None
            print(f"   有效: {'✅' if valid else '❌'}")
            
            if valid:
                # 计算分子式
                formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                print(f"   分子式: {formula}")
        except Exception as e:
            print(f"   验证失败: {e}")

def fix_failed_amino_acids():
    """修正失败的氨基酸SMILES"""
    
    print("🔧 修正失败的氨基酸SMILES:")
    print("=" * 50)
    
    # 已知失败的氨基酸ID（从您的输出中提取）
    failed_ids = [
        '0AF', '0BN', '200', '26P', '2AS', '2NP', '32T', '3GL', '4HT', 
        '6CL', '6CW', 'ABA', 'ADAM', 'AGM', 'ALN', 'ARO', 'AS2', 'AZDA',
        'BHD', 'BTH3', 'C2N', 'CCS', 'CIR', 'CPA3', 'CPG2', 'CSE', 'CTE',
        'DAB', 'DMK', 'DPP', 'FGL', 'GBUT', 'VAH'
    ]
    
    structures_dir = Path("data/structures")
    success_count = 0
    
    for amino_id in failed_ids[:10]:  # 测试前10个
        smi_file = structures_dir / amino_id / f"{amino_id}.smi"
        
        if not smi_file.exists():
            print(f"❌ {amino_id}: 文件不存在")
            continue
        
        try:
            # 读取原始SMILES
            with open(smi_file, 'r') as f:
                content = f.read().strip()
                if content:
                    original_smiles = content.split()[0]
                    
                    # 应用高级修正
                    fixed_smiles = advanced_smiles_fix(original_smiles)
                    
                    print(f"\n{amino_id}:")
                    print(f"  原始: {original_smiles}")
                    print(f"  修正: {fixed_smiles}")
                    
                    # 验证修正结果
                    try:
                        from rdkit import Chem
                        mol = Chem.MolFromSmiles(fixed_smiles)
                        if mol:
                            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                            print(f"  ✅ 有效 - {formula}")
                            success_count += 1
                        else:
                            print(f"  ❌ 仍然无效")
                    except Exception as e:
                        print(f"  ❌ 验证失败: {e}")
                        
        except Exception as e:
            print(f"❌ {amino_id}: 处理失败 - {e}")
    
    print(f"\n📊 修正结果: {success_count}/{min(10, len(failed_ids))} 成功")

def create_manual_fixes():
    """为特定的失败案例创建手动修正"""
    
    # 手动修正的SMILES（基于化学知识）
    manual_fixes = {
        'VAH': 'N[C@@H](C(=O)O)[C@H](CC)O',           # 羟基缬氨酸
        '0AF': 'N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2O',     # 5-羟基色氨酸
        '0BN': 'N[C@H](C(=O)O)Cc1ccc(cc1)C(=N)N',     # 4-胍基苯丙氨酸
        '200': 'N[C@H](C(=O)O)Cc1ccc(cc1)Cl',         # 4-氯苯丙氨酸
        '2AS': 'N[C@H](C(=O)O)[C@@H](C)C(=O)O',       # 苏氨酸衍生物
        'ARO': 'N[C@@H](CCCNC(=N)N)C(=O)O',           # 精氨酸衍生物
        'CSE': 'N[C@H](C(=O)O)C[Se]',                 # 硒代半胱氨酸
        'ABA': 'N[C@@H](CC)C(=O)O',                   # α-氨基丁酸
        'DAB': 'N[C@@H](CCN)C(=O)O',                  # 二氨基丁酸
        'DPP': 'N[C@@H](CN)C(=O)O',                   # 二氨基丙酸
    }
    
    print("📝 手动修正的SMILES:")
    print("=" * 40)
    
    for amino_id, fixed_smiles in manual_fixes.items():
        print(f"{amino_id}: {fixed_smiles}")
        
        # 验证手动修正
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(fixed_smiles)
            if mol:
                formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                print(f"  ✅ {formula}")
            else:
                print(f"  ❌ 无效")
        except:
            print(f"  ❌ 验证失败")
    
    return manual_fixes

def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(description="高级SMILES修正器")
    parser.add_argument("--test", action="store_true", help="测试有问题的SMILES")
    parser.add_argument("--fix", action="store_true", help="修正失败的氨基酸")
    parser.add_argument("--manual", action="store_true", help="显示手动修正")
    
    args = parser.parse_args()
    
    if args.test:
        test_problematic_smiles()
    elif args.fix:
        fix_failed_amino_acids()
    elif args.manual:
        create_manual_fixes()
    else:
        print("使用 --test, --fix, 或 --manual 选项")

if __name__ == "__main__":
    main()
