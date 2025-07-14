#!/usr/bin/env python3
"""
最终导入脚本 - 整合所有SMILES修正，导入所有可能的氨基酸
"""

import os
import sys
import sqlite3
import json
import re
from pathlib import Path
from typing import Dict, List, Optional

sys.path.append('.')

def get_manual_fixes() -> Dict[str, str]:
    """获取手动修正的SMILES"""
    return {
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
        '26P': 'N[C@H](C(=O)O)CCCC(=O)C(=O)O',       # 修正的26P
        '2NP': 'N[C@H](C(=O)O)CCCC(=C)C(=O)O',       # 修正的2NP
        '32T': 'N[C@H](C(=O)O)Cc1c[nH]c2c1scc2',     # 修正的32T
        '3GL': 'N[C@H](C(=O)O)C[C@@H](O)C(=O)O',     # 修正的3GL
        '4HT': 'N[C@H](C(=O)O)Cc1c[nH]c2c1c(O)ccc2', # 修正的4HT
        '6CL': 'N[C@H](C(=O)O)CCC[C@H](N)C(=O)O',    # 修正的6CL
        '6CW': 'N[C@H](C(=O)O)Cc1c[nH]c2c1ccc(Cl)c2', # 修正的6CW
        'ADAM': 'N[C@H](C(=O)O)C[C@]12C[C@H]3C[C@@H](C2)C[C@@H](C1)C3', # 金刚烷氨基酸
        'AGM': 'N[C@H](C(=O)O)CCC[C@H](C)NC(=N)N',    # 修正的AGM
        'ALN': 'N[C@H](C(=O)O)Cc1cccc2c1cccc2',       # 萘丙氨酸
        'AS2': 'N[C@H](C(=O)O)CCC(=O)O',             # 修正的AS2
        'AZDA': 'N[C@H](C(=O)O)CN=[N+]=[N-]',        # 叠氮氨基酸
        'BHD': 'N[C@H](C(=O)O)[C@@H](O)C(=O)O',      # β-羟基天冬氨酸
        'BTH3': 'N[C@H](C(=O)O)Cc1csc2c1cccc2',      # 苯并噻吩氨基酸
        'C2N': 'N[C@@H](CCl)C(=O)O',                 # 氯代丙氨酸
        'CCS': 'N[C@H](C(=O)O)CSC(=O)C',             # 修正的CCS
        'CIR': 'N[C@H](C(=O)O)CCCNC(=O)N',           # 瓜氨酸
        'CPA3': 'N[C@H](C(=O)O)CC1CCCC1',            # 环戊基氨基酸
        'CPG2': 'N[C@H](C(=O)O)c1ccccc1Cl',          # 邻氯苯甘氨酸
        'CTE': 'N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2Cl',  # 氯色氨酸
        'DMK': 'N[C@H](C(=O)O)C(C)(C)C(=O)O',        # 修正的DMK
        'FGL': 'N[C@H](C(=O)O)C(=O)O',               # 氨基丙二酸
        'GBUT': 'N[C@H](C(=O)O)CCNC(=N)N',           # γ-胍基丁氨酸
    }

def advanced_smiles_fix(smiles: str) -> str:
    """高级SMILES修正"""
    fixed = smiles
    
    # 基础修正
    basic_fixes = [
        ('[NH3]', 'N'), ('[NH2]', 'N'), ('[NH+]', 'N'), ('[N+]', 'N'),
        ('[OH2+]', 'O'), ('[O-]', 'O'), ('[COO-]', 'C(=O)O'),
    ]
    
    for old, new in basic_fixes:
        fixed = fixed.replace(old, new)
    
    # 价态修正
    fixed = re.sub(r'\[C\]\(=O\)=O', 'C(=O)O', fixed)
    fixed = re.sub(r'\[C\]\(=N\)=N', 'C(=N)N', fixed)
    fixed = re.sub(r'\[([CNOS])\]', r'\1', fixed)
    
    # 特殊修正
    special_fixes = [
        ('C(=O)=O', 'C(=O)O'), ('C(=N)=N', 'C(=N)N'),
        ('[SeH]', '[Se]'), ('ClC', 'CCl'),
    ]
    
    for old, new in special_fixes:
        fixed = fixed.replace(old, new)
    
    return fixed

def get_smiles_for_amino_acid(amino_id: str) -> Optional[str]:
    """获取氨基酸的SMILES（优先使用手动修正）"""
    
    # 首先检查手动修正
    manual_fixes = get_manual_fixes()
    if amino_id in manual_fixes:
        return manual_fixes[amino_id]
    
    # 然后尝试从文件读取并修正
    smi_file = Path("data/structures") / amino_id / f"{amino_id}.smi"
    if smi_file.exists():
        try:
            with open(smi_file, 'r') as f:
                content = f.read().strip()
                if content:
                    original_smiles = content.split()[0]
                    return advanced_smiles_fix(original_smiles)
        except:
            pass
    
    return None

def import_amino_acid_final(amino_id: str) -> bool:
    """最终导入单个氨基酸"""
    try:
        # 获取SMILES
        smiles = get_smiles_for_amino_acid(amino_id)
        if not smiles:
            return False
        
        # 验证SMILES
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return False
        
        # 计算属性
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        weight = round(Descriptors.MolWt(mol), 3)
        
        # 计算原子组成
        atom_counts = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        
        # 获取名称
        standard_names = {
            'ALA': '丙氨酸', 'ARG': '精氨酸', 'ASN': '天冬酰胺', 'ASP': '天冬氨酸',
            'CYS': '半胱氨酸', 'GLN': '谷氨酰胺', 'GLU': '谷氨酸', 'GLY': '甘氨酸',
            'HIS': '组氨酸', 'ILE': '异亮氨酸', 'LEU': '亮氨酸', 'LYS': '赖氨酸',
            'MET': '蛋氨酸', 'PHE': '苯丙氨酸', 'PRO': '脯氨酸', 'SER': '丝氨酸',
            'THR': '苏氨酸', 'TRP': '色氨酸', 'TYR': '酪氨酸', 'VAL': '缬氨酸'
        }
        name = standard_names.get(amino_id, f"氨基酸-{amino_id}")
        
        # 提取特征
        features = []
        if 'c1ccccc1' in smiles or 'c1cc' in smiles:
            features.append('aromatic')
        if 'O' in smiles and ('OH' in smiles or 'O)' in smiles):
            features.append('hydroxyl')
        if 'S' in smiles:
            features.append('sulfur_containing')
        if 'F' in smiles:
            features.append('fluorinated')
        if 'Cl' in smiles:
            features.append('chlorinated')
        if 'Br' in smiles:
            features.append('brominated')
        if 'Se' in smiles:
            features.append('selenium_containing')
        if 'N' in smiles and 'C(=O)O' in smiles:
            features.append('amino_acid')
        
        # 插入数据库
        conn = sqlite3.connect('amino_acids.db')
        cursor = conn.cursor()
        
        cursor.execute('''
            INSERT OR REPLACE INTO amino_acids 
            (id, name, molecular_formula, molecular_weight, smiles, atom_composition, key_features)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (
            amino_id, name, formula, weight, smiles,
            json.dumps(atom_counts), json.dumps(features)
        ))
        
        conn.commit()
        conn.close()
        
        print(f"✅ {amino_id}: {name} - {formula}")
        return True
        
    except Exception as e:
        print(f"❌ {amino_id}: {e}")
        return False

def import_all_amino_acids():
    """导入所有可能的氨基酸"""
    print("🚀 最终导入：所有氨基酸数据")
    print("=" * 60)
    
    # 检查RDKit
    try:
        from rdkit import Chem
        print("✅ RDKit可用")
    except ImportError:
        print("❌ RDKit不可用")
        return
    
    # 扫描所有氨基酸
    structures_dir = Path("data/structures")
    all_amino_acids = []
    
    for subdir in structures_dir.iterdir():
        if subdir.is_dir():
            all_amino_acids.append(subdir.name)
    
    all_amino_acids = sorted(all_amino_acids)
    print(f"📁 发现 {len(all_amino_acids)} 种氨基酸结构")
    
    # 导入统计
    success_count = 0
    failed_count = 0
    manual_count = 0
    
    manual_fixes = get_manual_fixes()
    
    for amino_id in all_amino_acids:
        if amino_id in manual_fixes:
            manual_count += 1
        
        if import_amino_acid_final(amino_id):
            success_count += 1
        else:
            failed_count += 1
    
    # 最终统计
    print("\n" + "=" * 60)
    print("📊 最终导入结果:")
    print(f"   成功导入: {success_count} 种")
    print(f"   导入失败: {failed_count} 种")
    print(f"   手动修正: {manual_count} 种")
    print(f"   成功率: {success_count/(success_count+failed_count)*100:.1f}%")
    
    # 验证数据库
    try:
        conn = sqlite3.connect('amino_acids.db')
        cursor = conn.cursor()
        cursor.execute("SELECT COUNT(*) FROM amino_acids")
        total_count = cursor.fetchone()[0]
        conn.close()
        print(f"\n🎉 数据库最终包含 {total_count} 种氨基酸！")
        
        if total_count >= 200:
            print("🏆 恭喜！成功突破200种氨基酸大关！")
        
    except Exception as e:
        print(f"❌ 数据库验证失败: {e}")

if __name__ == "__main__":
    import_all_amino_acids()
