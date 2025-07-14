#!/usr/bin/env python3
"""
为数据库中的氨基酸添加分子指纹数据
"""

import sqlite3
import json
import sys
from typing import Dict, List, Optional

def calculate_fingerprints(smiles: str) -> Optional[Dict]:
    """计算多种类型的分子指纹"""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        from rdkit.Chem import MACCSkeys
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        
        fingerprints = {}
        
        # 1. ECFP (Extended Connectivity Fingerprint) - 最常用
        ecfp2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        ecfp4 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 3, nBits=1024)
        
        # 2. MACCS Keys - 药物设计标准
        maccs = MACCSkeys.GenMACCSKeys(mol)
        
        # 3. Topological Fingerprint
        topo = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=1024)
        
        # 4. Atom Pair Fingerprint
        atom_pair = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=1024)
        
        # 转换为字符串存储
        fingerprints = {
            'ecfp2': ecfp2.ToBitString(),
            'ecfp4': ecfp4.ToBitString(), 
            'maccs': maccs.ToBitString(),
            'topological': topo.ToBitString(),
            'atom_pair': atom_pair.ToBitString(),
            'fingerprint_info': {
                'ecfp2_bits': 1024,
                'ecfp4_bits': 1024,
                'maccs_bits': 167,
                'topological_bits': 1024,
                'atom_pair_bits': 1024,
                'generated_by': 'RDKit',
                'version': '2023.09'
            }
        }
        
        return fingerprints
        
    except Exception as e:
        print(f"计算指纹失败: {e}")
        return None

def add_fingerprints_to_database():
    """为数据库中所有氨基酸添加指纹数据"""
    
    print("🧬 开始为氨基酸数据库添加分子指纹")
    print("=" * 60)
    
    # 检查RDKit
    try:
        from rdkit import Chem
        print("✅ RDKit可用")
    except ImportError:
        print("❌ RDKit不可用，请安装: conda install rdkit")
        return False
    
    conn = sqlite3.connect('amino_acids.db')
    cursor = conn.cursor()
    
    # 获取所有氨基酸的SMILES
    cursor.execute('SELECT id, smiles FROM amino_acids WHERE smiles IS NOT NULL')
    amino_acids = cursor.fetchall()
    
    print(f"📊 需要处理 {len(amino_acids)} 种氨基酸")
    
    success_count = 0
    failed_count = 0
    
    for amino_id, smiles in amino_acids:
        print(f"\n🧪 处理 {amino_id}: {smiles}")
        
        # 计算指纹
        fingerprints = calculate_fingerprints(smiles)
        
        if fingerprints:
            # 更新数据库
            fingerprints_json = json.dumps(fingerprints)
            cursor.execute(
                'UPDATE amino_acids SET fingerprints = ?, updated_at = CURRENT_TIMESTAMP WHERE id = ?',
                (fingerprints_json, amino_id)
            )
            
            print(f"   ✅ 成功计算 {len(fingerprints)-1} 种指纹")
            success_count += 1
        else:
            print(f"   ❌ 指纹计算失败")
            failed_count += 1
    
    conn.commit()
    conn.close()
    
    print(f"\n" + "=" * 60)
    print(f"📊 指纹添加结果:")
    print(f"   成功: {success_count} 种")
    print(f"   失败: {failed_count} 种")
    print(f"   成功率: {success_count/(success_count+failed_count)*100:.1f}%")
    
    return success_count == len(amino_acids)

def test_fingerprint_similarity():
    """测试指纹相似性计算"""
    
    print("\n🧪 测试指纹相似性计算:")
    print("=" * 40)
    
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        from rdkit import DataStructs
        
        # 测试两个相似的氨基酸
        smiles1 = "N[C@@H](Cc1ccccc1)C(=O)O"  # 苯丙氨酸
        smiles2 = "N[C@@H](Cc1ccc(O)cc1)C(=O)O"  # 酪氨酸
        
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        
        if mol1 and mol2:
            # 计算ECFP指纹
            fp1 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol1, 2)
            fp2 = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol2, 2)
            
            # 计算Tanimoto相似性
            similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
            
            print(f"苯丙氨酸 vs 酪氨酸:")
            print(f"  SMILES1: {smiles1}")
            print(f"  SMILES2: {smiles2}")
            print(f"  Tanimoto相似性: {similarity:.3f}")
            
            if similarity > 0.7:
                print("  ✅ 高相似性 - 指纹计算正常")
            else:
                print("  ⚠️ 相似性较低 - 需要检查")
        
    except Exception as e:
        print(f"❌ 相似性测试失败: {e}")

def verify_fingerprints_in_database():
    """验证数据库中的指纹数据"""
    
    print("\n🔍 验证数据库中的指纹数据:")
    print("=" * 40)
    
    conn = sqlite3.connect('amino_acids.db')
    cursor = conn.cursor()
    
    # 检查指纹数据
    cursor.execute('SELECT COUNT(*) FROM amino_acids WHERE fingerprints IS NOT NULL')
    fingerprint_count = cursor.fetchone()[0]
    
    cursor.execute('SELECT COUNT(*) FROM amino_acids')
    total_count = cursor.fetchone()[0]
    
    print(f"指纹数据覆盖率: {fingerprint_count}/{total_count} ({fingerprint_count/total_count*100:.1f}%)")
    
    # 查看示例指纹数据
    cursor.execute('SELECT id, fingerprints FROM amino_acids WHERE fingerprints IS NOT NULL LIMIT 2')
    samples = cursor.fetchall()
    
    for amino_id, fingerprints_json in samples:
        fingerprints = json.loads(fingerprints_json)
        info = fingerprints.get('fingerprint_info', {})
        
        print(f"\n{amino_id} 指纹信息:")
        print(f"  ECFP2长度: {len(fingerprints.get('ecfp2', ''))} bits")
        print(f"  MACCS长度: {len(fingerprints.get('maccs', ''))} bits")
        print(f"  生成工具: {info.get('generated_by', 'Unknown')}")
    
    conn.close()

def main():
    """主函数"""
    
    print("🚀 分子指纹数据添加工具")
    print("=" * 50)
    
    # 1. 添加指纹数据
    success = add_fingerprints_to_database()
    
    if success:
        print("\n🎉 所有氨基酸指纹添加成功！")
        
        # 2. 测试相似性计算
        test_fingerprint_similarity()
        
        # 3. 验证数据库数据
        verify_fingerprints_in_database()
        
        print(f"\n💡 下一步可以做:")
        print("1. 实现基于指纹的相似性搜索")
        print("2. 添加氨基酸聚类分析")
        print("3. 构建推荐系统")
        print("4. 优化搜索引擎性能")
        
    else:
        print("\n❌ 指纹添加过程中出现错误")

if __name__ == "__main__":
    main()
