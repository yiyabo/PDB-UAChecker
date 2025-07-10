#!/usr/bin/env python3
"""
调试同分异构体识别问题
"""

from isomer_identifier import IsomerIdentifier
from rdkit import Chem

def debug_isomer_identification():
    # 测试异亮氨酸 vs 亮氨酸（正确的SMILES）
    smiles1 = "CC[C@H](C)[C@@H](N)C(=O)O"  # 异亮氨酸（正确）
    smiles2 = "CC(C)C[C@@H](N)C(=O)O"      # 亮氨酸
    
    print("调试异亮氨酸 vs 亮氨酸识别问题")
    print(f"异亮氨酸 SMILES: {smiles1}")
    print(f"亮氨酸 SMILES: {smiles2}")
    
    # 检查RDKit标准化
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 and mol2:
        canonical1 = Chem.MolToSmiles(mol1)
        canonical2 = Chem.MolToSmiles(mol2)
        
        print(f"标准化后的异亮氨酸: {canonical1}")
        print(f"标准化后的亮氨酸: {canonical2}")
        print(f"标准化后是否相同: {canonical1 == canonical2}")
        
        # 检查分子式
        from rdkit.Chem import rdMolDescriptors
        formula1 = rdMolDescriptors.CalcMolFormula(mol1)
        formula2 = rdMolDescriptors.CalcMolFormula(mol2)
        
        print(f"异亮氨酸分子式: {formula1}")
        print(f"亮氨酸分子式: {formula2}")
    
    # 使用我们的算法分析
    identifier = IsomerIdentifier()
    result = identifier.analyze_isomers(smiles1, smiles2)
    
    print(f"\n我们的算法结果:")
    print(f"异构体类型: {result.isomer_type}")
    print(f"相似性分数: {result.similarity_score}")
    print(f"结构差异: {result.structural_differences}")
    
    # 生成指纹详细信息
    fp1 = identifier._get_fingerprint(smiles1)
    fp2 = identifier._get_fingerprint(smiles2)
    
    print(f"\n指纹分析:")
    print(f"拓扑哈希相同: {fp1.topological_hash == fp2.topological_hash}")
    print(f"立体化学相同: {fp1.stereochemistry == fp2.stereochemistry}")
    print(f"结构关键词1: {fp1.structural_keys}")
    print(f"结构关键词2: {fp2.structural_keys}")

if __name__ == "__main__":
    debug_isomer_identification()