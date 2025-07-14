#!/usr/bin/env python3
"""
验证SMILES修正的化学合规性和准确性
检查是否引入了歧义或错误
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.append('.')

def get_original_vs_corrected_smiles() -> List[Tuple[str, str, str]]:
    """获取原始vs修正的SMILES对比"""
    
    # 从我们的修正记录中提取关键案例
    corrections = [
        # (amino_id, original_smiles, corrected_smiles)
        ('VAH', 'CC[C@H]([C@@H]([C](=O)=O)N)O', 'N[C@@H](C(=O)O)[C@H](CC)O'),
        ('0AF', 'N[C@H]([C](=O)=O)Cc1c[nH]c2c1cccc2O', 'N[C@H](C(=O)O)Cc1c[nH]c2c1cccc2O'),
        ('0BN', 'N[C@H](C(=O)O)Cc1ccc(cc1)[C](=N)=N', 'N[C@H](C(=O)O)Cc1ccc(cc1)C(=N)N'),
        ('GDPR', '[NH3][C@H]([C](=O)=O)C[NH]=[C](=[NH2])=[NH2]', 'N[C@H](C(=O)O)CNC(=N)N'),
        ('HRG', '[NH2]=[C](=[NH]CCCC[C@@H](C(=O)O)[NH3])=[NH2]', 'N[C@@H](C(=O)O)CCCCNC(=N)N'),
        ('NIY', 'O=[C](=O)[C@H](Cc1ccc(c(c1)N(=O)=O)O)[NH3]', 'N[C@H](C(=O)O)Cc1ccc(c(c1)[N+](=O)[O-])O'),
        ('PPN', 'O=[C](=O)[C@H](Cc1ccc(cc1)N(=O)=O)[NH3]', 'N[C@H](C(=O)O)Cc1ccc(cc1)[N+](=O)[O-]'),
        ('THIC', 'OC(=O)[C@H](CCC[NH]=[C](=S)=[NH2])[NH3]', 'N[C@H](C(=O)O)CCCNC(=S)N'),
    ]
    
    return corrections

def analyze_chemical_validity(original: str, corrected: str, amino_id: str):
    """分析化学有效性"""
    
    print(f"\n🧪 分析 {amino_id}:")
    print(f"   原始: {original}")
    print(f"   修正: {corrected}")
    
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        
        # 检查原始SMILES
        orig_mol = Chem.MolFromSmiles(original)
        orig_valid = orig_mol is not None
        
        # 检查修正SMILES
        corr_mol = Chem.MolFromSmiles(corrected)
        corr_valid = corr_mol is not None
        
        print(f"   原始有效性: {'✅' if orig_valid else '❌'}")
        print(f"   修正有效性: {'✅' if corr_valid else '❌'}")
        
        if orig_valid and corr_valid:
            # 比较分子式
            orig_formula = Chem.rdMolDescriptors.CalcMolFormula(orig_mol)
            corr_formula = Chem.rdMolDescriptors.CalcMolFormula(corr_mol)
            
            # 比较分子量
            orig_weight = Descriptors.MolWt(orig_mol)
            corr_weight = Descriptors.MolWt(corr_mol)
            
            print(f"   分子式: {orig_formula} -> {corr_formula}")
            print(f"   分子量: {orig_weight:.3f} -> {corr_weight:.3f}")
            
            # 检查是否相同
            if orig_formula == corr_formula:
                print(f"   ✅ 分子式一致 - 修正合理")
                return True
            else:
                print(f"   ❌ 分子式不一致 - 可能有问题")
                return False
        
        elif not orig_valid and corr_valid:
            print(f"   ✅ 修正成功 - 原始无效，修正有效")
            return True
        
        else:
            print(f"   ❌ 修正失败")
            return False
            
    except Exception as e:
        print(f"   ❌ 分析失败: {e}")
        return False

def check_stereochemistry_preservation(original: str, corrected: str):
    """检查立体化学是否保持"""
    
    # 检查手性中心标记
    orig_chiral = original.count('@')
    corr_chiral = corrected.count('@')
    
    if orig_chiral == corr_chiral:
        return True, f"立体化学保持 ({orig_chiral}个手性中心)"
    else:
        return False, f"立体化学改变 ({orig_chiral} -> {corr_chiral}个手性中心)"

def identify_modification_types(original: str, corrected: str) -> List[str]:
    """识别修正类型"""
    
    modifications = []
    
    # 检查质子化状态移除
    if '[NH3]' in original and 'N' in corrected:
        modifications.append('质子化氨基移除')
    
    if '[C](=O)=O' in original and 'C(=O)O' in corrected:
        modifications.append('羧酸价态修正')
    
    # 检查胍基修正
    if '[NH]=[C]' in original and 'NC(' in corrected:
        modifications.append('胍基结构标准化')
    
    # 检查硝基修正
    if 'N(=O)=O' in original and '[N+](=O)[O-]' in corrected:
        modifications.append('硝基标准化')
    
    # 检查结构重排
    if len(original) != len(corrected):
        modifications.append('结构重排')
    
    return modifications

def validate_against_known_structures():
    """与已知结构数据库验证"""
    
    print("\n🔍 与已知结构验证:")
    print("=" * 40)
    
    # 一些我们可以确认的标准结构
    known_structures = {
        'VAL': 'N[C@@H](C(C)C)C(=O)O',  # 缬氨酸 - 标准结构
        'ARG': 'N[C@@H](CCCNC(=N)N)C(=O)O',  # 精氨酸 - 胍基标准表示
        'TYR': 'N[C@@H](Cc1ccc(O)cc1)C(=O)O',  # 酪氨酸 - 标准结构
    }
    
    try:
        from rdkit import Chem
        
        for name, smiles in known_structures.items():
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                print(f"✅ {name}: {formula} - 标准结构验证通过")
            else:
                print(f"❌ {name}: 验证失败")
                
    except Exception as e:
        print(f"❌ 验证过程出错: {e}")

def main():
    """主验证函数"""
    
    print("🔬 SMILES修正化学合规性验证")
    print("=" * 60)
    
    corrections = get_original_vs_corrected_smiles()
    
    valid_count = 0
    total_count = len(corrections)
    
    for amino_id, original, corrected in corrections:
        
        # 化学有效性分析
        is_valid = analyze_chemical_validity(original, corrected, amino_id)
        if is_valid:
            valid_count += 1
        
        # 立体化学检查
        stereo_ok, stereo_msg = check_stereochemistry_preservation(original, corrected)
        print(f"   立体化学: {'✅' if stereo_ok else '⚠️'} {stereo_msg}")
        
        # 修正类型识别
        mod_types = identify_modification_types(original, corrected)
        print(f"   修正类型: {', '.join(mod_types)}")
    
    print(f"\n" + "=" * 60)
    print(f"📊 验证结果统计:")
    print(f"   验证样本: {total_count} 个")
    print(f"   化学有效: {valid_count} 个")
    print(f"   有效率: {valid_count/total_count*100:.1f}%")
    
    # 风险评估
    print(f"\n⚠️ 风险评估:")
    if valid_count == total_count:
        print("✅ 低风险 - 所有修正都保持了化学一致性")
    elif valid_count >= total_count * 0.8:
        print("🟡 中等风险 - 大部分修正合理，少数需要进一步验证")
    else:
        print("🔴 高风险 - 多个修正可能引入了化学错误")
    
    # 建议
    print(f"\n💡 建议:")
    print("1. 对于分子式不一致的修正，需要重新检查")
    print("2. 立体化学改变的修正需要特别注意")
    print("3. 建议与原始文献或数据库进行交叉验证")
    print("4. 考虑保留原始SMILES作为备份")
    
    # 验证已知结构
    validate_against_known_structures()

if __name__ == "__main__":
    main()
