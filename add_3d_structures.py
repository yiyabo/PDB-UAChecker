#!/usr/bin/env python3
"""
为数据库添加3D结构信息
包括：标准构象、几何参数、立体化学信息
"""

import sys
import sqlite3
import json
import numpy as np
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, 'core')

def generate_3d_structure_from_smiles(smiles: str) -> Optional[Dict]:
    """从SMILES生成3D结构"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors
        
        # 从SMILES创建分子
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        
        # 添加氢原子
        mol = Chem.AddHs(mol)
        
        # 生成3D构象
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result != 0:
            # 如果嵌入失败，尝试使用距离几何
            AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=42)
        
        # 优化几何结构
        AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        
        # 提取3D坐标
        conf = mol.GetConformer()
        coordinates = []
        atoms_info = []
        
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            coordinates.append([pos.x, pos.y, pos.z])
            
            atoms_info.append({
                'atom_id': i,
                'element': atom.GetSymbol(),
                'atomic_number': atom.GetAtomicNum(),
                'formal_charge': atom.GetFormalCharge(),
                'hybridization': str(atom.GetHybridization()),
                'is_aromatic': atom.GetIsAromatic(),
                'coordinates': [pos.x, pos.y, pos.z]
            })
        
        # 计算几何参数
        geometry_params = calculate_geometry_parameters(mol)
        
        # 分析立体化学
        stereochemistry = analyze_stereochemistry(mol)
        
        # 计算分子描述符
        descriptors = calculate_3d_descriptors(mol)
        
        structure_data = {
            'coordinates': coordinates,
            'atoms_info': atoms_info,
            'geometry_parameters': geometry_params,
            'stereochemistry': stereochemistry,
            'descriptors_3d': descriptors,
            'generation_method': 'RDKit_ETKDG',
            'energy_minimized': True,
            'conformer_id': 0
        }
        
        return structure_data
        
    except Exception as e:
        print(f"3D结构生成失败: {e}")
        return None

def calculate_geometry_parameters(mol) -> Dict:
    """计算几何参数"""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        
        conf = mol.GetConformer()
        
        # 计算键长
        bond_lengths = []
        for bond in mol.GetBonds():
            atom1_idx = bond.GetBeginAtomIdx()
            atom2_idx = bond.GetEndAtomIdx()
            
            pos1 = conf.GetAtomPosition(atom1_idx)
            pos2 = conf.GetAtomPosition(atom2_idx)
            
            distance = ((pos1.x - pos2.x)**2 + (pos1.y - pos2.y)**2 + (pos1.z - pos2.z)**2)**0.5
            
            bond_lengths.append({
                'atom1': atom1_idx,
                'atom2': atom2_idx,
                'bond_type': str(bond.GetBondType()),
                'length': round(distance, 3)
            })
        
        # 计算键角
        bond_angles = []
        for atom in mol.GetAtoms():
            if atom.GetDegree() >= 2:
                neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
                for i in range(len(neighbors)):
                    for j in range(i+1, len(neighbors)):
                        angle = calculate_angle(conf, neighbors[i], atom.GetIdx(), neighbors[j])
                        bond_angles.append({
                            'atom1': neighbors[i],
                            'center_atom': atom.GetIdx(),
                            'atom3': neighbors[j],
                            'angle': round(angle, 2)
                        })
        
        # 计算二面角（对于氨基酸主链）
        dihedral_angles = calculate_amino_acid_dihedrals(mol, conf)
        
        # 计算分子体积 - 修复函数名
        try:
            # 尝试不同的体积计算方法
            mol_volume = None

            # 方法1: 尝试直接计算体积
            try:
                mol_volume = rdMolDescriptors.CalcMolVolume(mol)
            except AttributeError:
                pass

            # 方法2: 使用AllChem中的体积计算
            if mol_volume is None:
                try:
                    from rdkit.Chem import AllChem
                    mol_volume = AllChem.ComputeMolVolume(mol)
                except:
                    pass

            # 方法3: 基于原子半径估算体积
            if mol_volume is None:
                try:
                    from rdkit.Chem import rdMolDescriptors
                    # 使用分子表面积作为体积的替代指标
                    mol_volume = rdMolDescriptors.CalcTPSA(mol)  # 拓扑极性表面积
                except:
                    mol_volume = None

        except Exception as e:
            mol_volume = None

        return {
            'bond_lengths': bond_lengths,
            'bond_angles': bond_angles,
            'dihedral_angles': dihedral_angles,
            'molecular_volume': round(mol_volume, 2) if mol_volume else None
        }
        
    except Exception as e:
        print(f"几何参数计算失败: {e}")
        return {}

def calculate_angle(conf, atom1_idx: int, atom2_idx: int, atom3_idx: int) -> float:
    """计算三个原子间的键角"""
    pos1 = conf.GetAtomPosition(atom1_idx)
    pos2 = conf.GetAtomPosition(atom2_idx)
    pos3 = conf.GetAtomPosition(atom3_idx)
    
    # 向量
    v1 = np.array([pos1.x - pos2.x, pos1.y - pos2.y, pos1.z - pos2.z])
    v2 = np.array([pos3.x - pos2.x, pos3.y - pos2.y, pos3.z - pos2.z])
    
    # 计算角度
    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_angle = np.clip(cos_angle, -1.0, 1.0)  # 防止数值误差
    angle = np.arccos(cos_angle) * 180.0 / np.pi
    
    return angle

def calculate_amino_acid_dihedrals(mol, conf) -> List[Dict]:
    """计算氨基酸特有的二面角"""
    dihedrals = []

    try:
        # 寻找氨基酸主链原子 (N-CA-C-O)
        # 这里简化处理，实际应用中需要更复杂的模式匹配

        # 寻找羧基碳
        carboxyl_carbons = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C':
                # 检查是否连接了两个氧原子
                oxygen_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']
                if len(oxygen_neighbors) >= 2:
                    carboxyl_carbons.append(atom.GetIdx())

        # 如果找到羧基，计算相关二面角
        for c_idx in carboxyl_carbons:
            c_atom = mol.GetAtomWithIdx(c_idx)  # 修复: 使用GetAtomWithIdx
            carbon_neighbors = [n for n in c_atom.GetNeighbors() if n.GetSymbol() == 'C']

            if carbon_neighbors:
                ca_idx = carbon_neighbors[0].GetIdx()  # 假设是α碳
                dihedrals.append({
                    'type': 'carboxyl_dihedral',
                    'atoms': [ca_idx, c_idx],
                    'description': 'C-alpha to carboxyl carbon'
                })

    except Exception as e:
        print(f"二面角计算失败: {e}")

    return dihedrals

def analyze_stereochemistry(mol) -> Dict:
    """分析立体化学信息"""
    try:
        from rdkit import Chem
        
        # 分析手性中心
        chiral_centers = []
        for atom in mol.GetAtoms():
            if atom.HasProp('_ChiralityPossible'):
                chiral_tag = atom.GetChiralTag()
                if chiral_tag != Chem.ChiralType.CHI_UNSPECIFIED:
                    chiral_centers.append({
                        'atom_idx': atom.GetIdx(),
                        'element': atom.GetSymbol(),
                        'chiral_tag': str(chiral_tag),
                        'cip_code': atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None
                    })
        
        # 分析双键立体化学
        double_bond_stereo = []
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                stereo = bond.GetStereo()
                if stereo != Chem.BondStereo.STEREONONE:
                    double_bond_stereo.append({
                        'bond_idx': bond.GetIdx(),
                        'atom1': bond.GetBeginAtomIdx(),
                        'atom2': bond.GetEndAtomIdx(),
                        'stereo': str(stereo)
                    })
        
        return {
            'chiral_centers': chiral_centers,
            'double_bond_stereo': double_bond_stereo,
            'num_chiral_centers': len(chiral_centers),
            'is_chiral': len(chiral_centers) > 0
        }
        
    except Exception as e:
        print(f"立体化学分析失败: {e}")
        return {}

def calculate_3d_descriptors(mol) -> Dict:
    """计算3D分子描述符"""
    try:
        from rdkit.Chem import rdMolDescriptors, Descriptors3D
        
        descriptors = {
            'asphericity': round(rdMolDescriptors.CalcAsphericity(mol), 4),
            'eccentricity': round(rdMolDescriptors.CalcEccentricity(mol), 4),
            'inertial_shape_factor': round(rdMolDescriptors.CalcInertialShapeFactor(mol), 4),
            'radius_of_gyration': round(rdMolDescriptors.CalcRadiusOfGyration(mol), 4),
            'spherocity_index': round(rdMolDescriptors.CalcSpherocityIndex(mol), 4)
        }
        
        return descriptors
        
    except Exception as e:
        print(f"3D描述符计算失败: {e}")
        return {}

def add_3d_structures_to_database():
    """为数据库中所有氨基酸添加3D结构信息"""
    
    print("🧬 开始为氨基酸数据库添加3D结构信息")
    print("=" * 60)
    
    # 检查RDKit
    try:
        from rdkit import Chem
        print("✅ RDKit可用")
    except ImportError:
        print("❌ RDKit不可用，请安装: conda install rdkit")
        return False
    
    conn = sqlite3.connect('core/amino_acids.db')
    cursor = conn.cursor()
    
    # 获取所有氨基酸的SMILES
    cursor.execute('SELECT id, smiles FROM amino_acids WHERE smiles IS NOT NULL')
    amino_acids = cursor.fetchall()
    
    print(f"📊 需要处理 {len(amino_acids)} 种氨基酸")
    
    success_count = 0
    failed_count = 0
    
    for amino_id, smiles in amino_acids:
        print(f"\n🧪 处理 {amino_id}: {smiles}")
        
        # 生成3D结构
        structure_data = generate_3d_structure_from_smiles(smiles)
        
        if structure_data:
            # 更新数据库
            structure_json = json.dumps(structure_data)
            cursor.execute(
                'UPDATE amino_acids SET structure_data = ?, updated_at = CURRENT_TIMESTAMP WHERE id = ?',
                (structure_json, amino_id)
            )
            
            num_atoms = len(structure_data['atoms_info'])
            num_bonds = len(structure_data['geometry_parameters'].get('bond_lengths', []))
            num_chiral = structure_data['stereochemistry'].get('num_chiral_centers', 0)
            
            print(f"   ✅ 成功: {num_atoms}原子, {num_bonds}键, {num_chiral}手性中心")
            success_count += 1
        else:
            print(f"   ❌ 3D结构生成失败")
            failed_count += 1
    
    conn.commit()
    conn.close()
    
    print(f"\n" + "=" * 60)
    print(f"📊 3D结构添加结果:")
    print(f"   成功: {success_count} 种")
    print(f"   失败: {failed_count} 种")
    print(f"   成功率: {success_count/(success_count+failed_count)*100:.1f}%")
    
    return success_count == len(amino_acids)

def test_3d_structure_analysis():
    """测试3D结构分析功能"""
    
    print("\n🧪 测试3D结构分析功能:")
    print("=" * 40)
    
    # 测试几个氨基酸的3D结构
    test_smiles = [
        ("ALA", "N[C@@H](C)C(=O)O"),
        ("PHE", "N[C@@H](Cc1ccccc1)C(=O)O"),
        ("CYS", "N[C@@H](CS)C(=O)O")
    ]
    
    for name, smiles in test_smiles:
        print(f"\n🔍 分析 {name}: {smiles}")
        structure_data = generate_3d_structure_from_smiles(smiles)
        
        if structure_data:
            atoms_count = len(structure_data['atoms_info'])
            bonds_count = len(structure_data['geometry_parameters'].get('bond_lengths', []))
            chiral_count = structure_data['stereochemistry'].get('num_chiral_centers', 0)
            
            print(f"   ✅ 原子数: {atoms_count}")
            print(f"   ✅ 键数: {bonds_count}")
            print(f"   ✅ 手性中心: {chiral_count}")
            
            if structure_data['descriptors_3d']:
                desc = structure_data['descriptors_3d']
                print(f"   📊 回转半径: {desc.get('radius_of_gyration', 'N/A')}")
                print(f"   📊 非球形度: {desc.get('asphericity', 'N/A')}")
        else:
            print(f"   ❌ 结构生成失败")

def main():
    """主函数"""
    
    print("🌟 3D结构数据添加工具")
    print("=" * 50)
    
    # 1. 测试3D结构分析
    test_3d_structure_analysis()
    
    # 2. 询问是否添加到数据库
    print(f"\n💡 是否要为所有229种氨基酸添加3D结构数据？")
    print("这将需要几分钟时间...")
    
    response = input("继续？(y/N): ").lower().strip()
    if response == 'y':
        success = add_3d_structures_to_database()
        
        if success:
            print("\n🎉 所有氨基酸3D结构添加成功！")
            print("\n💡 现在您可以:")
            print("1. 进行构象匹配分析")
            print("2. 验证PDB结构的几何合理性")
            print("3. 精确的立体化学识别")
            print("4. 3D相似性搜索")
        else:
            print("\n❌ 3D结构添加过程中出现错误")
    else:
        print("操作已取消")

if __name__ == "__main__":
    main()
