#!/usr/bin/env python3
"""
测试分子式计算修复的效果
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.core.models import ResidueInfo, AtomInfo
from pdb_uachecker.core.parser import PDBParser


def test_residue_info_molecular_formula():
    """测试ResidueInfo的分子式计算"""
    print("🧪 测试ResidueInfo分子式计算...")
    
    # 创建测试原子（2AG的结构）
    atoms = [
        AtomInfo('N', 'N', 1.584, 0.063, -0.041, '2AG', 1, 'A'),
        AtomInfo('CA', 'C', 0.238, -0.494, -0.123, '2AG', 1, 'A'),
        AtomInfo('C', 'C', -0.010, -0.780, 1.266, '2AG', 1, 'A'),
        AtomInfo('O', 'O', 0.908, -0.875, 2.056, '2AG', 1, 'A'),
        AtomInfo('CB', 'C', 0.174, -1.720, -1.083, '2AG', 1, 'A'),
        AtomInfo('CG', 'C', 0.327, -1.324, -2.477, '2AG', 1, 'A'),
        AtomInfo('CD', 'C', 1.338, -1.711, -3.264, '2AG', 1, 'A'),
        AtomInfo('OXT', 'O', -1.144, -0.699, 1.708, '2AG', 1, 'A'),
        AtomInfo('HN1', 'H', 1.901, 0.316, -0.966, '2AG', 1, 'A'),
    ]
    
    # 创建ResidueInfo
    residue = ResidueInfo('2AG', 1, 'A', atoms)
    
    # 验证结果
    expected_formula = 'C5HNO2'  # 只有9个原子的简化版本
    expected_composition = {'N': 1, 'C': 5, 'O': 2, 'H': 1}
    
    print(f"   分子式: {residue.molecular_formula}")
    print(f"   原子组成: {residue.atom_composition}")
    
    assert residue.molecular_formula == expected_formula, f"分子式不匹配: 期望{expected_formula}, 实际{residue.molecular_formula}"
    assert residue.atom_composition == expected_composition, f"原子组成不匹配: 期望{expected_composition}, 实际{residue.atom_composition}"
    
    print("   ✅ ResidueInfo分子式计算测试通过")


def test_pdb_parser_molecular_formula():
    """测试PDB解析器的分子式计算修复"""
    print("🧪 测试PDB解析器分子式计算修复...")
    
    parser = PDBParser()
    residues = parser.parse_pdb_file('data/structures/2AG/2AG.pdb')
    
    assert len(residues) == 1, f"应该解析出1个残基，实际{len(residues)}个"
    
    residue = residues[0]
    expected_formula = 'C5H9NO2'
    expected_composition = {'N': 1, 'C': 5, 'O': 2, 'H': 9}
    
    print(f"   分子式: {residue.molecular_formula}")
    print(f"   原子组成: {residue.atom_composition}")
    
    assert residue.molecular_formula == expected_formula, f"分子式不匹配: 期望{expected_formula}, 实际{residue.molecular_formula}"
    assert residue.atom_composition == expected_composition, f"原子组成不匹配: 期望{expected_composition}, 实际{residue.atom_composition}"
    
    print("   ✅ PDB解析器分子式计算修复测试通过")


def test_heavy_atom_composition():
    """测试重原子组成计算"""
    print("🧪 测试重原子组成计算...")
    
    parser = PDBParser()
    residues = parser.parse_pdb_file('data/structures/2AG/2AG.pdb')
    residue = residues[0]
    
    expected_heavy_composition = {'N': 1, 'C': 5, 'O': 2}
    heavy_composition = residue.heavy_atom_composition
    
    print(f"   重原子组成: {heavy_composition}")
    
    assert heavy_composition == expected_heavy_composition, f"重原子组成不匹配: 期望{expected_heavy_composition}, 实际{heavy_composition}"
    
    print("   ✅ 重原子组成计算测试通过")


if __name__ == "__main__":
    print("🚀 开始分子式计算修复测试")
    print("=" * 50)
    
    try:
        test_residue_info_molecular_formula()
        test_pdb_parser_molecular_formula()
        test_heavy_atom_composition()
        
        print("=" * 50)
        print("🎉 所有测试通过！分子式计算修复成功！")
        
    except Exception as e:
        print(f"❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)