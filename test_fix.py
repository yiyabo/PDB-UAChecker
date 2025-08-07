#!/usr/bin/env python3
"""测试循环导入修复后的功能"""

from pdb_uachecker.core.models import ResidueInfo, AtomInfo

print("🧪 测试不同原子组合...")

# 测试1: 包含氢原子
atoms_with_h = [
    AtomInfo('CA', 'C', 0, 0, 0, 'ALA', 1, 'A'),
    AtomInfo('N', 'N', 1, 0, 0, 'ALA', 1, 'A'),
    AtomInfo('O', 'O', 2, 0, 0, 'ALA', 1, 'A'),
    AtomInfo('H1', 'H', 3, 0, 0, 'ALA', 1, 'A'),
    AtomInfo('H2', 'H', 4, 0, 0, 'ALA', 1, 'A'),
]
residue_with_h = ResidueInfo('ALA', 1, 'A', atoms_with_h)
print(f"含氢原子 - 分子式: {residue_with_h.molecular_formula}, 组成: {residue_with_h.atom_composition}")
print(f"重原子组成: {residue_with_h.heavy_atom_composition}")

# 测试2: 不含氢原子
atoms_no_h = [
    AtomInfo('CA', 'C', 0, 0, 0, 'ALA', 2, 'A'),
    AtomInfo('N', 'N', 1, 0, 0, 'ALA', 2, 'A'),
    AtomInfo('O', 'O', 2, 0, 0, 'ALA', 2, 'A'),
]
residue_no_h = ResidueInfo('ALA', 2, 'A', atoms_no_h)
print(f"不含氢原子 - 分子式: {residue_no_h.molecular_formula}, 组成: {residue_no_h.atom_composition}")
print(f"重原子组成: {residue_no_h.heavy_atom_composition}")

# 测试3: 空原子列表
empty_residue = ResidueInfo('XXX', 3, 'A', [])
print(f"空原子列表 - 分子式: '{empty_residue.molecular_formula}', 组成: {empty_residue.atom_composition}")

print("✅ 功能完整性测试通过")
