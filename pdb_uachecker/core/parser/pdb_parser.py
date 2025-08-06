"""
PDB文件解析器
解析PDB文件，提取残基和原子信息
"""

import re
from typing import List, Dict, Optional, Tuple
from pathlib import Path

from ..models import ResidueInfo, AtomInfo, ParsingError
from ...utils.chemistry import ChemistryUtils


class PDBParser:
    """PDB文件解析器"""
    
    def __init__(self):
        self.chemistry_utils = ChemistryUtils()
    
    def parse_pdb_file(self, pdb_file: str) -> List[ResidueInfo]:
        """
        解析PDB文件，提取残基信息
        
        Args:
            pdb_file: PDB文件路径
        
        Returns:
            残基信息列表
        
        Raises:
            ParsingError: 解析失败时抛出
        """
        if not Path(pdb_file).exists():
            raise ParsingError(f"PDB文件不存在: {pdb_file}")
        
        print(f"📄 解析PDB文件: {pdb_file}")
        
        residues = {}
        line_count = 0
        atom_count = 0
        
        try:
            with open(pdb_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line_count += 1
                    
                    # 解析ATOM和HETATM行
                    if line.startswith(('ATOM', 'HETATM')):
                        atom_info = self._parse_atom_line(line)
                        if atom_info:
                            atom_count += 1
                            residue_key = (
                                atom_info.chain_id,
                                atom_info.residue_name,
                                atom_info.residue_number
                            )
                            
                            if residue_key not in residues:
                                residues[residue_key] = ResidueInfo(
                                    residue_name=atom_info.residue_name,
                                    residue_number=atom_info.residue_number,
                                    chain_id=atom_info.chain_id,
                                    atoms=[]
                                )
                            
                            residues[residue_key].atoms.append(atom_info)
        
        except Exception as e:
            raise ParsingError(f"PDB文件解析失败 (行 {line_count}): {e}")
        
        residue_list = list(residues.values())
        
        # 验证解析结果
        if not residue_list:
            raise ParsingError("未找到有效的残基信息")
        
        print(f"✅ 解析完成: {len(residue_list)} 个残基, {atom_count} 个原子")
        
        # 后处理：计算分子式和原子组成
        self._post_process_residues(residue_list)
        
        return residue_list
    
    def _parse_atom_line(self, line: str) -> Optional[AtomInfo]:
        """
        解析PDB原子行
        
        Args:
            line: PDB原子行
        
        Returns:
            原子信息，解析失败返回None
        """
        if len(line) < 54:
            return None
        
        try:
            # 提取字段（按PDB格式标准）
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain_id = line[21:22].strip() or 'A'  # 默认链ID
            residue_number = int(line[22:26].strip())
            
            # 坐标
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            
            # 元素符号（优先使用77-78列，否则从原子名推断）
            element = ""
            if len(line) > 77:
                element = line[76:78].strip()
            
            if not element:
                element = self._extract_element_from_atom_name(atom_name)
            
            # 验证元素符号
            if not element or not self.chemistry_utils.is_valid_element(element):
                element = self.chemistry_utils.normalize_atom_name(atom_name)
            
            if not element:
                return None  # 无法确定元素类型
            
            return AtomInfo(
                atom_name=atom_name,
                element=element,
                x=x, y=y, z=z,
                residue_name=residue_name,
                residue_number=residue_number,
                chain_id=chain_id
            )
        
        except (ValueError, IndexError) as e:
            # 解析失败，跳过该原子
            return None
    
    def _extract_element_from_atom_name(self, atom_name: str) -> str:
        """
        从原子名称提取元素符号
        
        Args:
            atom_name: 原子名称，如"CA", "CB1", "N1"
        
        Returns:
            元素符号
        """
        if not atom_name:
            return ""
        
        atom_name = atom_name.strip()
        
        # 使用ChemistryUtils的标准化方法
        return self.chemistry_utils.normalize_atom_name(atom_name)
    
    def _post_process_residues(self, residues: List[ResidueInfo]):
        """
        后处理残基信息，计算分子式和原子组成
        
        Args:
            residues: 残基信息列表
        """
        for residue in residues:
            # 验证原子数量
            if len(residue.atoms) == 0:
                print(f"⚠️ 残基 {residue.residue_key} 没有原子")
                continue
            
            # 验证分子式
            if not residue.molecular_formula:
                print(f"⚠️ 残基 {residue.residue_key} 分子式计算失败，尝试修复")
                # 尝试手动计算并修复
                from ...utils.chemistry import ChemistryUtils
                atom_dicts = [
                    {'element': atom.element}
                    for atom in residue.atoms if atom.element
                ]
                if atom_dicts:
                    composition = ChemistryUtils.calculate_atom_composition(atom_dicts)
                    formula = ChemistryUtils.calculate_molecular_formula(composition)
                    # 修复分子式和原子组成
                    residue.molecular_formula = formula
                    residue.atom_composition = composition
                    print(f"   ✅ 修复成功: 组成={composition}, 分子式={formula}")
                else:
                    print(f"   ❌ 修复失败: 没有有效原子")
            
            # 统计信息（调试用）
            element_counts = {}
            for atom in residue.atoms:
                if atom.element:
                    element_counts[atom.element] = element_counts.get(atom.element, 0) + 1
    
    def parse_residue_from_atoms(self, atoms: List[Dict]) -> ResidueInfo:
        """
        从原子信息列表创建残基信息
        
        Args:
            atoms: 原子信息字典列表
        
        Returns:
            残基信息
        """
        if not atoms:
            raise ParsingError("原子列表为空")
        
        # 从第一个原子获取残基基本信息
        first_atom = atoms[0]
        residue_name = first_atom.get('residue_name', 'UNK')
        residue_number = first_atom.get('residue_number', 1)
        chain_id = first_atom.get('chain_id', 'A')
        
        # 转换为AtomInfo对象
        atom_infos = []
        for atom_dict in atoms:
            atom_info = AtomInfo(
                atom_name=atom_dict.get('atom_name', ''),
                element=atom_dict.get('element', ''),
                x=atom_dict.get('x', 0.0),
                y=atom_dict.get('y', 0.0),
                z=atom_dict.get('z', 0.0),
                residue_name=residue_name,
                residue_number=residue_number,
                chain_id=chain_id
            )
            atom_infos.append(atom_info)
        
        return ResidueInfo(
            residue_name=residue_name,
            residue_number=residue_number,
            chain_id=chain_id,
            atoms=atom_infos
        )
    
    def validate_pdb_file(self, pdb_file: str) -> Tuple[bool, List[str]]:
        """
        验证PDB文件格式
        
        Args:
            pdb_file: PDB文件路径
        
        Returns:
            (是否有效, 错误信息列表)
        """
        errors = []
        
        if not Path(pdb_file).exists():
            errors.append(f"文件不存在: {pdb_file}")
            return False, errors
        
        try:
            with open(pdb_file, 'r', encoding='utf-8') as f:
                lines = f.readlines()
            
            if not lines:
                errors.append("文件为空")
                return False, errors
            
            # 检查是否有ATOM或HETATM行
            has_atoms = any(line.startswith(('ATOM', 'HETATM')) for line in lines)
            if not has_atoms:
                errors.append("未找到ATOM或HETATM记录")
            
            # 检查行格式
            atom_lines = [line for line in lines if line.startswith(('ATOM', 'HETATM'))]
            if atom_lines:
                for i, line in enumerate(atom_lines[:10]):  # 检查前10行
                    if len(line) < 54:
                        errors.append(f"第{i+1}个原子行格式不正确（长度不足）")
                        break
        
        except Exception as e:
            errors.append(f"文件读取失败: {e}")
        
        return len(errors) == 0, errors
    
    def get_residue_statistics(self, residues: List[ResidueInfo]) -> Dict[str, int]:
        """
        获取残基统计信息
        
        Args:
            residues: 残基信息列表
        
        Returns:
            统计信息字典
        """
        stats = {
            'total_residues': len(residues),
            'total_atoms': sum(len(r.atoms) for r in residues),
            'residue_types': {},
            'chain_counts': {},
            'elements': {}
        }
        
        for residue in residues:
            # 残基类型统计
            residue_type = residue.residue_name
            stats['residue_types'][residue_type] = stats['residue_types'].get(residue_type, 0) + 1
            
            # 链统计
            chain = residue.chain_id
            stats['chain_counts'][chain] = stats['chain_counts'].get(chain, 0) + 1
            
            # 元素统计
            for atom in residue.atoms:
                element = atom.element
                stats['elements'][element] = stats['elements'].get(element, 0) + 1
        
        return stats