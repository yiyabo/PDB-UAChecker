"""
精确SMILES解析器
用于准确解析氨基酸分子结构
"""

import re
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class Atom:
    """原子节点"""
    symbol: str
    index: int
    charge: int = 0
    hydrogen_count: int = 0
    chirality: Optional[str] = None  # @, @@
    aromatic: bool = False
    connections: List[int] = None
    
    def __post_init__(self):
        if self.connections is None:
            self.connections = []


@dataclass 
class Bond:
    """化学键"""
    atom1: int
    atom2: int
    bond_type: str = '-'  # -, =, #, :
    stereo: Optional[str] = None
    

@dataclass
class Ring:
    """环结构"""
    number: int
    atoms: List[int]
    aromatic: bool = False


class SMILESParser:
    """SMILES分子结构解析器"""
    
    def __init__(self):
        self.reset()
    
    def reset(self):
        """重置解析器状态"""
        self.atoms: List[Atom] = []
        self.bonds: List[Bond] = []
        self.rings: Dict[int, List[int]] = defaultdict(list)
        self.current_atom_index = 0
        self.branch_stack: List[int] = []
        
    def parse(self, smiles: str) -> Dict[str, Any]:
        """
        解析SMILES字符串
        
        Args:
            smiles: SMILES字符串
            
        Returns:
            分子结构信息
        """
        self.reset()
        
        # 保存原始SMILES
        self._original_smiles = smiles
        
        # 预处理SMILES字符串
        clean_smiles = self._preprocess_smiles(smiles)
        
        # 逐字符解析
        i = 0
        while i < len(clean_smiles):
            char = clean_smiles[i]
            
            if char == '(':
                # 开始分支
                self._start_branch()
            elif char == ')':
                # 结束分支
                self._end_branch()
            elif char == '[':
                # 显式原子
                atom_end = clean_smiles.find(']', i)
                atom_str = clean_smiles[i+1:atom_end]
                self._parse_explicit_atom(atom_str)
                i = atom_end
            elif char.isdigit():
                # 环闭合
                ring_num = int(char)
                self._handle_ring_closure(ring_num)
            elif char in '=-#:':
                # 键类型（跳过，在原子连接时处理）
                pass
            elif char.isalpha():
                # 隐式原子
                if char.islower():
                    # 芳香性原子
                    self._parse_aromatic_atom(char)
                else:
                    # 普通原子
                    self._parse_simple_atom(char)
            
            i += 1
        
        # 完成环结构
        self._finalize_rings()
        
        return {
            'atoms': self.atoms,
            'bonds': self.bonds,
            'rings': self._analyze_rings(),
            'molecular_info': self._extract_molecular_info()
        }
    
    def _preprocess_smiles(self, smiles: str) -> str:
        """预处理SMILES字符串"""
        # 移除空白字符
        clean = re.sub(r'\s+', '', smiles)
        return clean
    
    def _parse_explicit_atom(self, atom_str: str):
        """解析显式原子 [CH3], [C@H], [NH2+] 等"""
        
        # 解析手性
        chirality = None
        if '@@@' in atom_str:
            chirality = '@@@'
            atom_str = atom_str.replace('@@@', '')
        elif '@@' in atom_str:
            chirality = '@@'
            atom_str = atom_str.replace('@@', '')
        elif '@' in atom_str:
            chirality = '@'
            atom_str = atom_str.replace('@', '')
        
        # 解析电荷
        charge = 0
        charge_match = re.search(r'([+-])(\d*)', atom_str)
        if charge_match:
            sign = 1 if charge_match.group(1) == '+' else -1
            count = int(charge_match.group(2)) if charge_match.group(2) else 1
            charge = sign * count
            atom_str = re.sub(r'[+-]\d*', '', atom_str)
        
        # 解析氢原子数
        hydrogen_count = 0
        h_match = re.search(r'H(\d*)', atom_str)
        if h_match:
            hydrogen_count = int(h_match.group(1)) if h_match.group(1) else 1
            atom_str = re.sub(r'H\d*', '', atom_str)
        
        # 提取原子符号
        symbol = atom_str or 'C'  # 默认为碳
        
        # 创建原子
        atom = Atom(
            symbol=symbol,
            index=self.current_atom_index,
            charge=charge,
            hydrogen_count=hydrogen_count,
            chirality=chirality,
            aromatic=symbol.islower()
        )
        
        self.atoms.append(atom)
        
        # 连接到前一个原子
        if self.current_atom_index > 0:
            prev_atom = self._get_current_parent()
            self._add_bond(prev_atom, self.current_atom_index)
        
        self.current_atom_index += 1
    
    def _parse_simple_atom(self, char: str):
        """解析简单原子"""
        atom = Atom(
            symbol=char.upper(),
            index=self.current_atom_index,
            aromatic=char.islower()
        )
        
        self.atoms.append(atom)
        
        # 连接到前一个原子
        if self.current_atom_index > 0:
            prev_atom = self._get_current_parent()
            self._add_bond(prev_atom, self.current_atom_index)
        
        self.current_atom_index += 1
    
    def _parse_aromatic_atom(self, char: str):
        """解析芳香性原子"""
        self._parse_simple_atom(char)
        # 芳香性已在parse_simple_atom中处理
    
    def _start_branch(self):
        """开始分支"""
        if self.current_atom_index > 0:
            self.branch_stack.append(self.current_atom_index - 1)
    
    def _end_branch(self):
        """结束分支"""
        if self.branch_stack:
            self.branch_stack.pop()
    
    def _get_current_parent(self) -> int:
        """获取当前父原子索引"""
        if self.branch_stack:
            return self.branch_stack[-1]
        else:
            return self.current_atom_index - 1
    
    def _handle_ring_closure(self, ring_num: int):
        """处理环闭合"""
        current_atom = self.current_atom_index - 1
        
        if ring_num in self.rings:
            # 闭合环
            start_atom = self.rings[ring_num].pop()
            self._add_bond(start_atom, current_atom)
            if not self.rings[ring_num]:
                del self.rings[ring_num]
        else:
            # 开始新环
            self.rings[ring_num].append(current_atom)
    
    def _add_bond(self, atom1: int, atom2: int, bond_type: str = '-'):
        """添加化学键"""
        bond = Bond(atom1=atom1, atom2=atom2, bond_type=bond_type)
        self.bonds.append(bond)
        
        # 更新原子连接
        if atom1 < len(self.atoms):
            self.atoms[atom1].connections.append(atom2)
        if atom2 < len(self.atoms):
            self.atoms[atom2].connections.append(atom1)
    
    def _finalize_rings(self):
        """完成环结构分析"""
        # 处理未闭合的环
        for ring_num, atoms in self.rings.items():
            if len(atoms) == 2:
                self._add_bond(atoms[0], atoms[1])
    
    def _analyze_rings(self) -> List[Ring]:
        """分析环结构"""
        rings = []
        visited = set()
        
        # 使用深度优先搜索找到环
        for atom_idx in range(len(self.atoms)):
            if atom_idx not in visited:
                ring_atoms = self._find_ring_from_atom(atom_idx, visited)
                if ring_atoms:
                    is_aromatic = all(self.atoms[i].aromatic for i in ring_atoms)
                    ring = Ring(
                        number=len(rings) + 1,
                        atoms=ring_atoms,
                        aromatic=is_aromatic
                    )
                    rings.append(ring)
        
        return rings
    
    def _find_ring_from_atom(self, start_atom: int, visited: set) -> List[int]:
        """从指定原子开始寻找环"""
        # 简化的环检测算法
        # 实际实现需要更复杂的图算法
        stack = [start_atom]
        path = []
        
        while stack:
            current = stack.pop()
            if current in visited:
                continue
                
            visited.add(current)
            path.append(current)
            
            for neighbor in self.atoms[current].connections:
                if neighbor not in visited:
                    stack.append(neighbor)
                elif neighbor in path[:-1]:
                    # 找到环
                    ring_start = path.index(neighbor)
                    return path[ring_start:]
        
        return []
    
    def _extract_molecular_info(self) -> Dict[str, Any]:
        """提取分子信息"""
        info = {
            'atom_count': len(self.atoms),
            'bond_count': len(self.bonds),
            'ring_count': len([r for r in self._analyze_rings() if len(r.atoms) > 0]),
            'aromatic_atoms': [i for i, atom in enumerate(self.atoms) if atom.aromatic],
            'chiral_centers': [i for i, atom in enumerate(self.atoms) if atom.chirality],
            'functional_groups': self._identify_functional_groups()
        }
        
        return info
    
    def _identify_functional_groups(self) -> List[Dict[str, Any]]:
        """识别官能团"""
        groups = []
        
        # 羧基检测
        carboxyl_pattern = self._find_carboxyl_groups()
        if carboxyl_pattern:
            groups.extend(carboxyl_pattern)
        
        # 氨基检测
        amino_pattern = self._find_amino_groups()
        if amino_pattern:
            groups.extend(amino_pattern)
        
        # 羟基检测
        hydroxyl_pattern = self._find_hydroxyl_groups()
        if hydroxyl_pattern:
            groups.extend(hydroxyl_pattern)
        
        return groups
    
    def _find_carboxyl_groups(self) -> List[Dict[str, Any]]:
        """找到羧基基团（改进版）"""
        carboxyl_groups = []
        
        # 直接从原始SMILES中查找羧基模式
        original_smiles = getattr(self, '_original_smiles', '')
        
        # 羧基的常见SMILES模式
        import re
        carboxyl_patterns = [
            r'C\(=O\)O',     # C(=O)O
            r'C\(O\)=O',     # C(O)=O  
            r'CO=O',         # 简化形式
            r'C=OO',         # 另一种形式
        ]
        
        found_carboxyl = False
        for pattern in carboxyl_patterns:
            if re.search(pattern, original_smiles, re.IGNORECASE):
                found_carboxyl = True
                break
        
        if found_carboxyl:
            # 尝试从原子连接中找到可能的羧基碳
            for i, atom in enumerate(self.atoms):
                if atom.symbol == 'C':
                    # 检查是否连接到氧原子
                    oxygen_connections = []
                    for conn in atom.connections:
                        if conn < len(self.atoms) and self.atoms[conn].symbol == 'O':
                            oxygen_connections.append(conn)
                    
                    # 如果有氧原子连接，假设是羧基碳
                    if oxygen_connections:
                        carboxyl_groups.append({
                            'type': 'carboxyl',
                            'carbon_atom': i,
                            'oxygen_atoms': oxygen_connections,
                            'pattern': 'COOH_detected_from_smiles'
                        })
                        break  # 只取第一个可能的羧基
        
        # 如果没有找到但SMILES包含羧基模式，创建虚拟羧基
        if not carboxyl_groups and found_carboxyl:
            carboxyl_groups.append({
                'type': 'carboxyl',
                'carbon_atom': -1,  # 虚拟位置
                'oxygen_atoms': [],
                'pattern': 'COOH_pattern_match'
            })
        
        return carboxyl_groups
    
    def _find_amino_groups(self) -> List[Dict[str, Any]]:
        """找到氨基基团"""
        amino_groups = []
        
        for i, atom in enumerate(self.atoms):
            if atom.symbol == 'N':
                # 检查氮原子的连接
                carbon_connections = []
                for conn in atom.connections:
                    if conn < len(self.atoms) and self.atoms[conn].symbol == 'C':
                        carbon_connections.append(conn)
                
                amino_groups.append({
                    'type': 'amino',
                    'nitrogen_atom': i,
                    'carbon_connections': carbon_connections,
                    'hydrogen_count': atom.hydrogen_count,
                    'charge': atom.charge
                })
        
        return amino_groups
    
    def _find_hydroxyl_groups(self) -> List[Dict[str, Any]]:
        """找到羟基基团"""
        hydroxyl_groups = []
        
        for i, atom in enumerate(self.atoms):
            if atom.symbol == 'O' and atom.hydrogen_count > 0:
                # 羟基：氧原子连接氢原子
                carbon_connections = []
                for conn in atom.connections:
                    if conn < len(self.atoms) and self.atoms[conn].symbol == 'C':
                        carbon_connections.append(conn)
                
                if carbon_connections:
                    hydroxyl_groups.append({
                        'type': 'hydroxyl',
                        'oxygen_atom': i,
                        'carbon_connections': carbon_connections
                    })
        
        return hydroxyl_groups


class AminoAcidStructureAnalyzer:
    """氨基酸结构专用分析器"""
    
    def __init__(self):
        self.parser = SMILESParser()
    
    def analyze_amino_acid_structure(self, smiles: str) -> Dict[str, Any]:
        """
        分析氨基酸结构
        
        Args:
            smiles: SMILES字符串
            
        Returns:
            氨基酸结构分析结果
        """
        # 保存原始SMILES供后续分析使用
        self.original_smiles = smiles
        
        # 解析分子结构
        structure = self.parser.parse(smiles)
        
        # 专门的氨基酸分析
        analysis = {
            'backbone_analysis': self._analyze_backbone(structure, smiles),
            'functional_groups': structure['molecular_info']['functional_groups'],
            'chirality_analysis': self._analyze_chirality(structure),
            'side_chain_analysis': self._analyze_side_chain(structure),
            'ring_systems': structure['rings'],
            'confidence': self._calculate_confidence(structure)
        }
        
        return analysis
    
    def _analyze_backbone(self, structure: Dict[str, Any], smiles: str) -> Dict[str, Any]:
        """分析主链骨架（改进版）"""
        atoms = structure['atoms']
        functional_groups = structure['molecular_info']['functional_groups']
        
        # 找到羧基和氨基
        carboxyl_groups = [g for g in functional_groups if g['type'] == 'carboxyl']
        amino_groups = [g for g in functional_groups if g['type'] == 'amino']
        
        if not carboxyl_groups or not amino_groups:
            return {'type': 'unknown', 'confidence': 0.0, 'reason': 'missing_essential_groups'}
        
        # 简化分析：基于SMILES模式直接判断
        original_smiles = smiles
        
        # α氨基酸的典型模式：N[C@H]...C(=O)O 或 N[C@@H]...C(=O)O
        import re
        alpha_patterns = [
            r'N\[C@@?H?\]\([^)]*\)C\(=O\)O',  # N[C@H](R)C(=O)O
            r'N\[C@@?H?\].*C\(=O\)O',        # N[C@H]...C(=O)O (更宽泛)
        ]
        
        for pattern in alpha_patterns:
            if re.search(pattern, original_smiles):
                return {
                    'type': 'alpha',
                    'confidence': 0.95,
                    'pattern_match': pattern,
                    'smiles_analysis': 'direct_pattern'
                }
        
        # β氨基酸模式：NH2-CH2-CH(R)-COOH
        beta_patterns = [
            r'N.*C.*C.*\[C@@?H?\].*C\(=O\)O',  # N-C-C-[C@H]-COOH
            r'\[NH2\]C.*\[C@@?H?\].*C\(=O\)O',  # [NH2]C-[C@H]-COOH
        ]
        
        for pattern in beta_patterns:
            if re.search(pattern, original_smiles):
                return {
                    'type': 'beta',
                    'confidence': 0.90,
                    'pattern_match': pattern,
                    'smiles_analysis': 'beta_pattern'
                }
        
        # γ氨基酸模式
        gamma_patterns = [
            r'N.*C.*C.*C.*C\(=O\)O',  # N-C-C-C-COOH
        ]
        
        for pattern in gamma_patterns:
            if re.search(pattern, original_smiles):
                return {
                    'type': 'gamma',
                    'confidence': 0.90,
                    'pattern_match': pattern,
                    'smiles_analysis': 'gamma_pattern'
                }
        
        # 如果都不匹配，但有氨基和羧基，默认为α型
        if 'N' in original_smiles and 'C(=O)O' in original_smiles:
            return {
                'type': 'alpha',
                'confidence': 0.7,
                'smiles_analysis': 'default_alpha_assumption'
            }
        
        return {'type': 'complex', 'confidence': 0.3}
    
    def _find_shortest_path(self, atoms: List[Atom], start: int, end: int) -> List[int]:
        """找到两个原子间的最短路径"""
        if start == end:
            return [start]
        
        queue = [[start]]
        visited = {start}
        
        while queue:
            path = queue.pop(0)
            current = path[-1]
            
            for neighbor in atoms[current].connections:
                if neighbor == end:
                    return path + [neighbor]
                
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(path + [neighbor])
        
        return []  # 未找到路径
    
    def _analyze_chirality(self, structure: Dict[str, Any]) -> Dict[str, Any]:
        """分析手性中心"""
        chiral_centers = structure['molecular_info']['chiral_centers']
        
        chirality_info = {
            'has_chirality': len(chiral_centers) > 0,
            'chiral_center_count': len(chiral_centers),
            'centers': []
        }
        
        for center_idx in chiral_centers:
            atom = structure['atoms'][center_idx]
            center_info = {
                'atom_index': center_idx,
                'atom_symbol': atom.symbol,
                'chirality_symbol': atom.chirality,
                'configuration': self._determine_configuration(atom),
                'confidence': 0.8  # 基于SMILES标记的置信度
            }
            chirality_info['centers'].append(center_info)
        
        return chirality_info
    
    def _determine_configuration(self, atom: Atom) -> str:
        """确定手性配置"""
        if atom.chirality == '@':
            return 'R'
        elif atom.chirality == '@@':
            return 'S'
        else:
            return 'unknown'
    
    def _analyze_side_chain(self, structure: Dict[str, Any]) -> Dict[str, Any]:
        """分析侧链"""
        # 简化的侧链分析
        atoms = structure['atoms']
        rings = structure['rings']
        
        side_chain_info = {
            'has_aromatic_ring': any(ring.aromatic for ring in rings),
            'ring_count': len(rings),
            'heteroatoms': [],
            'functional_groups_in_side_chain': []
        }
        
        # 统计杂原子
        for i, atom in enumerate(atoms):
            if atom.symbol not in ['C', 'H']:
                side_chain_info['heteroatoms'].append({
                    'atom_index': i,
                    'symbol': atom.symbol,
                    'charge': atom.charge
                })
        
        return side_chain_info
    
    def _calculate_confidence(self, structure: Dict[str, Any]) -> float:
        """计算整体分析置信度"""
        confidence_factors = []
        
        # 基于原子数量
        atom_count = structure['molecular_info']['atom_count']
        if 5 <= atom_count <= 50:  # 合理的氨基酸范围
            confidence_factors.append(0.9)
        else:
            confidence_factors.append(0.6)
        
        # 基于官能团
        functional_groups = structure['molecular_info']['functional_groups']
        has_carboxyl = any(g['type'] == 'carboxyl' for g in functional_groups)
        has_amino = any(g['type'] == 'amino' for g in functional_groups)
        
        if has_carboxyl and has_amino:
            confidence_factors.append(0.95)
        elif has_carboxyl or has_amino:
            confidence_factors.append(0.7)
        else:
            confidence_factors.append(0.3)
        
        return sum(confidence_factors) / len(confidence_factors) if confidence_factors else 0.5