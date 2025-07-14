#!/usr/bin/env python3
"""
原子级别分析器 - 精细化学结构分析
支持碳原子编号、连接性分析、功能基团识别等高级功能
"""

import re
import json
from typing import Dict, List, Tuple, Optional, Set, Any
from dataclasses import dataclass
from collections import defaultdict, deque

# 尝试导入RDKit，如果没有则使用内置解析器
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("警告: RDKit未安装，将使用内置SMILES解析器")

@dataclass
class AtomInfo:
    """原子信息数据类"""
    index: int
    element: str
    formal_charge: int = 0
    hybridization: str = "unknown"
    aromatic: bool = False
    in_ring: bool = False
    neighbors: List[int] = None
    bonds: List[Tuple[int, str]] = None  # (neighbor_index, bond_type)
    
    def __post_init__(self):
        if self.neighbors is None:
            self.neighbors = []
        if self.bonds is None:
            self.bonds = []

@dataclass
class CarbonAnalysis:
    """碳原子分析结果"""
    carbon_number: int
    carbon_index: int
    attached_atoms: List[str]
    attached_indices: List[int]
    bond_types: List[str]
    functional_groups: List[str]
    hybridization: str
    stereochemistry: Optional[str] = None
    in_ring: bool = False
    ring_size: Optional[int] = None
    chemical_environment: str = ""

class SMILESParser:
    """内置SMILES解析器"""
    
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.rings = defaultdict(list)
    
    def parse_smiles(self, smiles: str) -> Tuple[List[AtomInfo], List[Tuple[int, int, str]]]:
        """解析SMILES字符串"""
        self.atoms = []
        self.bonds = []
        self.rings = defaultdict(list)
        
        # 简化的SMILES解析（处理基本情况）
        atom_stack = []
        branch_stack = []
        current_atom = -1
        
        i = 0
        while i < len(smiles):
            char = smiles[i]
            
            if char == '(':
                # 开始分支
                branch_stack.append(current_atom)
            elif char == ')':
                # 结束分支
                if branch_stack:
                    branch_stack.pop()
            elif char.isalpha():
                # 原子
                element = char
                if i + 1 < len(smiles) and smiles[i + 1].islower():
                    element += smiles[i + 1]
                    i += 1
                
                atom_info = AtomInfo(
                    index=len(self.atoms),
                    element=element
                )
                self.atoms.append(atom_info)
                
                # 连接到前一个原子
                if current_atom >= 0:
                    self._add_bond(current_atom, len(self.atoms) - 1, 'single')
                
                current_atom = len(self.atoms) - 1
                
                # 处理分支连接
                if branch_stack:
                    parent = branch_stack[-1]
                    if parent != current_atom:
                        self._add_bond(parent, current_atom, 'single')
            
            elif char.isdigit():
                # 环连接
                ring_num = int(char)
                if ring_num in self.rings:
                    # 闭环
                    other_atom = self.rings[ring_num].pop()
                    self._add_bond(current_atom, other_atom, 'single')
                else:
                    # 开环
                    self.rings[ring_num].append(current_atom)
            
            elif char in ['=', '#']:
                # 双键或三键（简化处理）
                pass
            
            i += 1
        
        return self.atoms, self.bonds
    
    def _add_bond(self, atom1: int, atom2: int, bond_type: str):
        """添加化学键"""
        self.bonds.append((atom1, atom2, bond_type))
        if atom1 < len(self.atoms):
            self.atoms[atom1].neighbors.append(atom2)
            self.atoms[atom1].bonds.append((atom2, bond_type))
        if atom2 < len(self.atoms):
            self.atoms[atom2].neighbors.append(atom1)
            self.atoms[atom2].bonds.append((atom1, bond_type))

class AtomicAnalyzer:
    """原子级别分析器"""
    
    def __init__(self):
        self.parser = SMILESParser()
        self.functional_groups = {
            'hydroxyl': ['[OH]', 'O'],
            'amino': ['[NH2]', '[NH]', 'N'],
            'carboxyl': ['C(=O)O', 'COOH'],
            'carbonyl': ['C=O'],
            'methyl': ['CH3'],
            'phenyl': ['c1ccccc1'],
            'sulfhydryl': ['SH'],
            'phosphate': ['PO4'],
            'nitro': ['NO2'],
            'halogen': ['F', 'Cl', 'Br', 'I']
        }
    
    def analyze_molecule(self, smiles: str) -> Dict[str, Any]:
        """分析整个分子"""
        if RDKIT_AVAILABLE:
            return self._analyze_with_rdkit(smiles)
        else:
            return self._analyze_with_parser(smiles)
    
    def _analyze_with_rdkit(self, smiles: str) -> Dict[str, Any]:
        """使用RDKit进行分析"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("无效的SMILES字符串")
            
            # 添加氢原子
            mol = Chem.AddHs(mol)
            
            analysis = {
                'total_atoms': mol.GetNumAtoms(),
                'molecular_formula': Chem.rdMolDescriptors.CalcMolFormula(mol),
                'molecular_weight': Descriptors.MolWt(mol),
                'carbon_analysis': self._analyze_carbons_rdkit(mol),
                'functional_groups': self._identify_functional_groups_rdkit(mol),
                'ring_info': self._analyze_rings_rdkit(mol)
            }
            
            return analysis
            
        except Exception as e:
            print(f"RDKit分析失败: {e}")
            return self._analyze_with_parser(smiles)
    
    def _analyze_carbons_rdkit(self, mol) -> List[CarbonAnalysis]:
        """使用RDKit分析碳原子"""
        carbon_analyses = []
        carbon_count = 0
        
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C':
                carbon_count += 1
                
                # 获取邻接原子
                neighbors = atom.GetNeighbors()
                attached_atoms = [n.GetSymbol() for n in neighbors]
                attached_indices = [n.GetIdx() for n in neighbors]
                
                # 获取键类型
                bond_types = []
                for neighbor in neighbors:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    bond_types.append(str(bond.GetBondType()))
                
                # 杂化状态
                hybridization = str(atom.GetHybridization())
                
                # 环信息
                in_ring = atom.IsInRing()
                ring_size = None
                if in_ring:
                    ring_info = mol.GetRingInfo()
                    for ring in ring_info.AtomRings():
                        if atom.GetIdx() in ring:
                            ring_size = len(ring)
                            break
                
                # 功能基团识别
                functional_groups = self._identify_carbon_functional_groups(
                    atom, mol, attached_atoms
                )
                
                # 化学环境描述
                environment = self._describe_chemical_environment(
                    attached_atoms, bond_types, functional_groups
                )
                
                carbon_analysis = CarbonAnalysis(
                    carbon_number=carbon_count,
                    carbon_index=atom.GetIdx(),
                    attached_atoms=attached_atoms,
                    attached_indices=attached_indices,
                    bond_types=bond_types,
                    functional_groups=functional_groups,
                    hybridization=hybridization,
                    in_ring=in_ring,
                    ring_size=ring_size,
                    chemical_environment=environment
                )
                
                carbon_analyses.append(carbon_analysis)
        
        return carbon_analyses
    
    def _analyze_with_parser(self, smiles: str) -> Dict[str, Any]:
        """使用内置解析器进行分析"""
        atoms, bonds = self.parser.parse_smiles(smiles)
        
        # 统计原子
        atom_counts = defaultdict(int)
        for atom in atoms:
            atom_counts[atom.element] += 1
        
        # 分析碳原子
        carbon_analyses = []
        carbon_count = 0
        
        for atom in atoms:
            if atom.element == 'C':
                carbon_count += 1
                
                attached_atoms = [atoms[i].element for i in atom.neighbors]
                bond_types = [bond_type for _, bond_type in atom.bonds]
                
                functional_groups = self._identify_functional_groups_simple(
                    attached_atoms
                )
                
                environment = self._describe_chemical_environment(
                    attached_atoms, bond_types, functional_groups
                )
                
                carbon_analysis = CarbonAnalysis(
                    carbon_number=carbon_count,
                    carbon_index=atom.index,
                    attached_atoms=attached_atoms,
                    attached_indices=atom.neighbors,
                    bond_types=bond_types,
                    functional_groups=functional_groups,
                    hybridization="unknown",
                    chemical_environment=environment
                )
                
                carbon_analyses.append(carbon_analysis)
        
        return {
            'total_atoms': len(atoms),
            'atom_composition': dict(atom_counts),
            'carbon_analysis': carbon_analyses,
            'functional_groups': self._identify_functional_groups_simple_mol(atoms),
            'bonds_count': len(bonds)
        }

    def _identify_carbon_functional_groups(self, carbon_atom, mol, attached_atoms) -> List[str]:
        """识别碳原子相关的功能基团"""
        groups = []

        # 检查常见功能基团
        if 'O' in attached_atoms:
            # 可能是羟基、羰基等
            for neighbor in carbon_atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O':
                    # 检查氧原子的连接
                    oxygen_neighbors = [n.GetSymbol() for n in neighbor.GetNeighbors()]
                    if len(oxygen_neighbors) == 1:  # 只连一个原子，可能是羰基
                        bond = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), neighbor.GetIdx())
                        if str(bond.GetBondType()) == 'DOUBLE':
                            groups.append('carbonyl')
                    elif 'H' in oxygen_neighbors:  # 连接氢原子，是羟基
                        groups.append('hydroxyl')

        if 'N' in attached_atoms:
            groups.append('amino_related')

        if attached_atoms.count('H') == 3 and len(attached_atoms) == 4:
            groups.append('methyl')

        return groups

    def _identify_functional_groups_simple(self, attached_atoms) -> List[str]:
        """简单的功能基团识别"""
        groups = []

        if 'O' in attached_atoms:
            groups.append('oxygen_containing')
        if 'N' in attached_atoms:
            groups.append('nitrogen_containing')
        if 'S' in attached_atoms:
            groups.append('sulfur_containing')
        if attached_atoms.count('H') == 3:
            groups.append('methyl')

        return groups

    def _identify_functional_groups_simple_mol(self, atoms) -> List[str]:
        """识别整个分子的功能基团"""
        groups = []
        elements = [atom.element for atom in atoms]

        if 'O' in elements:
            groups.append('oxygen_containing')
        if 'N' in elements:
            groups.append('nitrogen_containing')
        if 'S' in elements:
            groups.append('sulfur_containing')

        return groups

    def _identify_functional_groups_rdkit(self, mol) -> List[str]:
        """使用RDKit识别功能基团"""
        groups = []

        # 定义SMARTS模式
        patterns = {
            'hydroxyl': '[OH]',
            'amino': '[NH2]',
            'carboxyl': 'C(=O)O',
            'carbonyl': 'C=O',
            'methyl': '[CH3]',
            'phenyl': 'c1ccccc1',
            'sulfhydryl': '[SH]',
            'nitro': '[N+](=O)[O-]'
        }

        for group_name, pattern in patterns.items():
            try:
                patt = Chem.MolFromSmarts(pattern)
                if mol.HasSubstructMatch(patt):
                    groups.append(group_name)
            except:
                continue

        return groups

    def _analyze_rings_rdkit(self, mol) -> Dict[str, Any]:
        """分析环结构"""
        ring_info = mol.GetRingInfo()

        return {
            'num_rings': ring_info.NumRings(),
            'ring_sizes': [len(ring) for ring in ring_info.AtomRings()],
            'aromatic_rings': len([ring for ring in ring_info.AtomRings()
                                 if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)])
        }

    def _describe_chemical_environment(self, attached_atoms: List[str],
                                     bond_types: List[str],
                                     functional_groups: List[str]) -> str:
        """描述化学环境"""
        descriptions = []

        # 连接的原子类型
        atom_counts = {}
        for atom in attached_atoms:
            atom_counts[atom] = atom_counts.get(atom, 0) + 1

        atom_desc = []
        for atom, count in atom_counts.items():
            if count == 1:
                atom_desc.append(atom)
            else:
                atom_desc.append(f"{count}{atom}")

        if atom_desc:
            descriptions.append(f"连接: {', '.join(atom_desc)}")

        # 功能基团
        if functional_groups:
            descriptions.append(f"功能基团: {', '.join(functional_groups)}")

        # 键类型
        if 'DOUBLE' in bond_types:
            descriptions.append("含双键")
        if 'TRIPLE' in bond_types:
            descriptions.append("含三键")

        return "; ".join(descriptions) if descriptions else "简单碳原子"

    def analyze_amino_acid(self, smiles: str, name: str = "") -> Dict[str, Any]:
        """专门分析氨基酸结构"""
        analysis = self.analyze_molecule(smiles)

        # 氨基酸特异性分析
        amino_acid_analysis = {
            'name': name,
            'smiles': smiles,
            'basic_analysis': analysis,
            'amino_acid_features': self._identify_amino_acid_features(analysis),
            'carbon_summary': self._summarize_carbon_analysis(analysis.get('carbon_analysis', []))
        }

        return amino_acid_analysis

    def _identify_amino_acid_features(self, analysis: Dict[str, Any]) -> Dict[str, Any]:
        """识别氨基酸特征"""
        features = {
            'has_amino_group': False,
            'has_carboxyl_group': False,
            'side_chain_carbons': 0,
            'special_features': []
        }

        functional_groups = analysis.get('functional_groups', [])

        if 'amino' in functional_groups or 'amino_related' in functional_groups:
            features['has_amino_group'] = True

        if 'carboxyl' in functional_groups:
            features['has_carboxyl_group'] = True

        # 计算侧链碳原子数（简化计算）
        carbon_analysis = analysis.get('carbon_analysis', [])
        if len(carbon_analysis) > 2:  # 除了α-碳和羧基碳
            features['side_chain_carbons'] = len(carbon_analysis) - 2

        return features

    def _summarize_carbon_analysis(self, carbon_analyses: List[CarbonAnalysis]) -> Dict[str, str]:
        """总结碳原子分析"""
        summary = {}

        for carbon in carbon_analyses:
            key = f"C{carbon.carbon_number}"
            summary[key] = {
                'attached_atoms': carbon.attached_atoms,
                'functional_groups': carbon.functional_groups,
                'environment': carbon.chemical_environment,
                'hybridization': carbon.hybridization
            }

        return summary

    def get_carbon_connectivity_report(self, smiles: str) -> str:
        """生成碳原子连接性报告"""
        analysis = self.analyze_molecule(smiles)
        carbon_analyses = analysis.get('carbon_analysis', [])

        report = []
        report.append("=== 碳原子连接性分析报告 ===\n")
        report.append(f"分子式: {analysis.get('molecular_formula', 'Unknown')}")
        report.append(f"总原子数: {analysis.get('total_atoms', 0)}")
        report.append(f"碳原子数: {len(carbon_analyses)}\n")

        for carbon in carbon_analyses:
            report.append(f"碳原子 C{carbon.carbon_number} (索引 {carbon.carbon_index}):")
            report.append(f"  连接原子: {', '.join(carbon.attached_atoms)}")
            report.append(f"  键类型: {', '.join(carbon.bond_types)}")
            report.append(f"  杂化状态: {carbon.hybridization}")

            if carbon.functional_groups:
                report.append(f"  功能基团: {', '.join(carbon.functional_groups)}")

            if carbon.in_ring:
                report.append(f"  环状结构: 是 (环大小: {carbon.ring_size})")

            report.append(f"  化学环境: {carbon.chemical_environment}")
            report.append("")

        return "\n".join(report)

def main():
    """演示原子级分析功能"""
    analyzer = AtomicAnalyzer()

    # 测试分子
    test_molecules = {
        '丝氨酸': 'N[C@@H](CO)C(=O)O',
        '苯丙氨酸': 'N[C@@H](Cc1ccccc1)C(=O)O',
        '酪氨酸': 'N[C@@H](Cc1ccc(O)cc1)C(=O)O'
    }

    print("原子级分析器演示")
    print("=" * 50)

    for name, smiles in test_molecules.items():
        print(f"\n分析 {name} ({smiles}):")
        print("-" * 30)

        try:
            # 生成详细报告
            report = analyzer.get_carbon_connectivity_report(smiles)
            print(report)

            # 氨基酸特异性分析
            aa_analysis = analyzer.analyze_amino_acid(smiles, name)
            print("氨基酸特征:")
            features = aa_analysis['amino_acid_features']
            for key, value in features.items():
                print(f"  {key}: {value}")

        except Exception as e:
            print(f"分析失败: {e}")

if __name__ == "__main__":
    main()
