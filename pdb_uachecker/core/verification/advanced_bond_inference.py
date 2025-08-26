"""
高级化学键推断模块
处理复杂分子结构的化学键推断和价态验证
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass
import logging
from collections import defaultdict

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, AllChem, rdDistGeom, rdDetermineBonds
    from rdkit.Geometry import Point3D
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    # 创建空的Chem类用于类型提示
    class Chem:
        class Mol:
            pass
    logging.warning("RDKit未安装，高级化学键推断功能受限")


@dataclass
class AtomProperties:
    """原子属性"""
    element: str
    atomic_number: int
    max_valence: int
    common_valences: List[int]
    electronegativity: float
    covalent_radius: float


@dataclass
class BondCandidate:
    """化学键候选"""
    atom1_idx: int
    atom2_idx: int
    distance: float
    bond_order: int = 1
    confidence: float = 0.0
    bond_type: str = "SINGLE"


class AdvancedBondInference:
    """高级化学键推断器"""
    
    # 原子属性表（扩展版）
    ATOM_PROPERTIES = {
        'H': AtomProperties('H', 1, 1, [1], 2.20, 0.31),
        'C': AtomProperties('C', 6, 4, [4], 2.55, 0.76),
        'N': AtomProperties('N', 7, 5, [3, 5], 3.04, 0.71),
        'O': AtomProperties('O', 8, 2, [2], 3.44, 0.66),
        'F': AtomProperties('F', 9, 1, [1], 3.98, 0.57),
        'P': AtomProperties('P', 15, 5, [3, 5], 2.19, 1.07),
        'S': AtomProperties('S', 16, 6, [2, 4, 6], 2.58, 1.05),
        'Cl': AtomProperties('Cl', 17, 1, [1], 3.16, 0.99),
        'Br': AtomProperties('Br', 35, 1, [1], 2.96, 1.14),
        'I': AtomProperties('I', 53, 1, [1], 2.66, 1.33),
        'Se': AtomProperties('Se', 34, 6, [2, 4, 6], 2.55, 1.20),
    }
    
    # 键长参考表（单键长度，Å）
    BOND_LENGTHS = {
        ('C', 'C'): 1.54, ('C', 'N'): 1.47, ('C', 'O'): 1.43, ('C', 'S'): 1.82,
        ('C', 'P'): 1.84, ('C', 'F'): 1.35, ('C', 'Cl'): 1.77, ('C', 'Br'): 1.94,
        ('C', 'I'): 2.14, ('C', 'H'): 1.09, ('N', 'N'): 1.45, ('N', 'O'): 1.36,
        ('N', 'S'): 1.74, ('N', 'P'): 1.77, ('N', 'H'): 1.01, ('O', 'O'): 1.48,
        ('O', 'S'): 1.70, ('O', 'P'): 1.63, ('O', 'H'): 0.96, ('S', 'S'): 2.05,
        ('S', 'P'): 2.10, ('S', 'H'): 1.34, ('P', 'P'): 2.21, ('P', 'H'): 1.42,
        ('F', 'F'): 1.42, ('Cl', 'Cl'): 1.99, ('Br', 'Br'): 2.28, ('I', 'I'): 2.67,
    }
    
    # 键长修正系数（双键、三键）
    BOND_ORDER_FACTORS = {
        1: 1.0,    # 单键
        2: 0.87,   # 双键
        3: 0.78    # 三键
    }
    
    def __init__(self):
        """初始化高级化学键推断器"""
        self.bond_tolerance = 1.4  # 键长容忍度
        self.angle_tolerance = 15.0  # 角度容忍度（度）
        
    def infer_bonds_advanced(self, mol, atoms_info: List[Dict]) -> Optional[Chem.Mol]:
        """
        高级化学键推断
        
        Args:
            mol: RDKit分子对象
            atoms_info: 原子信息列表
            
        Returns:
            推断化学键后的分子对象
        """
        try:
            # 方法1: 尝试RDKit自动推断（带参数优化）
            mol_copy = self._try_rdkit_inference_optimized(mol)
            if mol_copy is not None:
                return mol_copy
            
            # 方法2: 基于距离和化学规则的推断
            mol_copy = self._distance_based_inference(mol, atoms_info)
            if mol_copy is not None:
                return mol_copy
            
            # 方法3: 模板匹配推断
            mol_copy = self._template_based_inference(mol, atoms_info)
            if mol_copy is not None:
                return mol_copy
            
            # 方法4: 最保守的推断
            return self._conservative_inference(mol, atoms_info)
            
        except Exception as e:
            logging.error(f"高级化学键推断失败: {e}")
            return None
    
    def _try_rdkit_inference_optimized(self, mol) -> Optional[Chem.Mol]:
        """
        优化的RDKit自动推断
        
        Args:
            mol: 分子对象
            
        Returns:
            推断后的分子对象
        """
        try:
            mol_copy = Chem.RWMol(mol)
            
            # 设置更宽松的参数
            rdDetermineBonds.DetermineBonds(
                mol_copy, 
                charge=0,  # 假设电荷为0
                allowChargedFrags=True,  # 允许带电片段
                embedChiral=True,  # 保留手性信息
                useAtomMap=False
            )
            
            # 尝试清理分子
            try:
                Chem.SanitizeMol(mol_copy, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
                return mol_copy
            except Exception as sanitize_error:
                # 尝试部分清理
                try:
                    Chem.SanitizeMol(mol_copy, sanitizeOps=(
                        Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                        Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                        Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                    ))
                    return mol_copy
                except:
                    logging.warning(f"分子清理失败: {sanitize_error}")
                    return None
                    
        except Exception as e:
            logging.warning(f"优化RDKit推断失败: {e}")
            return None
    
    def _distance_based_inference(self, mol, atoms_info: List[Dict]) -> Optional[Chem.Mol]:
        """
        基于距离和化学规则的键推断
        
        Args:
            mol: 分子对象
            atoms_info: 原子信息列表
            
        Returns:
            推断后的分子对象
        """
        try:
            mol_copy = Chem.RWMol(mol)
            conf = mol.GetConformer()
            num_atoms = mol.GetNumAtoms()
            
            # 获取所有可能的键候选
            bond_candidates = self._generate_bond_candidates(mol_copy, conf, atoms_info)
            
            # 按置信度排序
            bond_candidates.sort(key=lambda x: x.confidence, reverse=True)
            
            # 跟踪每个原子的价态使用情况
            valence_usage = defaultdict(int)
            added_bonds = set()
            
            # 逐步添加化学键
            for candidate in bond_candidates:
                atom1_idx, atom2_idx = candidate.atom1_idx, candidate.atom2_idx
                bond_key = tuple(sorted([atom1_idx, atom2_idx]))
                
                # 避免重复添加键
                if bond_key in added_bonds:
                    continue
                
                # 检查价态限制
                atom1 = mol_copy.GetAtomWithIdx(atom1_idx)
                atom2 = mol_copy.GetAtomWithIdx(atom2_idx)
                
                atom1_props = self.ATOM_PROPERTIES.get(atom1.GetSymbol())
                atom2_props = self.ATOM_PROPERTIES.get(atom2.GetSymbol())
                
                if atom1_props and atom2_props:
                    # 检查是否违反价态限制
                    if (valence_usage[atom1_idx] + candidate.bond_order <= atom1_props.max_valence and
                        valence_usage[atom2_idx] + candidate.bond_order <= atom2_props.max_valence):
                        
                        try:
                            # 添加键
                            if candidate.bond_order == 1:
                                mol_copy.AddBond(atom1_idx, atom2_idx, Chem.BondType.SINGLE)
                            elif candidate.bond_order == 2:
                                mol_copy.AddBond(atom1_idx, atom2_idx, Chem.BondType.DOUBLE)
                            elif candidate.bond_order == 3:
                                mol_copy.AddBond(atom1_idx, atom2_idx, Chem.BondType.TRIPLE)
                            
                            # 更新价态使用
                            valence_usage[atom1_idx] += candidate.bond_order
                            valence_usage[atom2_idx] += candidate.bond_order
                            added_bonds.add(bond_key)
                            
                        except Exception as e:
                            logging.warning(f"添加键失败 {atom1_idx}-{atom2_idx}: {e}")
                            continue
            
            # 尝试清理分子
            try:
                # 使用更宽松的清理选项
                Chem.SanitizeMol(mol_copy, sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                    Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                    Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                    Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                ))
                return mol_copy
                
            except Exception as e:
                logging.warning(f"距离推断分子清理失败: {e}")
                # 返回未清理的分子，至少有键结构
                return mol_copy
                
        except Exception as e:
            logging.error(f"基于距离的键推断失败: {e}")
            return None
    
    def _generate_bond_candidates(self, mol, conf, atoms_info: List[Dict]) -> List[BondCandidate]:
        """
        生成化学键候选列表
        
        Args:
            mol: 分子对象
            conf: 构象对象
            atoms_info: 原子信息列表
            
        Returns:
            键候选列表
        """
        candidates = []
        num_atoms = mol.GetNumAtoms()
        
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                atom1 = mol.GetAtomWithIdx(i)
                atom2 = mol.GetAtomWithIdx(j)
                
                # 计算距离
                pos1 = conf.GetAtomPosition(i)
                pos2 = conf.GetAtomPosition(j)
                distance = pos1.Distance(pos2)
                
                # 判断是否可能形成键
                bond_info = self._analyze_potential_bond(atom1, atom2, distance)
                
                if bond_info['is_possible']:
                    candidate = BondCandidate(
                        atom1_idx=i,
                        atom2_idx=j,
                        distance=distance,
                        bond_order=bond_info['bond_order'],
                        confidence=bond_info['confidence'],
                        bond_type=bond_info['bond_type']
                    )
                    candidates.append(candidate)
        
        return candidates
    
    def _analyze_potential_bond(self, atom1, atom2, distance: float) -> Dict:
        """
        分析潜在的化学键
        
        Args:
            atom1: 第一个原子
            atom2: 第二个原子
            distance: 原子间距离
            
        Returns:
            键分析结果
        """
        elem1, elem2 = atom1.GetSymbol(), atom2.GetSymbol()
        bond_pair = tuple(sorted([elem1, elem2]))
        
        # 获取标准键长
        standard_length = self.BOND_LENGTHS.get(bond_pair)
        if standard_length is None:
            # 使用共价半径估算
            props1 = self.ATOM_PROPERTIES.get(elem1)
            props2 = self.ATOM_PROPERTIES.get(elem2)
            if props1 and props2:
                standard_length = props1.covalent_radius + props2.covalent_radius
            else:
                return {'is_possible': False}
        
        # 分析可能的键级
        bond_analysis = {'is_possible': False}
        best_confidence = 0.0
        
        for bond_order in [1, 2, 3]:
            expected_length = standard_length * self.BOND_ORDER_FACTORS[bond_order]
            length_tolerance = expected_length * (self.bond_tolerance - 1.0)
            
            if abs(distance - expected_length) <= length_tolerance:
                # 计算置信度
                deviation = abs(distance - expected_length) / expected_length
                confidence = max(0.0, 1.0 - deviation / 0.3)  # 30%偏差时置信度为0
                
                if confidence > best_confidence:
                    best_confidence = confidence
                    bond_analysis = {
                        'is_possible': True,
                        'bond_order': bond_order,
                        'confidence': confidence,
                        'bond_type': ['SINGLE', 'DOUBLE', 'TRIPLE'][bond_order - 1],
                        'expected_length': expected_length,
                        'actual_deviation': deviation
                    }
        
        return bond_analysis
    
    def _template_based_inference(self, mol, atoms_info: List[Dict]) -> Optional[Chem.Mol]:
        """
        基于分子模板的键推断
        
        适用于氨基酸等已知结构模式
        
        Args:
            mol: 分子对象
            atoms_info: 原子信息列表
            
        Returns:
            推断后的分子对象
        """
        try:
            mol_copy = Chem.RWMol(mol)
            num_atoms = mol.GetNumAtoms()
            
            # 识别氨基酸骨架
            backbone_indices = self._identify_amino_acid_backbone(mol_copy)
            
            if backbone_indices:
                # 添加骨架键
                self._add_backbone_bonds(mol_copy, backbone_indices)
                
                # 识别侧链
                sidechain_indices = self._identify_sidechain_atoms(mol_copy, backbone_indices)
                
                # 添加侧链键
                self._add_sidechain_bonds(mol_copy, sidechain_indices, backbone_indices)
                
                try:
                    # 尝试部分清理
                    Chem.SanitizeMol(mol_copy, sanitizeOps=(
                        Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                        Chem.SanitizeFlags.SANITIZE_KEKULIZE
                    ))
                    return mol_copy
                except:
                    return mol_copy  # 返回未完全清理的分子
            
            return None
            
        except Exception as e:
            logging.error(f"模板推断失败: {e}")
            return None
    
    def _identify_amino_acid_backbone(self, mol) -> Optional[Dict[str, int]]:
        """
        识别氨基酸主链原子
        
        Args:
            mol: 分子对象
            
        Returns:
            主链原子索引字典
        """
        try:
            backbone = {}
            conf = mol.GetConformer()
            
            # 寻找羧基碳（C=O, C-OH）
            carboxyl_c = None
            amino_n = None
            alpha_c = None
            
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                if atom.GetSymbol() == 'C':
                    pos = conf.GetAtomPosition(i)
                    
                    # 寻找附近的氧原子（羧基特征）
                    nearby_oxygens = []
                    for j in range(mol.GetNumAtoms()):
                        if j != i:
                            other_atom = mol.GetAtomWithIdx(j)
                            if other_atom.GetSymbol() == 'O':
                                other_pos = conf.GetAtomPosition(j)
                                distance = pos.Distance(other_pos)
                                if distance < 1.8:  # C-O键长范围
                                    nearby_oxygens.append(j)
                    
                    if len(nearby_oxygens) >= 2:  # 羧基通常有两个氧
                        carboxyl_c = i
                        break
            
            if carboxyl_c is not None:
                backbone['carboxyl_c'] = carboxyl_c
                
                # 寻找α碳（连接到羧基碳的碳）
                carboxyl_pos = conf.GetAtomPosition(carboxyl_c)
                for i in range(mol.GetNumAtoms()):
                    if i != carboxyl_c:
                        atom = mol.GetAtomWithIdx(i)
                        if atom.GetSymbol() == 'C':
                            pos = conf.GetAtomPosition(i)
                            distance = carboxyl_pos.Distance(pos)
                            if 1.3 < distance < 1.7:  # C-C键长范围
                                alpha_c = i
                                backbone['alpha_c'] = alpha_c
                                break
                
                # 寻找氨基氮（连接到α碳的氮）
                if alpha_c is not None:
                    alpha_pos = conf.GetAtomPosition(alpha_c)
                    for i in range(mol.GetNumAtoms()):
                        if i not in [carboxyl_c, alpha_c]:
                            atom = mol.GetAtomWithIdx(i)
                            if atom.GetSymbol() == 'N':
                                pos = conf.GetAtomPosition(i)
                                distance = alpha_pos.Distance(pos)
                                if 1.2 < distance < 1.8:  # C-N键长范围
                                    amino_n = i
                                    backbone['amino_n'] = amino_n
                                    break
            
            return backbone if len(backbone) >= 2 else None
            
        except Exception as e:
            logging.error(f"氨基酸主链识别失败: {e}")
            return None
    
    def _add_backbone_bonds(self, mol, backbone_indices: Dict[str, int]):
        """添加氨基酸主链化学键"""
        try:
            if 'carboxyl_c' in backbone_indices and 'alpha_c' in backbone_indices:
                mol.AddBond(
                    backbone_indices['carboxyl_c'], 
                    backbone_indices['alpha_c'], 
                    Chem.BondType.SINGLE
                )
            
            if 'alpha_c' in backbone_indices and 'amino_n' in backbone_indices:
                mol.AddBond(
                    backbone_indices['alpha_c'], 
                    backbone_indices['amino_n'], 
                    Chem.BondType.SINGLE
                )
                
        except Exception as e:
            logging.warning(f"主链键添加失败: {e}")
    
    def _identify_sidechain_atoms(self, mol, backbone_indices: Dict[str, int]) -> List[int]:
        """识别侧链原子"""
        sidechain = []
        backbone_atoms = set(backbone_indices.values())
        
        for i in range(mol.GetNumAtoms()):
            if i not in backbone_atoms:
                sidechain.append(i)
        
        return sidechain
    
    def _add_sidechain_bonds(self, mol, sidechain_indices: List[int], backbone_indices: Dict[str, int]):
        """添加侧链化学键"""
        try:
            conf = mol.GetConformer()
            alpha_c = backbone_indices.get('alpha_c')
            
            if alpha_c is not None:
                alpha_pos = conf.GetAtomPosition(alpha_c)
                
                # 寻找与α碳最近的侧链原子（通常是β碳）
                min_distance = float('inf')
                beta_c = None
                
                for idx in sidechain_indices:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetSymbol() == 'C':  # 寻找碳原子作为β碳
                        pos = conf.GetAtomPosition(idx)
                        distance = alpha_pos.Distance(pos)
                        if 1.3 < distance < 1.8 and distance < min_distance:
                            min_distance = distance
                            beta_c = idx
                
                if beta_c is not None:
                    mol.AddBond(alpha_c, beta_c, Chem.BondType.SINGLE)
                    
                    # 为侧链内部添加键
                    self._add_sidechain_internal_bonds(mol, sidechain_indices, beta_c)
                    
        except Exception as e:
            logging.warning(f"侧链键添加失败: {e}")
    
    def _add_sidechain_internal_bonds(self, mol, sidechain_indices: List[int], start_atom: int):
        """添加侧链内部的化学键"""
        try:
            conf = mol.GetConformer()
            bonded_atoms = {start_atom}
            
            # 使用贪婪算法连接侧链原子
            while len(bonded_atoms) < len(sidechain_indices):
                best_pair = None
                best_distance = float('inf')
                
                for bonded_idx in bonded_atoms:
                    bonded_pos = conf.GetAtomPosition(bonded_idx)
                    
                    for unbonded_idx in sidechain_indices:
                        if unbonded_idx not in bonded_atoms:
                            unbonded_pos = conf.GetAtomPosition(unbonded_idx)
                            distance = bonded_pos.Distance(unbonded_pos)
                            
                            # 检查距离是否合理
                            if 0.8 < distance < 2.2 and distance < best_distance:
                                best_distance = distance
                                best_pair = (bonded_idx, unbonded_idx)
                
                if best_pair:
                    mol.AddBond(best_pair[0], best_pair[1], Chem.BondType.SINGLE)
                    bonded_atoms.add(best_pair[1])
                else:
                    break  # 无法找到更多连接
                    
        except Exception as e:
            logging.warning(f"侧链内部键添加失败: {e}")
    
    def _conservative_inference(self, mol, atoms_info: List[Dict]) -> Optional[Chem.Mol]:
        """
        保守的化学键推断
        
        只添加最确定的化学键
        
        Args:
            mol: 分子对象
            atoms_info: 原子信息列表
            
        Returns:
            推断后的分子对象
        """
        try:
            mol_copy = Chem.RWMol(mol)
            conf = mol.GetConformer()
            num_atoms = mol.GetNumAtoms()
            
            # 只添加高置信度的键
            for i in range(num_atoms):
                for j in range(i + 1, num_atoms):
                    atom1 = mol_copy.GetAtomWithIdx(i)
                    atom2 = mol_copy.GetAtomWithIdx(j)
                    
                    pos1 = conf.GetAtomPosition(i)
                    pos2 = conf.GetAtomPosition(j)
                    distance = pos1.Distance(pos2)
                    
                    # 只连接非常接近的原子
                    elem1, elem2 = atom1.GetSymbol(), atom2.GetSymbol()
                    
                    # 使用非常严格的距离阈值
                    if self._is_very_likely_bond(elem1, elem2, distance):
                        try:
                            mol_copy.AddBond(i, j, Chem.BondType.SINGLE)
                        except:
                            continue  # 忽略添加失败的键
            
            return mol_copy  # 返回不经过清理的分子
            
        except Exception as e:
            logging.error(f"保守推断失败: {e}")
            return None
    
    def _is_very_likely_bond(self, elem1: str, elem2: str, distance: float) -> bool:
        """
        判断是否是非常可能的化学键
        
        使用严格的标准
        """
        bond_pair = tuple(sorted([elem1, elem2]))
        standard_length = self.BOND_LENGTHS.get(bond_pair)
        
        if standard_length is None:
            # 使用共价半径
            props1 = self.ATOM_PROPERTIES.get(elem1)
            props2 = self.ATOM_PROPERTIES.get(elem2)
            if props1 and props2:
                standard_length = props1.covalent_radius + props2.covalent_radius
            else:
                return False
        
        # 使用更严格的容忍度
        tolerance = standard_length * 0.2  # 20%容忍度
        return abs(distance - standard_length) <= tolerance
    
    def validate_molecule_valences(self, mol) -> Tuple[bool, List[str]]:
        """
        验证分子中原子的价态是否合理
        
        Args:
            mol: 分子对象
            
        Returns:
            (is_valid, error_messages)
        """
        is_valid = True
        errors = []
        
        try:
            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                element = atom.GetSymbol()
                
                # 计算实际价态
                actual_valence = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
                
                # 获取理论最大价态
                props = self.ATOM_PROPERTIES.get(element)
                if props:
                    if actual_valence > props.max_valence:
                        is_valid = False
                        errors.append(f"原子{i}({element})价态超限: {actual_valence} > {props.max_valence}")
                    elif actual_valence not in props.common_valences:
                        # 警告但不认为是错误
                        errors.append(f"原子{i}({element})价态不常见: {actual_valence}, 常见价态: {props.common_valences}")
        
        except Exception as e:
            is_valid = False
            errors.append(f"价态验证失败: {e}")
        
        return is_valid, errors