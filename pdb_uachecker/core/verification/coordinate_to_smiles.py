"""
从PDB坐标生成SMILES字符串的增强模块
支持立体化学信息保留和异构体识别
"""

import numpy as np
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass
import logging

from .unified_data_processor import UnifiedDataProcessor, StandardResidueInfo
from .advanced_bond_inference import AdvancedBondInference
from .valence_validator import SmartValenceValidator

try:
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds, rdMolAlign, AllChem
    from rdkit.Geometry import Point3D
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    logging.warning("RDKit未安装，坐标到SMILES转换功能受限")


@dataclass
class AtomInfo:
    """原子信息"""
    element: str
    x: float
    y: float  
    z: float
    formal_charge: int = 0
    atom_id: Optional[str] = None


@dataclass
class BondInfo:
    """化学键信息"""
    atom1_idx: int
    atom2_idx: int
    bond_type: str = "SINGLE"  # SINGLE, DOUBLE, TRIPLE, AROMATIC
    distance: float = 0.0


class CoordinateToSmilesConverter:
    """从PDB坐标生成SMILES的转换器"""
    
    # 标准原子半径（Å）- 用于化学键判断
    COVALENT_RADII = {
        'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05,
        'P': 1.07, 'F': 0.57, 'Cl': 0.99, 'Br': 1.14, 'I': 1.33,
        'Se': 1.20
    }
    
    # 化学键长度阈值（倍数）
    BOND_TOLERANCE = 1.3
    
    def __init__(self, use_rdkit: bool = True):
        """
        初始化转换器
        
        Args:
            use_rdkit: 是否使用RDKit进行转换
        """
        self.use_rdkit = use_rdkit and RDKIT_AVAILABLE
        if not self.use_rdkit:
            logging.warning("RDKit不可用，使用简化版SMILES生成")
        
        # 初始化辅助组件
        self.data_processor = UnifiedDataProcessor()
        self.bond_inference = AdvancedBondInference()
        self.valence_validator = SmartValenceValidator(tolerance_level="moderate")
    
    def convert_residue_to_smiles(self, residue_info, include_stereochemistry: bool = True) -> Optional[str]:
        """
        将残基信息转换为SMILES字符串（增强版）
        
        Args:
            residue_info: 残基信息对象（任意格式）
            include_stereochemistry: 是否包含立体化学信息
            
        Returns:
            SMILES字符串，失败返回None
        """
        try:
            # 标准化残基数据
            standard_residue = self.data_processor.standardize_residue(residue_info)
            if not standard_residue or not standard_residue.atoms:
                logging.warning("残基数据标准化失败或无原子")
                return None
            
            # 验证数据完整性
            is_valid, errors = self.data_processor.validate_standard_data(standard_residue)
            if not is_valid:
                logging.warning(f"残基数据验证失败: {errors}")
            
            if self.use_rdkit:
                return self._convert_with_rdkit_enhanced(standard_residue, include_stereochemistry)
            else:
                return self._convert_with_basic_algorithm_enhanced(standard_residue)
                
        except Exception as e:
            logging.error(f"增强坐标到SMILES转换失败: {e}")
            return None
    
    def _extract_atoms_from_residue(self, residue_info) -> List[AtomInfo]:
        """
        从残基信息中提取原子信息
        
        Args:
            residue_info: 残基信息对象
            
        Returns:
            原子信息列表
        """
        atoms = []
        
        for atom in residue_info.atoms:
            # 处理不同的原子数据格式
            if hasattr(atom, 'element'):
                # 原子对象格式
                element = atom.element.strip() if atom.element else ''
                if not element:
                    continue
                
                atom_info = AtomInfo(
                    element=element,
                    x=float(getattr(atom, 'x', 0)),
                    y=float(getattr(atom, 'y', 0)),
                    z=float(getattr(atom, 'z', 0)),
                    atom_id=getattr(atom, 'atom_name', '')
                )
            else:
                # 字典格式
                element = atom.get('element', '').strip()
                if not element:
                    continue
                    
                atom_info = AtomInfo(
                    element=element,
                    x=float(atom.get('x', 0)),
                    y=float(atom.get('y', 0)),
                    z=float(atom.get('z', 0)),
                    atom_id=atom.get('atom_name', '')
                )
            
            atoms.append(atom_info)
        
        return atoms
    
    def _convert_with_rdkit(self, atoms: List[AtomInfo], include_stereochemistry: bool = True) -> Optional[str]:
        """
        使用RDKit进行坐标到SMILES的转换
        
        Args:
            atoms: 原子信息列表
            include_stereochemistry: 是否包含立体化学
            
        Returns:
            SMILES字符串
        """
        try:
            # 创建RDKit分子对象
            mol = Chem.RWMol()
            
            # 添加原子
            atom_idx_map = {}
            for i, atom_info in enumerate(atoms):
                atom = Chem.Atom(atom_info.element)
                atom.SetFormalCharge(atom_info.formal_charge)
                idx = mol.AddAtom(atom)
                atom_idx_map[i] = idx
            
            # 设置3D坐标
            conf = Chem.Conformer(len(atoms))
            for i, atom_info in enumerate(atoms):
                point = Point3D(atom_info.x, atom_info.y, atom_info.z)
                conf.SetAtomPosition(atom_idx_map[i], point)
            
            mol.AddConformer(conf)
            
            # 推断化学键
            mol = self._infer_bonds_rdkit(mol)
            if mol is None:
                return None
            
            # 生成SMILES
            if include_stereochemistry:
                smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            else:
                smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
            
            return smiles
            
        except Exception as e:
            logging.error(f"RDKit转换失败: {e}")
            return None
    
    def _infer_bonds_rdkit(self, mol) -> Optional[Chem.Mol]:
        """
        使用RDKit推断化学键
        
        Args:
            mol: RDKit分子对象
            
        Returns:
            推断化学键后的分子对象
        """
        try:
            # 方法1: 使用RDKit的自动键推断
            rdDetermineBonds.DetermineBonds(mol)
            
            # 验证分子的有效性
            Chem.SanitizeMol(mol)
            
            return mol
            
        except Exception as e:
            logging.warning(f"RDKit自动键推断失败: {e}")
            
            # 方法2: 基于距离的化学键推断
            try:
                mol_copy = Chem.RWMol(mol)
                bonds = self._infer_bonds_by_distance(mol)
                
                for bond in bonds:
                    mol_copy.AddBond(bond.atom1_idx, bond.atom2_idx, 
                                   self._get_rdkit_bond_type(bond.bond_type))
                
                Chem.SanitizeMol(mol_copy)
                return mol_copy
                
            except Exception as e2:
                logging.error(f"基于距离的键推断也失败: {e2}")
                return None
    
    def _infer_bonds_by_distance(self, mol) -> List[BondInfo]:
        """
        基于原子间距离推断化学键
        
        Args:
            mol: RDKit分子对象
            
        Returns:
            推断的化学键列表
        """
        bonds = []
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                atom1 = mol.GetAtomWithIdx(i)
                atom2 = mol.GetAtomWithIdx(j)
                
                # 计算距离
                pos1 = conf.GetAtomPosition(i)
                pos2 = conf.GetAtomPosition(j)
                distance = pos1.Distance(pos2)
                
                # 判断是否形成化学键
                if self._is_bonded(atom1.GetSymbol(), atom2.GetSymbol(), distance):
                    bond_type = self._determine_bond_type(atom1, atom2, distance)
                    bonds.append(BondInfo(i, j, bond_type, distance))
        
        return bonds
    
    def _is_bonded(self, element1: str, element2: str, distance: float) -> bool:
        """
        判断两个原子是否形成化学键
        
        Args:
            element1: 第一个原子的元素
            element2: 第二个原子的元素
            distance: 原子间距离
            
        Returns:
            是否形成化学键
        """
        radius1 = self.COVALENT_RADII.get(element1, 1.0)
        radius2 = self.COVALENT_RADII.get(element2, 1.0)
        
        # 化学键长度阈值
        max_bond_length = (radius1 + radius2) * self.BOND_TOLERANCE
        
        return distance <= max_bond_length
    
    def _determine_bond_type(self, atom1, atom2, distance: float) -> str:
        """
        确定化学键类型
        
        Args:
            atom1: 第一个原子
            atom2: 第二个原子
            distance: 原子间距离
            
        Returns:
            化学键类型
        """
        # 简化的键类型判断逻辑
        # 实际应用中可能需要更复杂的判断
        
        element1 = atom1.GetSymbol()
        element2 = atom2.GetSymbol()
        
        # 基于距离的简单判断
        radius1 = self.COVALENT_RADII.get(element1, 1.0)
        radius2 = self.COVALENT_RADII.get(element2, 1.0)
        expected_single = radius1 + radius2
        
        if distance < expected_single * 0.9:
            # 可能是双键或三键
            if (element1 == 'C' and element2 == 'C') or \
               (element1 == 'C' and element2 == 'N') or \
               (element1 == 'C' and element2 == 'O'):
                return "DOUBLE" if distance < expected_single * 0.8 else "SINGLE"
        
        return "SINGLE"
    
    def _get_rdkit_bond_type(self, bond_type: str):
        """获取RDKit化学键类型"""
        bond_type_map = {
            "SINGLE": Chem.BondType.SINGLE,
            "DOUBLE": Chem.BondType.DOUBLE,
            "TRIPLE": Chem.BondType.TRIPLE,
            "AROMATIC": Chem.BondType.AROMATIC
        }
        return bond_type_map.get(bond_type, Chem.BondType.SINGLE)
    
    def _convert_with_basic_algorithm(self, atoms: List[AtomInfo]) -> Optional[str]:
        """
        使用基础算法生成简化SMILES
        
        Args:
            atoms: 原子信息列表
            
        Returns:
            简化的SMILES字符串
        """
        # 基础实现：仅支持简单的分子结构
        # 实际应用中建议使用RDKit
        
        try:
            # 统计原子组成
            composition = {}
            for atom in atoms:
                composition[atom.element] = composition.get(atom.element, 0) + 1
            
            # 生成简单的分子式表示（不是真正的SMILES）
            formula_parts = []
            for element in ['C', 'H', 'N', 'O', 'S', 'P']:
                if element in composition:
                    count = composition[element]
                    if count == 1:
                        formula_parts.append(element)
                    else:
                        formula_parts.append(f"{element}{count}")
            
            return ''.join(formula_parts) if formula_parts else None
            
        except Exception as e:
            logging.error(f"基础算法转换失败: {e}")
            return None
    
    def _convert_with_rdkit_enhanced(self, standard_residue: StandardResidueInfo, include_stereochemistry: bool = True) -> Optional[str]:
        """
        使用RDKit进行增强的坐标到SMILES转换
        
        Args:
            standard_residue: 标准化残基信息
            include_stereochemistry: 是否包含立体化学
            
        Returns:
            SMILES字符串
        """
        try:
            # 创建RDKit分子对象
            mol = self._create_rdkit_molecule_enhanced(standard_residue)
            if mol is None:
                logging.warning("RDKit分子对象创建失败")
                return None
            
            # 使用高级化学键推断
            mol = self.bond_inference.infer_bonds_advanced(mol, self._standard_atoms_to_dict_list(standard_residue.atoms))
            if mol is None:
                logging.warning("化学键推断失败")
                return None
            
            # 价态验证和修正
            is_valid, issues, report = self.valence_validator.validate_molecule_valences(mol)
            if not is_valid:
                logging.info(f"分子价态问题: {len(issues)}个问题")
                # 尝试自动修正
                mol, corrections = self.valence_validator.apply_automatic_corrections(mol, issues)
                if corrections:
                    logging.info(f"应用了{len(corrections)}个价态修正")
            
            # 生成SMILES
            try:
                if include_stereochemistry:
                    smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                else:
                    smiles = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
                
                # 验证生成的SMILES
                if self.validate_generated_smiles(smiles):
                    return smiles
                else:
                    logging.warning("生成的SMILES验证失败")
                    return None
                    
            except Exception as smiles_error:
                logging.error(f"SMILES生成失败: {smiles_error}")
                return None
            
        except Exception as e:
            logging.error(f"RDKit增强转换失败: {e}")
            return None
    
    def _create_rdkit_molecule_enhanced(self, standard_residue: StandardResidueInfo):
        """
        从标准化残基创建RDKit分子对象（增强版）
        
        Args:
            standard_residue: 标准化残基信息
            
        Returns:
            RDKit分子对象
        """
        try:
            mol = Chem.RWMol()
            atom_map = {}
            
            # 添加原子
            for i, std_atom in enumerate(standard_residue.atoms):
                try:
                    rd_atom = Chem.Atom(std_atom.element)
                    
                    # 设置形式电荷
                    if std_atom.formal_charge != 0:
                        rd_atom.SetFormalCharge(std_atom.formal_charge)
                    
                    atom_idx = mol.AddAtom(rd_atom)
                    atom_map[i] = atom_idx
                    
                except Exception as atom_error:
                    logging.warning(f"原子{i}({std_atom.element})添加失败: {atom_error}")
                    continue
            
            # 设置3D坐标
            if len(atom_map) > 0:
                conf = Chem.Conformer(len(atom_map))
                
                for i, std_atom in enumerate(standard_residue.atoms):
                    if i in atom_map:
                        point = Point3D(std_atom.x, std_atom.y, std_atom.z)
                        conf.SetAtomPosition(atom_map[i], point)
                
                mol.AddConformer(conf)
            
            return mol
            
        except Exception as e:
            logging.error(f"增强分子对象创建失败: {e}")
            return None
    
    def _standard_atoms_to_dict_list(self, atoms: List) -> List[Dict]:
        """将标准原子列表转换为字典列表"""
        return [
            {
                'element': atom.element,
                'x': atom.x,
                'y': atom.y,
                'z': atom.z,
                'atom_name': atom.atom_name,
                'formal_charge': atom.formal_charge
            }
            for atom in atoms
        ]
    
    def _convert_with_basic_algorithm_enhanced(self, standard_residue: StandardResidueInfo) -> Optional[str]:
        """
        使用基础算法生成增强SMILES（当RDKit不可用时）
        
        Args:
            standard_residue: 标准化残基信息
            
        Returns:
            简化的分子表示字符串
        """
        try:
            # 统计原子组成
            composition = {}
            for atom in standard_residue.atoms:
                element = atom.element
                composition[element] = composition.get(element, 0) + 1
            
            # 生成改进的分子式表示
            formula_parts = []
            
            # 按标准顺序排列元素
            element_order = ['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Se']
            
            for element in element_order:
                if element in composition:
                    count = composition[element]
                    if count == 1:
                        formula_parts.append(element)
                    else:
                        formula_parts.append(f"{element}{count}")
            
            # 添加其他元素
            other_elements = sorted(set(composition.keys()) - set(element_order))
            for element in other_elements:
                count = composition[element]
                if count == 1:
                    formula_parts.append(element)
                else:
                    formula_parts.append(f"{element}{count}")
            
            molecular_formula = ''.join(formula_parts)
            
            # 添加结构信息前缀，表明这不是真正的SMILES
            return f"FORMULA:{molecular_formula}" if molecular_formula else None
            
        except Exception as e:
            logging.error(f"基础增强算法转换失败: {e}")
            return None
    
    def generate_stereochemistry_info(self, atoms: List[AtomInfo]) -> Dict[str, str]:
        """
        生成立体化学信息
        
        Args:
            atoms: 原子信息列表
            
        Returns:
            立体化学信息字典
        """
        stereo_info = {}
        
        if not self.use_rdkit:
            return stereo_info
        
        try:
            # 创建分子并推断化学键
            mol = self._create_mol_from_atoms(atoms)
            if mol is None:
                return stereo_info
            
            # 检测手性中心
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            
            for center_idx, chirality in chiral_centers:
                atom = mol.GetAtomWithIdx(center_idx)
                stereo_info[f"chiral_center_{center_idx}"] = {
                    'atom_symbol': atom.GetSymbol(),
                    'chirality': chirality,
                    'coordinates': self._get_atom_coordinates(atoms, center_idx)
                }
            
            return stereo_info
            
        except Exception as e:
            logging.error(f"立体化学信息生成失败: {e}")
            return stereo_info
    
    def _create_mol_from_atoms(self, atoms: List[AtomInfo]):
        """从原子列表创建RDKit分子对象"""
        if not self.use_rdkit:
            return None
        
        try:
            mol = Chem.RWMol()
            
            # 添加原子
            for atom_info in atoms:
                atom = Chem.Atom(atom_info.element)
                mol.AddAtom(atom)
            
            # 设置坐标
            conf = Chem.Conformer(len(atoms))
            for i, atom_info in enumerate(atoms):
                point = Point3D(atom_info.x, atom_info.y, atom_info.z)
                conf.SetAtomPosition(i, point)
            
            mol.AddConformer(conf)
            
            # 推断化学键
            return self._infer_bonds_rdkit(mol)
            
        except Exception as e:
            logging.error(f"分子对象创建失败: {e}")
            return None
    
    def _get_atom_coordinates(self, atoms: List[AtomInfo], index: int) -> Tuple[float, float, float]:
        """获取指定原子的坐标"""
        if 0 <= index < len(atoms):
            atom = atoms[index]
            return (atom.x, atom.y, atom.z)
        return (0.0, 0.0, 0.0)
    
    def validate_generated_smiles(self, smiles: str) -> bool:
        """
        验证生成的SMILES字符串的有效性
        
        Args:
            smiles: SMILES字符串
            
        Returns:
            是否有效
        """
        if not smiles or not self.use_rdkit:
            return bool(smiles)
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
            
        except Exception:
            return False


class EnhancedFingerprintGenerator:
    """增强的分子指纹生成器"""
    
    def __init__(self, coordinate_converter: Optional[CoordinateToSmilesConverter] = None):
        """
        初始化指纹生成器
        
        Args:
            coordinate_converter: 坐标转换器
        """
        self.coordinate_converter = coordinate_converter or CoordinateToSmilesConverter()
    
    def generate_fingerprint_from_residue(self, residue_info, fingerprint_type: str = 'morgan') -> Optional[str]:
        """
        从残基信息生成分子指纹
        
        Args:
            residue_info: 残基信息
            fingerprint_type: 指纹类型 ('morgan', 'ecfp2', 'ecfp4')
            
        Returns:
            指纹字符串
        """
        try:
            # 首先尝试从坐标生成SMILES
            smiles = self.coordinate_converter.convert_residue_to_smiles(residue_info)
            
            if smiles:
                return self.generate_fingerprint_from_smiles(smiles, fingerprint_type)
            else:
                # 如果无法生成SMILES，使用原子组成生成简化指纹
                return self.generate_composition_fingerprint(residue_info)
                
        except Exception as e:
            logging.error(f"从残基生成指纹失败: {e}")
            return None
    
    def generate_fingerprint_from_smiles(self, smiles: str, fingerprint_type: str = 'morgan') -> Optional[str]:
        """
        从SMILES生成分子指纹
        
        Args:
            smiles: SMILES字符串
            fingerprint_type: 指纹类型
            
        Returns:
            指纹字符串
        """
        if not RDKIT_AVAILABLE:
            logging.warning("RDKit不可用，无法生成分子指纹")
            return None
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            if fingerprint_type.lower() == 'morgan' or fingerprint_type.lower() == 'ecfp2':
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=1, nBits=1024)
            elif fingerprint_type.lower() == 'ecfp4':
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)
            else:
                # 默认使用Morgan指纹
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=1, nBits=1024)
            
            return fp.ToBitString()
            
        except Exception as e:
            logging.error(f"从SMILES生成指纹失败: {e}")
            return None
    
    def generate_composition_fingerprint(self, residue_info) -> str:
        """
        基于原子组成生成简化指纹
        
        Args:
            residue_info: 残基信息
            
        Returns:
            组成指纹字符串
        """
        try:
            # 统计重原子组成
            composition = {}
            for atom in residue_info.atoms:
                element = atom.get('element', '').strip()
                if element and element != 'H':  # 排除氢原子
                    composition[element] = composition.get(element, 0) + 1
            
            # 生成组成指纹（排序的元素-数量对）
            sorted_elements = sorted(composition.items())
            fingerprint = '_'.join([f"{elem}:{count}" for elem, count in sorted_elements])
            
            return fingerprint
            
        except Exception as e:
            logging.error(f"生成组成指纹失败: {e}")
            return ""