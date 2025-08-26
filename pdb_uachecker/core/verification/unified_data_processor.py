"""
统一数据格式处理器
兼容不同来源的原子信息格式，提供标准化的数据接口
"""

import numpy as np
from typing import Dict, List, Optional, Union, Any, Tuple
from dataclasses import dataclass, asdict
import logging
import json
from abc import ABC, abstractmethod


@dataclass
class StandardAtomInfo:
    """标准化原子信息"""
    atom_id: str
    element: str
    x: float
    y: float
    z: float
    atom_name: Optional[str] = None
    residue_name: Optional[str] = None
    residue_number: Optional[int] = None
    chain_id: Optional[str] = None
    occupancy: float = 1.0
    temperature_factor: float = 0.0
    formal_charge: int = 0
    partial_charge: float = 0.0
    bond_partners: List[str] = None
    
    def __post_init__(self):
        if self.bond_partners is None:
            self.bond_partners = []
    
    def to_dict(self) -> Dict:
        """转换为字典格式"""
        return asdict(self)
    
    def distance_to(self, other: 'StandardAtomInfo') -> float:
        """计算与另一个原子的距离"""
        return np.sqrt(
            (self.x - other.x)**2 + 
            (self.y - other.y)**2 + 
            (self.z - other.z)**2
        )


@dataclass
class StandardResidueInfo:
    """标准化残基信息"""
    residue_name: str
    residue_number: int
    chain_id: str
    atoms: List[StandardAtomInfo]
    residue_type: Optional[str] = None  # "amino_acid", "nucleotide", "ligand", etc.
    secondary_structure: Optional[str] = None
    
    @property
    def atom_count(self) -> int:
        return len(self.atoms)
    
    @property
    def heavy_atom_count(self) -> int:
        return len([atom for atom in self.atoms if atom.element != 'H'])
    
    @property
    def elements(self) -> List[str]:
        return [atom.element for atom in self.atoms]
    
    @property
    def coordinates(self) -> np.ndarray:
        """获取所有原子坐标"""
        return np.array([[atom.x, atom.y, atom.z] for atom in self.atoms])
    
    def get_atoms_by_element(self, element: str) -> List[StandardAtomInfo]:
        """获取特定元素的原子"""
        return [atom for atom in self.atoms if atom.element == element]
    
    def get_atom_by_name(self, atom_name: str) -> Optional[StandardAtomInfo]:
        """根据原子名称获取原子"""
        for atom in self.atoms:
            if atom.atom_name == atom_name:
                return atom
        return None


class DataFormatConverter(ABC):
    """数据格式转换器基类"""
    
    @abstractmethod
    def can_handle(self, data: Any) -> bool:
        """检查是否能处理该数据格式"""
        pass
    
    @abstractmethod
    def convert_to_standard(self, data: Any) -> Union[StandardAtomInfo, StandardResidueInfo, List[StandardAtomInfo]]:
        """转换为标准格式"""
        pass
    
    @abstractmethod
    def get_format_name(self) -> str:
        """获取格式名称"""
        pass


class PDBAtomConverter(DataFormatConverter):
    """PDB原子格式转换器"""
    
    def can_handle(self, data: Any) -> bool:
        """检查是否是PDB原子格式"""
        if isinstance(data, dict):
            # 检查PDB ATOM记录的典型字段
            required_fields = {'element', 'x', 'y', 'z'}
            optional_fields = {'atom_name', 'residue_name', 'residue_number', 'chain_id'}
            
            data_keys = set(data.keys())
            return required_fields.issubset(data_keys)
        
        return False
    
    def convert_to_standard(self, data: Dict) -> StandardAtomInfo:
        """转换PDB原子格式到标准格式"""
        try:
            return StandardAtomInfo(
                atom_id=data.get('atom_id', f"{data.get('residue_name', 'UNK')}_{data.get('atom_name', 'X')}"),
                element=str(data['element']).strip(),
                x=float(data['x']),
                y=float(data['y']),
                z=float(data['z']),
                atom_name=data.get('atom_name'),
                residue_name=data.get('residue_name'),
                residue_number=data.get('residue_number'),
                chain_id=data.get('chain_id'),
                occupancy=float(data.get('occupancy', 1.0)),
                temperature_factor=float(data.get('temperature_factor', 0.0)),
                formal_charge=int(data.get('formal_charge', 0)),
                partial_charge=float(data.get('partial_charge', 0.0))
            )
        except (KeyError, ValueError, TypeError) as e:
            raise ValueError(f"PDB原子格式转换失败: {e}")
    
    def get_format_name(self) -> str:
        return "PDB_ATOM"


class RDKitAtomConverter(DataFormatConverter):
    """RDKit原子格式转换器"""
    
    def can_handle(self, data: Any) -> bool:
        """检查是否是RDKit原子对象"""
        try:
            # 检查是否有RDKit原子的典型方法
            return hasattr(data, 'GetSymbol') and hasattr(data, 'GetIdx')
        except:
            return False
    
    def convert_to_standard(self, data: Any) -> StandardAtomInfo:
        """转换RDKit原子到标准格式"""
        try:
            # 从RDKit原子获取基本信息
            atom_info = StandardAtomInfo(
                atom_id=str(data.GetIdx()),
                element=data.GetSymbol(),
                x=0.0, y=0.0, z=0.0,  # 坐标需要从构象获取
                formal_charge=data.GetFormalCharge()
            )
            
            # 尝试获取坐标信息（需要分子对象和构象）
            mol = data.GetOwningMol()
            if mol and mol.GetNumConformers() > 0:
                conf = mol.GetConformer()
                pos = conf.GetAtomPosition(data.GetIdx())
                atom_info.x = pos.x
                atom_info.y = pos.y
                atom_info.z = pos.z
            
            return atom_info
            
        except Exception as e:
            raise ValueError(f"RDKit原子格式转换失败: {e}")
    
    def get_format_name(self) -> str:
        return "RDKIT_ATOM"


class DictListConverter(DataFormatConverter):
    """字典格式残基数据转换器"""
    
    def can_handle(self, data: Any) -> bool:
        """检查是否是包含atoms列表的字典"""
        if isinstance(data, dict) and 'atoms' in data:
            atoms = data['atoms']
            if isinstance(atoms, list) and len(atoms) > 0:
                # 检查第一个原子是否是有效格式
                first_atom = atoms[0]
                return isinstance(first_atom, dict) and 'element' in first_atom
        return False
    
    def convert_to_standard(self, data: Dict) -> StandardResidueInfo:
        """转换字典格式残基到标准格式"""
        try:
            # 转换原子列表
            standard_atoms = []
            
            for atom_data in data['atoms']:
                try:
                    # 直接转换原子，避免循环依赖
                    if isinstance(atom_data, dict) and 'element' in atom_data:
                        standard_atom = StandardAtomInfo(
                            atom_id=atom_data.get('atom_name', f"ATOM_{len(standard_atoms)}"),
                            element=atom_data['element'].strip(),
                            x=float(atom_data.get('x', 0)),
                            y=float(atom_data.get('y', 0)),
                            z=float(atom_data.get('z', 0)),
                            atom_name=atom_data.get('atom_name'),
                            formal_charge=int(atom_data.get('formal_charge', 0))
                        )
                        standard_atoms.append(standard_atom)
                except Exception as e:
                    logging.warning(f"原子转换失败: {e}")
                    continue
            
            return StandardResidueInfo(
                residue_name=data.get('name', 'UNK'),
                residue_number=1,
                chain_id='A',
                atoms=standard_atoms,
                residue_type=data.get('type', 'unknown')
            )
            
        except Exception as e:
            raise ValueError(f"字典格式残基转换失败: {e}")
    
    def get_format_name(self) -> str:
        return "DICT_RESIDUE"


class ResidueInfoConverter(DataFormatConverter):
    """ResidueInfo对象转换器"""
    
    def can_handle(self, data: Any) -> bool:
        """检查是否是ResidueInfo对象"""
        return hasattr(data, 'atoms') and hasattr(data, 'residue_name')
    
    def convert_to_standard(self, data: Any) -> StandardResidueInfo:
        """转换ResidueInfo到标准格式"""
        try:
            # 转换原子列表
            standard_atoms = []
            processor = UnifiedDataProcessor()
            
            for atom in data.atoms:
                try:
                    standard_atom = processor.standardize_atom(atom)
                    if standard_atom:
                        standard_atoms.append(standard_atom)
                except Exception as e:
                    logging.warning(f"原子转换失败: {e}")
                    continue
            
            return StandardResidueInfo(
                residue_name=data.residue_name,
                residue_number=getattr(data, 'residue_number', 1),
                chain_id=getattr(data, 'chain_id', 'A'),
                atoms=standard_atoms,
                residue_type="amino_acid"  # 假设是氨基酸
            )
            
        except Exception as e:
            raise ValueError(f"ResidueInfo格式转换失败: {e}")
    
    def get_format_name(self) -> str:
        return "RESIDUE_INFO"


class AtomInfoObjectConverter(DataFormatConverter):
    """AtomInfo对象转换器"""
    
    def can_handle(self, data: Any) -> bool:
        """检查是否是AtomInfo对象"""
        return (hasattr(data, 'element') and 
                hasattr(data, 'x') and 
                hasattr(data, 'y') and 
                hasattr(data, 'z'))
    
    def convert_to_standard(self, data: Any) -> StandardAtomInfo:
        """转换AtomInfo对象到标准格式"""
        try:
            return StandardAtomInfo(
                atom_id=getattr(data, 'atom_id', f"{getattr(data, 'residue_name', 'UNK')}_{getattr(data, 'atom_name', 'X')}"),
                element=data.element.strip(),
                x=float(data.x),
                y=float(data.y),
                z=float(data.z),
                atom_name=getattr(data, 'atom_name', None),
                residue_name=getattr(data, 'residue_name', None),
                residue_number=getattr(data, 'residue_number', None),
                chain_id=getattr(data, 'chain_id', None)
            )
        except Exception as e:
            raise ValueError(f"AtomInfo对象格式转换失败: {e}")
    
    def get_format_name(self) -> str:
        return "ATOM_INFO_OBJECT"


class JSONAtomConverter(DataFormatConverter):
    """JSON格式原子转换器"""
    
    def can_handle(self, data: Any) -> bool:
        """检查是否是JSON格式的原子数据"""
        if isinstance(data, str):
            try:
                parsed = json.loads(data)
                return self._is_atom_json(parsed)
            except:
                return False
        elif isinstance(data, dict):
            return self._is_atom_json(data)
        return False
    
    def _is_atom_json(self, data: dict) -> bool:
        """检查字典是否包含原子信息"""
        required_fields = {'element', 'coordinates'}
        if required_fields.issubset(data.keys()):
            coords = data['coordinates']
            return (isinstance(coords, (list, tuple)) and 
                    len(coords) >= 3 and 
                    all(isinstance(c, (int, float)) for c in coords[:3]))
        return False
    
    def convert_to_standard(self, data: Union[str, Dict]) -> StandardAtomInfo:
        """转换JSON格式到标准格式"""
        try:
            if isinstance(data, str):
                data = json.loads(data)
            
            coords = data['coordinates']
            
            return StandardAtomInfo(
                atom_id=data.get('atom_id', data.get('id', 'unknown')),
                element=data['element'].strip(),
                x=float(coords[0]),
                y=float(coords[1]),
                z=float(coords[2]),
                atom_name=data.get('atom_name'),
                residue_name=data.get('residue_name'),
                residue_number=data.get('residue_number'),
                chain_id=data.get('chain_id'),
                formal_charge=data.get('formal_charge', 0),
                partial_charge=data.get('partial_charge', 0.0)
            )
            
        except Exception as e:
            raise ValueError(f"JSON格式转换失败: {e}")
    
    def get_format_name(self) -> str:
        return "JSON_ATOM"


class UnifiedDataProcessor:
    """统一数据处理器"""
    
    def __init__(self):
        """初始化处理器"""
        self.converters = [
            DictListConverter(),      # 字典格式残基（新增，优先级高）
            ResidueInfoConverter(),   # ResidueInfo对象
            PDBAtomConverter(),       # PDB原子格式
            RDKitAtomConverter(),     # RDKit原子
            AtomInfoObjectConverter(), # AtomInfo对象
            JSONAtomConverter()       # JSON格式
        ]
        self.format_stats = {}
    
    def register_converter(self, converter: DataFormatConverter):
        """注册新的格式转换器"""
        self.converters.insert(0, converter)  # 新转换器优先
    
    def standardize_atom(self, data: Any) -> Optional[StandardAtomInfo]:
        """
        标准化单个原子信息
        
        Args:
            data: 原子数据（任意格式）
            
        Returns:
            标准化的原子信息
        """
        for converter in self.converters:
            try:
                if converter.can_handle(data):
                    result = converter.convert_to_standard(data)
                    
                    # 更新格式统计
                    format_name = converter.get_format_name()
                    self.format_stats[format_name] = self.format_stats.get(format_name, 0) + 1
                    
                    if isinstance(result, StandardAtomInfo):
                        return result
                    elif isinstance(result, list) and result:
                        return result[0]  # 返回第一个原子
                    
            except Exception as e:
                logging.debug(f"转换器 {converter.get_format_name()} 处理失败: {e}")
                continue
        
        logging.warning(f"无法识别的原子数据格式: {type(data)}")
        return None
    
    def standardize_residue(self, data: Any) -> Optional[StandardResidueInfo]:
        """
        标准化残基信息
        
        Args:
            data: 残基数据（任意格式）
            
        Returns:
            标准化的残基信息
        """
        for converter in self.converters:
            try:
                if converter.can_handle(data):
                    result = converter.convert_to_standard(data)
                    
                    if isinstance(result, StandardResidueInfo):
                        format_name = converter.get_format_name()
                        self.format_stats[format_name] = self.format_stats.get(format_name, 0) + 1
                        return result
                    
            except Exception as e:
                logging.debug(f"残基转换器 {converter.get_format_name()} 处理失败: {e}")
                continue
        
        # 尝试将数据作为原子列表处理
        if isinstance(data, (list, tuple)):
            atoms = []
            for item in data:
                atom = self.standardize_atom(item)
                if atom:
                    atoms.append(atom)
            
            if atoms:
                # 从原子信息推断残基信息
                return self._infer_residue_from_atoms(atoms)
        
        logging.warning(f"无法识别的残基数据格式: {type(data)}")
        return None
    
    def _infer_residue_from_atoms(self, atoms: List[StandardAtomInfo]) -> StandardResidueInfo:
        """从原子列表推断残基信息"""
        # 尝试从原子信息中获取残基信息
        residue_name = "UNK"
        residue_number = 1
        chain_id = "A"
        
        for atom in atoms:
            if atom.residue_name:
                residue_name = atom.residue_name
            if atom.residue_number:
                residue_number = atom.residue_number
            if atom.chain_id:
                chain_id = atom.chain_id
        
        return StandardResidueInfo(
            residue_name=residue_name,
            residue_number=residue_number,
            chain_id=chain_id,
            atoms=atoms
        )
    
    def batch_standardize_atoms(self, data_list: List[Any]) -> List[StandardAtomInfo]:
        """
        批量标准化原子信息
        
        Args:
            data_list: 原子数据列表
            
        Returns:
            标准化原子列表
        """
        results = []
        
        for data in data_list:
            try:
                atom = self.standardize_atom(data)
                if atom:
                    results.append(atom)
            except Exception as e:
                logging.warning(f"批量处理中原子转换失败: {e}")
                continue
        
        return results
    
    def convert_to_legacy_format(self, standard_data: Union[StandardAtomInfo, StandardResidueInfo], 
                                target_format: str = "dict") -> Dict:
        """
        转换标准格式到传统格式
        
        Args:
            standard_data: 标准化数据
            target_format: 目标格式
            
        Returns:
            转换后的数据
        """
        if isinstance(standard_data, StandardAtomInfo):
            return self._atom_to_legacy_format(standard_data, target_format)
        elif isinstance(standard_data, StandardResidueInfo):
            return self._residue_to_legacy_format(standard_data, target_format)
        else:
            raise ValueError(f"不支持的数据类型: {type(standard_data)}")
    
    def _atom_to_legacy_format(self, atom: StandardAtomInfo, target_format: str) -> Dict:
        """将标准原子转换为传统格式"""
        if target_format == "dict":
            return {
                'element': atom.element,
                'x': atom.x,
                'y': atom.y,
                'z': atom.z,
                'atom_name': atom.atom_name,
                'residue_name': atom.residue_name,
                'residue_number': atom.residue_number,
                'chain_id': atom.chain_id
            }
        elif target_format == "pdb":
            return {
                'element': atom.element,
                'x': atom.x,
                'y': atom.y,
                'z': atom.z,
                'atom_name': atom.atom_name,
                'residue_name': atom.residue_name,
                'residue_number': atom.residue_number,
                'chain_id': atom.chain_id,
                'occupancy': atom.occupancy,
                'temperature_factor': atom.temperature_factor
            }
        else:
            return atom.to_dict()
    
    def _residue_to_legacy_format(self, residue: StandardResidueInfo, target_format: str) -> Dict:
        """将标准残基转换为传统格式"""
        atoms_legacy = [
            self._atom_to_legacy_format(atom, target_format) 
            for atom in residue.atoms
        ]
        
        return {
            'residue_name': residue.residue_name,
            'residue_number': residue.residue_number,
            'chain_id': residue.chain_id,
            'atoms': atoms_legacy,
            'atom_count': residue.atom_count,
            'heavy_atom_count': residue.heavy_atom_count
        }
    
    def get_format_statistics(self) -> Dict[str, int]:
        """获取格式使用统计"""
        return self.format_stats.copy()
    
    def validate_standard_data(self, data: Union[StandardAtomInfo, StandardResidueInfo]) -> Tuple[bool, List[str]]:
        """
        验证标准化数据的完整性
        
        Args:
            data: 标准化数据
            
        Returns:
            (is_valid, error_messages)
        """
        errors = []
        
        if isinstance(data, StandardAtomInfo):
            # 验证原子数据
            if not data.element or not data.element.strip():
                errors.append("元素符号为空")
            
            if not isinstance(data.x, (int, float)) or not isinstance(data.y, (int, float)) or not isinstance(data.z, (int, float)):
                errors.append("坐标数据类型错误")
            
            if abs(data.x) > 1000 or abs(data.y) > 1000 or abs(data.z) > 1000:
                errors.append("坐标数值异常（可能是单位错误）")
            
        elif isinstance(data, StandardResidueInfo):
            # 验证残基数据
            if not data.residue_name or not data.residue_name.strip():
                errors.append("残基名称为空")
            
            if not data.atoms:
                errors.append("残基不包含原子")
            
            # 验证每个原子
            for i, atom in enumerate(data.atoms):
                atom_valid, atom_errors = self.validate_standard_data(atom)
                if not atom_valid:
                    errors.extend([f"原子{i}: {err}" for err in atom_errors])
        
        return len(errors) == 0, errors
    
    def create_compatibility_layer(self, standard_residue: StandardResidueInfo):
        """
        为标准残基创建兼容层
        
        提供与旧接口兼容的属性访问
        """
        class CompatibilityResidueInfo:
            def __init__(self, standard_residue: StandardResidueInfo):
                self._standard = standard_residue
            
            @property
            def atoms(self):
                """返回兼容格式的原子列表"""
                return [self._create_compatible_atom(atom) for atom in self._standard.atoms]
            
            @property
            def residue_name(self):
                return self._standard.residue_name
            
            @property
            def residue_number(self):
                return self._standard.residue_number
            
            @property
            def chain_id(self):
                return self._standard.chain_id
            
            def _create_compatible_atom(self, atom: StandardAtomInfo):
                """创建兼容的原子对象"""
                class CompatibleAtom:
                    def __init__(self, standard_atom: StandardAtomInfo):
                        self.element = standard_atom.element
                        self.x = standard_atom.x
                        self.y = standard_atom.y
                        self.z = standard_atom.z
                        self.atom_name = standard_atom.atom_name
                        self.residue_name = standard_atom.residue_name
                        self.residue_number = standard_atom.residue_number
                        self.chain_id = standard_atom.chain_id
                    
                    def get(self, key, default=None):
                        """字典式访问"""
                        return getattr(self, key, default)
                
                return CompatibleAtom(atom)
        
        return CompatibilityResidueInfo(standard_residue)