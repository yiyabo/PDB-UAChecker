"""
化学计算工具模块
提供分子式计算、原子组成分析等基础化学功能
"""

import re
from typing import Dict, List, Optional, Tuple
from collections import Counter


class ChemistryUtils:
    """化学计算工具类"""
    
    # 标准元素顺序（用于分子式排列）
    ELEMENT_ORDER = ['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Se']
    
    # 元素权重（用于加权相似性计算）
    ELEMENT_WEIGHTS = {
        'C': 1.0,   # 碳骨架最重要
        'N': 0.9,   # 氮原子很重要（氨基酸特征）
        'O': 0.9,   # 氧原子很重要（羧基等）
        'S': 0.8,   # 硫原子重要
        'P': 0.8,   # 磷原子重要
        'F': 0.8,   # 氟原子重要（强电负性）
        'Cl': 0.7, 'Br': 0.6, 'I': 0.5,  # 卤素原子
        'Se': 0.7   # 硒原子
    }
    
    @staticmethod
    def calculate_atom_composition(atoms: List[Dict], include_hydrogen: Optional[bool] = None) -> Dict[str, int]:
        """
        计算原子组成
        
        Args:
            atoms: 原子信息列表，每个原子包含'element'字段
            include_hydrogen: 是否包含氢原子。None表示自动检测
        
        Returns:
            原子组成字典 {元素: 数量}
        """
        composition = {}
        
        for atom in atoms:
            element = atom.get('element', '').strip()
            if not element:
                continue
            
            # 智能氢原子处理
            if element == 'H' and include_hydrogen is False:
                continue
            
            composition[element] = composition.get(element, 0) + 1
        
        return composition
    
    @staticmethod
    def calculate_molecular_formula(atom_composition: Dict[str, int]) -> str:
        """
        从原子组成计算分子式
        
        Args:
            atom_composition: 原子组成字典
        
        Returns:
            标准格式的分子式字符串
        """
        if not atom_composition:
            return ""
        
        formula_parts = []
        
        # 按标准顺序排列元素
        for element in ChemistryUtils.ELEMENT_ORDER:
            if element in atom_composition:
                count = atom_composition[element]
                if count == 1:
                    formula_parts.append(element)
                else:
                    formula_parts.append(f"{element}{count}")
        
        # 添加其他元素（按字母顺序）
        other_elements = sorted(set(atom_composition.keys()) - set(ChemistryUtils.ELEMENT_ORDER))
        for element in other_elements:
            count = atom_composition[element]
            if count == 1:
                formula_parts.append(element)
            else:
                formula_parts.append(f"{element}{count}")
        
        return ''.join(formula_parts)
    
    @staticmethod
    def remove_hydrogen_from_formula(formula: str) -> str:
        """
        从分子式中移除氢原子
        
        Args:
            formula: 原始分子式
        
        Returns:
            移除氢原子后的分子式
        """
        if not formula:
            return ""
        
        # 移除H和H后面的数字
        formula_no_h = re.sub(r'H\d*', '', formula)
        return formula_no_h
    
    @staticmethod
    def parse_molecular_formula(formula: str) -> Dict[str, int]:
        """
        解析分子式为原子组成
        
        Args:
            formula: 分子式字符串，如"C3H7NO2"
        
        Returns:
            原子组成字典
        """
        if not formula:
            return {}
        
        # 正则表达式匹配元素和数量
        pattern = r'([A-Z][a-z]?)(\d*)'
        matches = re.findall(pattern, formula)
        
        composition = {}
        for element, count_str in matches:
            count = int(count_str) if count_str else 1
            composition[element] = composition.get(element, 0) + count
        
        return composition
    
    @staticmethod
    def calculate_jaccard_similarity(comp1: Dict[str, int], comp2: Dict[str, int]) -> float:
        """
        计算两个原子组成的Jaccard相似性
        
        Args:
            comp1, comp2: 原子组成字典
        
        Returns:
            Jaccard相似性 (0.0-1.0)
        """
        if not comp1 and not comp2:
            return 1.0
        if not comp1 or not comp2:
            return 0.0
        
        # 获取所有元素
        all_elements = set(comp1.keys()) | set(comp2.keys())
        
        intersection = 0
        union = 0
        
        for element in all_elements:
            count1 = comp1.get(element, 0)
            count2 = comp2.get(element, 0)
            
            intersection += min(count1, count2)
            union += max(count1, count2)
        
        return intersection / union if union > 0 else 0.0
    
    @staticmethod
    def calculate_weighted_jaccard_similarity(comp1: Dict[str, int], comp2: Dict[str, int], 
                                            weights: Optional[Dict[str, float]] = None) -> float:
        """
        计算加权Jaccard相似性
        
        Args:
            comp1, comp2: 原子组成字典
            weights: 元素权重字典，默认使用ELEMENT_WEIGHTS
        
        Returns:
            加权Jaccard相似性 (0.0-1.0)
        """
        if weights is None:
            weights = ChemistryUtils.ELEMENT_WEIGHTS
        
        if not comp1 and not comp2:
            return 1.0
        if not comp1 or not comp2:
            return 0.0
        
        all_elements = set(comp1.keys()) | set(comp2.keys())
        
        weighted_intersection = 0.0
        weighted_union = 0.0
        
        for element in all_elements:
            count1 = comp1.get(element, 0)
            count2 = comp2.get(element, 0)
            weight = weights.get(element, 0.5)  # 默认权重0.5
            
            weighted_intersection += min(count1, count2) * weight
            weighted_union += max(count1, count2) * weight
        
        return weighted_intersection / weighted_union if weighted_union > 0 else 0.0
    
    @staticmethod
    def extract_heavy_atoms(composition: Dict[str, int]) -> Dict[str, int]:
        """
        提取重原子组成（排除氢原子）
        
        Args:
            composition: 原子组成字典
        
        Returns:
            重原子组成字典
        """
        return {element: count for element, count in composition.items() if element != 'H'}
    
    @staticmethod
    def is_valid_element(element: str) -> bool:
        """
        检查是否为有效的化学元素符号
        
        Args:
            element: 元素符号
        
        Returns:
            是否有效
        """
        # 简单的元素符号验证
        if not element:
            return False
        
        # 第一个字符必须是大写字母
        if not element[0].isupper():
            return False
        
        # 如果有第二个字符，必须是小写字母
        if len(element) > 1 and not element[1].islower():
            return False
        
        # 长度不能超过2
        if len(element) > 2:
            return False
        
        return True
    
    @staticmethod
    def normalize_atom_name(atom_name: str) -> str:
        """
        标准化原子名称，提取元素符号
        
        Args:
            atom_name: 原子名称，如"CA", "CB1", "N1"
        
        Returns:
            标准化的元素符号
        """
        if not atom_name:
            return ""
        
        atom_name = atom_name.strip()
        
        # PDB原子名称映射表（基于化学知识）
        pdb_atom_mapping = {
            # 碳原子
            'CA': 'C', 'CB': 'C', 'CG': 'C', 'CG1': 'C', 'CG2': 'C',
            'CD': 'C', 'CD1': 'C', 'CD2': 'C', 'CE': 'C', 'CE1': 'C', 
            'CE2': 'C', 'CE3': 'C', 'CZ': 'C', 'CZ2': 'C', 'CZ3': 'C',
            'CH2': 'C', 'C': 'C',
            # 氮原子
            'N': 'N', 'ND1': 'N', 'ND2': 'N', 'NE': 'N', 'NE1': 'N', 
            'NE2': 'N', 'NH1': 'N', 'NH2': 'N', 'NZ': 'N',
            # 氧原子
            'O': 'O', 'OG': 'O', 'OG1': 'O', 'OH': 'O', 'OD1': 'O', 
            'OD2': 'O', 'OE1': 'O', 'OE2': 'O', 'OXT': 'O',
            # 硫原子
            'S': 'S', 'SG': 'S', 'SD': 'S',
            # 磷原子
            'P': 'P',
            # 氢原子
            'H': 'H', 'HA': 'H', 'HB': 'H', 'HG': 'H', 'HD': 'H', 
            'HE': 'H', 'HZ': 'H', 'HH': 'H',
        }
        
        # 首先尝试精确匹配
        if atom_name in pdb_atom_mapping:
            return pdb_atom_mapping[atom_name]
        
        # 处理带数字的原子名（如HB1, HB2, HB3）
        base_name = ''.join(c for c in atom_name if not c.isdigit())
        if base_name in pdb_atom_mapping:
            return pdb_atom_mapping[base_name]
        
        # 通用规则：提取元素符号
        if len(atom_name) >= 2 and atom_name[1].islower():
            # 双字符元素，如Cl, Br
            element = atom_name[:2]
        else:
            # 单字符元素，如C, N, O
            element = atom_name[0]
        
        return element if ChemistryUtils.is_valid_element(element) else ""


class MolecularFingerprint:
    """分子指纹计算工具"""
    
    @staticmethod
    def calculate_morgan_fingerprint(smiles: str, radius: int = 2, n_bits: int = 1024) -> Optional[str]:
        """
        计算Morgan指纹（ECFP）
        
        Args:
            smiles: SMILES字符串
            radius: 指纹半径
            n_bits: 指纹位数
        
        Returns:
            指纹位串，失败返回None
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors
            
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            return fp.ToBitString()
        
        except ImportError:
            print("⚠️ RDKit不可用，无法计算分子指纹")
            return None
        except Exception as e:
            print(f"⚠️ 分子指纹计算失败: {e}")
            return None
    
    @staticmethod
    def calculate_tanimoto_similarity(fp1: str, fp2: str) -> float:
        """
        计算Tanimoto相似性
        
        Args:
            fp1, fp2: 指纹位串
        
        Returns:
            Tanimoto相似性 (0.0-1.0)
        """
        if not fp1 or not fp2 or len(fp1) != len(fp2):
            return 0.0
        
        try:
            bits1 = [int(b) for b in fp1]
            bits2 = [int(b) for b in fp2]
            
            intersection = sum(a & b for a, b in zip(bits1, bits2))
            union = sum(a | b for a, b in zip(bits1, bits2))
            
            return intersection / union if union > 0 else 0.0
        
        except Exception:
            return 0.0
    
    @staticmethod
    def generate_molecular_formula_from_smiles(smiles: str) -> str:
        """
        从SMILES生成分子式
        
        Args:
            smiles: SMILES字符串
        
        Returns:
            分子式字符串
        """
        try:
            from rdkit import Chem
            
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return ""
            
            # 计算原子组成
            atom_counts = {}
            for atom in mol.GetAtoms():
                symbol = atom.GetSymbol()
                atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
            
            # 添加隐式氢原子
            total_h = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
            if total_h > 0:
                atom_counts['H'] = atom_counts.get('H', 0) + total_h
            
            return ChemistryUtils.calculate_molecular_formula(atom_counts)
        
        except ImportError:
            print("⚠️ RDKit不可用，无法从SMILES生成分子式")
            return ""
        except Exception as e:
            print(f"⚠️ 从SMILES生成分子式失败: {e}")
            return ""