"""
标准氨基酸注册表
集中管理20种标准氨基酸的定义和属性
"""

from typing import Dict, Optional, Set
from dataclasses import dataclass


@dataclass
class StandardAminoAcidInfo:
    """标准氨基酸信息"""
    code: str
    name: str
    full_name: str
    category: str
    properties: Dict[str, any]
    molecular_weight: float
    canonical_smiles: str


class StandardAminoAcids:
    """标准氨基酸注册表"""
    
    def __init__(self):
        """初始化标准氨基酸注册表"""
        self._registry = self._build_registry()
        self._codes = set(self._registry.keys())
    
    def _build_registry(self) -> Dict[str, StandardAminoAcidInfo]:
        """构建标准氨基酸注册表"""
        
        amino_acids = [
            StandardAminoAcidInfo(
                code="ALA", name="Alanine", full_name="L-Alanine",
                category="aliphatic", molecular_weight=89.09,
                canonical_smiles="N[C@@H](C)C(=O)O",
                properties={"hydrophobic": True, "polar": False, "charged": False, "small": True}
            ),
            StandardAminoAcidInfo(
                code="ARG", name="Arginine", full_name="L-Arginine", 
                category="basic", molecular_weight=174.20,
                canonical_smiles="N[C@@H](CCCNC(=N)N)C(=O)O",
                properties={"hydrophobic": False, "polar": True, "charged": True, "basic": True, "large": True}
            ),
            StandardAminoAcidInfo(
                code="ASN", name="Asparagine", full_name="L-Asparagine",
                category="polar", molecular_weight=132.12,
                canonical_smiles="N[C@@H](CC(=O)N)C(=O)O", 
                properties={"hydrophobic": False, "polar": True, "charged": False}
            ),
            StandardAminoAcidInfo(
                code="ASP", name="Aspartic acid", full_name="L-Aspartic acid",
                category="acidic", molecular_weight=133.10,
                canonical_smiles="N[C@@H](CC(=O)O)C(=O)O",
                properties={"hydrophobic": False, "polar": True, "charged": True, "acidic": True}
            ),
            StandardAminoAcidInfo(
                code="CYS", name="Cysteine", full_name="L-Cysteine",
                category="sulfur", molecular_weight=121.16,
                canonical_smiles="N[C@@H](CS)C(=O)O",
                properties={"hydrophobic": False, "polar": True, "charged": False, "sulfur_containing": True}
            ),
            StandardAminoAcidInfo(
                code="GLN", name="Glutamine", full_name="L-Glutamine",
                category="polar", molecular_weight=146.14,
                canonical_smiles="N[C@@H](CCC(=O)N)C(=O)O",
                properties={"hydrophobic": False, "polar": True, "charged": False}
            ),
            StandardAminoAcidInfo(
                code="GLU", name="Glutamic acid", full_name="L-Glutamic acid",
                category="acidic", molecular_weight=147.13,
                canonical_smiles="N[C@@H](CCC(=O)O)C(=O)O",
                properties={"hydrophobic": False, "polar": True, "charged": True, "acidic": True}
            ),
            StandardAminoAcidInfo(
                code="GLY", name="Glycine", full_name="Glycine",
                category="special", molecular_weight=75.07,
                canonical_smiles="NCC(=O)O",
                properties={"hydrophobic": False, "polar": False, "charged": False, "flexible": True, "small": True}
            ),
            StandardAminoAcidInfo(
                code="HIS", name="Histidine", full_name="L-Histidine",
                category="basic", molecular_weight=155.15,
                canonical_smiles="N[C@@H](Cc1c[nH]cn1)C(=O)O",
                properties={"hydrophobic": False, "polar": True, "charged": True, "aromatic": True, "basic": True}
            ),
            StandardAminoAcidInfo(
                code="ILE", name="Isoleucine", full_name="L-Isoleucine",
                category="aliphatic", molecular_weight=131.17,
                canonical_smiles="N[C@@H]([C@H](C)CC)C(=O)O",
                properties={"hydrophobic": True, "polar": False, "charged": False, "branched": True}
            ),
            StandardAminoAcidInfo(
                code="LEU", name="Leucine", full_name="L-Leucine",
                category="aliphatic", molecular_weight=131.17,
                canonical_smiles="N[C@@H](CC(C)C)C(=O)O",
                properties={"hydrophobic": True, "polar": False, "charged": False, "branched": True}
            ),
            StandardAminoAcidInfo(
                code="LYS", name="Lysine", full_name="L-Lysine", 
                category="basic", molecular_weight=146.19,
                canonical_smiles="N[C@@H](CCCCN)C(=O)O",
                properties={"hydrophobic": False, "polar": True, "charged": True, "basic": True}
            ),
            StandardAminoAcidInfo(
                code="MET", name="Methionine", full_name="L-Methionine",
                category="sulfur", molecular_weight=149.21,
                canonical_smiles="N[C@@H](CCSC)C(=O)O",
                properties={"hydrophobic": True, "polar": False, "charged": False, "sulfur_containing": True}
            ),
            StandardAminoAcidInfo(
                code="PHE", name="Phenylalanine", full_name="L-Phenylalanine",
                category="aromatic", molecular_weight=165.19,
                canonical_smiles="N[C@@H](Cc1ccccc1)C(=O)O",
                properties={"hydrophobic": True, "polar": False, "charged": False, "aromatic": True, "large": True}
            ),
            StandardAminoAcidInfo(
                code="PRO", name="Proline", full_name="L-Proline",
                category="cyclic", molecular_weight=115.13,
                canonical_smiles="N1[C@@H](CCC1)C(=O)O",
                properties={"hydrophobic": False, "polar": False, "charged": False, "cyclic": True, "rigid": True}
            ),
            StandardAminoAcidInfo(
                code="SER", name="Serine", full_name="L-Serine",
                category="polar", molecular_weight=105.09,
                canonical_smiles="N[C@@H](CO)C(=O)O",
                properties={"hydrophobic": False, "polar": True, "charged": False, "small": True}
            ),
            StandardAminoAcidInfo(
                code="THR", name="Threonine", full_name="L-Threonine",
                category="polar", molecular_weight=119.12,
                canonical_smiles="N[C@@H]([C@H](O)C)C(=O)O", 
                properties={"hydrophobic": False, "polar": True, "charged": False}
            ),
            StandardAminoAcidInfo(
                code="TRP", name="Tryptophan", full_name="L-Tryptophan",
                category="aromatic", molecular_weight=204.23,
                canonical_smiles="N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",
                properties={"hydrophobic": True, "polar": False, "charged": False, "aromatic": True, "large": True, "bulky": True}
            ),
            StandardAminoAcidInfo(
                code="TYR", name="Tyrosine", full_name="L-Tyrosine",
                category="aromatic", molecular_weight=181.19,
                canonical_smiles="N[C@@H](Cc1ccc(O)cc1)C(=O)O",
                properties={"hydrophobic": False, "polar": True, "charged": False, "aromatic": True, "large": True}
            ),
            StandardAminoAcidInfo(
                code="VAL", name="Valine", full_name="L-Valine",
                category="aliphatic", molecular_weight=117.15,
                canonical_smiles="N[C@@H](C(C)C)C(=O)O",
                properties={"hydrophobic": True, "polar": False, "charged": False, "branched": True}
            )
        ]
        
        return {aa.code: aa for aa in amino_acids}
    
    def is_standard(self, code: str) -> bool:
        """检查是否为标准氨基酸"""
        return code.upper() in self._codes
    
    def get_info(self, code: str) -> Optional[StandardAminoAcidInfo]:
        """获取标准氨基酸信息"""
        return self._registry.get(code.upper())
    
    def get_all_codes(self) -> Set[str]:
        """获取所有标准氨基酸代码"""
        return self._codes.copy()
    
    def get_by_category(self, category: str) -> Dict[str, StandardAminoAcidInfo]:
        """按类别获取氨基酸"""
        return {code: info for code, info in self._registry.items() if info.category == category}
    
    def get_by_property(self, property_name: str, property_value: any = True) -> Dict[str, StandardAminoAcidInfo]:
        """按属性获取氨基酸"""
        return {
            code: info for code, info in self._registry.items()
            if info.properties.get(property_name) == property_value
        }
    
    def search_by_name(self, name: str) -> Optional[StandardAminoAcidInfo]:
        """按名称搜索"""
        name_lower = name.lower()
        for info in self._registry.values():
            if (info.name.lower() == name_lower or 
                info.full_name.lower() == name_lower or
                info.code.lower() == name_lower):
                return info
        return None
    
    def get_statistics(self) -> Dict[str, any]:
        """获取统计信息"""
        categories = {}
        properties = {}
        
        for info in self._registry.values():
            # 统计类别
            categories[info.category] = categories.get(info.category, 0) + 1
            
            # 统计属性
            for prop, value in info.properties.items():
                if value:  # 只统计True值
                    properties[prop] = properties.get(prop, 0) + 1
        
        return {
            'total_standard_amino_acids': len(self._registry),
            'categories': categories,
            'properties': properties
        }
