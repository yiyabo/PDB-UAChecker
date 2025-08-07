"""
化学知识库
高精度氨基酸分类的化学知识数据库
"""

from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass
from .models import ChemicalKnowledge, AminoAcidInfo


@dataclass
class AminoAcidEntry:
    """氨基酸条目"""
    code: str
    name: str
    smiles: str
    category: str
    subcategory: str
    molecular_weight: float
    properties: Dict[str, any]


class ChemicalKnowledgeBase:
    """化学知识库"""
    
    def __init__(self):
        """初始化知识库"""
        self.amino_acids = self._load_amino_acid_data()
        self.smiles_index = self._build_smiles_index()
        self.name_index = self._build_name_index()
        self.code_index = self._build_code_index()
    
    def _load_amino_acid_data(self) -> Dict[str, AminoAcidEntry]:
        """加载氨基酸数据"""
        # 20种标准氨基酸
        standard_amino_acids = [
            AminoAcidEntry("ALA", "Alanine", "N[C@@H](C)C(=O)O", "standard", "aliphatic", 89.09, 
                          {"hydrophobic": True, "polar": False, "charged": False}),
            AminoAcidEntry("ARG", "Arginine", "N[C@@H](CCCNC(=N)N)C(=O)O", "standard", "basic", 174.20,
                          {"hydrophobic": False, "polar": True, "charged": True, "charge": "+1"}),
            AminoAcidEntry("ASN", "Asparagine", "N[C@@H](CC(=O)N)C(=O)O", "standard", "polar", 132.12,
                          {"hydrophobic": False, "polar": True, "charged": False}),
            AminoAcidEntry("ASP", "Aspartic acid", "N[C@@H](CC(=O)O)C(=O)O", "standard", "acidic", 133.10,
                          {"hydrophobic": False, "polar": True, "charged": True, "charge": "-1"}),
            AminoAcidEntry("CYS", "Cysteine", "N[C@@H](CS)C(=O)O", "standard", "sulfur", 121.16,
                          {"hydrophobic": False, "polar": True, "charged": False}),
            AminoAcidEntry("GLN", "Glutamine", "N[C@@H](CCC(=O)N)C(=O)O", "standard", "polar", 146.14,
                          {"hydrophobic": False, "polar": True, "charged": False}),
            AminoAcidEntry("GLU", "Glutamic acid", "N[C@@H](CCC(=O)O)C(=O)O", "standard", "acidic", 147.13,
                          {"hydrophobic": False, "polar": True, "charged": True, "charge": "-1"}),
            AminoAcidEntry("GLY", "Glycine", "NCC(=O)O", "standard", "special", 75.07,
                          {"hydrophobic": False, "polar": False, "charged": False, "flexible": True}),
            AminoAcidEntry("HIS", "Histidine", "N[C@@H](Cc1c[nH]cn1)C(=O)O", "standard", "basic", 155.15,
                          {"hydrophobic": False, "polar": True, "charged": True, "aromatic": True}),
            AminoAcidEntry("ILE", "Isoleucine", "N[C@@H]([C@H](C)CC)C(=O)O", "standard", "aliphatic", 131.17,
                          {"hydrophobic": True, "polar": False, "charged": False, "branched": True}),
            AminoAcidEntry("LEU", "Leucine", "N[C@@H](CC(C)C)C(=O)O", "standard", "aliphatic", 131.17,
                          {"hydrophobic": True, "polar": False, "charged": False, "branched": True}),
            AminoAcidEntry("LYS", "Lysine", "N[C@@H](CCCCN)C(=O)O", "standard", "basic", 146.19,
                          {"hydrophobic": False, "polar": True, "charged": True, "charge": "+1"}),
            AminoAcidEntry("MET", "Methionine", "N[C@@H](CCSC)C(=O)O", "standard", "sulfur", 149.21,
                          {"hydrophobic": True, "polar": False, "charged": False}),
            AminoAcidEntry("PHE", "Phenylalanine", "N[C@@H](Cc1ccccc1)C(=O)O", "standard", "aromatic", 165.19,
                          {"hydrophobic": True, "polar": False, "charged": False, "aromatic": True}),
            AminoAcidEntry("PRO", "Proline", "N1[C@@H](CCC1)C(=O)O", "standard", "cyclic", 115.13,
                          {"hydrophobic": False, "polar": False, "charged": False, "cyclic": True}),
            AminoAcidEntry("SER", "Serine", "N[C@@H](CO)C(=O)O", "standard", "polar", 105.09,
                          {"hydrophobic": False, "polar": True, "charged": False}),
            AminoAcidEntry("THR", "Threonine", "N[C@@H]([C@H](O)C)C(=O)O", "standard", "polar", 119.12,
                          {"hydrophobic": False, "polar": True, "charged": False}),
            AminoAcidEntry("TRP", "Tryptophan", "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", "standard", "aromatic", 204.23,
                          {"hydrophobic": True, "polar": False, "charged": False, "aromatic": True, "large": True}),
            AminoAcidEntry("TYR", "Tyrosine", "N[C@@H](Cc1ccc(O)cc1)C(=O)O", "standard", "aromatic", 181.19,
                          {"hydrophobic": False, "polar": True, "charged": False, "aromatic": True}),
            AminoAcidEntry("VAL", "Valine", "N[C@@H](C(C)C)C(=O)O", "standard", "aliphatic", 117.15,
                          {"hydrophobic": True, "polar": False, "charged": False, "branched": True}),
        ]
        
        # 常见非天然氨基酸
        non_standard_amino_acids = [
            # D-型氨基酸
            AminoAcidEntry("D-ALA", "D-Alanine", "N[C@H](C)C(=O)O", "d_amino_acid", "d_aliphatic", 89.09,
                          {"d_form": True, "hydrophobic": True, "polar": False, "charged": False}),
            AminoAcidEntry("D-PHE", "D-Phenylalanine", "N[C@H](Cc1ccccc1)C(=O)O", "d_amino_acid", "d_aromatic", 165.19,
                          {"d_form": True, "hydrophobic": True, "polar": False, "charged": False, "aromatic": True}),
            AminoAcidEntry("D-LEU", "D-Leucine", "N[C@H](CC(C)C)C(=O)O", "d_amino_acid", "d_aliphatic", 131.17,
                          {"d_form": True, "hydrophobic": True, "polar": False, "charged": False, "branched": True}),
            
            # β-氨基酸
            AminoAcidEntry("BETA-ALA", "β-Alanine", "NCCC(=O)O", "beta_amino_acid", "beta_aliphatic", 89.09,
                          {"beta_form": True, "hydrophobic": False, "polar": True, "charged": False}),
            AminoAcidEntry("BETA-PHE", "β-Phenylalanine", "NCCc1ccccc1C(=O)O", "beta_amino_acid", "beta_aromatic", 165.19,
                          {"beta_form": True, "aromatic": True}),
            
            # γ-氨基酸
            AminoAcidEntry("GAMMA-ABA", "γ-Aminobutyric acid", "NCCCC(=O)O", "gamma_amino_acid", "gamma_aliphatic", 103.12,
                          {"gamma_form": True, "neurotransmitter": True}),
            
            # N-甲基氨基酸
            AminoAcidEntry("N-ME-ALA", "N-Methylalanine", "CN[C@@H](C)C(=O)O", "n_methyl", "n_methyl_aliphatic", 103.12,
                          {"n_methyl": True, "hydrophobic": True}),
            AminoAcidEntry("N-ME-PHE", "N-Methylphenylalanine", "CN[C@@H](Cc1ccccc1)C(=O)O", "n_methyl", "n_methyl_aromatic", 179.22,
                          {"n_methyl": True, "aromatic": True}),
            
            # 其他常见修饰氨基酸
            AminoAcidEntry("HYP", "4-Hydroxyproline", "N1[C@@H](CC(O)C1)C(=O)O", "modified", "hydroxylated", 131.13,
                          {"hydroxylated": True, "cyclic": True}),
            AminoAcidEntry("ORN", "Ornithine", "N[C@@H](CCCN)C(=O)O", "modified", "basic", 132.16,
                          {"basic": True, "charged": True}),
        ]
        
        # 合并所有数据
        all_amino_acids = standard_amino_acids + non_standard_amino_acids
        
        return {entry.code: entry for entry in all_amino_acids}
    
    def _build_smiles_index(self) -> Dict[str, str]:
        """构建SMILES索引"""
        return {entry.smiles: code for code, entry in self.amino_acids.items()}
    
    def _build_name_index(self) -> Dict[str, str]:
        """构建名称索引"""
        index = {}
        for code, entry in self.amino_acids.items():
            index[entry.name.lower()] = code
            index[code.lower()] = code
        return index
    
    def _build_code_index(self) -> Dict[str, str]:
        """构建代码索引"""
        return {code.upper(): code for code in self.amino_acids.keys()}
    
    def lookup_by_smiles(self, smiles: str) -> Optional[AminoAcidEntry]:
        """根据SMILES查找氨基酸"""
        code = self.smiles_index.get(smiles)
        if code:
            return self.amino_acids[code]
        return None
    
    def lookup_by_name(self, name: str) -> Optional[AminoAcidEntry]:
        """根据名称查找氨基酸"""
        code = self.name_index.get(name.lower())
        if code:
            return self.amino_acids[code]
        return None
    
    def lookup_by_code(self, code: str) -> Optional[AminoAcidEntry]:
        """根据代码查找氨基酸"""
        normalized_code = self.code_index.get(code.upper())
        if normalized_code:
            return self.amino_acids[normalized_code]
        return None
    
    def find_similar_smiles(self, target_smiles: str, threshold: float = 0.8) -> List[Tuple[str, float]]:
        """查找相似的SMILES结构"""
        similar = []
        
        # 简单的字符串相似度（生产环境应使用分子指纹）
        for smiles, code in self.smiles_index.items():
            similarity = self._calculate_string_similarity(target_smiles, smiles)
            if similarity >= threshold:
                similar.append((code, similarity))
        
        return sorted(similar, key=lambda x: x[1], reverse=True)
    
    def _calculate_string_similarity(self, s1: str, s2: str) -> float:
        """计算字符串相似度（简化版）"""
        if not s1 or not s2:
            return 0.0
        
        # 使用Jaccard相似度
        set1 = set(s1.lower())
        set2 = set(s2.lower())
        
        intersection = len(set1 & set2)
        union = len(set1 | set2)
        
        return intersection / union if union > 0 else 0.0
    
    def get_category_members(self, category: str) -> List[str]:
        """获取特定类别的所有氨基酸"""
        return [code for code, entry in self.amino_acids.items() 
                if entry.category == category]
    
    def get_subcategory_members(self, subcategory: str) -> List[str]:
        """获取特定子类别的所有氨基酸"""
        return [code for code, entry in self.amino_acids.items() 
                if entry.subcategory == subcategory]
    
    def get_properties(self, code: str) -> Optional[Dict[str, any]]:
        """获取氨基酸属性"""
        entry = self.lookup_by_code(code)
        return entry.properties if entry else None
    
    def has_property(self, code: str, property_name: str) -> bool:
        """检查氨基酸是否具有特定属性"""
        properties = self.get_properties(code)
        if properties:
            return properties.get(property_name, False)
        return False
    
    def get_all_codes(self) -> Set[str]:
        """获取所有氨基酸代码"""
        return set(self.amino_acids.keys())
    
    def get_statistics(self) -> Dict[str, int]:
        """获取知识库统计信息"""
        stats = {
            'total_amino_acids': len(self.amino_acids),
            'standard': 0,
            'd_amino_acid': 0,
            'beta_amino_acid': 0,
            'gamma_amino_acid': 0,
            'n_methyl': 0,
            'modified': 0
        }
        
        for entry in self.amino_acids.values():
            if entry.category in stats:
                stats[entry.category] += 1
        
        return stats


# 全局知识库实例
knowledge_base = ChemicalKnowledgeBase()
