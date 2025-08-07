"""
化学知识库
存储和管理氨基酸的化学知识，提供高精度查询功能
"""

from typing import Dict, List, Optional, Set, Any, Tuple
import json
from pathlib import Path

from ..core.models import ChemicalKnowledge, AminoAcidInfo


class ChemicalKnowledgeBase:
    """化学知识库"""
    
    def __init__(self):
        """初始化化学知识库"""
        self.knowledge_db: Dict[str, ChemicalKnowledge] = {}
        self.smiles_index: Dict[str, str] = {}  # SMILES -> amino_acid_id
        self.category_index: Dict[str, Set[str]] = {}  # category -> set of amino_acid_ids
        
        # 初始化已知的高质量数据
        self._initialize_curated_data()
        print(f"📚 化学知识库初始化完成，包含 {len(self.knowledge_db)} 条记录")
    
    def _initialize_curated_data(self):
        """初始化精心整理的化学数据"""
        
        # 标准氨基酸 - 100%确定
        standard_amino_acids = {
            'ALA': ChemicalKnowledge(
                amino_acid_id='ALA',
                canonical_smiles='N[C@@H](C)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['aliphatic'],
                alternative_names=['L-丙氨酸', 'Alanine']
            ),
            'ARG': ChemicalKnowledge(
                amino_acid_id='ARG',
                canonical_smiles='N[C@@H](CCCNC(=N)N)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl', 'guanidinium'],
                structural_features=['basic', 'polar'],
                alternative_names=['L-精氨酸', 'Arginine']
            ),
            'ASN': ChemicalKnowledge(
                amino_acid_id='ASN',
                canonical_smiles='N[C@@H](CC(=O)N)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl', 'amide'],
                structural_features=['polar'],
                alternative_names=['L-天冬酰胺', 'Asparagine']
            ),
            'ASP': ChemicalKnowledge(
                amino_acid_id='ASP',
                canonical_smiles='N[C@@H](CC(=O)O)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['acidic', 'polar'],
                alternative_names=['L-天冬氨酸', 'Aspartic acid']
            ),
            'CYS': ChemicalKnowledge(
                amino_acid_id='CYS',
                canonical_smiles='N[C@@H](CS)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl', 'thiol'],
                structural_features=['sulfur_containing'],
                alternative_names=['L-半胱氨酸', 'Cysteine']
            ),
            'GLN': ChemicalKnowledge(
                amino_acid_id='GLN',
                canonical_smiles='N[C@@H](CCC(=O)N)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl', 'amide'],
                structural_features=['polar'],
                alternative_names=['L-谷氨酰胺', 'Glutamine']
            ),
            'GLU': ChemicalKnowledge(
                amino_acid_id='GLU',
                canonical_smiles='N[C@@H](CCC(=O)O)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['acidic', 'polar'],
                alternative_names=['L-谷氨酸', 'Glutamic acid']
            ),
            'GLY': ChemicalKnowledge(
                amino_acid_id='GLY',
                canonical_smiles='NCC(=O)O',
                stereochemistry=None,  # 甘氨酸没有手性中心
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['smallest'],
                alternative_names=['甘氨酸', 'Glycine']
            ),
            'HIS': ChemicalKnowledge(
                amino_acid_id='HIS',
                canonical_smiles='N[C@@H](Cc1c[nH]cn1)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl', 'imidazole'],
                structural_features=['aromatic', 'basic', 'polar'],
                alternative_names=['L-组氨酸', 'Histidine']
            ),
            'ILE': ChemicalKnowledge(
                amino_acid_id='ILE',
                canonical_smiles='N[C@@H]([C@H](C)CC)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['aliphatic', 'branched'],
                alternative_names=['L-异亮氨酸', 'Isoleucine']
            ),
            'LEU': ChemicalKnowledge(
                amino_acid_id='LEU',
                canonical_smiles='N[C@@H](CC(C)C)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['aliphatic', 'branched'],
                alternative_names=['L-亮氨酸', 'Leucine']
            ),
            'LYS': ChemicalKnowledge(
                amino_acid_id='LYS',
                canonical_smiles='N[C@@H](CCCCN)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['basic', 'polar'],
                alternative_names=['L-赖氨酸', 'Lysine']
            ),
            'MET': ChemicalKnowledge(
                amino_acid_id='MET',
                canonical_smiles='N[C@@H](CCSC)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl', 'thioether'],
                structural_features=['sulfur_containing'],
                alternative_names=['L-蛋氨酸', 'Methionine']
            ),
            'PHE': ChemicalKnowledge(
                amino_acid_id='PHE',
                canonical_smiles='N[C@@H](Cc1ccccc1)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['aromatic', 'hydrophobic'],
                alternative_names=['L-苯丙氨酸', 'Phenylalanine']
            ),
            'PRO': ChemicalKnowledge(
                amino_acid_id='PRO',
                canonical_smiles='N1[C@@H](CCC1)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['cyclic', 'imino_acid'],
                alternative_names=['L-脯氨酸', 'Proline']
            ),
            'SER': ChemicalKnowledge(
                amino_acid_id='SER',
                canonical_smiles='N[C@@H](CO)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl', 'hydroxyl'],
                structural_features=['polar'],
                alternative_names=['L-丝氨酸', 'Serine']
            ),
            'THR': ChemicalKnowledge(
                amino_acid_id='THR',
                canonical_smiles='N[C@@H]([C@H](C)O)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl', 'hydroxyl'],
                structural_features=['polar'],
                alternative_names=['L-苏氨酸', 'Threonine']
            ),
            'TRP': ChemicalKnowledge(
                amino_acid_id='TRP',
                canonical_smiles='N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl', 'indole'],
                structural_features=['aromatic', 'hydrophobic'],
                alternative_names=['L-色氨酸', 'Tryptophan']
            ),
            'TYR': ChemicalKnowledge(
                amino_acid_id='TYR',
                canonical_smiles='N[C@@H](Cc1ccc(O)cc1)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl', 'phenol'],
                structural_features=['aromatic', 'polar'],
                alternative_names=['L-酪氨酸', 'Tyrosine']
            ),
            'VAL': ChemicalKnowledge(
                amino_acid_id='VAL',
                canonical_smiles='N[C@@H](C(C)C)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['aliphatic', 'branched'],
                alternative_names=['L-缬氨酸', 'Valine']
            )
        }
        
        # 已知D氨基酸
        d_amino_acids = {
            'DAL': ChemicalKnowledge(
                amino_acid_id='DAL',
                canonical_smiles='N[C@H](C)C(=O)O',
                stereochemistry='D',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['aliphatic'],
                alternative_names=['D-丙氨酸', 'D-Alanine']
            ),
            'DLE': ChemicalKnowledge(
                amino_acid_id='DLE',
                canonical_smiles='N[C@H](CC(C)C)C(=O)O',
                stereochemistry='D',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['aliphatic', 'branched'],
                alternative_names=['D-亮氨酸', 'D-Leucine']
            ),
            'DPR': ChemicalKnowledge(
                amino_acid_id='DPR',
                canonical_smiles='N1[C@H](CCC1)C(=O)O',
                stereochemistry='D',
                backbone_type='alpha',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['cyclic', 'imino_acid'],
                alternative_names=['D-脯氨酸', 'D-Proline']
            ),
        }
        
        # 已知Beta氨基酸
        beta_amino_acids = {
            'DPP': ChemicalKnowledge(
                amino_acid_id='DPP',
                canonical_smiles='NCCC(=O)O',
                stereochemistry=None,
                backbone_type='beta',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['linear'],
                alternative_names=['β-丙氨酸', 'beta-Alanine', '3-氨基丙酸']
            ),
        }
        
        # 已知Gamma氨基酸
        gamma_amino_acids = {
            'DAB': ChemicalKnowledge(
                amino_acid_id='DAB',
                canonical_smiles='NCCCC(=O)O',
                stereochemistry=None,
                backbone_type='gamma',
                functional_groups=['amino', 'carboxyl'],
                structural_features=['linear'],
                alternative_names=['γ-氨基丁酸', 'GABA', '4-氨基丁酸']
            ),
        }
        
        # 已知N-甲基氨基酸
        n_methyl_amino_acids = {
            'MEN': ChemicalKnowledge(
                amino_acid_id='MEN',
                canonical_smiles='CN[C@@H](CC(=O)N)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['n_methyl_amino', 'carboxyl', 'amide'],
                structural_features=['n_methylated'],
                alternative_names=['N-甲基-L-天冬酰胺', 'N-Methyl-L-asparagine']
            ),
            'MEQ': ChemicalKnowledge(
                amino_acid_id='MEQ',
                canonical_smiles='CN[C@@H](CCC(=O)N)C(=O)O',
                stereochemistry='L',
                backbone_type='alpha',
                functional_groups=['n_methyl_amino', 'carboxyl', 'amide'],
                structural_features=['n_methylated'],
                alternative_names=['N-甲基-L-谷氨酰胺', 'N-Methyl-L-glutamine']
            ),
        }
        
        # 合并所有数据
        all_knowledge = {
            **standard_amino_acids,
            **d_amino_acids,
            **beta_amino_acids,
            **gamma_amino_acids,
            **n_methyl_amino_acids
        }
        
        # 添加到知识库并建立索引
        for knowledge in all_knowledge.values():
            self.add_knowledge(knowledge)
    
    def add_knowledge(self, knowledge: ChemicalKnowledge):
        """添加化学知识条目"""
        self.knowledge_db[knowledge.amino_acid_id] = knowledge
        
        # 建立SMILES索引
        if knowledge.canonical_smiles:
            self.smiles_index[knowledge.canonical_smiles] = knowledge.amino_acid_id
        
        # 建立分类索引
        categories = []
        if knowledge.stereochemistry:
            categories.append(f"{knowledge.stereochemistry.lower()}_amino_acid")
        if knowledge.backbone_type:
            categories.append(f"{knowledge.backbone_type}_amino_acid")
        
        categories.extend(knowledge.structural_features)
        
        for category in categories:
            if category not in self.category_index:
                self.category_index[category] = set()
            self.category_index[category].add(knowledge.amino_acid_id)
    
    def lookup_by_id(self, amino_acid_id: str) -> Optional[ChemicalKnowledge]:
        """根据氨基酸ID查询"""
        return self.knowledge_db.get(amino_acid_id)
    
    def lookup_by_smiles(self, smiles: str) -> Optional[ChemicalKnowledge]:
        """根据SMILES查询"""
        amino_acid_id = self.smiles_index.get(smiles)
        if amino_acid_id:
            return self.knowledge_db.get(amino_acid_id)
        return None
    
    def get_by_category(self, category: str) -> List[ChemicalKnowledge]:
        """获取特定分类的所有氨基酸"""
        amino_acid_ids = self.category_index.get(category, set())
        return [self.knowledge_db[aid] for aid in amino_acid_ids if aid in self.knowledge_db]
    
    def is_known_amino_acid(self, amino_acid_id: str) -> bool:
        """检查是否为已知氨基酸"""
        return amino_acid_id in self.knowledge_db
    
    def is_standard_amino_acid(self, amino_acid_id: str) -> bool:
        """检查是否为标准氨基酸"""
        knowledge = self.lookup_by_id(amino_acid_id)
        if knowledge:
            return (knowledge.stereochemistry == 'L' and 
                   knowledge.backbone_type == 'alpha' and
                   amino_acid_id in {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 
                                   'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 
                                   'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'})
        return False
    
    def get_statistics(self) -> Dict[str, Any]:
        """获取知识库统计信息"""
        stats = {
            'total_entries': len(self.knowledge_db),
            'standard_amino_acids': len([k for k in self.knowledge_db.values() if self.is_standard_amino_acid(k.amino_acid_id)]),
            'stereochemistry_distribution': {},
            'backbone_distribution': {},
            'category_distribution': {cat: len(aids) for cat, aids in self.category_index.items()}
        }
        
        # 统计立体化学分布
        for knowledge in self.knowledge_db.values():
            stereo = knowledge.stereochemistry or 'none'
            stats['stereochemistry_distribution'][stereo] = stats['stereochemistry_distribution'].get(stereo, 0) + 1
            
            backbone = knowledge.backbone_type or 'unknown'
            stats['backbone_distribution'][backbone] = stats['backbone_distribution'].get(backbone, 0) + 1
        
        return stats
    
    def export_to_file(self, filepath: str):
        """导出知识库到文件"""
        data = []
        for knowledge in self.knowledge_db.values():
            data.append({
                'amino_acid_id': knowledge.amino_acid_id,
                'canonical_smiles': knowledge.canonical_smiles,
                'stereochemistry': knowledge.stereochemistry,
                'backbone_type': knowledge.backbone_type,
                'functional_groups': knowledge.functional_groups,
                'structural_features': knowledge.structural_features,
                'alternative_names': knowledge.alternative_names,
                'confidence': knowledge.confidence,
                'source': knowledge.source
            })
        
        with open(filepath, 'w', encoding='utf-8') as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        
        print(f"📁 知识库已导出到: {filepath}")
    
    def load_from_file(self, filepath: str):
        """从文件加载知识库"""
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            for entry in data:
                knowledge = ChemicalKnowledge(**entry)
                self.add_knowledge(knowledge)
            
            print(f"📁 从文件加载了 {len(data)} 条知识库记录")
        except Exception as e:
            print(f"⚠️ 加载知识库文件失败: {e}")
    
    def search_similar(self, query_features: List[str]) -> List[Tuple[ChemicalKnowledge, float]]:
        """基于特征搜索相似氨基酸"""
        results = []
        
        for knowledge in self.knowledge_db.values():
            # 计算特征相似度
            all_features = knowledge.structural_features + knowledge.functional_groups
            if knowledge.stereochemistry:
                all_features.append(f"{knowledge.stereochemistry.lower()}_amino_acid")
            if knowledge.backbone_type:
                all_features.append(f"{knowledge.backbone_type}_amino_acid")
            
            # Jaccard相似度
            intersection = set(query_features) & set(all_features)
            union = set(query_features) | set(all_features)
            
            if union:
                similarity = len(intersection) / len(union)
                if similarity > 0:
                    results.append((knowledge, similarity))
        
        # 按相似度排序
        results.sort(key=lambda x: x[1], reverse=True)
        return results
