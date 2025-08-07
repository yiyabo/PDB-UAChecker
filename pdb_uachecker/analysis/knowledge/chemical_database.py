"""
重构的化学数据库
整合所有化学知识，提供统一的数据接口
"""

from typing import Dict, List, Optional, Set, Tuple
from dataclasses import dataclass

# 重用现有的化学知识库数据，但提供更清晰的接口
from ...core.chemical_knowledge_base import knowledge_base, AminoAcidEntry


class ChemicalDatabase:
    """
    重构的化学数据库
    
    提供统一的化学知识接口，整合：
    1. 标准氨基酸数据
    2. 非天然氨基酸数据  
    3. 结构相似性搜索
    4. 属性查询
    """
    
    def __init__(self):
        """初始化化学数据库"""
        # 使用现有的知识库数据
        self._knowledge_base = knowledge_base
        
        # 构建快速查询索引
        self._smiles_index = self._build_smiles_index()
        self._name_index = self._build_name_index()
        self._code_index = self._build_code_index()
        self._category_index = self._build_category_index()
        self._property_index = self._build_property_index()
    
    def _build_smiles_index(self) -> Dict[str, str]:
        """构建SMILES索引"""
        return {entry.smiles: code for code, entry in self._knowledge_base.amino_acids.items()}
    
    def _build_name_index(self) -> Dict[str, str]:
        """构建名称索引"""
        index = {}
        for code, entry in self._knowledge_base.amino_acids.items():
            # 名称的多种形式
            index[entry.name.lower()] = code
            index[code.lower()] = code
            
            # 去掉特殊字符的版本
            clean_name = entry.name.lower().replace('-', '').replace(' ', '')
            index[clean_name] = code
        
        return index
    
    def _build_code_index(self) -> Dict[str, str]:
        """构建代码索引"""
        return {code.upper(): code for code in self._knowledge_base.amino_acids.keys()}
    
    def _build_category_index(self) -> Dict[str, List[str]]:
        """构建类别索引"""
        index = {}
        for code, entry in self._knowledge_base.amino_acids.items():
            if entry.category not in index:
                index[entry.category] = []
            index[entry.category].append(code)
        return index
    
    def _build_property_index(self) -> Dict[str, List[str]]:
        """构建属性索引"""
        index = {}
        for code, entry in self._knowledge_base.amino_acids.items():
            for prop, value in entry.properties.items():
                if value:  # 只索引True值
                    if prop not in index:
                        index[prop] = []
                    index[prop].append(code)
        return index
    
    def lookup_by_smiles(self, smiles: str) -> Optional[AminoAcidEntry]:
        """根据SMILES查找氨基酸"""
        code = self._smiles_index.get(smiles)
        if code:
            return self._knowledge_base.amino_acids[code]
        return None
    
    def lookup_by_name(self, name: str) -> Optional[AminoAcidEntry]:
        """根据名称查找氨基酸"""
        code = self._name_index.get(name.lower())
        if code:
            return self._knowledge_base.amino_acids[code]
        return None
    
    def lookup_by_code(self, code: str) -> Optional[AminoAcidEntry]:
        """根据代码查找氨基酸"""
        normalized_code = self._code_index.get(code.upper())
        if normalized_code:
            return self._knowledge_base.amino_acids[normalized_code]
        return None
    
    def find_similar_smiles(self, target_smiles: str, threshold: float = 0.8) -> List[Tuple[str, float]]:
        """查找相似的SMILES结构"""
        return self._knowledge_base.find_similar_smiles(target_smiles, threshold)
    
    def get_by_category(self, category: str) -> List[AminoAcidEntry]:
        """按类别获取氨基酸"""
        codes = self._category_index.get(category, [])
        return [self._knowledge_base.amino_acids[code] for code in codes]
    
    def get_by_property(self, property_name: str) -> List[AminoAcidEntry]:
        """按属性获取氨基酸"""
        codes = self._property_index.get(property_name, [])
        return [self._knowledge_base.amino_acids[code] for code in codes]
    
    def search_multi_criteria(self, **criteria) -> List[AminoAcidEntry]:
        """多条件搜索"""
        results = []
        
        for code, entry in self._knowledge_base.amino_acids.items():
            match = True
            
            # 检查每个条件
            for key, value in criteria.items():
                if key == 'category':
                    if entry.category != value:
                        match = False
                        break
                elif key == 'subcategory':
                    if entry.subcategory != value:
                        match = False
                        break
                elif key in entry.properties:
                    if entry.properties.get(key) != value:
                        match = False
                        break
                elif key == 'molecular_weight_range':
                    min_mw, max_mw = value
                    if not (min_mw <= entry.molecular_weight <= max_mw):
                        match = False
                        break
            
            if match:
                results.append(entry)
        
        return results
    
    def get_all_categories(self) -> Set[str]:
        """获取所有类别"""
        return set(self._category_index.keys())
    
    def get_all_properties(self) -> Set[str]:
        """获取所有属性"""
        return set(self._property_index.keys())
    
    def get_statistics(self) -> Dict[str, any]:
        """获取数据库统计信息"""
        total_entries = len(self._knowledge_base.amino_acids)
        
        # 类别统计
        category_stats = {cat: len(codes) for cat, codes in self._category_index.items()}
        
        # 属性统计
        property_stats = {prop: len(codes) for prop, codes in self._property_index.items()}
        
        # 分子量统计
        molecular_weights = [entry.molecular_weight for entry in self._knowledge_base.amino_acids.values()]
        mw_stats = {
            'min': min(molecular_weights),
            'max': max(molecular_weights),
            'average': sum(molecular_weights) / len(molecular_weights)
        }
        
        return {
            'total_entries': total_entries,
            'categories': category_stats,
            'properties': property_stats,
            'molecular_weight_stats': mw_stats,
            'smiles_coverage': len([e for e in self._knowledge_base.amino_acids.values() if e.smiles]),
        }
    
    def validate_entry(self, entry: AminoAcidEntry) -> Dict[str, any]:
        """验证数据库条目的完整性"""
        issues = []
        warnings = []
        
        # 必需字段检查
        if not entry.code:
            issues.append("Missing amino acid code")
        if not entry.name:
            issues.append("Missing amino acid name") 
        if not entry.smiles:
            warnings.append("Missing SMILES structure")
        
        # 数据合理性检查
        if entry.molecular_weight <= 0:
            issues.append("Invalid molecular weight")
        
        if entry.molecular_weight > 1000:
            warnings.append("Unusually high molecular weight")
        
        # 类别一致性检查
        if entry.category not in self.get_all_categories():
            warnings.append(f"Unknown category: {entry.category}")
        
        return {
            'valid': len(issues) == 0,
            'issues': issues,
            'warnings': warnings
        }
    
    def export_summary(self) -> Dict[str, any]:
        """导出数据库摘要"""
        return {
            'database_info': {
                'name': 'PDB-UAChecker Chemical Database',
                'version': '1.0',
                'entries': len(self._knowledge_base.amino_acids)
            },
            'statistics': self.get_statistics(),
            'available_categories': list(self.get_all_categories()),
            'available_properties': list(self.get_all_properties())
        }
