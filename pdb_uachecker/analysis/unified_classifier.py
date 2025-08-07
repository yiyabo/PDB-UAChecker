"""
统一氨基酸分类器
整合所有分类方法的单一入口，提供高精度分类
"""

from typing import Dict, List, Optional, Any
from dataclasses import dataclass
from enum import Enum

from ..core.models import AminoAcidInfo, ClassificationResult, ChemicalKnowledge
from .knowledge.chemical_database import ChemicalDatabase
from .knowledge.amino_acid_registry import StandardAminoAcids
from .molecular_structure_analyzer import MolecularStructureAnalyzer
from .analyzers.backbone_analyzer import BackboneAnalyzer
from .analyzers.stereochemistry_analyzer import StereochemistryAnalyzer


class ClassificationTier(Enum):
    """分类层级"""
    DEFINITIVE = "definitive"           # 100%确定
    HIGH_CONFIDENCE = "high_confidence"  # >90%置信度
    INTELLIGENT = "intelligent"         # >80%置信度  
    NEEDS_REVIEW = "needs_review"       # <80%置信度


@dataclass
class ClassificationStrategy:
    """分类策略"""
    tier: ClassificationTier
    method: str
    confidence_threshold: float


class UnifiedClassifier:
    """
    统一氨基酸分类器
    
    整合了所有分类方法：
    1. 标准氨基酸直接匹配
    2. 化学知识库查询
    3. 分子结构分析
    4. 智能推断
    """
    
    def __init__(self):
        """初始化统一分类器"""
        # 核心组件
        self.chemical_database = ChemicalDatabase()
        self.standard_registry = StandardAminoAcids()
        self.molecular_analyzer = MolecularStructureAnalyzer()
        self.backbone_analyzer = BackboneAnalyzer()
        self.stereochemistry_analyzer = StereochemistryAnalyzer()
        
        # 分类策略
        self.strategies = [
            ClassificationStrategy(ClassificationTier.DEFINITIVE, "standard_lookup", 1.0),
            ClassificationStrategy(ClassificationTier.DEFINITIVE, "database_lookup", 1.0),
            ClassificationStrategy(ClassificationTier.HIGH_CONFIDENCE, "molecular_analysis", 0.95),
            ClassificationStrategy(ClassificationTier.HIGH_CONFIDENCE, "structure_matching", 0.90),
            ClassificationStrategy(ClassificationTier.INTELLIGENT, "pattern_analysis", 0.80),
            ClassificationStrategy(ClassificationTier.NEEDS_REVIEW, "fallback_rules", 0.70)
        ]
        
        # 统计信息
        self.stats = {
            'total_classified': 0,
            'definitive_matches': 0,
            'high_confidence_matches': 0,
            'intelligent_matches': 0,
            'needs_review': 0
        }
    
    def classify(self, amino_acid: AminoAcidInfo) -> ClassificationResult:
        """
        统一分类入口
        
        Args:
            amino_acid: 氨基酸信息
            
        Returns:
            分类结果
        """
        self.stats['total_classified'] += 1
        
        # 按策略优先级依次尝试分类
        for strategy in self.strategies:
            result = self._try_classification_strategy(amino_acid, strategy)
            if result and result.confidence >= strategy.confidence_threshold:
                self._update_stats(strategy.tier)
                return result
        
        # 如果所有策略都失败，返回未分类结果
        return self._create_unclassified_result(amino_acid)
    
    def _try_classification_strategy(
        self, 
        amino_acid: AminoAcidInfo, 
        strategy: ClassificationStrategy
    ) -> Optional[ClassificationResult]:
        """尝试特定分类策略"""
        
        try:
            if strategy.method == "standard_lookup":
                return self._classify_standard(amino_acid)
            elif strategy.method == "database_lookup":
                return self._classify_from_database(amino_acid)
            elif strategy.method == "molecular_analysis":
                return self._classify_by_molecular_analysis(amino_acid)
            elif strategy.method == "structure_matching":
                return self._classify_by_structure_matching(amino_acid)
            elif strategy.method == "pattern_analysis":
                return self._classify_by_pattern_analysis(amino_acid)
            elif strategy.method == "fallback_rules":
                return self._classify_by_fallback_rules(amino_acid)
        except Exception as e:
            print(f"⚠️ 分类策略 {strategy.method} 失败: {e}")
            return None
        
        return None
    
    def _classify_standard(self, amino_acid: AminoAcidInfo) -> Optional[ClassificationResult]:
        """标准氨基酸分类"""
        if self.standard_registry.is_standard(amino_acid.id):
            standard_info = self.standard_registry.get_info(amino_acid.id)
            return ClassificationResult(
                amino_acid_id=amino_acid.id,
                amino_acid_name=amino_acid.name,
                categories=["standard", standard_info.category],
                confidence=1.0,
                classification_method="standard_lookup",
                evidence=[f"标准氨基酸: {amino_acid.id}"],
                details={"standard_info": standard_info}
            )
        return None
    
    def _classify_from_database(self, amino_acid: AminoAcidInfo) -> Optional[ClassificationResult]:
        """从化学数据库分类"""
        # 按SMILES查找
        if amino_acid.smiles:
            entry = self.chemical_database.lookup_by_smiles(amino_acid.smiles)
            if entry:
                return ClassificationResult(
                    amino_acid_id=amino_acid.id,
                    amino_acid_name=amino_acid.name,
                    categories=[entry.category, entry.subcategory],
                    confidence=1.0,
                    classification_method="database_lookup",
                    evidence=[f"SMILES匹配: {amino_acid.smiles}"],
                    details={"database_entry": entry}
                )
        
        # 按名称查找
        entry = self.chemical_database.lookup_by_name(amino_acid.name)
        if entry:
            return ClassificationResult(
                amino_acid_id=amino_acid.id,
                amino_acid_name=amino_acid.name,
                categories=[entry.category, entry.subcategory],
                confidence=0.98,
                classification_method="database_lookup",
                evidence=[f"名称匹配: {amino_acid.name}"],
                details={"database_entry": entry}
            )
        
        return None
    
    def _classify_by_molecular_analysis(self, amino_acid: AminoAcidInfo) -> Optional[ClassificationResult]:
        """基于分子结构分析分类"""
        if not amino_acid.smiles:
            return None
        
        analysis = self.molecular_analyzer.analyze(amino_acid.smiles)
        if not analysis.is_valid:
            return None
        
        categories = []
        evidence = []
        confidence = 0.95
        
        # 芳香性分析
        if analysis.aromatic_atoms:
            categories.append("aromatic")
            evidence.append(f"检测到芳香原子: {len(analysis.aromatic_atoms)}个")
        
        # 环系统分析
        if analysis.ring_systems:
            categories.append("cyclic")
            evidence.append(f"检测到环系统: {len(analysis.ring_systems)}个")
        
        # 骨架分析
        if analysis.backbone_analysis:
            backbone_type = analysis.backbone_analysis.get('backbone_type')
            if backbone_type != 'unknown':
                categories.append(f"{backbone_type}_amino_acid")
                evidence.append(f"骨架类型: {backbone_type}")
        
        # 立体化学分析
        stereochemistry = self.stereochemistry_analyzer.analyze(amino_acid.smiles)
        if stereochemistry.get('stereochemistry') not in ['unknown', 'achiral_or_unknown']:
            stereo_type = stereochemistry['stereochemistry']
            if 'D_' in stereo_type:
                categories.append("d_amino_acid")
                evidence.append("D型立体化学")
            elif 'L_' in stereo_type:
                categories.append("l_amino_acid") 
                evidence.append("L型立体化学")
        
        if categories:
            return ClassificationResult(
                amino_acid_id=amino_acid.id,
                amino_acid_name=amino_acid.name,
                categories=categories,
                confidence=confidence,
                classification_method="molecular_analysis",
                evidence=evidence,
                details={"molecular_analysis": analysis}
            )
        
        return None
    
    def _classify_by_structure_matching(self, amino_acid: AminoAcidInfo) -> Optional[ClassificationResult]:
        """基于结构匹配分类"""
        if not amino_acid.smiles:
            return None
        
        # 查找相似结构
        similar_entries = self.chemical_database.find_similar_smiles(amino_acid.smiles, threshold=0.85)
        
        if similar_entries:
            best_match = similar_entries[0]  # 最高相似度
            entry_code, similarity = best_match
            entry = self.chemical_database.lookup_by_code(entry_code)
            
            if entry and similarity >= 0.85:
                return ClassificationResult(
                    amino_acid_id=amino_acid.id,
                    amino_acid_name=amino_acid.name,
                    categories=[entry.category, "structurally_similar"],
                    confidence=similarity * 0.9,  # 稍微降低置信度
                    classification_method="structure_matching",
                    evidence=[f"结构相似于 {entry_code} (相似度: {similarity:.2f})"],
                    details={"similar_entry": entry, "similarity": similarity}
                )
        
        return None
    
    def _classify_by_pattern_analysis(self, amino_acid: AminoAcidInfo) -> Optional[ClassificationResult]:
        """基于模式分析分类（智能规则）"""
        categories = []
        evidence = []
        confidence = 0.80
        
        # 骨架分析
        backbone_result = self.backbone_analyzer.analyze(amino_acid.smiles, amino_acid.id)
        if backbone_result:
            categories.extend(backbone_result.get('categories', []))
            evidence.extend(backbone_result.get('evidence', []))
        
        # 基于名称的智能推断
        name_lower = amino_acid.name.lower()
        if 'd-' in name_lower or name_lower.startswith('d'):
            categories.append("d_amino_acid")
            evidence.append("名称表明D型")
        
        if 'beta' in name_lower or 'β' in name_lower:
            categories.append("beta_amino_acid")
            evidence.append("名称表明β型")
        
        if 'gamma' in name_lower or 'γ' in name_lower:
            categories.append("gamma_amino_acid")
            evidence.append("名称表明γ型")
        
        if 'methyl' in name_lower and 'n-' in name_lower:
            categories.append("n_methyl")
            evidence.append("名称表明N-甲基化")
        
        if categories:
            return ClassificationResult(
                amino_acid_id=amino_acid.id,
                amino_acid_name=amino_acid.name,
                categories=categories,
                confidence=confidence,
                classification_method="pattern_analysis",
                evidence=evidence,
                details={"pattern_based": True}
            )
        
        return None
    
    def _classify_by_fallback_rules(self, amino_acid: AminoAcidInfo) -> Optional[ClassificationResult]:
        """后备规则分类"""
        # 基本的非标准分类
        return ClassificationResult(
            amino_acid_id=amino_acid.id,
            amino_acid_name=amino_acid.name,
            categories=["non_standard", "unspecified"],
            confidence=0.70,
            classification_method="fallback_rules",
            evidence=["基于后备规则的基本分类"],
            details={"fallback": True}
        )
    
    def _create_unclassified_result(self, amino_acid: AminoAcidInfo) -> ClassificationResult:
        """创建未分类结果"""
        self.stats['needs_review'] += 1
        
        return ClassificationResult(
            amino_acid_id=amino_acid.id,
            amino_acid_name=amino_acid.name,
            categories=["unclassified"],
            confidence=0.0,
            classification_method="none",
            evidence=["所有分类方法均失败"],
            details={"requires_manual_review": True}
        )
    
    def _update_stats(self, tier: ClassificationTier):
        """更新统计信息"""
        if tier == ClassificationTier.DEFINITIVE:
            self.stats['definitive_matches'] += 1
        elif tier == ClassificationTier.HIGH_CONFIDENCE:
            self.stats['high_confidence_matches'] += 1
        elif tier == ClassificationTier.INTELLIGENT:
            self.stats['intelligent_matches'] += 1
        else:
            self.stats['needs_review'] += 1
    
    def get_statistics(self) -> Dict[str, Any]:
        """获取分类统计信息"""
        total = self.stats['total_classified']
        if total == 0:
            return self.stats
        
        return {
            **self.stats,
            'accuracy_rate': (self.stats['definitive_matches'] + self.stats['high_confidence_matches']) / total,
            'coverage_rate': (total - self.stats['needs_review']) / total,
            'precision_distribution': {
                'definitive': self.stats['definitive_matches'] / total,
                'high_confidence': self.stats['high_confidence_matches'] / total,
                'intelligent': self.stats['intelligent_matches'] / total,
                'needs_review': self.stats['needs_review'] / total
            }
        }
    
    def batch_classify(self, amino_acids: List[AminoAcidInfo]) -> List[ClassificationResult]:
        """批量分类"""
        results = []
        for amino_acid in amino_acids:
            result = self.classify(amino_acid)
            results.append(result)
        return results
