"""
重构后架构的使用示例
展示新的统一分类系统如何使用
"""

from pdb_uachecker.analysis import (
    UnifiedClassifier,
    ChemicalDatabase, 
    StandardAminoAcids,
    MolecularAnalyzer,
    BackboneAnalyzer,
    StereochemistryAnalyzer
)
from pdb_uachecker.core.models import AminoAcidInfo


def demo_unified_classifier():
    """演示统一分类器"""
    print("🎯 统一分类器演示")
    print("=" * 50)
    
    classifier = UnifiedClassifier()
    
    # 测试数据
    test_amino_acids = [
        AminoAcidInfo(
            id="ALA", 
            name="Alanine", 
            molecular_formula="C3H7NO2", 
            molecular_weight=89.09,
            smiles="N[C@@H](C)C(=O)O", 
            atom_composition={'C': 3, 'H': 7, 'N': 1, 'O': 2}
        ),
        AminoAcidInfo(
            id="D-ALA",
            name="D-Alanine",
            molecular_formula="C3H7NO2",
            molecular_weight=89.09,
            smiles="N[C@H](C)C(=O)O",
            atom_composition={'C': 3, 'H': 7, 'N': 1, 'O': 2}
        ),
        AminoAcidInfo(
            id="BETA-ALA",
            name="β-Alanine", 
            molecular_formula="C3H7NO2",
            molecular_weight=89.09,
            smiles="NCCC(=O)O",
            atom_composition={'C': 3, 'H': 7, 'N': 1, 'O': 2}
        )
    ]
    
    # 批量分类
    for amino_acid in test_amino_acids:
        print(f"\n分析氨基酸: {amino_acid.name} ({amino_acid.id})")
        result = classifier.classify(amino_acid)
        
        print(f"  类别: {result.categories}")
        print(f"  置信度: {result.confidence:.2f}")
        print(f"  方法: {result.classification_method}")
        print(f"  证据: {result.evidence}")
        
        if result.is_high_confidence:
            print("  ✅ 高置信度分类")
        elif result.requires_review:
            print("  ⚠️ 需要人工审查")
    
    # 显示统计信息
    print(f"\n📊 分类统计:")
    stats = classifier.get_statistics()
    print(f"  总计: {stats['total_classified']}")
    print(f"  确定性匹配: {stats['definitive_matches']}")
    print(f"  高置信度: {stats['high_confidence_matches']}")
    print(f"  准确率: {stats.get('accuracy_rate', 0):.2%}")


def demo_individual_analyzers():
    """演示各个分析器"""
    print("\n🧪 分析器组件演示")
    print("=" * 50)
    
    # 分子分析器
    molecular_analyzer = MolecularAnalyzer()
    smiles = "N[C@@H](Cc1ccccc1)C(=O)O"  # 苯丙氨酸
    
    print(f"\n分子分析 - SMILES: {smiles}")
    analysis = molecular_analyzer.get_analysis_summary(smiles)
    if analysis['valid']:
        features = analysis['basic_features']
        print(f"  芳香性: {features['is_aromatic']}")
        print(f"  环数: {features['ring_count']}")
        print(f"  手性中心: {features['chiral_centers_count']}")
    
    # 骨架分析器
    backbone_analyzer = BackboneAnalyzer()
    
    test_cases = [
        ("N[C@@H](C)C(=O)O", "ALA", "α-丙氨酸"),
        ("NCCC(=O)O", "BETA-ALA", "β-丙氨酸"),
        ("NCCCC(=O)O", "GABA", "γ-氨基丁酸")
    ]
    
    print(f"\n骨架分析:")
    for smiles, code, name in test_cases:
        result = backbone_analyzer.analyze(smiles, code)
        print(f"  {name}: {result['backbone_type']} (置信度: {result['confidence']:.2f})")
    
    # 立体化学分析器
    stereo_analyzer = StereochemistryAnalyzer()
    
    print(f"\n立体化学分析:")
    d_ala = "N[C@H](C)C(=O)O"
    l_ala = "N[C@@H](C)C(=O)O"
    
    d_result = stereo_analyzer.analyze(d_ala)
    l_result = stereo_analyzer.analyze(l_ala)
    
    print(f"  D-丙氨酸: {d_result['stereochemistry']} (置信度: {d_result['confidence']:.2f})")
    print(f"  L-丙氨酸: {l_result['stereochemistry']} (置信度: {l_result['confidence']:.2f})")


def demo_knowledge_base():
    """演示知识库"""
    print("\n📚 知识库演示") 
    print("=" * 50)
    
    # 化学数据库
    chemical_db = ChemicalDatabase()
    print(f"化学数据库统计:")
    stats = chemical_db.get_statistics()
    print(f"  总条目: {stats['total_entries']}")
    print(f"  类别: {len(stats['categories'])} 种")
    print(f"  SMILES覆盖率: {stats['smiles_coverage']}")
    
    # 标准氨基酸注册表
    standard_aa = StandardAminoAcids()
    print(f"\n标准氨基酸:")
    print(f"  总数: {len(standard_aa.get_all_codes())}")
    
    # 按类别统计
    aromatic_aa = standard_aa.get_by_property("aromatic")
    basic_aa = standard_aa.get_by_property("basic")
    print(f"  芳香性: {len(aromatic_aa)} 种")
    print(f"  碱性: {len(basic_aa)} 种")


def performance_comparison():
    """性能对比"""
    print("\n⚡ 架构重构效果")
    print("=" * 50)
    
    print("✅ 重构优势:")
    print("  1. 单一入口 - UnifiedClassifier统一所有分类逻辑")
    print("  2. 模块清晰 - 知识库、分析器、分类器职责分离") 
    print("  3. 消除重复 - 不再有4个重复的分类器")
    print("  4. 易于维护 - 修改只需要在一个地方进行")
    print("  5. 高精度 - 基于化学知识而非regex模式匹配")
    print("  6. 可扩展 - 新的分类方法容易添加")
    print("  7. 向后兼容 - Legacy代码保留用于参考")
    
    print("\n📁 新的架构结构:")
    print("""
    analysis/
    ├── unified_classifier.py      # 🎯 统一分类入口
    ├── knowledge/                  # 📚 知识库
    │   ├── chemical_database.py
    │   └── amino_acid_registry.py
    ├── analyzers/                  # 🧪 专门分析器
    │   ├── molecular_analyzer.py
    │   ├── backbone_analyzer.py
    │   └── stereochemistry_analyzer.py
    └── legacy/                     # 📜 历史代码
        ├── classifier.py
        ├── expert_classifier.py
        └── intelligent_classifier.py
    """)


if __name__ == "__main__":
    print("🏗️ PDB-UAChecker 架构重构演示")
    print("=" * 60)
    
    try:
        demo_unified_classifier()
        demo_individual_analyzers()
        demo_knowledge_base()
        performance_comparison()
        
        print(f"\n🎉 演示完成！")
        print(f"新架构已就绪，建议使用 UnifiedClassifier 进行氨基酸分类。")
        
    except Exception as e:
        print(f"❌ 演示出错: {e}")
        print("请检查导入和依赖是否正确。")
