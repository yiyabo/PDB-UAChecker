#!/usr/bin/env python3
"""
第三阶段高级功能演示脚本
展示3D结构分析、机器学习分类、立体化学分析和高级分析工具
"""

import os
import time
import json
import numpy as np
from typing import Dict, List, Any

# 导入第三阶段组件
from advanced_features_engine import (
    Structure3DAnalyzer,
    MLClassifier,
    StereochemistryAnalyzer,
    Structure3D
)
from advanced_analytics import (
    MolecularVisualizer,
    DrugLikenessAnalyzer,
    ChemicalSpaceAnalyzer
)
from performance_optimized_engine import PerformanceOptimizedSearchEngine, PerformanceConfig
from scalable_search_engine import AminoAcidRecord

def demo_3d_structure_analysis():
    """演示3D结构分析功能"""
    print("=" * 60)
    print("3D结构分析演示")
    print("=" * 60)
    
    analyzer = Structure3DAnalyzer()
    
    # 创建示例PDB内容
    sample_pdb = """
ATOM      1  N   0A1 A   1      20.154  16.967  14.421  1.00 20.00           N  
ATOM      2  CA  0A1 A   1      19.030  16.045  14.421  1.00 20.00           C  
ATOM      3  C   0A1 A   1      17.618  16.696  14.421  1.00 20.00           C  
ATOM      4  O   0A1 A   1      17.618  17.928  14.421  1.00 20.00           O  
ATOM      5  CB  0A1 A   1      19.030  15.215  13.135  1.00 20.00           C  
ATOM      6  CG  0A1 A   1      19.030  15.215  11.849  1.00 20.00           C  
ATOM      7  CD1 0A1 A   1      17.849  15.215  11.135  1.00 20.00           C  
ATOM      8  CD2 0A1 A   1      20.211  15.215  11.135  1.00 20.00           C  
ATOM      9  CE1 0A1 A   1      17.849  15.215   9.849  1.00 20.00           C  
ATOM     10  CE2 0A1 A   1      20.211  15.215   9.849  1.00 20.00           C  
ATOM     11  CZ  0A1 A   1      19.030  15.215   9.135  1.00 20.00           C  
ATOM     12  OH  0A1 A   1      19.030  15.215   7.849  1.00 20.00           O  
CONECT    1    2
CONECT    2    1    3    5
CONECT    3    2    4
CONECT    4    3
CONECT    5    2    6
CONECT    6    5    7    8
CONECT    7    6    9
CONECT    8    6   10
CONECT    9    7   11
CONECT   10    8   11
CONECT   11    9   10   12
CONECT   12   11
"""
    
    # 解析3D结构
    print("解析PDB结构...")
    structure = analyzer.parse_pdb_structure(sample_pdb, "0A1")
    
    if structure:
        print(f"成功解析结构: {structure.amino_acid_id}")
        print(f"原子数: {len(structure.atoms)}")
        print(f"键数: {len(structure.bonds)}")
        print(f"分子体积: {structure.molecular_volume:.2f} Ų")
        print(f"表面积: {structure.surface_area:.2f} Ų")
        print(f"质心: ({structure.center_of_mass[0]:.2f}, {structure.center_of_mass[1]:.2f}, {structure.center_of_mass[2]:.2f})")
        
        # 提取几何特征
        print("\n提取几何特征...")
        geometric_features = analyzer.extract_geometric_features(structure)
        
        print(f"键长数量: {len(geometric_features.bond_lengths)}")
        if geometric_features.bond_lengths:
            print(f"平均键长: {np.mean(geometric_features.bond_lengths):.3f} Å")
            print(f"键长范围: {min(geometric_features.bond_lengths):.3f} - {max(geometric_features.bond_lengths):.3f} Å")
        
        print(f"键角数量: {len(geometric_features.bond_angles)}")
        if geometric_features.bond_angles:
            print(f"平均键角: {np.mean(geometric_features.bond_angles):.1f}°")
        
        print(f"环系统数量: {len(geometric_features.ring_systems)}")
        for i, ring in enumerate(geometric_features.ring_systems):
            print(f"  环 {i+1}: {ring['size']}元环, 类型: {ring['type']}")
        
        print(f"手性中心数量: {len(geometric_features.chiral_centers)}")
        for chiral in geometric_features.chiral_centers:
            print(f"  手性中心: 原子{chiral['atom_index']} ({chiral['atom_name']}), 构型: {chiral['chirality']}")
        
        # 几何验证
        print("\n几何验证...")
        validation = analyzer.validate_geometry(structure)
        print(f"整体几何有效性: {'通过' if validation['overall_valid'] else '未通过'}")
        print(f"有效键长: {len(validation['valid_bond_lengths'])}")
        print(f"无效键长: {len(validation['invalid_bond_lengths'])}")
        print(f"有效键角: {len(validation['valid_bond_angles'])}")
        print(f"无效键角: {len(validation['invalid_bond_angles'])}")
        
        return structure
    else:
        print("结构解析失败")
        return None

def demo_stereochemistry_analysis(structure: Structure3D):
    """演示立体化学分析"""
    print("\n" + "=" * 60)
    print("立体化学分析演示")
    print("=" * 60)
    
    analyzer = StereochemistryAnalyzer()
    
    # 手性分析
    print("手性分析...")
    chirality_analysis = analyzer.analyze_chirality(structure)
    
    print(f"是否为手性分子: {'是' if chirality_analysis['is_chiral'] else '否'}")
    print(f"手性中心数量: {chirality_analysis['num_chiral_centers']}")
    
    for descriptor in chirality_analysis['stereochemical_descriptors']:
        print(f"  原子 {descriptor['atom_index']} ({descriptor['atom_name']}): {descriptor['configuration']}")
        print(f"    优先级顺序: {descriptor['priority_order']}")
    
    # E/Z异构体检测
    print("\nE/Z异构体检测...")
    ez_isomers = analyzer.detect_ez_isomers(structure)
    
    if ez_isomers:
        print(f"发现 {len(ez_isomers)} 个潜在的E/Z异构体位点:")
        for i, ez_info in enumerate(ez_isomers):
            print(f"  位点 {i+1}: 双键 {ez_info['double_bond']}")
            print(f"    构型: {ez_info['configuration']}")
    else:
        print("未发现E/Z异构体位点")

def demo_machine_learning_classification():
    """演示机器学习分类"""
    print("\n" + "=" * 60)
    print("机器学习分类演示")
    print("=" * 60)
    
    # 初始化搜索引擎获取数据
    config = PerformanceConfig()
    search_engine = PerformanceOptimizedSearchEngine(config=config)
    
    # 获取所有氨基酸记录
    all_records = search_engine.database.get_all_amino_acids()
    print(f"数据库包含 {len(all_records)} 种氨基酸")
    
    if len(all_records) < 2:
        print("数据量不足，无法进行机器学习演示")
        return
    
    # 初始化ML分类器
    ml_classifier = MLClassifier()
    
    # 准备训练数据
    print("\n准备训练数据...")
    X, y = ml_classifier.prepare_training_data(all_records)
    print(f"特征维度: {X.shape}")
    print(f"样本数量: {len(y)}")
    print(f"类别数量: {len(set(y))}")
    
    # 训练模型
    print("\n训练机器学习模型...")
    training_results = ml_classifier.train_models(X, y, ['random_forest', 'svm'])
    
    for model_name, results in training_results.items():
        print(f"\n{model_name} 模型结果:")
        print(f"  准确率: {results['accuracy']:.3f}")
        print(f"  交叉验证: {results['cv_mean']:.3f} ± {results['cv_std']:.3f}")
        
        if results['feature_importance']:
            print("  重要特征 (前5个):")
            for feature_name, importance in results['feature_importance'][:5]:
                print(f"    {feature_name}: {importance:.3f}")
    
    # 预测示例
    if ml_classifier.is_trained and all_records:
        print("\n预测示例...")
        test_record = all_records[0]
        prediction = ml_classifier.predict(test_record, model_name='random_forest')
        
        print(f"测试氨基酸: {test_record.id}")
        print(f"预测类别: {prediction['predicted_class']}")
        print(f"置信度: {prediction['confidence']:.3f}")
        
        if prediction['class_probabilities']:
            print("类别概率:")
            for class_name, prob in sorted(prediction['class_probabilities'].items(), 
                                         key=lambda x: x[1], reverse=True)[:3]:
                print(f"  {class_name}: {prob:.3f}")
    
    return ml_classifier

def demo_drug_likeness_analysis():
    """演示药物相似性分析"""
    print("\n" + "=" * 60)
    print("药物相似性分析演示")
    print("=" * 60)
    
    # 初始化分析器
    drug_analyzer = DrugLikenessAnalyzer()
    
    # 获取氨基酸数据
    config = PerformanceConfig()
    search_engine = PerformanceOptimizedSearchEngine(config=config)
    all_records = search_engine.database.get_all_amino_acids()
    
    if not all_records:
        print("无氨基酸数据进行分析")
        return
    
    print(f"分析 {len(all_records)} 种氨基酸的药物相似性...")
    
    # 批量分析
    drug_scores = drug_analyzer.batch_analyze_drug_likeness(all_records)
    
    # 统计结果
    drug_like_count = sum(1 for score in drug_scores if score.drug_like)
    avg_score = np.mean([score.overall_score for score in drug_scores])
    
    print(f"\n分析结果:")
    print(f"  具有药物相似性的氨基酸: {drug_like_count}/{len(all_records)} ({drug_like_count/len(all_records)*100:.1f}%)")
    print(f"  平均药物相似性评分: {avg_score:.3f}")
    
    # 显示详细结果
    print(f"\n详细结果:")
    for i, (record, score) in enumerate(zip(all_records, drug_scores)):
        print(f"{record.id}:")
        print(f"  分子量: {score.molecular_weight:.1f} Da")
        print(f"  LogP: {score.logp:.2f}")
        print(f"  氢键供体: {score.hbd}")
        print(f"  氢键受体: {score.hba}")
        print(f"  TPSA: {score.tpsa:.1f} Ų")
        print(f"  Lipinski违规: {score.lipinski_violations}")
        print(f"  药物相似性: {'是' if score.drug_like else '否'} (评分: {score.overall_score:.3f})")
        print()

def demo_chemical_space_analysis():
    """演示化学空间分析"""
    print("\n" + "=" * 60)
    print("化学空间分析演示")
    print("=" * 60)
    
    try:
        # 初始化分析器
        space_analyzer = ChemicalSpaceAnalyzer()
        
        # 获取氨基酸数据
        config = PerformanceConfig()
        search_engine = PerformanceOptimizedSearchEngine(config=config)
        all_records = search_engine.database.get_all_amino_acids()
        
        if len(all_records) < 3:
            print("数据量不足，无法进行化学空间分析")
            return
        
        print(f"分析 {len(all_records)} 种氨基酸的化学空间...")
        
        # 化学空间分析
        space_analysis = space_analyzer.analyze_chemical_space(all_records, feature_type='basic')
        
        print(f"\nPCA分析结果:")
        print(f"  特征维度: {space_analysis['features'].shape}")
        print(f"  主成分解释方差比:")
        for i, ratio in enumerate(space_analysis['pca']['explained_variance_ratio'][:5]):
            print(f"    PC{i+1}: {ratio:.3f}")
        print(f"  前3个主成分累计解释方差: {space_analysis['pca']['cumulative_variance'][2]:.3f}")
        
        print(f"\nK-means聚类结果:")
        kmeans_result = space_analysis['clustering']['kmeans']
        print(f"  最优聚类数: {kmeans_result['n_clusters']}")
        print(f"  轮廓系数: {kmeans_result['silhouette_score']:.3f}")
        
        # 显示聚类分布
        cluster_distribution = {}
        for label, cluster in zip(space_analysis['labels'], kmeans_result['labels']):
            if cluster not in cluster_distribution:
                cluster_distribution[cluster] = []
            cluster_distribution[cluster].append(label)
        
        print(f"  聚类分布:")
        for cluster_id, members in cluster_distribution.items():
            print(f"    聚类 {cluster_id}: {', '.join(members)}")
        
        print(f"\nDBSCAN聚类结果:")
        dbscan_result = space_analysis['clustering']['dbscan']
        print(f"  聚类数: {dbscan_result['n_clusters']}")
        print(f"  噪声点数: {dbscan_result['n_noise']}")
        
    except ImportError as e:
        print(f"化学空间分析需要scikit-learn库: {e}")

def demo_molecular_visualization():
    """演示分子可视化"""
    print("\n" + "=" * 60)
    print("分子可视化演示")
    print("=" * 60)
    
    try:
        visualizer = MolecularVisualizer()
        
        # 获取氨基酸数据
        config = PerformanceConfig()
        search_engine = PerformanceOptimizedSearchEngine(config=config)
        all_records = search_engine.database.get_all_amino_acids()
        
        if not all_records:
            print("无氨基酸数据进行可视化")
            return
        
        print("创建分子性质分布图...")
        
        # 分子量分布
        mw_fig = visualizer.create_property_distribution(all_records, 'molecular_weight')
        if mw_fig:
            print("  分子量分布图已创建")
        
        # 原子数分布
        atom_fig = visualizer.create_property_distribution(all_records, 'atom_count')
        if atom_fig:
            print("  原子数分布图已创建")
        
        print("可视化功能演示完成")
        
    except ImportError as e:
        print(f"分子可视化需要plotly库: {e}")

def main():
    """主程序"""
    print("可扩展非天然氨基酸PDB搜索引擎 - 第三阶段高级功能演示")
    print("=" * 80)
    
    # 演示3D结构分析
    structure = demo_3d_structure_analysis()
    
    # 演示立体化学分析
    if structure:
        demo_stereochemistry_analysis(structure)
    
    # 演示机器学习分类
    ml_classifier = demo_machine_learning_classification()
    
    # 演示药物相似性分析
    demo_drug_likeness_analysis()
    
    # 演示化学空间分析
    demo_chemical_space_analysis()
    
    # 演示分子可视化
    demo_molecular_visualization()
    
    print("\n" + "=" * 80)
    print("第三阶段高级功能演示完成！")
    print("✓ 3D结构分析和几何验证")
    print("✓ 立体化学和异构体识别")
    print("✓ 机器学习分类和预测")
    print("✓ 药物相似性分析")
    print("✓ 化学空间分析和聚类")
    print("✓ 分子可视化功能")
    print("=" * 80)

if __name__ == "__main__":
    main()
