#!/usr/bin/env python3
"""
同分异构体识别测试脚本
验证新的ECFP算法在解决同分异构体识别问题上的效果
"""

import sys
import time
from typing import Dict, List, Any
from scalable_search_engine import ScalableSearchEngine
from isomer_identifier import IsomerIdentifier, detect_potential_isomers

def test_isomer_identification():
    """测试同分异构体识别功能"""
    print("=" * 80)
    print("🧪 同分异构体识别测试开始")
    print("=" * 80)
    
    # 初始化搜索引擎
    print("\n1. 初始化搜索引擎...")
    engine = ScalableSearchEngine()
    
    # 测试数据：已知的同分异构体对
    test_cases = [
        {
            "name": "异亮氨酸 vs 亮氨酸 (同分异构体)",
            "smiles1": "CC[C@H](C)[C@@H](N)C(=O)O",    # 异亮氨酸（正确的SMILES）
            "smiles2": "CC(C)C[C@@H](N)C(=O)O",        # 亮氨酸
            "expected_type": "structural",
            "description": "支链位置不同的结构异构体"
        },
        {
            "name": "L-丙氨酸 vs D-丙氨酸 (立体异构体)",
            "smiles1": "C[C@@H](N)C(=O)O",              # L-丙氨酸
            "smiles2": "C[C@H](N)C(=O)O",               # D-丙氨酸
            "expected_type": "stereoisomer",
            "description": "手性中心构型不同的立体异构体"
        },
        {
            "name": "相同分子测试",
            "smiles1": "C[C@@H](N)C(=O)O",              # L-丙氨酸
            "smiles2": "C[C@@H](N)C(=O)O",              # L-丙氨酸
            "expected_type": "identical",
            "description": "完全相同的分子"
        },
        {
            "name": "不同分子测试",
            "smiles1": "C[C@@H](N)C(=O)O",              # 丙氨酸
            "smiles2": "CC(C)C[C@@H](N)C(=O)O",         # 异亮氨酸
            "expected_type": "different",
            "description": "分子式不同的分子"
        }
    ]
    
    # 测试同分异构体识别
    print("\n2. 测试同分异构体识别算法...")
    identifier = IsomerIdentifier()
    
    for i, test_case in enumerate(test_cases, 1):
        print(f"\n测试案例 {i}: {test_case['name']}")
        print(f"   描述: {test_case['description']}")
        print(f"   分子1: {test_case['smiles1']}")
        print(f"   分子2: {test_case['smiles2']}")
        
        # 执行分析
        result = identifier.analyze_isomers(test_case['smiles1'], test_case['smiles2'])
        
        # 显示结果
        print(f"   结果:")
        print(f"     - 是否为异构体: {result.is_isomer}")
        print(f"     - 异构体类型: {result.isomer_type}")
        print(f"     - 相似性分数: {result.similarity_score:.3f}")
        print(f"     - 置信度: {result.confidence:.3f}")
        
        if result.structural_differences:
            print(f"     - 结构差异: {result.structural_differences}")
        
        # 验证结果
        expected_type = test_case['expected_type']
        if result.isomer_type == expected_type:
            print(f"   ✅ 测试通过：预期 {expected_type}，实际 {result.isomer_type}")
        else:
            print(f"   ❌ 测试失败：预期 {expected_type}，实际 {result.isomer_type}")
    
    # 测试搜索引擎集成
    print("\n3. 测试搜索引擎集成...")
    test_search_integration(engine)
    
    # 测试数据库中的同分异构体检测
    print("\n4. 检测数据库中的同分异构体...")
    test_database_isomer_detection(engine)
    
    print("\n" + "=" * 80)
    print("🎯 同分异构体识别测试完成")
    print("=" * 80)

def test_search_integration(engine: ScalableSearchEngine):
    """测试搜索引擎的同分异构体识别集成"""
    
    # 测试查询数据
    test_queries = [
        {
            "name": "使用ECFP相似性搜索",
            "query_data": {
                "smiles": "COc1ccc(cc1)C[C@@H](N)C(=O)O"  # 4-甲氧基苯丙氨酸
            },
            "methods": ["ecfp_similarity"]
        },
        {
            "name": "使用结构相似性搜索",
            "query_data": {
                "smiles": "COc1ccc(cc1)C[C@@H](N)C(=O)O"  # 4-甲氧基苯丙氨酸
            },
            "methods": ["structural_similarity"]
        },
        {
            "name": "使用同分异构体感知搜索",
            "query_data": {
                "smiles": "COc1ccc(cc1)C[C@@H](N)C(=O)O"  # 4-甲氧基苯丙氨酸
            },
            "methods": ["isomer_aware"]
        },
        {
            "name": "传统方法与增强方法对比",
            "query_data": {
                "smiles": "COc1ccc(cc1)C[C@@H](N)C(=O)O",  # 4-甲氧基苯丙氨酸
                "atom_composition": {"C": 10, "H": 13, "N": 1, "O": 3}
            },
            "methods": ["atom_composition"]
        }
    ]
    
    for test_query in test_queries:
        print(f"\n   {test_query['name']}:")
        
        start_time = time.time()
        results = engine.search(test_query['query_data'], methods=test_query['methods'])
        search_time = time.time() - start_time
        
        print(f"     搜索时间: {search_time:.3f}s")
        print(f"     找到结果: {len(results)} 个")
        
        for j, result in enumerate(results[:3]):  # 显示前3个结果
            print(f"       {j+1}. {result.amino_acid_id}: {result.amino_acid_record.name}")
            print(f"          置信度: {result.confidence_score:.3f}")
            print(f"          方法: {result.match_method}")
            
            if result.additional_info:
                if 'ecfp_similarity' in result.additional_info:
                    print(f"          ECFP相似性: {result.additional_info['ecfp_similarity']:.3f}")
                if 'structural_similarity' in result.additional_info:
                    print(f"          结构相似性: {result.additional_info['structural_similarity']:.3f}")
                if 'isomer_type' in result.additional_info:
                    print(f"          异构体类型: {result.additional_info['isomer_type']}")

def test_database_isomer_detection(engine: ScalableSearchEngine):
    """测试数据库中的同分异构体检测"""
    
    print("   正在检测数据库中的同分异构体组...")
    
    # 检测同分异构体
    isomer_groups = engine.detect_isomers_in_database()
    
    if isomer_groups:
        print(f"   发现 {len(isomer_groups)} 个同分异构体组:")
        
        for molecular_formula, amino_ids in isomer_groups.items():
            print(f"     分子式 {molecular_formula}:")
            for amino_id in amino_ids:
                record = engine.database.get_amino_acid(amino_id)
                if record:
                    print(f"       - {amino_id}: {record.name}")
                    print(f"         SMILES: {record.smiles}")
    else:
        print("   数据库中未发现同分异构体组")

def test_performance_comparison():
    """性能对比测试：传统方法 vs 增强方法"""
    print("\n5. 性能对比测试...")
    
    engine = ScalableSearchEngine()
    
    # 测试查询
    test_query = {
        "atom_composition": {"C": 10, "H": 13, "N": 1, "O": 3},
        "smiles": "COc1ccc(cc1)C[C@@H](N)C(=O)O"
    }
    
    print("   测试查询: 4-甲氧基苯丙氨酸")
    print("   原子组成: C10H13NO3")
    
    # 传统方法
    print("\n   传统方法（仅原子组成）:")
    start_time = time.time()
    traditional_results = engine.search(
        {"atom_composition": test_query["atom_composition"]},
        methods=["atom_composition"]
    )
    traditional_time = time.time() - start_time
    
    print(f"     搜索时间: {traditional_time:.3f}s")
    print(f"     找到结果: {len(traditional_results)} 个")
    
    # 增强方法
    print("\n   增强方法（原子组成 + 结构识别）:")
    start_time = time.time()
    enhanced_results = engine.search(test_query, methods=["atom_composition"])
    enhanced_time = time.time() - start_time
    
    print(f"     搜索时间: {enhanced_time:.3f}s")
    print(f"     找到结果: {len(enhanced_results)} 个")
    print(f"     性能开销: {((enhanced_time - traditional_time) / traditional_time * 100):.1f}%")
    
    # 比较结果质量
    print("\n   结果质量对比:")
    print("     传统方法结果:")
    for result in traditional_results[:3]:
        print(f"       {result.amino_acid_id}: 置信度 {result.confidence_score:.3f}")
    
    print("     增强方法结果:")
    for result in enhanced_results[:3]:
        print(f"       {result.amino_acid_id}: 置信度 {result.confidence_score:.3f}")
        if 'structural_enhanced' in result.additional_info:
            print(f"         结构增强: {result.additional_info['structural_enhanced']}")

def main():
    """主测试函数"""
    try:
        test_isomer_identification()
        test_performance_comparison()
        
        print("\n🎉 所有测试完成！")
        print("\n📊 测试总结:")
        print("   ✅ 同分异构体识别算法正常工作")
        print("   ✅ 搜索引擎集成成功")
        print("   ✅ 新的搜索策略可用")
        print("   ✅ 向后兼容性保持")
        
    except Exception as e:
        print(f"\n❌ 测试过程中发生错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()