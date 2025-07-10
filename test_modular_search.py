#!/usr/bin/env python3
"""
模块化搜索引擎测试脚本
验证重构后的搜索引擎功能
"""

import sys
import os

# 添加项目路径到 Python 路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.search import ScalableSearchEngine
from src.core.models import AminoAcidRecord

def test_modular_search_engine():
    """测试模块化搜索引擎"""
    print("=== 模块化搜索引擎测试 ===")
    print()
    
    try:
        # 初始化搜索引擎
        print("1. 初始化搜索引擎...")
        engine = ScalableSearchEngine()
        print("   ✓ 搜索引擎初始化成功")
        print()
        
        # 获取统计信息
        print("2. 获取统计信息...")
        stats = engine.get_search_statistics()
        print(f"   数据库大小: {stats['database_size']}")
        print(f"   可用搜索策略: {', '.join(stats['available_strategies'])}")
        isomer_symbol = '✓' if stats['isomer_identifier_enabled'] else '✗'
        atomic_symbol = '✓' if stats['atomic_analyzer_enabled'] else '✗'
        print(f"   同分异构体识别: {isomer_symbol}")
        print(f"   原子级分析器: {atomic_symbol}")
        print()
        
        # 测试残基名搜索
        print("3. 测试残基名搜索...")
        query_data = {'residue_name': 'TYR'}
        results = engine.search(query_data, methods=['residue_name'])
        
        if results:
            print(f"   找到 {len(results)} 个结果")
            for result in results[:2]:
                print(f"   - {result.amino_acid_id}: {result.amino_acid_record.name} (置信度: {result.confidence_score:.2f})")
        else:
            print("   未找到结果")
        print()
        
        # 测试分子式搜索
        print("4. 测试分子式搜索...")
        query_data = {'molecular_formula': 'C9H11NO3'}
        results = engine.search(query_data, methods=['molecular_formula'])
        
        if results:
            print(f"   找到 {len(results)} 个结果")
            for result in results[:2]:
                print(f"   - {result.amino_acid_id}: {result.amino_acid_record.name} (置信度: {result.confidence_score:.2f})")
        else:
            print("   未找到结果")
        print()
        
        # 测试分子量搜索
        print("5. 测试分子量搜索...")
        query_data = {'molecular_weight': 181.19, 'weight_tolerance': 2.0}
        results = engine.search(query_data, methods=['molecular_weight'])
        
        if results:
            print(f"   找到 {len(results)} 个结果")
            for result in results[:2]:
                actual_weight = result.amino_acid_record.molecular_weight
                print(f"   - {result.amino_acid_id}: {result.amino_acid_record.name} ({actual_weight:.2f} Da, 置信度: {result.confidence_score:.2f})")
        else:
            print("   未找到结果")
        print()
        
        # 测试多策略搜索
        print("6. 测试多策略搜索...")
        query_data = {
            'residue_name': 'TYR',
            'molecular_formula': 'C9H11NO3',
            'molecular_weight': 181.19
        }
        results = engine.search(query_data, methods=['residue_name', 'molecular_formula', 'molecular_weight'])
        
        if results:
            print(f"   找到 {len(results)} 个结果")
            for result in results[:3]:
                print(f"   - {result.amino_acid_id}: {result.amino_acid_record.name} (方法: {result.match_method}, 置信度: {result.confidence_score:.2f})")
        else:
            print("   未找到结果")
        print()
        
        # 测试缓存统计
        print("7. 缓存统计信息...")
        residue_stats = engine.residue_matcher.get_match_statistics()
        fingerprint_stats = engine.fingerprint_matcher.get_fingerprint_statistics()
        print(f"   残基名缓存: {residue_stats['exact_cache_size']} 个精确匹配, {residue_stats['fuzzy_cache_size']} 个模糊匹配")
        print(f"   指纹缓存: {fingerprint_stats['cache_size']} 个指纹")
        print()
        
        print("✓ 模块化搜索引擎测试完成！")
        return True
        
    except Exception as e:
        print(f"\u2717 测试失败: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_modular_search_engine()
    sys.exit(0 if success else 1)