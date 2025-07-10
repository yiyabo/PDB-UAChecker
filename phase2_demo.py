#!/usr/bin/env python3
"""
第二阶段性能优化演示脚本
展示LSH索引、分层缓存、并行处理和高级指纹算法
"""

import os
import time
import json
import numpy as np
from typing import Dict, List, Any

from performance_optimized_engine import (
    PerformanceOptimizedSearchEngine,
    PerformanceConfig,
    OptimizedCompatibilityWrapper,
    AminoAcidRecord
)

def demo_lsh_index():
    """演示LSH索引性能"""
    print("=" * 60)
    print("LSH索引性能演示")
    print("=" * 60)
    
    # 创建配置
    config = PerformanceConfig(
        lsh_num_tables=8,
        lsh_hash_size=12,
        memory_cache_size=5000,
        max_workers=4
    )
    
    # 初始化引擎
    engine = PerformanceOptimizedSearchEngine(config=config)
    
    # 获取所有氨基酸记录
    all_records = engine.database.get_all_amino_acids()
    print(f"数据库包含 {len(all_records)} 种氨基酸")
    
    # 选择一个记录作为查询
    query_record = all_records[0]
    print(f"查询氨基酸: {query_record.id} - {query_record.name}")
    
    # 生成ECFP指纹
    print("\n生成ECFP指纹...")
    start_time = time.time()
    query_fingerprint = engine.fingerprint_generator.generate_fingerprint(
        query_record, 'ecfp'
    )
    fingerprint_time = time.time() - start_time
    print(f"指纹生成时间: {fingerprint_time*1000:.2f}ms")
    print(f"指纹维度: {len(query_fingerprint)}")
    print(f"设置位数: {np.sum(query_fingerprint)}")
    
    # 执行LSH相似性搜索
    print("\n执行LSH相似性搜索...")
    start_time = time.time()
    similar_ids = engine.lsh_index.query_similar(query_fingerprint, top_k=5)
    lsh_time = time.time() - start_time
    print(f"LSH搜索时间: {lsh_time*1000:.2f}ms")
    print(f"找到 {len(similar_ids)} 个相似氨基酸")
    
    # 计算精确相似度
    print("\n计算精确相似度...")
    similarities = []
    for amino_id in similar_ids:
        record = engine.database.get_amino_acid(amino_id)
        if record:
            record_fingerprint = engine.fingerprint_generator.generate_fingerprint(
                record, 'ecfp'
            )
            similarity = engine.fingerprint_generator.get_fingerprint_similarity(
                query_fingerprint, record_fingerprint
            )
            similarities.append((amino_id, similarity))
    
    # 显示结果
    print("\nLSH相似性搜索结果:")
    for amino_id, similarity in sorted(similarities, key=lambda x: x[1], reverse=True):
        record = engine.database.get_amino_acid(amino_id)
        print(f"  {amino_id} - {record.name}: 相似度 {similarity:.3f}")
    
    # 显示LSH索引统计
    lsh_stats = engine.lsh_index.get_statistics()
    print(f"\nLSH索引统计:")
    print(f"  哈希表数量: {lsh_stats['num_tables']}")
    print(f"  哈希大小: {lsh_stats['hash_size']}")
    print(f"  总桶数: {lsh_stats['total_buckets']}")
    print(f"  平均桶大小: {lsh_stats['average_bucket_size']}")
    
    return engine

def demo_advanced_fingerprints(engine):
    """演示高级分子指纹"""
    print("\n" + "=" * 60)
    print("高级分子指纹演示")
    print("=" * 60)
    
    # 获取所有氨基酸记录
    all_records = engine.database.get_all_amino_acids()
    
    # 选择一个记录
    record = all_records[0]
    print(f"氨基酸: {record.id} - {record.name}")
    print(f"分子式: {record.molecular_formula}")
    print(f"SMILES: {record.smiles}")
    
    # 生成不同类型的指纹
    print("\n生成不同类型的指纹...")
    fingerprint_types = ['basic', 'ecfp', 'maccs', 'topological', 'pharmacophore']
    
    for fp_type in fingerprint_types:
        start_time = time.time()
        fingerprint = engine.fingerprint_generator.generate_fingerprint(record, fp_type)
        gen_time = time.time() - start_time
        
        print(f"\n{fp_type.upper()}指纹:")
        print(f"  维度: {len(fingerprint)}")
        print(f"  设置位数: {np.sum(fingerprint)}")
        print(f"  生成时间: {gen_time*1000:.2f}ms")
        
        # 显示部分指纹
        if len(fingerprint) <= 64:
            print(f"  指纹: {fingerprint}")
        else:
            print(f"  指纹片段: {fingerprint[:20]}...")
    
    # 指纹相似性比较
    print("\n指纹相似性比较:")
    
    # 选择两个记录
    if len(all_records) >= 2:
        record1 = all_records[0]
        record2 = all_records[1]
        
        print(f"比较 {record1.id} 和 {record2.id}:")
        
        for fp_type in fingerprint_types:
            fp1 = engine.fingerprint_generator.generate_fingerprint(record1, fp_type)
            fp2 = engine.fingerprint_generator.generate_fingerprint(record2, fp_type)
            
            similarity = engine.fingerprint_generator.get_fingerprint_similarity(fp1, fp2)
            print(f"  {fp_type.upper()} 相似度: {similarity:.3f}")
    
    return engine

def demo_cache_system(engine):
    """演示分层缓存系统"""
    print("\n" + "=" * 60)
    print("分层缓存系统演示")
    print("=" * 60)
    
    # 准备测试查询
    test_queries = [
        {'residue_name': '0A1'},
        {'residue_name': '2AG'},
        {'molecular_formula': 'C10H13NO3'}
    ]
    
    # 第一次查询（无缓存）
    print("\n第一次查询（无缓存）:")
    first_times = []
    
    for i, query in enumerate(test_queries):
        start_time = time.time()
        results = engine.optimized_search(query)
        query_time = time.time() - start_time
        
        first_times.append(query_time)
        print(f"  查询 {i+1}: {query_time*1000:.2f}ms, 找到 {len(results)} 个结果")
    
    # 第二次查询（有缓存）
    print("\n第二次查询（有缓存）:")
    second_times = []
    
    for i, query in enumerate(test_queries):
        start_time = time.time()
        results = engine.optimized_search(query)
        query_time = time.time() - start_time
        
        second_times.append(query_time)
        print(f"  查询 {i+1}: {query_time*1000:.2f}ms, 找到 {len(results)} 个结果")
    
    # 计算加速比
    avg_first = sum(first_times) / len(first_times)
    avg_second = sum(second_times) / len(second_times)
    speedup = avg_first / avg_second if avg_second > 0 else 0
    
    print(f"\n缓存性能:")
    print(f"  无缓存平均时间: {avg_first*1000:.2f}ms")
    print(f"  有缓存平均时间: {avg_second*1000:.2f}ms")
    print(f"  加速比: {speedup:.2f}x")
    
    # 显示缓存统计
    cache_stats = engine.cache_manager.get_statistics()
    print(f"\n缓存统计:")
    print(f"  内存缓存大小: {cache_stats['memory_cache_size']}")
    print(f"  缓存命中率: {cache_stats['hit_rate']:.3f}")
    print(f"  内存命中: {cache_stats['cache_stats']['memory_hits']}")
    print(f"  磁盘命中: {cache_stats['cache_stats']['disk_hits']}")
    print(f"  缓存未命中: {cache_stats['cache_stats']['misses']}")
    
    return engine

def demo_parallel_processing(engine):
    """演示并行处理"""
    print("\n" + "=" * 60)
    print("并行处理演示")
    print("=" * 60)
    
    # 查找测试PDB文件
    test_dir = "test"
    pdb_files = []
    
    if os.path.exists(test_dir):
        pdb_files = [os.path.join(test_dir, f) for f in os.listdir(test_dir) if f.endswith('.pdb')]
    
    if not pdb_files:
        print("未找到PDB文件进行测试")
        return engine
    
    print(f"找到 {len(pdb_files)} 个PDB文件:")
    for pdb_file in pdb_files:
        print(f"  - {os.path.basename(pdb_file)}")
    
    # 串行处理
    print("\n串行处理:")
    start_time = time.time()
    
    # 使用兼容包装器
    from scalable_search_engine import CompatibilityWrapper
    serial_wrapper = CompatibilityWrapper()
    serial_results = serial_wrapper.search_pdb_files(pdb_files)
    
    serial_time = time.time() - start_time
    print(f"  总时间: {serial_time:.3f}s")
    print(f"  找到 {len(serial_results)} 个匹配")
    print(f"  平均每个文件: {serial_time/len(pdb_files)*1000:.2f}ms")
    
    # 并行处理
    print("\n并行处理:")
    start_time = time.time()
    
    parallel_wrapper = OptimizedCompatibilityWrapper()
    parallel_results = parallel_wrapper.search_pdb_files(pdb_files)
    
    parallel_time = time.time() - start_time
    print(f"  总时间: {parallel_time:.3f}s")
    print(f"  找到 {len(parallel_results)} 个匹配")
    print(f"  平均每个文件: {parallel_time/len(pdb_files)*1000:.2f}ms")
    
    # 计算加速比
    speedup = serial_time / parallel_time if parallel_time > 0 else 0
    print(f"\n并行加速比: {speedup:.2f}x")
    
    # 显示并行统计
    parallel_stats = engine.parallel_engine.get_statistics()
    print(f"\n并行处理统计:")
    print(f"  总并行查询: {parallel_stats['total_parallel_queries']}")
    print(f"  平均加速比: {parallel_stats['avg_parallel_speedup']:.2f}x")
    
    return engine

def main():
    """主程序"""
    print("可扩展非天然氨基酸PDB搜索引擎 - 第二阶段性能优化演示")
    print("=" * 80)
    
    # 演示LSH索引
    engine = demo_lsh_index()
    
    # 演示高级分子指纹
    engine = demo_advanced_fingerprints(engine)
    
    # 演示分层缓存系统
    engine = demo_cache_system(engine)
    
    # 演示并行处理
    engine = demo_parallel_processing(engine)
    
    # 显示综合统计
    print("\n" + "=" * 60)
    print("综合性能统计")
    print("=" * 60)
    
    stats = engine.get_comprehensive_statistics()
    
    print("性能指标:")
    print(f"  平均查询时间: {stats['performance_stats']['avg_query_time']*1000:.2f}ms")
    print(f"  缓存命中率: {stats['cache_stats']['hit_rate']:.3f}")
    print(f"  内存使用: {stats['memory_stats']['memory_usage']['total']/1024/1024:.1f}MB")
    
    print("\n优化建议:")
    for rec in stats['optimization_recommendations']:
        print(f"  - {rec}")
    
    print("\n" + "=" * 80)
    print("第二阶段性能优化演示完成！")
    print("=" * 80)
    
    # 关闭引擎
    engine.shutdown()

if __name__ == "__main__":
    main()
