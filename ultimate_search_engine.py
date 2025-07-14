#!/usr/bin/env python3
"""
终极搜索引擎 - 整合所有搜索功能
包括：精确搜索、相似性搜索、聚类分析、推荐系统
"""

import sys
import time
from typing import List, Dict, Optional, Union

sys.path.append('.')

# 导入现有的搜索引擎
from scalable_search_engine import ScalableSearchEngine
from enhanced_search_engine import FingerprintSearchEngine, SimilarityResult

class UltimateSearchEngine:
    """终极氨基酸搜索引擎"""
    
    def __init__(self):
        print("🚀 初始化终极搜索引擎...")
        
        # 初始化子引擎
        self.exact_engine = ScalableSearchEngine()
        self.similarity_engine = FingerprintSearchEngine()
        
        print("✅ 终极搜索引擎就绪！")
        print(f"   支持 {len(self.similarity_engine.amino_acids_data)} 种氨基酸")
        print(f"   功能: 精确搜索 + 相似性搜索 + 推荐系统")
    
    def search_comprehensive(self, query: Union[str, Dict], 
                           include_similar: bool = True,
                           similarity_threshold: float = 0.6,
                           max_results: int = 10) -> Dict:
        """综合搜索 - 同时进行精确搜索和相似性搜索"""
        
        start_time = time.time()
        results = {
            'exact_matches': [],
            'similar_matches': [],
            'search_time': 0,
            'total_results': 0
        }
        
        # 1. 精确搜索
        if isinstance(query, dict):
            exact_results = self.exact_engine.search(query)
            results['exact_matches'] = exact_results
        
        # 2. 相似性搜索
        if include_similar and isinstance(query, str):
            # 假设字符串是SMILES或氨基酸ID
            if len(query) <= 5 and query.isupper():
                # 可能是氨基酸ID
                similar_results = self.similarity_engine.search_similar_by_id(
                    query, similarity_threshold, max_results
                )
            else:
                # 可能是SMILES
                similar_results = self.similarity_engine.search_similar_by_smiles(
                    query, similarity_threshold, max_results
                )
            results['similar_matches'] = similar_results
        
        results['search_time'] = time.time() - start_time
        results['total_results'] = len(results['exact_matches']) + len(results['similar_matches'])
        
        return results
    
    def recommend_similar(self, amino_acid_id: str, 
                         count: int = 5) -> List[SimilarityResult]:
        """推荐系统 - 基于给定氨基酸推荐相似的"""
        
        print(f"🎯 为 {amino_acid_id} 推荐相似的氨基酸...")
        
        return self.similarity_engine.search_similar_by_id(
            amino_acid_id, 
            similarity_threshold=0.5, 
            max_results=count
        )
    
    def cluster_amino_acids(self, similarity_threshold: float = 0.7) -> Dict[str, List[str]]:
        """氨基酸聚类分析"""
        
        print(f"🔬 进行氨基酸聚类分析 (相似性阈值: {similarity_threshold})...")
        
        clusters = {}
        processed = set()
        
        for amino_id in self.similarity_engine.amino_acids_data.keys():
            if amino_id in processed:
                continue
            
            # 找到与当前氨基酸相似的所有氨基酸
            similar = self.similarity_engine.search_similar_by_id(
                amino_id, similarity_threshold, max_results=50
            )
            
            if similar:
                cluster_name = f"cluster_{amino_id}"
                cluster_members = [amino_id]
                
                for result in similar:
                    if result.amino_acid_id not in processed:
                        cluster_members.append(result.amino_acid_id)
                        processed.add(result.amino_acid_id)
                
                if len(cluster_members) > 1:
                    clusters[cluster_name] = cluster_members
                
                processed.add(amino_id)
        
        return clusters
    
    def analyze_chemical_space(self) -> Dict:
        """化学空间分析"""
        
        print("🧬 分析化学空间...")
        
        # 找到最相似的氨基酸对
        similar_pairs = self.similarity_engine.find_most_similar_pairs(top_n=10)
        
        # 统计相似性分布
        all_similarities = [pair[2] for pair in similar_pairs]
        
        analysis = {
            'most_similar_pairs': similar_pairs[:5],
            'avg_similarity': sum(all_similarities) / len(all_similarities),
            'max_similarity': max(all_similarities),
            'min_similarity': min(all_similarities),
            'total_amino_acids': len(self.similarity_engine.amino_acids_data)
        }
        
        return analysis
    
    def search_by_features(self, features: List[str], 
                          similarity_search: bool = False) -> List:
        """基于化学特征搜索"""
        
        print(f"🔍 基于特征搜索: {features}")
        
        # 这里可以扩展为更复杂的特征搜索
        # 目前使用现有的搜索引擎
        results = []
        
        # 可以添加特征组合搜索逻辑
        
        return results

def demo_ultimate_search():
    """演示终极搜索引擎的功能"""
    
    print("🎪 终极搜索引擎功能演示")
    print("=" * 70)
    
    engine = UltimateSearchEngine()
    
    # 演示1: 综合搜索
    print("\n📋 演示1: 综合搜索 - 苯丙氨酸")
    results = engine.search_comprehensive("PHE", include_similar=True)
    
    print(f"   搜索时间: {results['search_time']*1000:.1f}ms")
    print(f"   总结果数: {results['total_results']}")
    
    if results['similar_matches']:
        print("   相似氨基酸:")
        for i, result in enumerate(results['similar_matches'][:3], 1):
            print(f"     {i}. {result.amino_acid_id}: {result.name} (相似性: {result.similarity_score:.3f})")
    
    # 演示2: 推荐系统
    print("\n📋 演示2: 推荐系统 - 基于色氨酸推荐")
    recommendations = engine.recommend_similar("TRP", count=3)
    
    for i, rec in enumerate(recommendations, 1):
        print(f"   {i}. {rec.amino_acid_id}: {rec.name}")
        print(f"      相似性: {rec.similarity_score:.3f}")
        print(f"      分子式: {rec.molecular_formula}")
    
    # 演示3: 化学空间分析
    print("\n📋 演示3: 化学空间分析")
    analysis = engine.analyze_chemical_space()
    
    print(f"   数据库规模: {analysis['total_amino_acids']} 种氨基酸")
    print(f"   平均相似性: {analysis['avg_similarity']:.3f}")
    print(f"   最高相似性: {analysis['max_similarity']:.3f}")
    
    print("   最相似的氨基酸对:")
    for i, (id1, id2, sim) in enumerate(analysis['most_similar_pairs'], 1):
        name1 = engine.similarity_engine.amino_acids_data[id1]['name']
        name2 = engine.similarity_engine.amino_acids_data[id2]['name']
        print(f"     {i}. {id1}({name1}) vs {id2}({name2}) - {sim:.3f}")
    
    # 演示4: 聚类分析
    print("\n📋 演示4: 聚类分析 (高相似性阈值)")
    clusters = engine.cluster_amino_acids(similarity_threshold=0.9)
    
    print(f"   发现 {len(clusters)} 个高相似性聚类:")
    for cluster_name, members in list(clusters.items())[:3]:
        print(f"     {cluster_name}: {len(members)} 个成员")
        print(f"       成员: {', '.join(members[:5])}")

def performance_benchmark():
    """性能基准测试"""
    
    print("\n⚡ 性能基准测试")
    print("=" * 40)
    
    engine = UltimateSearchEngine()
    
    # 测试搜索性能
    test_queries = ["PHE", "TRP", "VAL", "ADAM", "CSE"]
    
    total_time = 0
    for query in test_queries:
        start_time = time.time()
        results = engine.search_comprehensive(query, include_similar=True)
        search_time = time.time() - start_time
        total_time += search_time
        
        print(f"   {query}: {search_time*1000:.1f}ms ({results['total_results']} 结果)")
    
    avg_time = total_time / len(test_queries)
    print(f"\n   平均搜索时间: {avg_time*1000:.1f}ms")
    print(f"   性能评级: {'🚀 优秀' if avg_time < 0.1 else '✅ 良好' if avg_time < 0.5 else '⚠️ 需优化'}")

def main():
    """主函数"""
    
    print("🌟 终极氨基酸搜索引擎")
    print("=" * 50)
    
    # 运行演示
    demo_ultimate_search()
    
    # 性能测试
    performance_benchmark()
    
    print(f"\n🎉 终极搜索引擎演示完成！")
    print(f"\n💡 现在您拥有:")
    print("   ✅ 精确搜索 (分子式、残基名等)")
    print("   ✅ 相似性搜索 (基于分子指纹)")
    print("   ✅ 推荐系统 (智能推荐相似氨基酸)")
    print("   ✅ 聚类分析 (化学结构分组)")
    print("   ✅ 化学空间分析 (整体数据洞察)")
    print("   ✅ 高性能 (毫秒级响应)")

if __name__ == "__main__":
    main()
