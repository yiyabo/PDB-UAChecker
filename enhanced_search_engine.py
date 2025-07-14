#!/usr/bin/env python3
"""
增强版搜索引擎 - 支持基于分子指纹的相似性搜索
"""

import sqlite3
import json
import sys
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass

sys.path.append('.')

@dataclass
class SimilarityResult:
    """相似性搜索结果"""
    amino_acid_id: str
    name: str
    molecular_formula: str
    smiles: str
    similarity_score: float
    similarity_type: str

class FingerprintSearchEngine:
    """基于分子指纹的搜索引擎"""
    
    def __init__(self):
        self.db_path = 'amino_acids.db'
        self._load_fingerprints()
    
    def _load_fingerprints(self):
        """加载所有氨基酸的指纹数据"""
        print("🔄 加载分子指纹数据...")
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute('''
            SELECT id, name, molecular_formula, smiles, fingerprints 
            FROM amino_acids 
            WHERE fingerprints IS NOT NULL
        ''')
        
        self.amino_acids_data = {}
        self.fingerprints_cache = {}
        
        for row in cursor.fetchall():
            amino_id, name, formula, smiles, fingerprints_json = row
            
            self.amino_acids_data[amino_id] = {
                'name': name,
                'molecular_formula': formula,
                'smiles': smiles
            }
            
            # 解析指纹数据
            fingerprints = json.loads(fingerprints_json)
            self.fingerprints_cache[amino_id] = fingerprints
        
        conn.close()
        print(f"✅ 加载了 {len(self.amino_acids_data)} 种氨基酸的指纹数据")
    
    def calculate_tanimoto_similarity(self, fp1_str: str, fp2_str: str) -> float:
        """计算Tanimoto相似性"""
        try:
            # 转换为bit向量
            fp1_bits = [int(b) for b in fp1_str]
            fp2_bits = [int(b) for b in fp2_str]
            
            # 计算交集和并集
            intersection = sum(a & b for a, b in zip(fp1_bits, fp2_bits))
            union = sum(a | b for a, b in zip(fp1_bits, fp2_bits))
            
            if union == 0:
                return 0.0
            
            return intersection / union
            
        except Exception as e:
            print(f"相似性计算错误: {e}")
            return 0.0
    
    def search_similar_by_smiles(self, query_smiles: str, 
                                similarity_threshold: float = 0.6,
                                max_results: int = 10,
                                fingerprint_type: str = 'ecfp2') -> List[SimilarityResult]:
        """基于SMILES搜索相似的氨基酸"""
        
        print(f"🔍 搜索与 {query_smiles} 相似的氨基酸")
        print(f"   相似性阈值: {similarity_threshold}")
        print(f"   指纹类型: {fingerprint_type}")
        
        # 计算查询分子的指纹
        query_fingerprint = self._calculate_fingerprint(query_smiles, fingerprint_type)
        if not query_fingerprint:
            print("❌ 无法计算查询分子的指纹")
            return []
        
        results = []
        
        # 与数据库中每个氨基酸比较
        for amino_id, fingerprints in self.fingerprints_cache.items():
            if fingerprint_type not in fingerprints:
                continue
            
            target_fingerprint = fingerprints[fingerprint_type]
            similarity = self.calculate_tanimoto_similarity(query_fingerprint, target_fingerprint)
            
            if similarity >= similarity_threshold:
                amino_data = self.amino_acids_data[amino_id]
                
                result = SimilarityResult(
                    amino_acid_id=amino_id,
                    name=amino_data['name'],
                    molecular_formula=amino_data['molecular_formula'],
                    smiles=amino_data['smiles'],
                    similarity_score=similarity,
                    similarity_type=fingerprint_type
                )
                results.append(result)
        
        # 按相似性排序
        results.sort(key=lambda x: x.similarity_score, reverse=True)
        
        return results[:max_results]
    
    def search_similar_by_id(self, amino_acid_id: str,
                           similarity_threshold: float = 0.6,
                           max_results: int = 10,
                           fingerprint_type: str = 'ecfp2') -> List[SimilarityResult]:
        """基于氨基酸ID搜索相似的氨基酸"""
        
        if amino_acid_id not in self.amino_acids_data:
            print(f"❌ 未找到氨基酸: {amino_acid_id}")
            return []
        
        query_smiles = self.amino_acids_data[amino_acid_id]['smiles']
        print(f"🔍 搜索与 {amino_acid_id} ({query_smiles}) 相似的氨基酸")
        
        if amino_acid_id not in self.fingerprints_cache:
            print(f"❌ 未找到 {amino_acid_id} 的指纹数据")
            return []
        
        query_fingerprint = self.fingerprints_cache[amino_acid_id][fingerprint_type]
        results = []
        
        # 与数据库中每个氨基酸比较（排除自己）
        for target_id, fingerprints in self.fingerprints_cache.items():
            if target_id == amino_acid_id:  # 跳过自己
                continue
                
            if fingerprint_type not in fingerprints:
                continue
            
            target_fingerprint = fingerprints[fingerprint_type]
            similarity = self.calculate_tanimoto_similarity(query_fingerprint, target_fingerprint)
            
            if similarity >= similarity_threshold:
                amino_data = self.amino_acids_data[target_id]
                
                result = SimilarityResult(
                    amino_acid_id=target_id,
                    name=amino_data['name'],
                    molecular_formula=amino_data['molecular_formula'],
                    smiles=amino_data['smiles'],
                    similarity_score=similarity,
                    similarity_type=fingerprint_type
                )
                results.append(result)
        
        # 按相似性排序
        results.sort(key=lambda x: x.similarity_score, reverse=True)
        
        return results[:max_results]
    
    def _calculate_fingerprint(self, smiles: str, fingerprint_type: str) -> Optional[str]:
        """计算单个分子的指纹"""
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors, MACCSkeys
            
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            
            if fingerprint_type == 'ecfp2':
                fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            elif fingerprint_type == 'ecfp4':
                fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 3, nBits=1024)
            elif fingerprint_type == 'maccs':
                fp = MACCSkeys.GenMACCSKeys(mol)
            elif fingerprint_type == 'topological':
                fp = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=1024)
            elif fingerprint_type == 'atom_pair':
                fp = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=1024)
            else:
                return None
            
            return fp.ToBitString()
            
        except Exception as e:
            print(f"指纹计算错误: {e}")
            return None
    
    def find_most_similar_pairs(self, top_n: int = 10) -> List[Tuple[str, str, float]]:
        """找到数据库中最相似的氨基酸对"""
        
        print(f"🔍 寻找最相似的 {top_n} 对氨基酸...")
        
        similarities = []
        amino_ids = list(self.fingerprints_cache.keys())
        
        for i, id1 in enumerate(amino_ids):
            for j, id2 in enumerate(amino_ids[i+1:], i+1):
                fp1 = self.fingerprints_cache[id1]['ecfp2']
                fp2 = self.fingerprints_cache[id2]['ecfp2']
                
                similarity = self.calculate_tanimoto_similarity(fp1, fp2)
                similarities.append((id1, id2, similarity))
        
        # 按相似性排序
        similarities.sort(key=lambda x: x[2], reverse=True)
        
        return similarities[:top_n]

def test_similarity_search():
    """测试相似性搜索功能"""
    
    print("🧪 测试相似性搜索功能")
    print("=" * 60)
    
    engine = FingerprintSearchEngine()
    
    # 测试1: 基于ID搜索相似氨基酸
    print("\n📋 测试1: 寻找与苯丙氨酸(PHE)相似的氨基酸")
    results = engine.search_similar_by_id('PHE', similarity_threshold=0.5, max_results=5)
    
    for i, result in enumerate(results, 1):
        print(f"  {i}. {result.amino_acid_id}: {result.name}")
        print(f"     相似性: {result.similarity_score:.3f}")
        print(f"     分子式: {result.molecular_formula}")
    
    # 测试2: 基于SMILES搜索
    print("\n📋 测试2: 基于SMILES搜索相似结构")
    query_smiles = "N[C@@H](Cc1ccc(O)cc1)C(=O)O"  # 酪氨酸
    results = engine.search_similar_by_smiles(query_smiles, similarity_threshold=0.6, max_results=5)
    
    for i, result in enumerate(results, 1):
        print(f"  {i}. {result.amino_acid_id}: {result.name}")
        print(f"     相似性: {result.similarity_score:.3f}")
        print(f"     SMILES: {result.smiles}")
    
    # 测试3: 找到最相似的氨基酸对
    print("\n📋 测试3: 数据库中最相似的氨基酸对")
    similar_pairs = engine.find_most_similar_pairs(top_n=5)
    
    for i, (id1, id2, similarity) in enumerate(similar_pairs, 1):
        name1 = engine.amino_acids_data[id1]['name']
        name2 = engine.amino_acids_data[id2]['name']
        print(f"  {i}. {id1}({name1}) vs {id2}({name2})")
        print(f"     相似性: {similarity:.3f}")

def main():
    """主函数"""
    
    print("🚀 增强版氨基酸相似性搜索引擎")
    print("=" * 60)
    
    # 运行测试
    test_similarity_search()
    
    print(f"\n💡 使用示例:")
    print("from enhanced_search_engine import FingerprintSearchEngine")
    print("engine = FingerprintSearchEngine()")
    print("results = engine.search_similar_by_id('PHE', similarity_threshold=0.7)")

if __name__ == "__main__":
    main()
