#!/usr/bin/env python3
"""
统一PDB分析器 - 专注于PDB文件中非天然氨基酸的识别
整合所有搜索功能，提供统一的分析入口
"""

import sys
import time
import json
import sqlite3
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass

sys.path.append('.')

@dataclass
class ResidueInfo:
    """PDB残基信息"""
    residue_name: str
    residue_number: int
    chain_id: str
    atoms: List[Dict]  # 原子信息列表
    molecular_formula: Optional[str] = None
    atom_composition: Optional[Dict] = None

@dataclass
class MatchResult:
    """匹配结果"""
    amino_acid_id: str
    amino_acid_name: str
    confidence_score: float
    match_method: str
    molecular_formula: str
    smiles: str
    residue_info: ResidueInfo

class UnifiedSearchEngine:
    """统一搜索引擎 - 整合所有搜索功能"""
    
    def __init__(self):
        # 获取数据库的正确路径
        import os
        current_dir = os.path.dirname(os.path.abspath(__file__))
        self.db_path = os.path.join(current_dir, 'amino_acids.db')
        self._load_database()
    
    def _load_database(self):
        """加载数据库数据"""
        print("🔄 加载氨基酸数据库...")
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute('''
            SELECT id, name, molecular_formula, molecular_weight, smiles, 
                   atom_composition, key_features, fingerprints
            FROM amino_acids
        ''')
        
        self.amino_acids = {}
        self.fingerprints = {}
        
        for row in cursor.fetchall():
            amino_id, name, formula, weight, smiles, atom_comp, features, fingerprints = row
            
            self.amino_acids[amino_id] = {
                'name': name,
                'molecular_formula': formula,
                'molecular_weight': weight,
                'smiles': smiles,
                'atom_composition': json.loads(atom_comp) if atom_comp else {},
                'key_features': json.loads(features) if features else []
            }
            
            if fingerprints:
                self.fingerprints[amino_id] = json.loads(fingerprints)
        
        conn.close()
        print(f"✅ 加载了 {len(self.amino_acids)} 种氨基酸数据")
    
    def search_by_residue_name(self, residue_name: str) -> List[MatchResult]:
        """基于残基名精确搜索"""
        results = []
        
        if residue_name in self.amino_acids:
            aa_data = self.amino_acids[residue_name]
            result = MatchResult(
                amino_acid_id=residue_name,
                amino_acid_name=aa_data['name'],
                confidence_score=1.0,
                match_method="residue_name_exact",
                molecular_formula=aa_data['molecular_formula'],
                smiles=aa_data['smiles'],
                residue_info=None  # 将在外部设置
            )
            results.append(result)
        
        return results
    
    def search_by_molecular_formula(self, formula: str) -> List[MatchResult]:
        """基于分子式搜索 - 优化版：降低置信度，作为筛选条件"""
        results = []

        for amino_id, aa_data in self.amino_acids.items():
            if aa_data['molecular_formula'] == formula:
                result = MatchResult(
                    amino_acid_id=amino_id,
                    amino_acid_name=aa_data['name'],
                    confidence_score=0.60,  # 降低置信度：作为筛选条件而非最终判断
                    match_method="molecular_formula_screening",  # 更新方法名
                    molecular_formula=aa_data['molecular_formula'],
                    smiles=aa_data['smiles'],
                    residue_info=None
                )
                results.append(result)

        return results
    
    def search_by_atom_composition(self, atom_comp: Dict) -> List[MatchResult]:
        """基于原子组成搜索 - 优化版：降低置信度，作为辅助验证条件"""
        results = []

        for amino_id, aa_data in self.amino_acids.items():
            db_atom_comp = aa_data['atom_composition']

            if db_atom_comp == atom_comp:
                result = MatchResult(
                    amino_acid_id=amino_id,
                    amino_acid_name=aa_data['name'],
                    confidence_score=0.65,  # 降低置信度：作为辅助验证条件
                    match_method="atom_composition_verification",  # 更新方法名
                    molecular_formula=aa_data['molecular_formula'],
                    smiles=aa_data['smiles'],
                    residue_info=None
                )
                results.append(result)

        return results
    
    def search_by_fingerprint_similarity(self, target_smiles: str, 
                                       threshold: float = 0.7) -> List[MatchResult]:
        """基于指纹相似性搜索"""
        results = []
        
        # 计算目标分子的指纹
        target_fingerprint = self._calculate_fingerprint(target_smiles)
        if not target_fingerprint:
            return results
        
        for amino_id, fingerprints in self.fingerprints.items():
            if 'ecfp2' not in fingerprints:
                continue
            
            similarity = self._calculate_tanimoto_similarity(
                target_fingerprint, fingerprints['ecfp2']
            )
            
            if similarity >= threshold:
                aa_data = self.amino_acids[amino_id]
                result = MatchResult(
                    amino_acid_id=amino_id,
                    amino_acid_name=aa_data['name'],
                    confidence_score=similarity * 0.8,  # 相似性匹配置信度稍低
                    match_method="fingerprint_similarity",
                    molecular_formula=aa_data['molecular_formula'],
                    smiles=aa_data['smiles'],
                    residue_info=None
                )
                results.append(result)
        
        # 按相似性排序
        results.sort(key=lambda x: x.confidence_score, reverse=True)
        return results[:5]  # 返回前5个最相似的
    
    def _calculate_fingerprint(self, smiles: str) -> Optional[str]:
        """计算SMILES的ECFP2指纹"""
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors
            
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return None
            
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            return fp.ToBitString()
        except:
            return None
    
    def _calculate_tanimoto_similarity(self, fp1_str: str, fp2_str: str) -> float:
        """计算Tanimoto相似性"""
        try:
            fp1_bits = [int(b) for b in fp1_str]
            fp2_bits = [int(b) for b in fp2_str]
            
            intersection = sum(a & b for a, b in zip(fp1_bits, fp2_bits))
            union = sum(a | b for a, b in zip(fp1_bits, fp2_bits))
            
            return intersection / union if union > 0 else 0.0
        except:
            return 0.0

class PDBParser:
    """PDB文件解析器"""
    
    def parse_pdb_file(self, pdb_file: str) -> List[ResidueInfo]:
        """解析PDB文件，提取残基信息"""
        print(f"📄 解析PDB文件: {pdb_file}")
        
        residues = {}
        
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        # 解析原子行
                        atom_info = self._parse_atom_line(line)
                        if atom_info:
                            residue_key = (atom_info['chain_id'], 
                                         atom_info['residue_name'], 
                                         atom_info['residue_number'])
                            
                            if residue_key not in residues:
                                residues[residue_key] = ResidueInfo(
                                    residue_name=atom_info['residue_name'],
                                    residue_number=atom_info['residue_number'],
                                    chain_id=atom_info['chain_id'],
                                    atoms=[]
                                )
                            
                            residues[residue_key].atoms.append(atom_info)
        
        except Exception as e:
            print(f"❌ PDB文件解析失败: {e}")
            return []
        
        # 计算每个残基的分子式和原子组成
        residue_list = list(residues.values())
        for residue in residue_list:
            residue.atom_composition = self._calculate_atom_composition(residue.atoms)
            residue.molecular_formula = self._calculate_molecular_formula(residue.atom_composition)
        
        print(f"✅ 解析完成，发现 {len(residue_list)} 个残基")
        return residue_list
    
    def _parse_atom_line(self, line: str) -> Optional[Dict]:
        """解析PDB原子行"""
        try:
            return {
                'atom_name': line[12:16].strip(),
                'residue_name': line[17:20].strip(),
                'chain_id': line[21:22].strip(),
                'residue_number': int(line[22:26].strip()),
                'x': float(line[30:38].strip()),
                'y': float(line[38:46].strip()),
                'z': float(line[46:54].strip()),
                'element': line[76:78].strip() or line[12:14].strip()
            }
        except:
            return None
    
    def _calculate_atom_composition(self, atoms: List[Dict]) -> Dict:
        """计算原子组成"""
        composition = {}
        for atom in atoms:
            element = atom['element']
            if element:
                composition[element] = composition.get(element, 0) + 1
        return composition
    
    def _calculate_molecular_formula(self, atom_composition: Dict) -> str:
        """计算分子式"""
        if not atom_composition:
            return ""
        
        # 按标准顺序排列元素
        element_order = ['C', 'H', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Se']
        formula_parts = []
        
        for element in element_order:
            if element in atom_composition:
                count = atom_composition[element]
                if count == 1:
                    formula_parts.append(element)
                else:
                    formula_parts.append(f"{element}{count}")
        
        # 添加其他元素
        for element, count in sorted(atom_composition.items()):
            if element not in element_order:
                if count == 1:
                    formula_parts.append(element)
                else:
                    formula_parts.append(f"{element}{count}")
        
        return ''.join(formula_parts)

class PDBAnalyzer:
    """统一PDB分析器 - 主入口"""
    
    def __init__(self):
        print("🚀 初始化统一PDB分析器...")
        self.parser = PDBParser()
        self.search_engine = UnifiedSearchEngine()
        print("✅ PDB分析器就绪！")
    
    def analyze_pdb(self, pdb_file: str) -> Dict:
        """分析PDB文件中的非天然氨基酸"""
        start_time = time.time()
        
        print(f"\n🔬 开始分析PDB文件: {pdb_file}")
        print("=" * 60)
        
        # 1. 解析PDB文件
        residues = self.parser.parse_pdb_file(pdb_file)
        if not residues:
            return {'error': 'PDB文件解析失败'}
        
        # 2. 分析每个残基
        analysis_results = []
        
        for residue in residues:
            print(f"\n🧪 分析残基: {residue.chain_id}:{residue.residue_name}{residue.residue_number}")
            
            # 多策略搜索
            matches = self._search_residue(residue)
            
            if matches:
                best_match = matches[0]  # 取置信度最高的
                best_match.residue_info = residue
                analysis_results.append(best_match)
                
                print(f"   ✅ 匹配: {best_match.amino_acid_name}")
                print(f"   📊 置信度: {best_match.confidence_score:.3f}")
                print(f"   🔧 方法: {best_match.match_method}")
            else:
                print(f"   ❌ 未找到匹配")
        
        analysis_time = time.time() - start_time
        
        # 3. 生成分析报告
        report = {
            'pdb_file': pdb_file,
            'analysis_time': analysis_time,
            'total_residues': len(residues),
            'identified_residues': len(analysis_results),
            'identification_rate': len(analysis_results) / len(residues) * 100,
            'results': analysis_results
        }
        
        self._print_summary(report)
        return report
    
    def _search_residue(self, residue: ResidueInfo) -> List[MatchResult]:
        """对单个残基进行多策略搜索"""
        all_matches = []
        
        # 策略1: 残基名精确匹配
        matches = self.search_engine.search_by_residue_name(residue.residue_name)
        all_matches.extend(matches)
        
        # 策略2: 分子式匹配
        if residue.molecular_formula:
            matches = self.search_engine.search_by_molecular_formula(residue.molecular_formula)
            all_matches.extend(matches)
        
        # 策略3: 原子组成匹配
        if residue.atom_composition:
            matches = self.search_engine.search_by_atom_composition(residue.atom_composition)
            all_matches.extend(matches)
        
        # 去重并按置信度排序
        unique_matches = {}
        for match in all_matches:
            if match.amino_acid_id not in unique_matches:
                unique_matches[match.amino_acid_id] = match
            else:
                # 保留置信度更高的
                if match.confidence_score > unique_matches[match.amino_acid_id].confidence_score:
                    unique_matches[match.amino_acid_id] = match
        
        results = list(unique_matches.values())
        results.sort(key=lambda x: x.confidence_score, reverse=True)
        
        return results
    
    def _print_summary(self, report: Dict):
        """打印分析摘要"""
        print(f"\n" + "=" * 60)
        print(f"📊 PDB分析摘要:")
        print(f"   文件: {report['pdb_file']}")
        print(f"   分析时间: {report['analysis_time']:.2f}s")
        print(f"   总残基数: {report['total_residues']}")
        print(f"   识别残基数: {report['identified_residues']}")
        print(f"   识别率: {report['identification_rate']:.1f}%")
        
        if report['results']:
            print(f"\n🎯 识别的非天然氨基酸:")
            for i, result in enumerate(report['results'][:5], 1):
                residue = result.residue_info
                print(f"   {i}. {residue.chain_id}:{residue.residue_name}{residue.residue_number}")
                print(f"      → {result.amino_acid_name} ({result.confidence_score:.3f})")

def main():
    """主函数 - 演示统一PDB分析器"""
    
    print("🌟 统一PDB分析器演示")
    print("=" * 50)
    
    analyzer = PDBAnalyzer()
    
    # 这里可以测试实际的PDB文件
    print(f"\n💡 使用方法:")
    print("analyzer = PDBAnalyzer()")
    print("results = analyzer.analyze_pdb('your_protein.pdb')")
    print("\n🎯 核心功能:")
    print("✅ 统一入口 - 一个函数完成所有分析")
    print("✅ 多策略搜索 - 残基名/分子式/原子组成/指纹相似性")
    print("✅ 智能匹配 - 自动选择最佳匹配结果")
    print("✅ 详细报告 - 完整的分析结果和统计信息")

if __name__ == "__main__":
    main()
