#!/usr/bin/env python3
"""
优化的PDB分析器
集成多重验证并行搜索策略
"""

import sys
import time
import argparse
from pathlib import Path
from typing import List, Dict, Optional

# 添加路径
sys.path.insert(0, 'core')
sys.path.insert(0, 'legacy')

from unified_pdb_analyzer import PDBParser, ResidueInfo
from parallel_search_engine import ParallelSearchEngine, ParallelMatchResult

class OptimizedPDBAnalyzer:
    """优化的PDB分析器"""
    
    def __init__(self, search_strategy: str = "parallel"):
        """
        初始化分析器
        
        Args:
            search_strategy: 搜索策略 ("progressive" 或 "parallel")
        """
        print(f"🚀 初始化优化PDB分析器...")
        
        self.parser = PDBParser()
        self.search_strategy = search_strategy
        
        if search_strategy == "parallel":
            self.search_engine = ParallelSearchEngine()
            print("✅ 使用多重验证并行搜索策略")
        else:
            # 保持向后兼容性
            from unified_pdb_analyzer import UnifiedSearchEngine
            self.search_engine = UnifiedSearchEngine()
            print("⚠️ 使用传统渐进式搜索策略")
        
        print("✅ 优化PDB分析器就绪！")
    
    def analyze_pdb(self, pdb_file: str, 
                   enable_fingerprint: bool = True,
                   enable_3d: bool = True,
                   save_report: bool = False,
                   output_file: Optional[str] = None) -> Dict:
        """分析PDB文件"""
        
        start_time = time.time()
        
        print(f"\n🔬 开始PDB分析: {pdb_file}")
        print(f"📊 搜索策略: {self.search_strategy}")
        print("=" * 70)
        
        # 1. 解析PDB文件
        residues = self.parser.parse_pdb_file(pdb_file)
        if not residues:
            return {'error': 'PDB文件解析失败'}
        
        print(f"📄 解析完成，发现 {len(residues)} 个残基")
        
        # 2. 分析每个残基
        analysis_results = []
        high_confidence_matches = []
        
        for residue in residues:
            print(f"\n🧪 分析残基: {residue.chain_id}:{residue.residue_name}{residue.residue_number}")
            
            if self.search_strategy == "parallel":
                # 使用并行验证策略
                matches = self.search_engine.parallel_search(
                    residue, enable_fingerprint, enable_3d
                )
            else:
                # 使用传统策略
                matches = self._traditional_search(residue)
            
            if matches:
                best_match = matches[0]
                analysis_results.append(best_match)
                
                print(f"   ✅ 匹配: {best_match.amino_acid_name}")
                print(f"   📊 置信度: {best_match.confidence_score:.3f}")
                print(f"   🔧 方法: {best_match.match_method}")
                
                # 显示验证详情（如果是并行策略）
                if hasattr(best_match, 'verification_details') and best_match.verification_details:
                    details = best_match.verification_details
                    if not details.get('bypass_parallel', False):
                        print(f"   🔍 验证详情:")
                        print(f"     分子式: {'✅' if details.get('formula_match') else '❌'}")
                        print(f"     原子组成: {'✅' if details.get('atom_match') else '❌'}")
                        if enable_fingerprint:
                            print(f"     指纹相似性: {'✅' if details.get('fingerprint_match') else '❌'}")
                        if enable_3d:
                            print(f"     3D结构: {'✅' if details.get('structure_match') else '❌'}")
                
                # 收集高置信度匹配
                if best_match.confidence_score >= 0.8:
                    high_confidence_matches.append(best_match)
            else:
                print(f"   ❌ 未找到匹配")
        
        analysis_time = time.time() - start_time
        
        # 3. 生成分析报告
        report = self._generate_report(
            pdb_file, residues, analysis_results, high_confidence_matches,
            analysis_time, enable_fingerprint, enable_3d
        )
        
        # 4. 保存报告（如果需要）
        if save_report and output_file:
            self._save_report(report, output_file)
        
        self._print_summary(report)
        return report
    
    def _traditional_search(self, residue: ResidueInfo) -> List:
        """传统搜索方法（向后兼容）"""
        all_matches = []
        
        # 残基名匹配
        matches = self.search_engine.search_by_residue_name(residue.residue_name)
        all_matches.extend(matches)
        
        # 如果没有精确匹配，尝试其他方法
        if not matches:
            if residue.molecular_formula:
                matches = self.search_engine.search_by_molecular_formula(residue.molecular_formula)
                all_matches.extend(matches)
            
            if residue.atom_composition:
                matches = self.search_engine.search_by_atom_composition(residue.atom_composition)
                all_matches.extend(matches)
        
        # 去重并排序
        unique_matches = {}
        for match in all_matches:
            if match.amino_acid_id not in unique_matches:
                unique_matches[match.amino_acid_id] = match
            elif match.confidence_score > unique_matches[match.amino_acid_id].confidence_score:
                unique_matches[match.amino_acid_id] = match
        
        sorted_matches = list(unique_matches.values())
        sorted_matches.sort(key=lambda x: x.confidence_score, reverse=True)
        
        return sorted_matches
    
    def _generate_report(self, pdb_file: str, residues: List[ResidueInfo],
                        results: List, high_confidence: List,
                        analysis_time: float, enable_fingerprint: bool,
                        enable_3d: bool) -> Dict:
        """生成分析报告"""
        
        # 统计不同类型的匹配
        exact_matches = [r for r in results if r.confidence_score >= 0.95]
        parallel_matches = [r for r in results if hasattr(r, 'verification_details') 
                           and r.verification_details and not r.verification_details.get('bypass_parallel', False)]
        
        report = {
            'pdb_file': pdb_file,
            'analysis_time': analysis_time,
            'search_strategy': self.search_strategy,
            'total_residues': len(residues),
            'identified_residues': len(results),
            'identification_rate': len(results) / len(residues) * 100 if residues else 0,
            'results': results,
            'statistics': {
                'exact_matches': len(exact_matches),
                'parallel_verified_matches': len(parallel_matches),
                'high_confidence_matches': len(high_confidence),
            },
            'capabilities': {
                'search_strategy': self.search_strategy,
                'fingerprint_similarity': enable_fingerprint,
                '3d_structure_verification': enable_3d,
                'parallel_verification': self.search_strategy == "parallel"
            }
        }
        
        return report
    
    def _save_report(self, report: Dict, output_file: str):
        """保存分析报告"""
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("PDB分析报告 - 优化版\n")
            f.write("=" * 50 + "\n")
            f.write(f"文件: {report['pdb_file']}\n")
            f.write(f"搜索策略: {report['search_strategy']}\n")
            f.write(f"分析时间: {report['analysis_time']:.2f}s\n")
            f.write(f"识别率: {report['identification_rate']:.1f}%\n")
            f.write(f"识别的氨基酸数量: {report['identified_residues']}\n\n")
            
            # 统计信息
            stats = report['statistics']
            f.write("匹配统计:\n")
            f.write(f"  精确匹配: {stats['exact_matches']}\n")
            if self.search_strategy == "parallel":
                f.write(f"  并行验证匹配: {stats['parallel_verified_matches']}\n")
            f.write(f"  高置信度匹配: {stats['high_confidence_matches']}\n\n")
            
            # 详细结果
            f.write("详细结果:\n")
            for result in report['results']:
                residue = result.residue_info
                if residue:  # 检查residue_info是否存在
                    f.write(f"  {residue.chain_id}:{residue.residue_name}{residue.residue_number} -> ")
                else:
                    f.write(f"  {result.amino_acid_id} -> ")
                f.write(f"{result.amino_acid_name} (置信度: {result.confidence_score:.3f})\n")

                if hasattr(result, 'verification_details') and result.verification_details:
                    details = result.verification_details
                    if not details.get('bypass_parallel', False):
                        f.write(f"    验证: 分子式{'✅' if details.get('formula_match') else '❌'} ")
                        f.write(f"原子组成{'✅' if details.get('atom_match') else '❌'} ")
                        f.write(f"指纹{'✅' if details.get('fingerprint_match') else '❌'} ")
                        f.write(f"3D结构{'✅' if details.get('structure_match') else '❌'}\n")
            
            # 系统能力
            caps = report['capabilities']
            f.write(f"\n系统能力:\n")
            f.write(f"  搜索策略: {caps['search_strategy']}\n")
            f.write(f"  指纹相似性: {'✅' if caps['fingerprint_similarity'] else '❌'}\n")
            f.write(f"  3D结构验证: {'✅' if caps['3d_structure_verification'] else '❌'}\n")
            f.write(f"  并行验证: {'✅' if caps['parallel_verification'] else '❌'}\n")
        
        print(f"📄 报告已保存到: {output_file}")
    
    def _print_summary(self, report: Dict):
        """打印分析摘要"""
        print(f"\n" + "=" * 70)
        print(f"📊 PDB分析摘要 ({report['search_strategy']} 策略):")
        print(f"   文件: {report['pdb_file']}")
        print(f"   分析时间: {report['analysis_time']:.2f}s")
        print(f"   总残基数: {report['total_residues']}")
        print(f"   识别残基数: {report['identified_residues']}")
        print(f"   识别率: {report['identification_rate']:.1f}%")
        
        stats = report['statistics']
        print(f"\n📈 匹配统计:")
        print(f"   精确匹配: {stats['exact_matches']}")
        if self.search_strategy == "parallel":
            print(f"   并行验证匹配: {stats['parallel_verified_matches']}")
        print(f"   高置信度匹配: {stats['high_confidence_matches']}")
        
        capabilities = report['capabilities']
        print(f"\n🔧 系统能力:")
        print(f"   搜索策略: {capabilities['search_strategy']}")
        print(f"   指纹相似性: {'✅' if capabilities['fingerprint_similarity'] else '❌'}")
        print(f"   3D结构验证: {'✅' if capabilities['3d_structure_verification'] else '❌'}")
        print(f"   并行验证: {'✅' if capabilities['parallel_verification'] else '❌'}")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description="优化的非天然氨基酸识别系统")
    parser.add_argument("pdb_file", help="PDB文件路径")
    parser.add_argument("--output", "-o", help="输出报告文件")
    parser.add_argument("--strategy", choices=["progressive", "parallel"], 
                       default="parallel", help="搜索策略 (默认: parallel)")
    parser.add_argument("--no-fingerprint", action="store_true", 
                       help="禁用指纹相似性匹配")
    parser.add_argument("--no-3d", action="store_true", 
                       help="禁用3D结构验证")
    
    args = parser.parse_args()
    
    # 检查PDB文件
    if not Path(args.pdb_file).exists():
        print(f"❌ PDB文件不存在: {args.pdb_file}")
        return
    
    # 初始化分析器
    analyzer = OptimizedPDBAnalyzer(search_strategy=args.strategy)
    
    # 分析PDB文件
    results = analyzer.analyze_pdb(
        args.pdb_file,
        enable_fingerprint=not args.no_fingerprint,
        enable_3d=not args.no_3d,
        save_report=bool(args.output),
        output_file=args.output
    )
    
    print(f"\n🎉 分析完成！")

if __name__ == "__main__":
    main()
