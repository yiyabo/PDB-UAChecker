#!/usr/bin/env python3
"""
可扩展搜索引擎命令行接口
保持与原有系统的兼容性，同时提供新功能
"""

import os
import json
import argparse
import time
from typing import List, Dict, Any

from scalable_search_engine import ScalableSearchEngine, CompatibilityWrapper

class SearchResultAnalyzer:
    """搜索结果分析器（简化版）"""
    
    def __init__(self, results: List[Dict[str, Any]]):
        self.results = results
    
    def generate_summary_report(self) -> str:
        """生成摘要报告"""
        if not self.results:
            return "未找到任何匹配结果。"
        
        report = []
        report.append("=" * 80)
        report.append("非天然氨基酸搜索结果摘要")
        report.append("=" * 80)
        
        # 总体统计
        total_matches = len(self.results)
        unique_pdbs = len(set(r.get('pdb_id', 'unknown') for r in self.results))
        unique_amino_acids = len(set(r.get('amino_acid_id', 'unknown') for r in self.results))
        
        report.append(f"总匹配数: {total_matches}")
        report.append(f"涉及PDB文件: {unique_pdbs}")
        report.append(f"发现的氨基酸类型: {unique_amino_acids}")
        report.append("")
        
        # 按氨基酸类型统计
        from collections import Counter
        amino_acid_counts = Counter(r.get('amino_acid_id', 'unknown') for r in self.results)
        report.append("按氨基酸类型统计:")
        for amino_id, count in amino_acid_counts.most_common():
            report.append(f"  {amino_id}: {count} 个匹配")
        report.append("")
        
        # 按匹配方法统计
        method_counts = Counter(r.get('match_method', 'unknown') for r in self.results)
        report.append("按匹配方法统计:")
        for method, count in method_counts.most_common():
            report.append(f"  {method}: {count} 个匹配")
        report.append("")
        
        # 置信度分布
        high_conf = len([r for r in self.results if r.get('confidence_score', 0) >= 0.9])
        med_conf = len([r for r in self.results if 0.7 <= r.get('confidence_score', 0) < 0.9])
        low_conf = len([r for r in self.results if r.get('confidence_score', 0) < 0.7])
        
        report.append("置信度分布:")
        report.append(f"  高置信度 (≥0.9): {high_conf}")
        report.append(f"  中等置信度 (0.7-0.9): {med_conf}")
        report.append(f"  低置信度 (<0.7): {low_conf}")
        report.append("")
        
        # 详细结果
        report.append("详细匹配结果:")
        report.append("-" * 80)
        
        for result in sorted(self.results, key=lambda x: x.get('confidence_score', 0), reverse=True):
            pdb_id = result.get('pdb_id', 'unknown')
            amino_acid_id = result.get('amino_acid_id', 'unknown')
            residue_name = result.get('residue_name', 'unknown')
            chain_id = result.get('chain_id', '')
            residue_number = result.get('residue_number', '')
            match_method = result.get('match_method', 'unknown')
            confidence = result.get('confidence_score', 0)
            atom_composition = result.get('atom_composition', {})
            
            report.append(f"PDB: {pdb_id} | 氨基酸: {amino_acid_id} | "
                         f"残基: {residue_name} | 链: {chain_id} | "
                         f"位置: {residue_number}")
            report.append(f"  方法: {match_method} | 置信度: {confidence:.3f}")
            report.append(f"  原子组成: {atom_composition}")
            report.append("")
        
        return "\n".join(report)
    
    def export_to_json(self, filename: str):
        """导出结果为JSON格式"""
        data = {
            "search_timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            "total_matches": len(self.results),
            "results": self.results
        }
        
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        
        print(f"结果已导出到: {filename}")

def search_pdb_files(args):
    """搜索PDB文件"""
    print("使用可扩展搜索引擎搜索PDB文件...")
    
    # 使用兼容性包装器保持向后兼容
    wrapper = CompatibilityWrapper()
    
    pdb_files = []
    
    # 处理输入文件
    if args.input:
        if os.path.isfile(args.input) and args.input.endswith('.pdb'):
            pdb_files.append(args.input)
        elif os.path.isdir(args.input):
            for file in os.listdir(args.input):
                if file.endswith('.pdb'):
                    pdb_files.append(os.path.join(args.input, file))
    
    # 如果没有指定输入，搜索当前目录
    if not pdb_files:
        current_dir_pdbs = [f for f in os.listdir('.') if f.endswith('.pdb')]
        pdb_files.extend(current_dir_pdbs)
    
    if not pdb_files:
        print("未找到PDB文件进行搜索")
        return
    
    print(f"将搜索 {len(pdb_files)} 个PDB文件")
    
    # 执行搜索
    start_time = time.time()
    results = wrapper.search_pdb_files(pdb_files, args.methods)
    search_time = time.time() - start_time
    
    print(f"搜索完成，耗时 {search_time:.3f}s")
    
    # 分析结果
    analyzer = SearchResultAnalyzer(results)
    
    # 输出结果
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # 生成报告
    summary = analyzer.generate_summary_report()
    print("\n" + summary)
    
    if args.report:
        report_file = os.path.join(args.output, "search_report.txt")
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write(summary)
        print(f"详细报告已保存到: {report_file}")
    
    # 导出JSON结果
    json_file = os.path.join(args.output, "search_results.json")
    analyzer.export_to_json(json_file)
    
    # 显示最佳匹配
    if results:
        print("\n" + "=" * 50)
        print("所有匹配结果:")
        print("=" * 50)
        for i, match in enumerate(results, 1):
            print(f"{i}. PDB: {match.get('pdb_id', 'unknown')} | "
                  f"氨基酸: {match.get('amino_acid_id', 'unknown')} | "
                  f"置信度: {match.get('confidence_score', 0):.3f}")

def interactive_search(args):
    """交互式搜索"""
    print("可扩展搜索引擎 - 交互式模式")
    print("=" * 50)
    
    engine = ScalableSearchEngine()
    
    # 显示系统信息
    stats = engine.get_system_statistics()
    print(f"数据库包含 {stats['database_stats']['total_amino_acids']} 种氨基酸")
    print(f"支持的搜索方法: {', '.join(stats['supported_search_methods'])}")
    print()
    
    while True:
        print("请选择搜索方式:")
        print("1. 残基名搜索")
        print("2. 分子式搜索")
        print("3. 原子组成搜索")
        print("4. 分子量搜索")
        print("5. 组合搜索")
        print("6. 系统统计")
        print("0. 退出")
        
        choice = input("\n请输入选择 (0-6): ").strip()
        
        if choice == '0':
            print("退出交互式搜索")
            break
        elif choice == '1':
            residue_name = input("请输入残基名: ").strip()
            if residue_name:
                results = engine.search({'residue_name': residue_name})
                _display_search_results(results)
        elif choice == '2':
            formula = input("请输入分子式 (如 C10H13NO3): ").strip()
            if formula:
                results = engine.search({'molecular_formula': formula})
                _display_search_results(results)
        elif choice == '3':
            print("请输入原子组成 (格式: C:10,H:13,N:1,O:3)")
            comp_str = input("原子组成: ").strip()
            if comp_str:
                try:
                    atom_comp = {}
                    for pair in comp_str.split(','):
                        atom, count = pair.split(':')
                        atom_comp[atom.strip()] = int(count.strip())
                    results = engine.search({'atom_composition': atom_comp})
                    _display_search_results(results)
                except ValueError:
                    print("格式错误，请使用 C:10,H:13,N:1,O:3 格式")
        elif choice == '4':
            try:
                weight = float(input("请输入分子量: ").strip())
                tolerance = float(input("请输入容差 (默认5.0): ").strip() or "5.0")
                results = engine.search({
                    'molecular_weight': weight,
                    'weight_tolerance': tolerance
                })
                _display_search_results(results)
            except ValueError:
                print("请输入有效的数值")
        elif choice == '5':
            print("组合搜索 - 可以输入多个条件")
            query = {}
            
            residue_name = input("残基名 (可选): ").strip()
            if residue_name:
                query['residue_name'] = residue_name
            
            formula = input("分子式 (可选): ").strip()
            if formula:
                query['molecular_formula'] = formula
            
            if query:
                results = engine.search(query)
                _display_search_results(results)
            else:
                print("请至少输入一个搜索条件")
        elif choice == '6':
            _display_system_stats(engine)
        else:
            print("无效选择，请重新输入")
        
        print()

def _display_search_results(results):
    """显示搜索结果"""
    if not results:
        print("未找到匹配结果")
        return
    
    print(f"\n找到 {len(results)} 个结果:")
    print("-" * 60)
    
    for i, result in enumerate(results, 1):
        record = result.amino_acid_record
        print(f"{i}. {result.amino_acid_id} - {record.name}")
        print(f"   分子式: {record.molecular_formula}")
        print(f"   分子量: {record.molecular_weight}")
        print(f"   匹配方法: {result.match_method}")
        print(f"   置信度: {result.confidence_score:.3f}")
        if result.additional_info:
            print(f"   附加信息: {result.additional_info}")
        print()

def _display_system_stats(engine):
    """显示系统统计信息"""
    stats = engine.get_system_statistics()
    
    print("\n系统统计信息:")
    print("-" * 40)
    print(f"数据库统计:")
    db_stats = stats['database_stats']
    print(f"  总氨基酸数: {db_stats['total_amino_acids']}")
    print(f"  平均分子量: {db_stats['average_molecular_weight']}")
    print(f"  分子量范围: {db_stats['molecular_weight_range']}")
    print(f"  缓存大小: {db_stats['cache_size']}")
    
    print(f"\n索引统计:")
    idx_stats = stats['index_stats']
    print(f"  残基名索引: {idx_stats['residue_name_index_size']}")
    print(f"  分子式索引: {idx_stats['molecular_formula_index_size']}")
    print(f"  分子量索引: {idx_stats['molecular_weight_index_size']}")
    print(f"  特征索引: {idx_stats['feature_index_size']}")

def main():
    """主程序"""
    parser = argparse.ArgumentParser(description="可扩展非天然氨基酸PDB搜索引擎")
    
    subparsers = parser.add_subparsers(dest='command', help='可用命令')
    
    # PDB文件搜索命令（向后兼容）
    search_parser = subparsers.add_parser('search', help='搜索PDB文件')
    search_parser.add_argument("--input", "-i", help="输入PDB文件或目录")
    search_parser.add_argument("--output", "-o", default="search_results", help="输出目录")
    search_parser.add_argument("--methods", "-m", nargs="+", 
                              default=["residue_name", "molecular_formula", "atom_composition"],
                              help="搜索方法")
    search_parser.add_argument("--report", "-r", action="store_true", help="生成详细报告")
    
    # 交互式搜索命令
    interactive_parser = subparsers.add_parser('interactive', help='交互式搜索')
    
    # 演示命令
    demo_parser = subparsers.add_parser('demo', help='运行演示')
    
    args = parser.parse_args()
    
    if args.command == 'search':
        search_pdb_files(args)
    elif args.command == 'interactive':
        interactive_search(args)
    elif args.command == 'demo':
        # 运行演示脚本
        import scalable_search_demo
        scalable_search_demo.main()
    else:
        # 默认行为：如果没有子命令，执行PDB搜索（向后兼容）
        search_pdb_files(args)

if __name__ == "__main__":
    main()
