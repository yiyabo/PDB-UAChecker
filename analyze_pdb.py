#!/usr/bin/env python3
"""
非天然氨基酸识别系统 - 主入口
输入PDB文件，识别其中的非天然氨基酸
支持精确匹配、指纹相似性、立体化学识别
"""

import sys
from pathlib import Path

# 添加路径
sys.path.insert(0, str(Path(__file__).parent / "core"))
sys.path.insert(0, str(Path(__file__).parent / "legacy"))

# 尝试导入最高级的分析器
try:
    from pdb_analyzer_with_3d import PDBAnalyzerWith3D as PDBAnalyzer
    ANALYSIS_MODE = "3D_ENHANCED"
    print("🚀 启用3D增强模式: 精确匹配 + 指纹相似性 + 立体化学识别 + 3D结构匹配")
except ImportError:
    try:
        from enhanced_pdb_analyzer import EnhancedPDBAnalyzer as PDBAnalyzer
        ANALYSIS_MODE = "ENHANCED"
        print("🚀 启用增强模式: 精确匹配 + 指纹相似性 + 立体化学识别")
    except ImportError:
        from unified_pdb_analyzer import PDBAnalyzer
        ANALYSIS_MODE = "BASIC"
        print("⚠️ 基础模式: 仅精确匹配功能")

def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(description="非天然氨基酸识别系统")
    parser.add_argument("pdb_file", help="PDB文件路径")
    parser.add_argument("--output", "-o", help="输出报告文件")
    parser.add_argument("--format", choices=["json", "txt"], default="txt", help="输出格式")

    # 高级模式选项
    if ANALYSIS_MODE in ["ENHANCED", "3D_ENHANCED"]:
        parser.add_argument("--fingerprint", action="store_true", default=True,
                          help="启用指纹相似性搜索 (默认启用)")
        parser.add_argument("--threshold", type=float, default=0.6,
                          help="指纹相似性阈值 (默认0.6)")
        parser.add_argument("--basic-only", action="store_true",
                          help="仅使用基础匹配，禁用高级功能")

    # 3D增强模式选项
    if ANALYSIS_MODE == "3D_ENHANCED":
        parser.add_argument("--enable-3d", action="store_true", default=True,
                          help="启用3D结构匹配 (默认启用)")
        parser.add_argument("--geometry-threshold", type=float, default=0.7,
                          help="3D几何匹配阈值 (默认0.7)")
        parser.add_argument("--no-3d", action="store_true",
                          help="禁用3D结构匹配")

    args = parser.parse_args()

    # 检查PDB文件是否存在
    if not Path(args.pdb_file).exists():
        print(f"❌ PDB文件不存在: {args.pdb_file}")
        return

    # 初始化分析器
    analyzer = PDBAnalyzer()

    # 分析PDB文件
    if ANALYSIS_MODE == "3D_ENHANCED" and not getattr(args, 'basic_only', False):
        # 使用3D增强模式
        enable_3d = getattr(args, 'enable_3d', True) and not getattr(args, 'no_3d', False)
        results = analyzer.analyze_pdb_with_3d(
            args.pdb_file,
            enable_fingerprint=getattr(args, 'fingerprint', True),
            fingerprint_threshold=getattr(args, 'threshold', 0.6),
            enable_3d_matching=enable_3d,
            geometry_threshold=getattr(args, 'geometry_threshold', 0.7)
        )
    elif ANALYSIS_MODE == "ENHANCED" and not getattr(args, 'basic_only', False):
        # 使用增强模式
        results = analyzer.analyze_pdb_enhanced(
            args.pdb_file,
            enable_fingerprint=getattr(args, 'fingerprint', True),
            fingerprint_threshold=getattr(args, 'threshold', 0.6)
        )
    else:
        # 使用基础模式
        if hasattr(analyzer, 'analyze_pdb_enhanced'):
            results = analyzer.analyze_pdb_enhanced(args.pdb_file, enable_fingerprint=False)
        elif hasattr(analyzer, 'analyze_pdb'):
            results = analyzer.analyze_pdb(args.pdb_file)
        else:
            print("❌ 分析器不可用")
            return

    # 输出结果
    if args.output:
        save_results(results, args.output, args.format)

    print(f"\n🎉 分析完成！")

def save_results(results, output_file, format_type):
    """保存分析结果"""
    if format_type == "json":
        import json
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"📄 结果已保存到: {output_file}")
    else:
        with open(output_file, 'w') as f:
            f.write(f"PDB分析报告\n")
            f.write(f"=" * 50 + "\n")
            f.write(f"文件: {results['pdb_file']}\n")
            f.write(f"识别率: {results['identification_rate']:.1f}%\n")
            f.write(f"识别的氨基酸数量: {results['identified_residues']}\n")

            # 增强模式的额外信息
            if 'isomer_findings' in results and results['isomer_findings']:
                f.write(f"\n同分异构体发现:\n")
                for finding in results['isomer_findings']:
                    residue = finding.residue_info
                    f.write(f"  {residue.chain_id}:{residue.residue_name}{residue.residue_number} -> ")
                    f.write(f"{finding.amino_acid_name} ({finding.isomer_type})\n")

            # 3D结构匹配信息
            if 'structure_matches' in results and results['structure_matches']:
                f.write(f"\n3D结构验证通过:\n")
                for match in results['structure_matches']:
                    residue = match.residue_info
                    f.write(f"  {residue.chain_id}:{residue.residue_name}{residue.residue_number} -> ")
                    f.write(f"{match.amino_acid_name}\n")
                    f.write(f"    3D分数: {match.conformation_similarity:.3f}\n")
                    if hasattr(match, 'rmsd') and match.rmsd is not None:
                        f.write(f"    RMSD: {match.rmsd:.3f} Å\n")

            if 'capabilities' in results:
                f.write(f"\n系统能力:\n")
                caps = results['capabilities']
                f.write(f"  精确匹配: {'✅' if caps.get('exact_matching') else '❌'}\n")
                f.write(f"  指纹相似性: {'✅' if caps.get('fingerprint_similarity') else '❌'}\n")
                f.write(f"  立体化学分析: {'✅' if caps.get('stereochemistry_analysis') else '❌'}\n")
                f.write(f"  3D结构匹配: {'✅' if caps.get('3d_structure_matching') else '❌'}\n")

        print(f"📄 报告已保存到: {output_file}")

if __name__ == "__main__":
    main()
