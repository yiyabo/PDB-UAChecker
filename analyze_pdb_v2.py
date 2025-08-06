#!/usr/bin/env python3
"""
PDB-UAChecker v2.0 - 重构版主入口
非天然氨基酸识别系统，基于四重验证策略
"""

import sys
import argparse
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker import PDBAnalyzer, Config
from pdb_uachecker.core.models import VerificationMethod
from pdb_uachecker.api.cli import main as cli_main


def create_simple_parser():
    """创建简化的命令行解析器（兼容旧版本）"""
    parser = argparse.ArgumentParser(description="PDB-UAChecker v2.0: 非天然氨基酸识别系统")
    parser.add_argument("pdb_file", help="PDB文件路径")
    parser.add_argument("--output", "-o", help="输出报告文件")
    parser.add_argument("--format", choices=["json", "txt"], default="txt", help="输出格式")
    parser.add_argument("--config", help="配置文件路径")
    
    # 验证方法控制
    parser.add_argument("--disable-formula", action="store_true", help="禁用分子式验证")
    parser.add_argument("--disable-composition", action="store_true", help="禁用原子组成验证")
    parser.add_argument("--disable-fingerprint", action="store_true", help="禁用指纹相似性验证")
    parser.add_argument("--disable-3d", action="store_true", help="禁用3D结构验证")
    
    # 阈值设置
    parser.add_argument("--formula-threshold", type=float, help="分子式验证阈值 (默认: 1.0)")
    parser.add_argument("--composition-threshold", type=float, help="原子组成验证阈值 (默认: 0.9)")
    parser.add_argument("--fingerprint-threshold", type=float, help="指纹相似性阈值 (默认: 0.7)")
    parser.add_argument("--structure-threshold", type=float, help="3D结构验证阈值 (默认: 0.5)")
    
    # 高级选项
    parser.add_argument("--parallel", action="store_true", default=True, help="启用并行处理")
    parser.add_argument("--no-parallel", action="store_false", dest="parallel", help="禁用并行处理")
    parser.add_argument("--max-workers", type=int, help="最大并行工作线程数")
    
    return parser


def load_config_from_args(args):
    """从命令行参数加载配置"""
    # 加载基础配置
    if args.config and Path(args.config).exists():
        config = Config(args.config)
    else:
        config = Config()
    
    # 更新阈值
    if args.formula_threshold is not None:
        config.thresholds.molecular_formula = args.formula_threshold
    if args.composition_threshold is not None:
        config.thresholds.atom_composition = args.composition_threshold
    if args.fingerprint_threshold is not None:
        config.thresholds.fingerprint_similarity = args.fingerprint_threshold
    if args.structure_threshold is not None:
        config.thresholds.structure_3d = args.structure_threshold
    
    # 更新性能配置
    config.performance.enable_parallel = args.parallel
    if args.max_workers is not None:
        config.performance.max_workers = args.max_workers
    
    return config


def get_enabled_methods_from_args(args):
    """从命令行参数获取启用的验证方法"""
    enabled_methods = []
    
    if not args.disable_formula:
        enabled_methods.append(VerificationMethod.MOLECULAR_FORMULA)
    
    if not args.disable_composition:
        enabled_methods.append(VerificationMethod.ATOM_COMPOSITION)
    
    if not args.disable_fingerprint:
        enabled_methods.append(VerificationMethod.FINGERPRINT_SIMILARITY)
    
    if not args.disable_3d:
        enabled_methods.append(VerificationMethod.STRUCTURE_3D)
    
    return enabled_methods


def print_welcome():
    """打印欢迎信息"""
    print("🌟 PDB-UAChecker v2.0 - 非天然氨基酸识别系统")
    print("=" * 60)
    print("🎯 基于四重验证策略的高精度氨基酸识别")
    print("🚀 重构版本 - 模块化架构，性能优化")
    print()


def main():
    """主函数"""
    # 如果没有参数，显示帮助
    if len(sys.argv) == 1:
        print_welcome()
        print("💡 使用方法:")
        print("  python analyze_pdb_v2.py protein.pdb")
        print("  python analyze_pdb_v2.py protein.pdb --output report.txt")
        print("  python analyze_pdb_v2.py protein.pdb --disable-3d --fingerprint-threshold 0.8")
        print()
        print("🔧 高级功能:")
        print("  python -m pdb_uachecker.api.cli analyze protein.pdb")
        print("  python -m pdb_uachecker.api.cli status")
        print("  python -m pdb_uachecker.api.cli database --stats")
        print()
        return
    
    # 解析命令行参数
    parser = create_simple_parser()
    args = parser.parse_args()
    
    print_welcome()
    
    try:
        # 检查PDB文件
        if not Path(args.pdb_file).exists():
            print(f"❌ PDB文件不存在: {args.pdb_file}")
            return 1
        
        # 加载配置
        config = load_config_from_args(args)
        
        # 获取启用的验证方法
        enabled_methods = get_enabled_methods_from_args(args)
        
        if not enabled_methods:
            print("❌ 至少需要启用一种验证方法")
            return 1
        
        print(f"🔧 启用的验证方法: {[m.value for m in enabled_methods]}")
        print(f"📊 验证阈值: 分子式={config.thresholds.molecular_formula}, "
              f"原子组成={config.thresholds.atom_composition}, "
              f"指纹相似性={config.thresholds.fingerprint_similarity}, "
              f"3D结构={config.thresholds.structure_3d}")
        
        # 初始化分析器
        analyzer = PDBAnalyzer(config)
        
        # 执行分析
        result = analyzer.analyze_pdb(
            args.pdb_file,
            enabled_methods=enabled_methods,
            save_report=bool(args.output),
            output_file=args.output
        )
        
        # 输出JSON格式结果（如果需要）
        if args.output and args.format == "json":
            import json
            
            # 转换为可序列化的格式
            serializable_result = {
                'pdb_file': result.pdb_file,
                'total_residues': result.total_residues,
                'identified_residues': result.identified_residues,
                'identification_rate': result.identification_rate,
                'analysis_time': result.analysis_time,
                'capabilities': result.capabilities,
                'matches': [
                    {
                        'residue_key': match.residue_info.residue_key,
                        'amino_acid_id': match.amino_acid_id,
                        'amino_acid_name': match.amino_acid_name,
                        'confidence_score': match.confidence_score,
                        'match_method': match.match_method,
                        'verification_summary': match.verification_result.verification_summary
                    }
                    for match in result.matches
                ]
            }
            
            json_file = args.output.replace('.txt', '.json') if args.output.endswith('.txt') else args.output
            with open(json_file, 'w', encoding='utf-8') as f:
                json.dump(serializable_result, f, indent=2, ensure_ascii=False)
            print(f"📄 JSON结果已保存到: {json_file}")
        
        print(f"\n🎉 分析完成！识别率: {result.identification_rate:.1f}%")
        
        return 0
    
    except KeyboardInterrupt:
        print("\n⚠️ 用户中断操作")
        return 1
    except Exception as e:
        print(f"❌ 分析失败: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())