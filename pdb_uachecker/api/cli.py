"""
命令行接口
提供用户友好的命令行工具
"""

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional

from ..analysis import PDBAnalyzer, AminoAcidClassifier
from ..core.models import VerificationMethod
from ..utils.config import Config


def create_parser() -> argparse.ArgumentParser:
    """创建命令行参数解析器"""
    parser = argparse.ArgumentParser(
        description="PDB-UAChecker: 非天然氨基酸识别系统",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 基本分析
  pdb-uachecker analyze protein.pdb
  
  # 保存报告
  pdb-uachecker analyze protein.pdb -o report.txt
  
  # JSON格式输出
  pdb-uachecker analyze protein.pdb -o results.json --format json
  
  # 禁用某些验证方法
  pdb-uachecker analyze protein.pdb --disable-3d --disable-fingerprint
  
  # 系统状态检查
  pdb-uachecker status
  
  # 数据库统计
  pdb-uachecker database --stats
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='可用命令')
    
    # 分析命令
    analyze_parser = subparsers.add_parser('analyze', help='分析PDB文件')
    analyze_parser.add_argument('pdb_file', help='PDB文件路径')
    analyze_parser.add_argument('-o', '--output', help='输出文件路径')
    analyze_parser.add_argument('--format', choices=['txt', 'json'], default='txt', help='输出格式')
    analyze_parser.add_argument('--config', help='配置文件路径')
    
    # 验证方法控制
    analyze_parser.add_argument('--disable-formula', action='store_true', help='禁用分子式验证')
    analyze_parser.add_argument('--disable-composition', action='store_true', help='禁用原子组成验证')
    analyze_parser.add_argument('--disable-fingerprint', action='store_true', help='禁用指纹相似性验证')
    analyze_parser.add_argument('--disable-3d', action='store_true', help='禁用3D结构验证')
    
    # 阈值设置
    analyze_parser.add_argument('--formula-threshold', type=float, help='分子式验证阈值')
    analyze_parser.add_argument('--composition-threshold', type=float, help='原子组成验证阈值')
    analyze_parser.add_argument('--fingerprint-threshold', type=float, help='指纹相似性阈值')
    analyze_parser.add_argument('--structure-threshold', type=float, help='3D结构验证阈值')
    
    # 系统状态命令
    status_parser = subparsers.add_parser('status', help='检查系统状态')
    status_parser.add_argument('--config', help='配置文件路径')
    
    # 数据库命令
    db_parser = subparsers.add_parser('database', help='数据库操作')
    db_parser.add_argument('--stats', action='store_true', help='显示数据库统计信息')
    db_parser.add_argument('--migrate', help='从旧版数据库迁移数据')
    db_parser.add_argument('--config', help='配置文件路径')
    
    # 分类命令
    classify_parser = subparsers.add_parser('classify', help='氨基酸分类')
    classify_parser.add_argument('--output', help='输出文件路径')
    classify_parser.add_argument('--format', choices=['txt', 'json'], default='txt', help='输出格式')
    classify_parser.add_argument('--config', help='配置文件路径')
    
    return parser


def load_config(config_file: Optional[str]) -> Config:
    """加载配置"""
    if config_file and Path(config_file).exists():
        return Config(config_file)
    else:
        from ..utils.config import default_config
        return default_config


def get_enabled_methods(args) -> List[VerificationMethod]:
    """根据命令行参数获取启用的验证方法"""
    all_methods = [
        VerificationMethod.MOLECULAR_FORMULA,
        VerificationMethod.ATOM_COMPOSITION,
        VerificationMethod.FINGERPRINT_SIMILARITY,
        VerificationMethod.STRUCTURE_3D
    ]
    
    enabled_methods = []
    
    if not getattr(args, 'disable_formula', False):
        enabled_methods.append(VerificationMethod.MOLECULAR_FORMULA)
    
    if not getattr(args, 'disable_composition', False):
        enabled_methods.append(VerificationMethod.ATOM_COMPOSITION)
    
    if not getattr(args, 'disable_fingerprint', False):
        enabled_methods.append(VerificationMethod.FINGERPRINT_SIMILARITY)
    
    if not getattr(args, 'disable_3d', False):
        enabled_methods.append(VerificationMethod.STRUCTURE_3D)
    
    return enabled_methods


def update_config_from_args(config: Config, args):
    """根据命令行参数更新配置"""
    threshold_updates = {}
    
    if hasattr(args, 'formula_threshold') and args.formula_threshold is not None:
        threshold_updates['molecular_formula'] = args.formula_threshold
    
    if hasattr(args, 'composition_threshold') and args.composition_threshold is not None:
        threshold_updates['atom_composition'] = args.composition_threshold
    
    if hasattr(args, 'fingerprint_threshold') and args.fingerprint_threshold is not None:
        threshold_updates['fingerprint_similarity'] = args.fingerprint_threshold
    
    if hasattr(args, 'structure_threshold') and args.structure_threshold is not None:
        threshold_updates['structure_3d'] = args.structure_threshold
    
    if threshold_updates:
        for key, value in threshold_updates.items():
            setattr(config.thresholds, key, value)
        print(f"✅ 已更新阈值: {threshold_updates}")


def handle_analyze_command(args):
    """处理分析命令"""
    # 检查PDB文件
    if not Path(args.pdb_file).exists():
        print(f"❌ PDB文件不存在: {args.pdb_file}")
        return 1
    
    # 加载配置
    config = load_config(getattr(args, 'config', None))
    update_config_from_args(config, args)
    
    # 获取启用的验证方法
    enabled_methods = get_enabled_methods(args)
    
    if not enabled_methods:
        print("❌ 至少需要启用一种验证方法")
        return 1
    
    print(f"🔧 启用的验证方法: {[m.value for m in enabled_methods]}")
    
    try:
        # 初始化分析器
        analyzer = PDBAnalyzer(config)
        
        # 执行分析
        result = analyzer.analyze_pdb(
            args.pdb_file,
            enabled_methods=enabled_methods,
            save_report=bool(args.output),
            output_file=args.output
        )
        
        # 输出结果
        if args.output:
            save_results(result, args.output, args.format)
        
        return 0
    
    except Exception as e:
        print(f"❌ 分析失败: {e}")
        return 1


def handle_status_command(args):
    """处理状态检查命令"""
    config = load_config(getattr(args, 'config', None))
    
    try:
        analyzer = PDBAnalyzer(config)
        validation = analyzer.validate_system()
        
        print("🔍 系统状态检查")
        print("=" * 50)
        print(f"整体状态: {validation['overall_status']}")
        
        print("\n组件状态:")
        for component, status in validation['components'].items():
            if isinstance(status, dict):
                if status.get('status') == 'healthy':
                    print(f"  ✅ {component}: 正常")
                    if 'amino_acids_count' in status:
                        print(f"     氨基酸数量: {status['amino_acids_count']}")
                else:
                    print(f"  ❌ {component}: 错误")
                    if 'error' in status:
                        print(f"     错误: {status['error']}")
        
        if validation['warnings']:
            print("\n⚠️ 警告:")
            for warning in validation['warnings']:
                print(f"  - {warning}")
        
        if validation['errors']:
            print("\n❌ 错误:")
            for error in validation['errors']:
                print(f"  - {error}")
            return 1
        
        return 0
    
    except Exception as e:
        print(f"❌ 状态检查失败: {e}")
        return 1


def handle_database_command(args):
    """处理数据库命令"""
    config = load_config(getattr(args, 'config', None))
    
    try:
        analyzer = PDBAnalyzer(config)
        
        if args.stats:
            stats = analyzer.get_database_stats()
            print("📊 数据库统计信息")
            print("=" * 30)
            for key, value in stats.items():
                print(f"{key.replace('_', ' ').title()}: {value}")
        
        if args.migrate:
            success = analyzer.database.migrate_from_legacy_database(args.migrate)
            if success:
                print("✅ 数据迁移完成")
                return 0
            else:
                print("❌ 数据迁移失败")
                return 1
        
        return 0
    
    except Exception as e:
        print(f"❌ 数据库操作失败: {e}")
        return 1


def handle_classify_command(args):
    """处理分类命令"""
    config = load_config(getattr(args, 'config', None))
    
    try:
        analyzer = PDBAnalyzer(config)
        classifier = AminoAcidClassifier()
        
        # 获取所有氨基酸
        amino_acids = analyzer.database.get_all_amino_acids()
        
        # 执行分类
        classification_results = classifier.batch_classify(amino_acids)
        
        # 生成摘要
        summary = classifier.get_classification_summary(classification_results)
        print(summary)
        
        # 保存结果
        if args.output:
            save_classification_results(classification_results, args.output, args.format)
        
        return 0
    
    except Exception as e:
        print(f"❌ 分类失败: {e}")
        return 1


def save_results(result, output_file: str, format_type: str):
    """保存分析结果"""
    try:
        if format_type == 'json':
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
            
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(serializable_result, f, indent=2, ensure_ascii=False)
        
        print(f"📄 结果已保存到: {output_file}")
    
    except Exception as e:
        print(f"⚠️ 结果保存失败: {e}")


def save_classification_results(results, output_file: str, format_type: str):
    """保存分类结果"""
    try:
        if format_type == 'json':
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(results, f, indent=2, ensure_ascii=False)
        else:
            summary = AminoAcidClassifier().get_classification_summary(results)
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(summary)
        
        print(f"📄 分类结果已保存到: {output_file}")
    
    except Exception as e:
        print(f"⚠️ 分类结果保存失败: {e}")


def main():
    """主函数"""
    parser = create_parser()
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 1
    
    # 根据命令分发处理
    if args.command == 'analyze':
        return handle_analyze_command(args)
    elif args.command == 'status':
        return handle_status_command(args)
    elif args.command == 'database':
        return handle_database_command(args)
    elif args.command == 'classify':
        return handle_classify_command(args)
    else:
        print(f"❌ 未知命令: {args.command}")
        return 1


if __name__ == '__main__':
    sys.exit(main())