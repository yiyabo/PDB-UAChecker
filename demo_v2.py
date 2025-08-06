#!/usr/bin/env python3
"""
PDB-UAChecker v2.0 演示脚本
展示重构后系统的主要功能
"""

import sys
import tempfile
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker import PDBAnalyzer, Config
from pdb_uachecker.core.models import VerificationMethod, AtomInfo, ResidueInfo, AminoAcidInfo
from pdb_uachecker.analysis import AminoAcidClassifier


def create_demo_pdb_file():
    """创建演示用的PDB文件"""
    pdb_content = """HEADER    DEMO PROTEIN                            01-JAN-24   DEMO            
ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.101  14.618  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.664  16.849  14.897  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.764  18.067  15.086  1.00 20.00           O  
ATOM      5  CB  ALA A   1      18.756  15.178  13.425  1.00 20.00           C  
ATOM      6  N   GLY A   2      16.498  16.189  14.932  1.00 20.00           N  
ATOM      7  CA  GLY A   2      15.168  16.759  15.200  1.00 20.00           C  
ATOM      8  C   GLY A   2      14.021  15.759  15.394  1.00 20.00           C  
ATOM      9  O   GLY A   2      13.956  14.668  14.825  1.00 20.00           O  
ATOM     10  N   VAL A   3      13.067  16.032  16.284  1.00 20.00           N  
ATOM     11  CA  VAL A   3      11.890  15.174  16.553  1.00 20.00           C  
ATOM     12  C   VAL A   3      10.594  15.973  16.789  1.00 20.00           C  
ATOM     13  O   VAL A   3      10.648  17.201  16.895  1.00 20.00           O  
ATOM     14  CB  VAL A   3      11.651  14.217  15.378  1.00 20.00           C  
ATOM     15  CG1 VAL A   3      10.434  13.334  15.598  1.00 20.00           C  
ATOM     16  CG2 VAL A   3      12.890  13.364  15.178  1.00 20.00           C  
END
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
        f.write(pdb_content)
        return f.name


def demo_basic_analysis():
    """演示基础分析功能"""
    print("🔬 演示1: 基础PDB分析")
    print("-" * 40)
    
    # 创建演示PDB文件
    pdb_file = create_demo_pdb_file()
    print(f"📁 创建演示PDB文件: {pdb_file}")
    
    try:
        # 初始化分析器
        config = Config()
        analyzer = PDBAnalyzer(config)
        
        # 执行分析
        result = analyzer.analyze_pdb(pdb_file)
        
        # 显示结果
        print(f"✅ 分析完成:")
        print(f"   总残基数: {result.total_residues}")
        print(f"   识别残基数: {result.identified_residues}")
        print(f"   识别率: {result.identification_rate:.1f}%")
        print(f"   分析时间: {result.analysis_time:.3f}s")
        
        if result.matches:
            print(f"   识别结果:")
            for match in result.matches:
                print(f"     {match.residue_info.residue_key} -> {match.amino_acid_name} (置信度: {match.confidence_score:.3f})")
    
    finally:
        # 清理临时文件
        import os
        os.unlink(pdb_file)


def demo_custom_config():
    """演示自定义配置"""
    print("\n🔧 演示2: 自定义配置")
    print("-" * 40)
    
    # 创建自定义配置
    config = Config()
    config.thresholds.fingerprint_similarity = 0.8
    config.performance.enable_parallel = False
    
    print(f"📊 自定义阈值:")
    print(f"   分子式验证: {config.thresholds.molecular_formula}")
    print(f"   原子组成验证: {config.thresholds.atom_composition}")
    print(f"   指纹相似性: {config.thresholds.fingerprint_similarity}")
    print(f"   3D结构验证: {config.thresholds.structure_3d}")
    print(f"🚀 并行处理: {config.performance.enable_parallel}")


def demo_verification_methods():
    """演示验证方法控制"""
    print("\n🧪 演示3: 验证方法控制")
    print("-" * 40)
    
    # 创建演示残基
    atoms = [
        AtomInfo('N', 'N', 0.0, 0.0, 0.0, 'ALA', 1, 'A'),
        AtomInfo('CA', 'C', 1.5, 0.0, 0.0, 'ALA', 1, 'A'),
        AtomInfo('C', 'C', 1.5, 1.5, 0.0, 'ALA', 1, 'A'),
        AtomInfo('O', 'O', 1.5, 2.5, 0.0, 'ALA', 1, 'A'),
        AtomInfo('CB', 'C', 2.5, 0.0, 0.0, 'ALA', 1, 'A'),
    ]
    
    residue = ResidueInfo('ALA', 1, 'A', atoms)
    
    print(f"🧬 演示残基: {residue.residue_key}")
    print(f"   分子式: {residue.molecular_formula}")
    print(f"   原子组成: {residue.atom_composition}")
    print(f"   重原子组成: {residue.heavy_atom_composition}")
    
    # 演示不同验证方法组合
    method_combinations = [
        [VerificationMethod.MOLECULAR_FORMULA],
        [VerificationMethod.MOLECULAR_FORMULA, VerificationMethod.ATOM_COMPOSITION],
        [VerificationMethod.MOLECULAR_FORMULA, VerificationMethod.ATOM_COMPOSITION, VerificationMethod.FINGERPRINT_SIMILARITY]
    ]
    
    analyzer = PDBAnalyzer()
    
    for i, methods in enumerate(method_combinations, 1):
        print(f"\n   组合{i}: {[m.value for m in methods]}")
        matches = analyzer.analyze_residue(residue, enabled_methods=methods)
        print(f"   匹配数量: {len(matches)}")


def demo_classification():
    """演示分类功能"""
    print("\n🏷️ 演示4: 氨基酸分类")
    print("-" * 40)
    
    # 创建演示氨基酸
    demo_amino_acids = [
        AminoAcidInfo(
            id='ALA',
            name='L-Alanine',
            molecular_formula='C3H7NO2',
            molecular_weight=89.09,
            smiles='N[C@@H](C)C(=O)O',
            atom_composition={'C': 3, 'H': 7, 'N': 1, 'O': 2}
        ),
        AminoAcidInfo(
            id='PHE',
            name='L-Phenylalanine',
            molecular_formula='C9H11NO2',
            molecular_weight=165.19,
            smiles='N[C@@H](Cc1ccccc1)C(=O)O',
            atom_composition={'C': 9, 'H': 11, 'N': 1, 'O': 2}
        ),
        AminoAcidInfo(
            id='PRO',
            name='L-Proline',
            molecular_formula='C5H9NO2',
            molecular_weight=115.13,
            smiles='N1CCC[C@H]1C(=O)O',
            atom_composition={'C': 5, 'H': 9, 'N': 1, 'O': 2}
        )
    ]
    
    # 初始化分类器
    classifier = AminoAcidClassifier()
    
    # 分类演示
    for amino_acid in demo_amino_acids:
        result = classifier.classify_amino_acid(amino_acid)
        print(f"🧬 {amino_acid.name} ({amino_acid.id}):")
        print(f"   SMILES: {amino_acid.smiles}")
        print(f"   分类: {result['classifications']}")
        print(f"   骨架类型: {result['backbone_type']}")


def demo_system_status():
    """演示系统状态检查"""
    print("\n🔍 演示5: 系统状态检查")
    print("-" * 40)
    
    analyzer = PDBAnalyzer()
    validation = analyzer.validate_system()
    
    print(f"🎯 整体状态: {validation['overall_status']}")
    
    print(f"📊 组件状态:")
    for component, status in validation['components'].items():
        if isinstance(status, dict):
            if status.get('status') == 'healthy':
                print(f"   ✅ {component}: 正常")
            else:
                print(f"   ❌ {component}: 错误")
        else:
            print(f"   ℹ️ {component}: {status}")
    
    if validation['warnings']:
        print(f"⚠️ 警告:")
        for warning in validation['warnings']:
            print(f"   - {warning}")
    
    if validation['errors']:
        print(f"❌ 错误:")
        for error in validation['errors']:
            print(f"   - {error}")


def demo_database_stats():
    """演示数据库统计"""
    print("\n📊 演示6: 数据库统计")
    print("-" * 40)
    
    try:
        analyzer = PDBAnalyzer()
        stats = analyzer.get_database_stats()
        
        print(f"📚 数据库统计:")
        for key, value in stats.items():
            print(f"   {key.replace('_', ' ').title()}: {value}")
    
    except Exception as e:
        print(f"⚠️ 数据库统计获取失败: {e}")
        print("   可能需要先运行数据库迁移: python migrate_database.py")


def main():
    """主演示函数"""
    print("🌟 PDB-UAChecker v2.0 功能演示")
    print("=" * 60)
    print("🎯 展示重构后系统的核心功能和特性")
    print()
    
    try:
        # 运行各个演示
        demo_basic_analysis()
        demo_custom_config()
        demo_verification_methods()
        demo_classification()
        demo_system_status()
        demo_database_stats()
        
        print("\n🎉 演示完成！")
        print("\n💡 更多功能:")
        print("   - 运行 python analyze_pdb_v2.py --help 查看完整用法")
        print("   - 运行 python -m pdb_uachecker.api.cli --help 查看高级功能")
        print("   - 查看 README_v2.md 了解详细文档")
    
    except KeyboardInterrupt:
        print("\n⚠️ 演示被用户中断")
    except Exception as e:
        print(f"\n❌ 演示过程中出错: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()