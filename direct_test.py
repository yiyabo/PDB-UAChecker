#!/usr/bin/env python3
"""
直接测试脚本 - 避免复杂的导入问题
直接调用分类验证器测试几个样本
"""

import sys
import os
import random
import time
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

try:
    from pdb_uachecker.analysis.classification_validator import ClassificationValidator, ValidationLevel
    print("✅ 成功导入 ClassificationValidator")
except ImportError as e:
    print(f"❌ 导入失败: {e}")
    sys.exit(1)


def load_sample_smiles(sample_id: str, data_dir: Path) -> str:
    """加载样本SMILES"""
    sample_dir = data_dir / sample_id
    smiles_file = sample_dir / f"{sample_id}.smi"
    
    if not smiles_file.exists():
        raise FileNotFoundError(f"SMILES文件不存在: {smiles_file}")
    
    return smiles_file.read_text().strip()


def simple_classification_test():
    """简单分类测试"""
    print("🧬 PDB-UAChecker 直接分类测试")
    print("测试 classification_validator 模块")
    print("=" * 50)
    
    # 数据目录
    data_dir = Path("/Users/apple/AIBD/non-natural_amino_acids/amino_acids_data_demo/data/structures")
    
    if not data_dir.exists():
        print(f"❌ 数据目录不存在: {data_dir}")
        return
    
    # 获取所有样本
    all_samples = []
    for item in data_dir.iterdir():
        if item.is_dir() and (item / f"{item.name}.smi").exists():
            all_samples.append(item.name)
    
    print(f"📁 发现 {len(all_samples)} 个氨基酸样本")
    
    # 随机选择10个样本
    test_samples = random.sample(all_samples, min(10, len(all_samples)))
    print(f"🎲 随机选择: {', '.join(test_samples)}")
    print()
    
    # 创建验证器
    validator = ClassificationValidator(ValidationLevel.MODERATE)
    
    # 测试结果
    success_count = 0
    results = []
    
    for i, sample_id in enumerate(test_samples, 1):
        print(f"🔍 [{i:2d}/10] 测试 {sample_id:8s}", end=" ... ")
        
        start_time = time.time()
        
        try:
            # 加载SMILES
            smiles = load_sample_smiles(sample_id, data_dir)
            
            # 运行分类验证
            result = validator.validate_classification(
                smiles=smiles,
                amino_acid_code=sample_id,
                amino_acid_name=sample_id
            )
            
            processing_time = time.time() - start_time
            
            if result.is_consistent:
                categories_str = ', '.join(result.final_categories[:2]) if result.final_categories else "无"
                confidence = result.confidence_score
                print(f"✅ {categories_str:15s} (置信度: {confidence:.3f}) [{processing_time:.3f}s]")
                success_count += 1
                
                results.append({
                    'sample': sample_id,
                    'success': True,
                    'categories': result.final_categories,
                    'confidence': confidence,
                    'time': processing_time
                })
            else:
                inconsistencies = ', '.join(result.inconsistencies[:1]) if result.inconsistencies else "未知"
                print(f"⚠️ 不一致: {inconsistencies[:20]:20s} [{processing_time:.3f}s]")
                
                results.append({
                    'sample': sample_id,
                    'success': False,
                    'error': f"不一致: {inconsistencies}",
                    'time': processing_time
                })
                
        except Exception as e:
            processing_time = time.time() - start_time
            error_msg = str(e)[:30]
            print(f"❌ 错误: {error_msg:30s} [{processing_time:.3f}s]")
            
            results.append({
                'sample': sample_id,
                'success': False,
                'error': str(e),
                'time': processing_time
            })
    
    print()
    print("=" * 50)
    print("📊 测试结果统计:")
    
    failed_count = 10 - success_count
    print(f"✅ 成功: {success_count}/10 ({success_count*10}%)")
    print(f"❌ 失败: {failed_count}/10 ({failed_count*10}%)")
    
    if success_count > 0:
        successful_results = [r for r in results if r['success']]
        avg_confidence = sum(r['confidence'] for r in successful_results) / len(successful_results)
        avg_time = sum(r['time'] for r in successful_results) / len(successful_results)
        
        print(f"⚡ 平均处理时间: {avg_time:.3f}秒")
        print(f"🎯 平均置信度: {avg_confidence:.3f}")
        
        # 分类统计
        all_categories = []
        for r in successful_results:
            if r['categories']:
                all_categories.extend(r['categories'])
        
        if all_categories:
            category_count = {}
            for cat in all_categories:
                category_count[cat] = category_count.get(cat, 0) + 1
            
            print("\n📋 分类结果分布:")
            for category, count in sorted(category_count.items(), key=lambda x: x[1], reverse=True):
                print(f"   {category:20s}: {count} 次")
    
    if failed_count > 0:
        failed_results = [r for r in results if not r['success']]
        print(f"\n❌ 主要失败原因:")
        error_types = {}
        for r in failed_results:
            error = r['error'][:40]
            error_types[error] = error_types.get(error, 0) + 1
        
        for error, count in sorted(error_types.items(), key=lambda x: x[1], reverse=True)[:3]:
            print(f"   {error:40s}: {count} 次")
    
    print("\n" + "=" * 50)
    total_time = sum(r['time'] for r in results)
    print(f"⏱️ 总测试时间: {total_time:.2f}秒")
    
    # 系统评估
    if success_count >= 8:
        print("🎉 系统运行良好！大部分分类成功")
    elif success_count >= 6:
        print("⚠️ 系统基本正常，但需要优化")
    elif success_count >= 3:
        print("🔧 系统部分功能正常，需要调试")
    else:
        print("🚨 系统存在严重问题，需要全面检查")


if __name__ == "__main__":
    simple_classification_test()