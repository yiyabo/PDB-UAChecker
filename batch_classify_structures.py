#!/usr/bin/env python3
"""
批量分类 data/structures 中的所有氨基酸结构
使用核心分类验证系统进行高精度分类并生成CSV报告
"""

import sys
import os
import csv
import json
import time
from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import asdict

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

try:
    from pdb_uachecker.analysis.classification_validator import ClassificationValidator, ValidationLevel
    from pdb_uachecker.core.models import AminoAcidInfo
    print("✅ 成功导入分类验证系统")
except ImportError as e:
    print(f"❌ 导入失败: {e}")
    sys.exit(1)


class StructureBatchClassifier:
    """批量结构分类器"""
    
    def __init__(self, validation_level: ValidationLevel = ValidationLevel.MODERATE):
        """初始化批量分类器"""
        self.validator = ClassificationValidator(validation_level)
        self.results = []
        self.total_processed = 0
        self.successful_classifications = 0
        
    def load_amino_acid_data(self, structure_dir: Path) -> Optional[Dict[str, str]]:
        """加载氨基酸数据"""
        amino_acid_code = structure_dir.name
        smiles_file = structure_dir / f"{amino_acid_code}.smi"
        
        # 读取SMILES
        smiles = None
        if smiles_file.exists():
            try:
                with open(smiles_file, 'r', encoding='utf-8') as f:
                    content = f.read().strip()
                    # 提取SMILES (通常是第一个非空行或第一个词)
                    for line in content.split('\n'):
                        line = line.strip()
                        if line and not line.startswith('#'):
                            smiles = line.split()[0]  # 取第一个词作为SMILES
                            break
            except Exception as e:
                print(f"⚠️ 读取{smiles_file}失败: {e}")
        
        if not smiles:
            return None
            
        return {
            'amino_acid_code': amino_acid_code,
            'smiles': smiles,
            'structure_dir': str(structure_dir)
        }
    
    def classify_single_structure(self, structure_data: Dict[str, str]) -> Dict[str, Any]:
        """分类单个结构"""
        code = structure_data['amino_acid_code']
        smiles = structure_data['smiles']
        
        print(f"🧪 正在分类: {code} - {smiles}")
        
        try:
            # 创建AminoAcidInfo对象
            amino_acid_info = AminoAcidInfo(
                id=code,
                name=code,  # 使用代码作为名称
                molecular_formula="",  # 将在验证过程中计算
                molecular_weight=0.0,  # 将在验证过程中计算
                smiles=smiles,
                atom_composition={}  # 将在验证过程中计算
            )
            
            # 执行分类验证
            start_time = time.time()
            validation_result = self.validator.validate_classification(
                smiles=smiles,
                amino_acid_code=code,
                amino_acid_name=code
            )
            processing_time = time.time() - start_time
            
            # 整理结果
            result = {
                'amino_acid_code': code,
                'smiles': smiles,
                'classification_successful': validation_result.is_consistent,
                'final_categories': validation_result.final_categories,
                'confidence_score': round(validation_result.confidence_score, 3),
                'processing_time_sec': round(processing_time, 3),
                'inconsistencies_count': len(validation_result.inconsistencies),
                'recommendations': validation_result.recommendations,
                'status': 'SUCCESS' if validation_result.is_consistent else 'WARNING',
                'main_category': validation_result.final_categories[0] if validation_result.final_categories else 'unknown',
                'all_categories': ', '.join(validation_result.final_categories) if validation_result.final_categories else 'none'
            }
            
            if validation_result.is_consistent:
                self.successful_classifications += 1
                print(f"✅ {code}: {result['all_categories']} (置信度: {result['confidence_score']})")
            else:
                print(f"⚠️ {code}: 分类存在不确定性 (置信度: {result['confidence_score']})")
            
            return result
            
        except Exception as e:
            print(f"❌ {code} 分类失败: {e}")
            return {
                'amino_acid_code': code,
                'smiles': smiles,
                'classification_successful': False,
                'final_categories': [],
                'confidence_score': 0.0,
                'processing_time_sec': 0.0,
                'inconsistencies_count': 999,
                'recommendations': [f"分类失败: {str(e)}"],
                'status': 'FAILED',
                'main_category': 'error',
                'all_categories': 'processing_failed'
            }
    
    def classify_all_structures(self, structures_dir: Path) -> List[Dict[str, Any]]:
        """批量分类所有结构"""
        print(f"🚀 开始批量分类 {structures_dir}")
        print("=" * 60)
        
        # 收集所有结构目录
        structure_dirs = [d for d in structures_dir.iterdir() if d.is_dir()]
        total_structures = len(structure_dirs)
        
        print(f"📊 发现 {total_structures} 个氨基酸结构")
        print()
        
        results = []
        failed_loads = 0
        
        for i, structure_dir in enumerate(structure_dirs, 1):
            print(f"[{i}/{total_structures}] ", end="")
            
            # 加载数据
            structure_data = self.load_amino_acid_data(structure_dir)
            if not structure_data:
                print(f"❌ {structure_dir.name}: 无法加载SMILES数据")
                failed_loads += 1
                continue
            
            # 分类
            result = self.classify_single_structure(structure_data)
            results.append(result)
            self.total_processed += 1
            
            # 进度报告
            if i % 10 == 0:
                success_rate = (self.successful_classifications / self.total_processed) * 100 if self.total_processed > 0 else 0
                print(f"\n📈 进度: {i}/{total_structures} | 成功率: {success_rate:.1f}%\n")
        
        self.results = results
        
        # 最终统计
        print("\n" + "=" * 60)
        print("📊 批量分类完成")
        print(f"总计处理: {self.total_processed}")
        print(f"成功分类: {self.successful_classifications}")
        print(f"数据加载失败: {failed_loads}")
        print(f"整体成功率: {(self.successful_classifications / self.total_processed) * 100:.1f}%" if self.total_processed > 0 else "0%")
        
        return results
    
    def save_results_csv(self, output_file: Path):
        """保存结果为CSV格式"""
        if not self.results:
            print("❌ 没有结果可保存")
            return
        
        fieldnames = [
            '氨基酸代码', 'SMILES', '分类状态', '主要类别', '所有分类', 
            '置信度', '处理时间(秒)', '不一致数', '推荐建议'
        ]
        
        with open(output_file, 'w', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for result in self.results:
                writer.writerow({
                    '氨基酸代码': result['amino_acid_code'],
                    'SMILES': result['smiles'],
                    '分类状态': result['status'],
                    '主要类别': result['main_category'],
                    '所有分类': result['all_categories'],
                    '置信度': result['confidence_score'],
                    '处理时间(秒)': result['processing_time_sec'],
                    '不一致数': result['inconsistencies_count'],
                    '推荐建议': '; '.join(result['recommendations']) if result['recommendations'] else ''
                })
        
        print(f"📄 CSV结果已保存到: {output_file}")
    
    def save_results_json(self, output_file: Path):
        """保存详细结果为JSON格式"""
        if not self.results:
            print("❌ 没有结果可保存")
            return
        
        # 添加汇总统计
        summary = {
            'batch_summary': {
                'total_processed': self.total_processed,
                'successful_classifications': self.successful_classifications,
                'success_rate': (self.successful_classifications / self.total_processed) * 100 if self.total_processed > 0 else 0,
                'processing_timestamp': time.strftime("%Y-%m-%d %H:%M:%S")
            },
            'classification_results': self.results
        }
        
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2, ensure_ascii=False)
        
        print(f"📋 详细结果已保存到: {output_file}")


def main():
    """主函数"""
    print("🔬 氨基酸结构批量分类系统")
    print("=" * 60)
    
    # 检查structures目录
    structures_dir = Path("data/structures")
    if not structures_dir.exists():
        print(f"❌ 找不到目录: {structures_dir}")
        return
    
    # 初始化分类器
    classifier = StructureBatchClassifier(ValidationLevel.MODERATE)
    
    # 执行批量分类
    start_time = time.time()
    results = classifier.classify_all_structures(structures_dir)
    total_time = time.time() - start_time
    
    print(f"\n⏱️ 总耗时: {total_time:.2f} 秒")
    print(f"⚡ 平均每个: {total_time / len(results):.2f} 秒" if results else "")
    
    # 保存结果
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    csv_file = Path(f"classification_results_{timestamp}.csv")
    json_file = Path(f"classification_results_{timestamp}.json")
    
    classifier.save_results_csv(csv_file)
    classifier.save_results_json(json_file)
    
    print(f"\n🎉 批量分类完成！共处理 {len(results)} 个氨基酸结构")


if __name__ == "__main__":
    main()