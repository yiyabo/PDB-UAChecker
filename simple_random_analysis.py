#!/usr/bin/env python3
"""
简化版随机氨基酸分析脚本
直接从SMILES文件读取数据并使用分类验证器进行分析
"""

import os
import json
import random
from typing import Dict, List, Any
from pathlib import Path

def read_smiles_file(file_path: str) -> str:
    """读取SMILES文件内容"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            smiles = f.read().strip()
            return smiles
    except Exception as e:
        print(f"读取SMILES文件失败 {file_path}: {e}")
        return ""

def load_amino_acid_data_from_smiles(data_dir: str = "data/structures") -> Dict[str, Dict]:
    """从SMILES文件加载氨基酸数据"""
    all_data = {}
    
    if not os.path.exists(data_dir):
        print(f"❌ 数据目录不存在: {data_dir}")
        return {}
    
    # 遍历所有氨基酸目录
    for amino_code in os.listdir(data_dir):
        amino_dir = os.path.join(data_dir, amino_code)
        if not os.path.isdir(amino_dir):
            continue
            
        smi_file = os.path.join(amino_dir, f"{amino_code}.smi")
        if os.path.exists(smi_file):
            smiles = read_smiles_file(smi_file)
            if smiles:
                all_data[amino_code] = {
                    'smiles': smiles,
                    'name': amino_code,  # 使用代码作为名称
                    'code': amino_code
                }
    
    print(f"📂 成功加载 {len(all_data)} 个氨基酸SMILES数据")
    return all_data

def analyze_with_classification_validator(code: str, data: Dict) -> Dict[str, Any]:
    """使用分类验证器分析氨基酸"""
    smiles = data['smiles']
    name = data['name']
    
    print(f"\n🔬 分析: {code} - {name}")
    print(f"   SMILES: {smiles}")
    
    try:
        # 导入分类验证器
        from pdb_uachecker.analysis.classification_validator import ClassificationValidator
        
        validator = ClassificationValidator()
        result = validator.validate_classification(smiles, code, name)
        
        analysis = {
            'code': code,
            'name': name,
            'smiles': smiles,
            'classification': {
                'categories': result.final_categories,
                'confidence': result.confidence_score,
                'is_consistent': result.is_consistent
            },
            'details': {
                'inconsistencies': result.inconsistencies,
                'recommendations': result.recommendations[:3]  # 只保留前3个建议
            }
        }
        
        print(f"   ✅ 分类: {', '.join(result.final_categories) if result.final_categories else '无分类'}")
        print(f"   📊 置信度: {result.confidence_score:.3f}")
        print(f"   🎯 一致性: {'是' if result.is_consistent else '否'}")
        
        return analysis
        
    except Exception as e:
        print(f"   ❌ 分析失败: {e}")
        return {
            'code': code,
            'name': name, 
            'smiles': smiles,
            'error': str(e)
        }

def generate_llm_insights(results: List[Dict[str, Any]]) -> str:
    """生成LLM风格的化学见解分析"""
    
    insights = []
    insights.append("🧪 LLM化学结构深度分析报告")
    insights.append("=" * 60)
    
    # 数据概览
    total_samples = len(results)
    successful_analyses = [r for r in results if 'error' not in r]
    success_rate = (len(successful_analyses) / total_samples) * 100 if total_samples > 0 else 0
    
    insights.append(f"📈 数据概览:")
    insights.append(f"   • 总样本数: {total_samples}")
    insights.append(f"   • 成功分析: {len(successful_analyses)}")
    insights.append(f"   • 成功率: {success_rate:.1f}%")
    
    if not successful_analyses:
        insights.append("\n❌ 无成功分析的样本，无法生成进一步分析")
        return "\n".join(insights)
    
    # 分类类型统计
    category_counts = {}
    confidence_scores = []
    consistency_count = 0
    
    for result in successful_analyses:
        categories = result['classification']['categories']
        confidence = result['classification']['confidence']
        
        confidence_scores.append(confidence)
        
        if result['classification']['is_consistent']:
            consistency_count += 1
            
        for category in categories:
            category_counts[category] = category_counts.get(category, 0) + 1
    
    # 分类分布分析
    insights.append(f"\n🔍 分类类型分布:")
    if category_counts:
        sorted_categories = sorted(category_counts.items(), key=lambda x: x[1], reverse=True)
        for category, count in sorted_categories:
            percentage = (count / len(successful_analyses)) * 100
            insights.append(f"   • {category}: {count}次 ({percentage:.1f}%)")
    else:
        insights.append("   • 未检测到明确的分类类型")
    
    # 置信度统计
    if confidence_scores:
        avg_confidence = sum(confidence_scores) / len(confidence_scores)
        max_confidence = max(confidence_scores)
        min_confidence = min(confidence_scores)
        
        insights.append(f"\n📊 分类置信度分析:")
        insights.append(f"   • 平均置信度: {avg_confidence:.3f}")
        insights.append(f"   • 最高置信度: {max_confidence:.3f}")
        insights.append(f"   • 最低置信度: {min_confidence:.3f}")
        
        # 置信度评级
        if avg_confidence >= 0.8:
            confidence_rating = "极高可信度"
        elif avg_confidence >= 0.7:
            confidence_rating = "高可信度"
        elif avg_confidence >= 0.6:
            confidence_rating = "中等可信度"
        else:
            confidence_rating = "需要进一步验证"
        
        insights.append(f"   • 整体评级: {confidence_rating}")
    
    # 一致性分析
    consistency_rate = (consistency_count / len(successful_analyses)) * 100
    insights.append(f"\n🎯 分类一致性分析:")
    insights.append(f"   • 一致性率: {consistency_rate:.1f}% ({consistency_count}/{len(successful_analyses)})")
    
    if consistency_rate >= 90:
        consistency_rating = "极高一致性"
    elif consistency_rate >= 70:
        consistency_rating = "高一致性"
    elif consistency_rate >= 50:
        consistency_rating = "中等一致性"
    else:
        consistency_rating = "一致性待改进"
    
    insights.append(f"   • 一致性评级: {consistency_rating}")
    
    # 化学特征深度分析
    insights.append(f"\n🧬 化学特征深度洞察:")
    
    # D型氨基酸分析
    d_type_count = category_counts.get('d_type_amino_acid', 0)
    if d_type_count > 0:
        d_percentage = (d_type_count / len(successful_analyses)) * 100
        insights.append(f"   • D型立体化学: {d_type_count}个样本 ({d_percentage:.1f}%)")
        insights.append(f"     → 表明显著的非天然立体构型特征")
    
    # 芳香族分析
    aromatic_count = category_counts.get('aromatic_amino_acid', 0)
    if aromatic_count > 0:
        aromatic_percentage = (aromatic_count / len(successful_analyses)) * 100
        insights.append(f"   • 芳香族结构: {aromatic_count}个样本 ({aromatic_percentage:.1f}%)")
        insights.append(f"     → 丰富的π电子共轭系统，影响蛋白质稳定性")
    
    # 环状结构分析
    cyclic_count = category_counts.get('cyclic_amino_acid', 0)
    if cyclic_count > 0:
        cyclic_percentage = (cyclic_count / len(successful_analyses)) * 100
        insights.append(f"   • 环状结构: {cyclic_count}个样本 ({cyclic_percentage:.1f}%)")
        insights.append(f"     → 结构刚性增强，构象受限")
    
    # Alpha氨基酸分析
    alpha_count = category_counts.get('alpha_amino_acid', 0)
    if alpha_count > 0:
        alpha_percentage = (alpha_count / len(successful_analyses)) * 100
        insights.append(f"   • Alpha型骨架: {alpha_count}个样本 ({alpha_percentage:.1f}%)")
        insights.append(f"     → 符合标准蛋白质氨基酸骨架结构")
    
    # N-甲基化分析
    n_methyl_count = category_counts.get('n_methyl_amino_acid', 0)
    if n_methyl_count > 0:
        n_methyl_percentage = (n_methyl_count / len(successful_analyses)) * 100
        insights.append(f"   • N-甲基化: {n_methyl_count}个样本 ({n_methyl_percentage:.1f}%)")
        insights.append(f"     → 氢键能力受损，可能影响二级结构")
    
    # 样本亮点
    insights.append(f"\n⭐ 样本亮点分析:")
    
    # 最高置信度样本
    if confidence_scores:
        highest_conf_result = max(successful_analyses, key=lambda x: x['classification']['confidence'])
        insights.append(f"   • 最可信样本: {highest_conf_result['code']} (置信度: {highest_conf_result['classification']['confidence']:.3f})")
        insights.append(f"     → 分类: {', '.join(highest_conf_result['classification']['categories'])}")
    
    # 最复杂样本（分类最多）
    most_complex = max(successful_analyses, key=lambda x: len(x['classification']['categories']))
    if len(most_complex['classification']['categories']) > 1:
        insights.append(f"   • 最复杂样本: {most_complex['code']} ({len(most_complex['classification']['categories'])}个分类)")
        insights.append(f"     → 分类: {', '.join(most_complex['classification']['categories'])}")
    
    return "\n".join(insights)

def main():
    """主函数"""
    print("🎲 随机氨基酸LLM分析系统")
    print("=" * 50)
    
    # 加载数据
    all_data = load_amino_acid_data_from_smiles()
    
    if len(all_data) < 15:
        print(f"⚠️  可用数据不足15个 (共{len(all_data)}个)，将分析全部")
        sample_size = len(all_data)
    else:
        sample_size = 15
    
    if sample_size == 0:
        print("❌ 没有可用的氨基酸数据")
        return
    
    # 随机选择样本
    selected_codes = random.sample(list(all_data.keys()), sample_size)
    print(f"🎯 随机选择 {sample_size} 个样本进行LLM分析:")
    print(f"   {', '.join(selected_codes)}")
    
    # 分析每个样本
    analysis_results = []
    
    for code in selected_codes:
        data = all_data[code]
        result = analyze_with_classification_validator(code, data)
        analysis_results.append(result)
    
    # 生成LLM风格的深度分析
    llm_insights = generate_llm_insights(analysis_results)
    print(f"\n{llm_insights}")
    
    # 保存结果
    output_file = f"llm_analysis_{sample_size}samples.json"
    output_data = {
        'metadata': {
            'total_available': len(all_data),
            'analyzed_samples': sample_size,
            'selected_codes': selected_codes,
            'analysis_type': 'LLM深度分析',
            'timestamp': str(Path().cwd())
        },
        'results': analysis_results,
        'llm_insights': llm_insights
    }
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(output_data, f, indent=2, ensure_ascii=False)
    
    print(f"\n💾 LLM分析结果已保存: {output_file}")
    print("✅ LLM随机分析完成!")

if __name__ == "__main__":
    # 设置随机种子
    random.seed(42)
    main()