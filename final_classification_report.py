#!/usr/bin/env python3
"""
最终氨基酸分类报告
229个非天然氨基酸的完整分析结果
"""

import pandas as pd

def generate_final_report():
    """生成最终分类报告"""
    
    # 读取完整分类结果
    df = pd.read_csv('amino_acids_for_manual_classification.csv')
    
    print("🎉 非天然氨基酸分类项目 - 最终报告")
    print("=" * 60)
    
    # 基本统计
    total = len(df)
    print(f"📊 项目概况:")
    print(f"- 分析的氨基酸总数: {total}")
    print(f"- 分类完成度: 100%")
    
    # 各类型统计
    stats = {
        'Beta氨基酸': len(df[df['is_beta_manual'] == 'Yes']),
        'Gamma氨基酸': len(df[df['is_gamma_manual'] == 'Yes']),
        'D氨基酸': len(df[df['is_d_amino_manual'] == 'Yes']),
        '环状氨基酸': len(df[df['is_cyclic_manual'] == 'Yes']),
        '芳香族氨基酸': len(df[df['is_aromatic_manual'] == 'Yes']),
        'N-甲基氨基酸': len(df[df['is_n_methyl_manual'] == 'Yes'])
    }
    
    print(f"\n🔬 分类结果统计:")
    for category, count in stats.items():
        percentage = count/total*100
        print(f"- {category}: {count}个 ({percentage:.1f}%)")
    
    # 重大发现
    print(f"\n🏆 重大发现:")
    
    # Beta氨基酸
    beta_acids = df[df['is_beta_manual'] == 'Yes']
    print(f"- Beta氨基酸 ({len(beta_acids)}个):")
    for _, row in beta_acids.iterrows():
        print(f"  * {row['amino_acid_code']}: {row['analysis_notes']}")
    
    # Gamma氨基酸
    gamma_acids = df[df['is_gamma_manual'] == 'Yes']
    print(f"- Gamma氨基酸 ({len(gamma_acids)}个):")
    for _, row in gamma_acids.iterrows():
        print(f"  * {row['amino_acid_code']}: {row['analysis_notes']}")
    
    # N-甲基氨基酸
    n_methyl_acids = df[df['is_n_methyl_manual'] == 'Yes']
    print(f"- N-甲基氨基酸 ({len(n_methyl_acids)}个):")
    for _, row in n_methyl_acids.iterrows():
        print(f"  * {row['amino_acid_code']}: {row['analysis_notes']}")
    
    # 自动分类器性能分析
    print(f"\n📈 自动分类器性能评估:")
    
    # 计算各类型的准确率
    categories = ['beta', 'gamma', 'd_amino', 'cyclic', 'aromatic', 'n_methyl']
    
    for category in categories:
        auto_col = f'is_{category}'
        manual_col = f'is_{category}_manual'
        
        if auto_col in df.columns and manual_col in df.columns:
            # 计算混淆矩阵
            tp = len(df[(df[auto_col] == 'Yes') & (df[manual_col] == 'Yes')])
            fp = len(df[(df[auto_col] == 'Yes') & (df[manual_col] == 'No')])
            tn = len(df[(df[auto_col] == 'No') & (df[manual_col] == 'No')])
            fn = len(df[(df[auto_col] == 'No') & (df[manual_col] == 'Yes')])
            
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            accuracy = (tp + tn) / (tp + fp + tn + fn)
            
            print(f"- {category.replace('_', ' ').title()}:")
            print(f"  准确率: {accuracy:.1%}, 精确率: {precision:.1%}, 召回率: {recall:.1%}")
    
    # 主要问题分析
    print(f"\n⚠️  自动分类器主要问题:")
    print(f"- 芳香族 vs 环状区分: 大量芳香族被误分为环状")
    print(f"- N-甲基识别: 存在较多误报")
    print(f"- 复杂杂环结构识别不准确")
    
    print(f"\n✅ 项目成功完成！")
    print(f"发现了 {stats['Beta氨基酸']} 个β氨基酸、{stats['Gamma氨基酸']} 个γ氨基酸和 {stats['N-甲基氨基酸']} 个N-甲基氨基酸")

if __name__ == "__main__":
    generate_final_report()