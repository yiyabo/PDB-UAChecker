#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import json
import csv
from collections import defaultdict

def analyze_amino_acid_structure(smiles):
    """
    分析氨基酸结构，识别氨基和羧基的位置关系
    返回骨架类型和详细信息
    """
    if not smiles:
        return None, None
    
    # 查找羧基模式
    carboxyl_patterns = [
        r'C\(=O\)O',  # C(=O)O
        r'\[C\]\(=O\)=O',  # [C](=O)=O
        r'C\(=O\)=O',  # C(=O)=O
    ]
    
    # 查找氨基模式
    amino_patterns = [
        r'\[NH3\]',  # [NH3]
        r'\[NH2\]',  # [NH2]
        r'N(?![=\(])',  # N (不跟=或()
    ]
    
    # 特殊情况处理
    # 1. DAB: [NH3]CC[C@@H]([C](=O)=O)[NH3] - γ氨基酸
    if re.match(r'\[NH3\]CC\[C@@?H\]\(\[C\]\(=O\)=O\)\[NH3\]', smiles):
        return "gamma", "DAB型γ氨基酸"
    
    # 2. ORN: [NH3]CCC[C@@H]([C](=O)=O)[NH3] - α氨基酸（带额外氨基）
    if re.match(r'\[NH3\]CCC\[C@@?H\]\(\[C\]\(=O\)=O\)\[NH3\]', smiles):
        return "alpha", "鸟氨酸型（α氨基酸带侧链氨基）"
    
    # 3. HHK: [NH3]CCCCCC[C@@H](C(=O)O)[NH3] - α氨基酸（带额外氨基）
    if re.match(r'\[NH3\]CCCCCC\[C@@?H\]\(C\(=O\)O\)\[NH3\]', smiles):
        return "alpha", "赖氨酸衍生物（α氨基酸带侧链氨基）"
    
    # β氨基酸检测 - 更灵活的模式
    # β氨基酸：NH2-CH2-CH(R)-COOH 或类似结构
    beta_patterns = [
        # [NH3]后跟一个碳，然后是含羧基的碳
        r'\[NH3\]C(?!C)[^C]*\[C@@?H\].*C\(=O\)O',
        r'\[NH3\]C(?!C).*\[C\].*\(=O\)=O',
        # 反向：羧基碳后跟一个碳，然后是氨基
        r'C\(=O\)O.*\[C@@?H\].*C\[NH3\]',
        # 简单的β氨基酸模式
        r'NCC.*C\(=O\)O',
        r'\[NH3\]C[^C\[].*C\(=O\)O',
    ]
    
    # γ氨基酸检测 - 更灵活的模式
    # γ氨基酸：NH2-CH2-CH2-CH(R)-COOH 或类似结构
    gamma_patterns = [
        # [NH3]后跟两个碳，然后是含羧基的碳
        r'\[NH3\]CC(?!C)[^C]*\[C@@?H\].*C\(=O\)O',
        r'\[NH3\]CC.*\[C\].*\(=O\)=O',
        # DAB型模式
        r'\[NH3\]CC\[C@@?H\].*\[C\]\(=O\)=O',
        # 反向：羧基碳后跟两个碳，然后是氨基
        r'C\(=O\)O.*\[C@@?H\].*CC\[NH3\]',
        # 简单的γ氨基酸模式
        r'NCCC.*C\(=O\)O',
    ]
    
    # 检测γ氨基酸
    for pattern in gamma_patterns:
        if re.search(pattern, smiles):
            return "gamma", f"匹配模式: {pattern}"
    
    # 检测β氨基酸
    for pattern in beta_patterns:
        if re.search(pattern, smiles):
            # 确保不是γ氨基酸被误判
            if not any(re.search(p, smiles) for p in gamma_patterns):
                return "beta", f"匹配模式: {pattern}"
    
    # α氨基酸检测 - 氨基直接连在含羧基的碳上
    alpha_patterns = [
        r'\[NH3\]\[C@@?H\].*C\(=O\)O',
        r'\[NH3\]\[C@@?H\].*\[C\]\(=O\)=O',
        r'C\(=O\)O.*\[C@@?H\].*\[NH3\]',
        r'NCC\(=O\)O',  # 甘氨酸
    ]
    
    for pattern in alpha_patterns:
        if re.search(pattern, smiles):
            return "alpha", f"匹配模式: {pattern}"
    
    return "alpha", "默认分类"

def is_aromatic_amino_acid(smiles):
    """检测芳香族氨基酸"""
    if not smiles:
        return False
    
    aromatic_patterns = [
        r'c1ccccc1',  # 苯环
        r'c1cccc[nH]1',  # 吡咯环
        r'c1cccnc1',  # 吡啶环
        r'c1cnc[nH]1',  # 咪唑环
        r'c1cccs1',   # 噻吩环
        r'c1ccco1',   # 呋喃环
        r'c1ccc2c\(c1\)cccc2',  # 萘环
    ]
    
    return any(re.search(pattern, smiles) for pattern in aromatic_patterns)

def is_cyclic_amino_acid(smiles):
    """检测环状氨基酸"""
    if not smiles:
        return False
    
    # 检测SMILES中的环闭合数字对
    ring_numbers = re.findall(r'\d', smiles)
    if not ring_numbers:
        return False
    
    # 计算每个数字出现的次数
    from collections import Counter
    ring_counts = Counter(ring_numbers)
    
    # 如果有数字出现2次，说明形成了环
    return any(count >= 2 for count in ring_counts.values())

def is_d_amino_acid(smiles):
    """检测D型氨基酸"""
    if not smiles:
        return False
    
    # D型：C@H，L型：C@@H
    return bool(re.search(r'C@H(?!@)', smiles))

def is_n_methyl_amino_acid(smiles):
    """检测N-甲基氨基酸"""
    if not smiles:
        return False
    
    n_methyl_patterns = [
        r'NC\(=O\)',  # N-C=O键（酰胺）
        r'N\[C@@?H\]',  # N直接连手性碳
        r'CN\[C@@?H\]',  # 甲基-N-手性碳
        r'\[NH\]C',   # 仲胺
    ]
    
    return any(re.search(pattern, smiles) for pattern in n_methyl_patterns)

def classify_amino_acid(smiles):
    """改进的氨基酸分类函数"""
    if not smiles:
        return [], None, None
    
    classifications = []
    
    # 检测各种类型
    if is_aromatic_amino_acid(smiles):
        classifications.append("aromatic")
    
    if is_cyclic_amino_acid(smiles):
        classifications.append("cyclic")
    
    if is_d_amino_acid(smiles):
        classifications.append("d_amino")
    
    # 分析骨架类型
    backbone_type, backbone_info = analyze_amino_acid_structure(smiles)
    
    if backbone_type == "beta":
        classifications.append("beta_amino")
    elif backbone_type == "gamma":
        classifications.append("gamma_amino")
    
    if is_n_methyl_amino_acid(smiles):
        classifications.append("n_methyl")
    
    # 如果没有分类，标记为未分类
    if not classifications:
        classifications.append("unclassified")
    
    return classifications, backbone_type, backbone_info

def main():
    """主函数：使用改进的分类逻辑重新分类所有氨基酸"""
    data_dir = "data/structures"
    
    if not os.path.exists(data_dir):
        print(f"错误：找不到数据目录 {data_dir}")
        return
    
    results = []
    classification_counts = defaultdict(int)
    detailed_classifications = defaultdict(list)
    
    print("开始改进版氨基酸分类分析...")
    
    # 遍历所有氨基酸目录
    amino_acid_dirs = [d for d in os.listdir(data_dir) 
                      if os.path.isdir(os.path.join(data_dir, d))]
    
    # 特别关注的氨基酸
    special_cases = ["DAB", "ORN", "HHK", "ABA", "GBUT"]
    
    for amino_acid_code in sorted(amino_acid_dirs):
        amino_acid_path = os.path.join(data_dir, amino_acid_code)
        smi_file = os.path.join(amino_acid_path, f"{amino_acid_code}.smi")
        
        smiles = None
        if os.path.exists(smi_file):
            try:
                with open(smi_file, 'r', encoding='utf-8') as f:
                    content = f.read().strip()
                    if content:
                        smiles = content.split()[0]
            except Exception as e:
                print(f"读取 {smi_file} 时出错: {e}")
        
        # 分类氨基酸
        classifications, backbone_type, backbone_info = classify_amino_acid(smiles)
        
        # 记录结果
        result = {
            "code": amino_acid_code,
            "smiles": smiles,
            "classifications": classifications,
            "backbone_type": backbone_type,
            "backbone_info": backbone_info
        }
        results.append(result)
        
        # 统计分类
        for classification in classifications:
            classification_counts[classification] += 1
            detailed_classifications[classification].append(amino_acid_code)
        
        # 打印特殊案例
        if amino_acid_code in special_cases:
            print(f"\n特殊案例 - {amino_acid_code}:")
            print(f"  SMILES: {smiles}")
            print(f"  分类: {classifications}")
            print(f"  骨架类型: {backbone_type}")
            print(f"  详情: {backbone_info}")
    
    # 生成报告
    report = {
        "total_amino_acids": len(results),
        "classification_counts": dict(classification_counts),
        "detailed_classifications": dict(detailed_classifications),
        "examples": {},
        "special_cases": {}
    }
    
    # 添加示例
    for classification in classification_counts:
        for result in results:
            if classification in result["classifications"]:
                report["examples"][classification] = {
                    "code": result["code"],
                    "smiles": result["smiles"],
                    "backbone": result["backbone_type"],
                    "info": result["backbone_info"]
                }
                break
    
    # 添加特殊案例详情
    for code in special_cases:
        for result in results:
            if result["code"] == code:
                report["special_cases"][code] = result
                break
    
    # 保存结果
    output_data = {
        "results": results,
        "report": report
    }
    
    # 保存JSON文件
    with open("amino_acid_classification_improved.json", "w", encoding="utf-8") as f:
        json.dump(output_data, f, ensure_ascii=False, indent=2)
    
    # 保存CSV文件
    with open("amino_acid_classification_improved.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Code", "SMILES", "Classifications", "Backbone_Type", "Backbone_Info"])
        for result in results:
            writer.writerow([
                result["code"],
                result["smiles"],
                ";".join(result["classifications"]),
                result["backbone_type"],
                result["backbone_info"]
            ])
    
    # 生成摘要文件
    with open("classification_improved_summary.txt", "w", encoding="utf-8") as f:
        f.write("改进版非天然氨基酸分类摘要\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"总计氨基酸数量: {report['total_amino_acids']}\n\n")
        
        # 分类统计
        f.write("分类统计:\n")
        for classification, count in sorted(report['classification_counts'].items()):
            if classification == "unclassified":
                continue
            
            f.write(f"\n{classification.replace('_', ' ').title()}氨基酸: {count}个\n")
            
            # 写入示例
            if classification in report["examples"]:
                example = report["examples"][classification]
                f.write(f"  示例: {example['code']} - {example['smiles']}\n")
                if example.get('info'):
                    f.write(f"  详情: {example['info']}\n")
            
            # 写入前10个
            codes = detailed_classifications[classification][:10]
            f.write(f"  包含: {', '.join(codes)}")
            if len(detailed_classifications[classification]) > 10:
                f.write(f" 等{len(detailed_classifications[classification])}个")
            f.write("\n")
        
        # 未分类
        if "unclassified" in classification_counts:
            f.write(f"\n未分类: {classification_counts['unclassified']}个\n")
            if "unclassified" in report["examples"]:
                example = report["examples"]["unclassified"]
                f.write(f"  示例: {example['code']} - {example['smiles']}\n")
            codes = detailed_classifications["unclassified"][:10]
            f.write(f"  包含: {', '.join(codes)}")
            if len(detailed_classifications["unclassified"]) > 10:
                f.write(f" 等{len(detailed_classifications['unclassified'])}个")
            f.write("\n")
        
        # 特殊案例分析
        f.write("\n\n特殊案例分析:\n")
        f.write("-" * 50 + "\n")
        for code, details in report.get("special_cases", {}).items():
            f.write(f"\n{code}:\n")
            f.write(f"  SMILES: {details['smiles']}\n")
            f.write(f"  分类: {', '.join(details['classifications'])}\n")
            f.write(f"  骨架类型: {details['backbone_type']}\n")
            if details.get('backbone_info'):
                f.write(f"  详情: {details['backbone_info']}\n")
    
    print(f"\n改进版分析完成！共分析了 {len(results)} 个氨基酸")
    print("\n改进后分类统计:")
    for classification, count in sorted(classification_counts.items()):
        print(f"  {classification}: {count}个")
    
    # 显示β和γ氨基酸的详细列表
    print("\n\nβ氨基酸列表:")
    if "beta_amino" in detailed_classifications:
        for code in detailed_classifications["beta_amino"]:
            for result in results:
                if result["code"] == code:
                    print(f"  {code}: {result['smiles']}")
                    break
    else:
        print("  无")
    
    print("\nγ氨基酸列表:")
    if "gamma_amino" in detailed_classifications:
        for code in detailed_classifications["gamma_amino"]:
            for result in results:
                if result["code"] == code:
                    print(f"  {code}: {result['smiles']}")
                    break
    else:
        print("  无")

if __name__ == "__main__":
    main()