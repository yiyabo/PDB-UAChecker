#!/usr/bin/env python3
"""
智能分类器准确性验证脚本
随机抽取氨基酸样本，详细验证分类准确性
"""

import sys
import random
import re
from pathlib import Path
from typing import Dict, List, Any
import json

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.analysis.intelligent_classifier import IntelligentAminoAcidClassifier
from pdb_uachecker.core.database.enhanced_manager import EnhancedDatabaseManager
from pdb_uachecker.utils.config import Config


class ClassifierAccuracyValidator:
    """分类器准确性验证器"""
    
    def __init__(self):
        """初始化验证器"""
        self.config = Config()
        self.db = EnhancedDatabaseManager(self.config)
        self.classifier = IntelligentAminoAcidClassifier()
        
        # SMILES结构分析规则
        self.structure_patterns = {
            'modified_aromatic': {
                'pattern': r'c1c.*[FNOS][^H]',
                'description': '修饰的芳香环氨基酸'
            },
            'modified_aliphatic': {
                'pattern': r'C[CNOS].*[^c1]',
                'description': '修饰的脂肪族氨基酸'
            },
            'phosphorylated': {
                'pattern': r'P\([=O]\)\([OH]\)',
                'description': '磷酸化氨基酸'
            },
            'methylated': {
                'pattern': r'N\(C\)',
                'description': '甲基化氨基酸'
            },
            'halogenated': {
                'pattern': r'[FClBrI]',
                'description': '卤化氨基酸'
            },
            'hydroxylated': {
                'pattern': r'O[^=P]',
                'description': '羟基化氨基酸'
            }
        }
    
    def analyze_smiles_structure(self, smiles: str) -> Dict[str, Any]:
        """分析SMILES结构特征"""
        analysis = {
            'structural_features': [],
            'complexity_score': 0,
            'predicted_category': None,
            'confidence': 0.0
        }
        
        if not smiles:
            return analysis
        
        # 检测结构特征
        for feature, info in self.structure_patterns.items():
            if re.search(info['pattern'], smiles):
                analysis['structural_features'].append({
                    'feature': feature,
                    'description': info['description']
                })
        
        # 复杂度评分
        analysis['complexity_score'] = self._calculate_complexity_score(smiles)
        
        # 预测分类
        analysis['predicted_category'] = self._predict_category_from_smiles(smiles)
        analysis['confidence'] = self._calculate_structure_confidence(smiles)
        
        return analysis
    
    def _calculate_complexity_score(self, smiles: str) -> float:
        """计算结构复杂度评分"""
        if not smiles:
            return 0.0
        
        score = 0.0
        
        # 原子数量
        atoms = len([c for c in smiles if c.isupper()])
        score += atoms * 0.1
        
        # 环结构
        rings = smiles.count('1') + smiles.count('2') + smiles.count('3')
        score += rings * 0.5
        
        # 双键/三键
        double_bonds = smiles.count('=')
        triple_bonds = smiles.count('#')
        score += double_bonds * 0.3 + triple_bonds * 0.5
        
        # 杂原子
        heteroatoms = len([c for c in smiles if c in 'NOSPF'])
        score += heteroatoms * 0.2
        
        return min(score, 10.0)  # 最大评分10
    
    def _predict_category_from_smiles(self, smiles: str) -> str:
        """基于SMILES预测氨基酸分类"""
        if not smiles:
            return 'unknown'
        
        # 标准氨基酸SMILES片段
        standard_patterns = {
            'ALA': 'C[C@@H](N)C(=O)O',
            'GLY': 'C(C(=O)O)N',
            'VAL': 'CC(C)[C@@H](N)C(=O)O',
            'LEU': 'CC(C)C[C@@H](N)C(=O)O',
            'ILE': 'CC[C@H](C)[C@@H](N)C(=O)O',
            'PRO': 'C1C[C@H](NC1)C(=O)O',
            'PHE': 'c1ccc(cc1)C[C@@H](N)C(=O)O',
            'TRP': 'c1ccc2c(c1)c(c[nH]2)C[C@@H](N)C(=O)O',
            'TYR': 'c1cc(ccc1C[C@@H](N)C(=O)O)O'
        }
        
        # 检查是否为标准氨基酸
        for aa_code, pattern in standard_patterns.items():
            if pattern in smiles:
                return 'standard'
        
        # 检查修饰类型
        if 'P(=O)' in smiles:
            return 'modified_phosphorylated'
        elif 'c1' in smiles and any(x in smiles for x in 'FClBrI'):
            return 'modified_aromatic'
        elif 'c1' in smiles:
            return 'aromatic_modified'
        elif any(x in smiles for x in 'FClBrI'):
            return 'halogenated'
        elif smiles.count('O') > 2:
            return 'hydroxylated'
        else:
            return 'modified_other'
    
    def _calculate_structure_confidence(self, smiles: str) -> float:
        """计算结构分析置信度"""
        if not smiles:
            return 0.0
        
        confidence = 0.7  # 基础置信度
        
        # 结构完整性
        if smiles.count('(') == smiles.count(')'):
            confidence += 0.1
        
        # 标准氨基酸骨架
        if 'C(=O)O' in smiles and 'N' in smiles:
            confidence += 0.1
        
        # 合理的原子比例
        c_count = smiles.count('C')
        n_count = smiles.count('N')
        o_count = smiles.count('O')
        
        if c_count > 0 and n_count > 0 and o_count >= 2:
            confidence += 0.1
        
        return min(confidence, 1.0)
    
    def validate_classification(self, amino_acid: Dict[str, Any], classification: Dict[str, Any]) -> Dict[str, Any]:
        """验证单个氨基酸的分类结果"""
        validation = {
            'amino_acid_id': amino_acid.get('three_letter_code', 'Unknown'),
            'amino_acid_name': amino_acid.get('name', 'Unknown'),
            'smiles': amino_acid.get('smiles', ''),
            'classification_result': classification,
            'structure_analysis': None,
            'validation_result': {
                'is_reasonable': True,
                'confidence_assessment': 'good',
                'issues': [],
                'manual_verification_needed': False
            }
        }
        
        # 分析SMILES结构
        if validation['smiles']:
            validation['structure_analysis'] = self.analyze_smiles_structure(validation['smiles'])
        
        # 验证分类合理性
        validation['validation_result'] = self._assess_classification_reasonableness(
            validation['smiles'],
            classification,
            validation['structure_analysis']
        )
        
        return validation
    
    def _assess_classification_reasonableness(self, smiles: str, classification: Dict[str, Any], 
                                           structure_analysis: Dict[str, Any]) -> Dict[str, Any]:
        """评估分类合理性"""
        assessment = {
            'is_reasonable': True,
            'confidence_assessment': 'good',
            'issues': [],
            'manual_verification_needed': False
        }
        
        classified_category = classification.get('standard_category', '')
        detection_method = classification.get('detection_method', '')
        confidence = classification.get('confidence', 0.0)
        
        # 置信度评估
        if confidence < 0.5:
            assessment['confidence_assessment'] = 'low'
            assessment['issues'].append(f"置信度较低: {confidence:.2f}")
        elif confidence < 0.7:
            assessment['confidence_assessment'] = 'medium'
        
        # 结构与分类一致性检查
        if structure_analysis:
            predicted_category = structure_analysis.get('predicted_category', '')
            structural_features = structure_analysis.get('structural_features', [])
            
            # 检查分类一致性
            if predicted_category and classified_category:
                if not self._categories_compatible(predicted_category, classified_category):
                    assessment['is_reasonable'] = False
                    assessment['issues'].append(
                        f"结构分析预测 '{predicted_category}' 与分类结果 '{classified_category}' 不一致"
                    )
            
            # 检查特殊结构特征
            for feature_info in structural_features:
                feature = feature_info['feature']
                if feature == 'phosphorylated' and 'phosph' not in classified_category.lower():
                    assessment['issues'].append("检测到磷酸化特征但分类中未体现")
                elif feature == 'halogenated' and 'halogen' not in classified_category.lower():
                    assessment['issues'].append("检测到卤化特征但分类中未体现")
        
        # 检测方法合理性
        if detection_method == 'structure_analysis' and not smiles:
            assessment['is_reasonable'] = False
            assessment['issues'].append("使用结构分析方法但无SMILES数据")
        
        # 决定是否需要人工验证
        if (not assessment['is_reasonable'] or 
            assessment['confidence_assessment'] == 'low' or 
            len(assessment['issues']) >= 2):
            assessment['manual_verification_needed'] = True
        
        return assessment
    
    def _categories_compatible(self, predicted: str, classified: str) -> bool:
        """检查预测分类和实际分类是否兼容"""
        compatibility_map = {
            'standard': ['standard'],
            'modified_aromatic': ['modified', 'aromatic', 'non_standard'],
            'aromatic_modified': ['modified', 'aromatic', 'non_standard'],
            'modified_phosphorylated': ['modified', 'phosphorylated', 'non_standard'],
            'halogenated': ['modified', 'halogenated', 'non_standard'],
            'hydroxylated': ['modified', 'hydroxylated', 'non_standard'],
            'modified_other': ['modified', 'non_standard', 'unusual']
        }
        
        if predicted not in compatibility_map:
            return True  # 未知预测，不判断不兼容
        
        compatible_terms = compatibility_map[predicted]
        classified_lower = classified.lower()
        
        return any(term in classified_lower for term in compatible_terms)
    
    def run_validation(self, sample_size: int = 15) -> Dict[str, Any]:
        """运行完整的验证流程"""
        print("🔍 智能分类器准确性验证")
        print("=" * 80)
        
        # 获取所有氨基酸
        all_amino_acids = self.db.get_all_amino_acids()
        print(f"📊 数据库中共有 {len(all_amino_acids)} 个氨基酸")
        
        # 随机抽取样本
        if len(all_amino_acids) < sample_size:
            sample_size = len(all_amino_acids)
            print(f"⚠️ 调整样本大小为 {sample_size}")
        
        random.seed(42)  # 设置随机种子以便重现结果
        sample_amino_acids = random.sample(all_amino_acids, sample_size)
        print(f"🎯 随机抽取 {len(sample_amino_acids)} 个氨基酸进行测试")
        
        # 对样本进行分类
        print("\n📋 执行分类...")
        classification_results = []
        for amino_acid in sample_amino_acids:
            result = self.classifier.classify(amino_acid)
            classification_results.append(result)
        
        # 验证每个分类结果
        print("🔎 验证分类结果...")
        validations = []
        for i, (amino_acid, classification) in enumerate(zip(sample_amino_acids, classification_results)):
            print(f"  验证中... ({i+1}/{len(sample_amino_acids)})", end='\r')
            validation = self.validate_classification(amino_acid, classification)
            validations.append(validation)
        
        print("\n")
        
        # 生成验证报告
        return self._generate_validation_report(validations)
    
    def _generate_validation_report(self, validations: List[Dict[str, Any]]) -> Dict[str, Any]:
        """生成验证报告"""
        report = {
            'summary': {
                'total_samples': len(validations),
                'reasonable_classifications': 0,
                'needs_manual_verification': 0,
                'low_confidence': 0,
                'issues_found': 0
            },
            'detailed_results': validations,
            'statistics': {
                'classification_tiers': {},
                'confidence_distribution': {'high': 0, 'medium': 0, 'low': 0},
                'detection_methods': {},
                'issue_types': {}
            },
            'problematic_cases': [],
            'recommendations': []
        }
        
        # 统计分析
        for validation in validations:
            assessment = validation['validation_result']
            classification = validation['classification_result']
            
            # 合理性统计
            if assessment['is_reasonable']:
                report['summary']['reasonable_classifications'] += 1
            
            if assessment['manual_verification_needed']:
                report['summary']['needs_manual_verification'] += 1
            
            if assessment['confidence_assessment'] == 'low':
                report['summary']['low_confidence'] += 1
            
            if assessment['issues']:
                report['summary']['issues_found'] += 1
                # 记录问题案例
                if len(assessment['issues']) >= 2 or not assessment['is_reasonable']:
                    report['problematic_cases'].append({
                        'amino_acid_id': validation['amino_acid_id'],
                        'issues': assessment['issues'],
                        'classification': classification.get('standard_category', 'Unknown'),
                        'confidence': classification.get('confidence', 0.0)
                    })
            
            # 分层统计
            tier = classification.get('classification_tier', 'unknown')
            report['statistics']['classification_tiers'][tier] = \
                report['statistics']['classification_tiers'].get(tier, 0) + 1
            
            # 置信度分布
            conf_level = assessment['confidence_assessment']
            report['statistics']['confidence_distribution'][conf_level] += 1
            
            # 检测方法统计
            method = classification.get('detection_method', 'unknown')
            report['statistics']['detection_methods'][method] = \
                report['statistics']['detection_methods'].get(method, 0) + 1
            
            # 问题类型统计
            for issue in assessment['issues']:
                issue_type = issue.split(':')[0] if ':' in issue else issue[:20]
                report['statistics']['issue_types'][issue_type] = \
                    report['statistics']['issue_types'].get(issue_type, 0) + 1
        
        # 生成建议
        report['recommendations'] = self._generate_recommendations(report)
        
        return report
    
    def _generate_recommendations(self, report: Dict[str, Any]) -> List[str]:
        """生成改进建议"""
        recommendations = []
        
        total = report['summary']['total_samples']
        reasonable_rate = report['summary']['reasonable_classifications'] / total
        
        if reasonable_rate < 0.8:
            recommendations.append("分类合理性偏低，建议检查和优化分类规则")
        
        if report['summary']['low_confidence'] > total * 0.3:
            recommendations.append("低置信度案例较多，建议增强结构分析能力")
        
        if report['summary']['needs_manual_verification'] > total * 0.2:
            recommendations.append("需要人工验证的案例较多，考虑添加更多自动化规则")
        
        # 问题类型分析
        issue_types = report['statistics']['issue_types']
        if '结构分析预测' in str(issue_types):
            recommendations.append("存在结构分析与分类结果不一致的情况，需要调优一致性检查")
        
        if '检测到' in str(issue_types):
            recommendations.append("某些结构特征未在分类中体现，建议完善特征识别规则")
        
        return recommendations
    
    def print_detailed_report(self, report: Dict[str, Any]):
        """打印详细验证报告"""
        print("\n" + "=" * 80)
        print("📈 验证结果总结")
        print("=" * 80)
        
        # 总体统计
        summary = report['summary']
        total = summary['total_samples']
        
        print(f"🎯 样本总数: {total}")
        print(f"✅ 合理分类: {summary['reasonable_classifications']} ({summary['reasonable_classifications']/total*100:.1f}%)")
        print(f"⚠️ 需要人工验证: {summary['needs_manual_verification']} ({summary['needs_manual_verification']/total*100:.1f}%)")
        print(f"📉 低置信度: {summary['low_confidence']} ({summary['low_confidence']/total*100:.1f}%)")
        print(f"🚨 发现问题: {summary['issues_found']} ({summary['issues_found']/total*100:.1f}%)")
        
        # 分类统计
        print(f"\n📊 分类层级分布:")
        tier_names = {
            'tier1_definitive': 'Tier 1 - 100%确定',
            'tier2_high_confidence': 'Tier 2 - 高置信度',
            'tier3_intelligent': 'Tier 3 - 智能推断',
            'tier4_needs_review': 'Tier 4 - 需要审查'
        }
        
        for tier, count in report['statistics']['classification_tiers'].items():
            tier_name = tier_names.get(tier, tier)
            print(f"  {tier_name}: {count} ({count/total*100:.1f}%)")
        
        # 置信度分布
        print(f"\n🎯 置信度评估分布:")
        conf_dist = report['statistics']['confidence_distribution']
        for level, count in conf_dist.items():
            print(f"  {level.capitalize()}: {count} ({count/total*100:.1f}%)")
        
        # 检测方法统计
        print(f"\n🔬 检测方法使用:")
        for method, count in report['statistics']['detection_methods'].items():
            print(f"  {method}: {count} ({count/total*100:.1f}%)")
        
        # 详细结果展示
        print(f"\n📋 详细分类结果:")
        print("-" * 80)
        
        for validation in report['detailed_results']:
            aa_id = validation['amino_acid_id']
            aa_name = validation['amino_acid_name']
            classification = validation['classification_result']
            assessment = validation['validation_result']
            structure_analysis = validation['structure_analysis']
            
            # 状态图标
            if assessment['is_reasonable']:
                status_icon = "✅" if not assessment['manual_verification_needed'] else "⚠️"
            else:
                status_icon = "❌"
            
            print(f"\n{status_icon} {aa_id} - {aa_name}")
            print(f"   分类: {classification.get('standard_category', 'Unknown')}")
            print(f"   置信度: {classification.get('confidence', 0.0):.3f}")
            print(f"   方法: {classification.get('detection_method', 'Unknown')}")
            print(f"   层级: {classification.get('classification_tier', 'Unknown')}")
            
            # SMILES信息
            smiles = validation.get('smiles', '')
            if smiles:
                if len(smiles) > 60:
                    print(f"   SMILES: {smiles[:60]}...")
                else:
                    print(f"   SMILES: {smiles}")
            
            # 结构分析结果
            if structure_analysis:
                features = structure_analysis.get('structural_features', [])
                if features:
                    feature_names = [f['feature'] for f in features]
                    print(f"   结构特征: {', '.join(feature_names)}")
                
                complexity = structure_analysis.get('complexity_score', 0)
                print(f"   复杂度: {complexity:.1f}")
                
                predicted = structure_analysis.get('predicted_category')
                if predicted:
                    print(f"   结构预测: {predicted}")
            
            # 问题报告
            issues = assessment.get('issues', [])
            if issues:
                print(f"   ⚠️ 问题:")
                for issue in issues:
                    print(f"      - {issue}")
        
        # 问题案例总结
        if report['problematic_cases']:
            print(f"\n🚨 问题案例总结 ({len(report['problematic_cases'])}个):")
            print("-" * 50)
            
            for i, case in enumerate(report['problematic_cases'], 1):
                print(f"\n{i}. {case['amino_acid_id']} - {case['classification']} (置信度: {case['confidence']:.2f})")
                for issue in case['issues']:
                    print(f"   • {issue}")
        
        # 改进建议
        if report['recommendations']:
            print(f"\n💡 改进建议:")
            print("-" * 30)
            for i, rec in enumerate(report['recommendations'], 1):
                print(f"{i}. {rec}")
        
        print(f"\n🎉 验证完成！")
        
        # 总体评估
        reasonable_rate = summary['reasonable_classifications'] / total
        if reasonable_rate >= 0.9:
            print("🏆 分类器表现优秀！")
        elif reasonable_rate >= 0.8:
            print("👍 分类器表现良好")
        elif reasonable_rate >= 0.7:
            print("👌 分类器表现中等，有改进空间")
        else:
            print("⚠️ 分类器需要重要改进")


def main():
    """主函数"""
    print("🚀 启动智能分类器准确性验证")
    
    try:
        validator = ClassifierAccuracyValidator()
        
        # 运行验证
        report = validator.run_validation(sample_size=15)
        
        # 打印详细报告
        validator.print_detailed_report(report)
        
        # 保存报告到文件
        output_file = Path(__file__).parent / "validation_report.json"
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(report, f, ensure_ascii=False, indent=2)
        
        print(f"\n📄 详细报告已保存至: {output_file}")
        
        return report
        
    except Exception as e:
        print(f"❌ 验证过程中发生错误: {e}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == "__main__":
    result = main()