"""
氨基酸分类器
基于SMILES结构分析对氨基酸进行分类
"""

import re
from typing import List, Dict, Any, Optional, Tuple
from collections import defaultdict

from ..core.models import AminoAcidInfo


class AminoAcidClassifier:
    """氨基酸分类器"""
    
    def __init__(self):
        self.classification_rules = self._initialize_classification_rules()
    
    def _initialize_classification_rules(self) -> Dict[str, List[str]]:
        """初始化分类规则"""
        return {
            'aromatic': [
                r'c1ccccc1',      # 苯环
                r'c1cccc[nH]1',   # 吡咯环
                r'c1cccnc1',      # 吡啶环
                r'c1cnc[nH]1',    # 咪唑环
                r'c1cccs1',       # 噻吩环
                r'c1ccco1',       # 呋喃环
                r'c1ccc2c\(c1\)cccc2',  # 萘环
            ],
            'cyclic': [
                r'\d',  # 任何数字表示环闭合
            ],
            'd_amino': [
                r'C@H(?!@)',  # D型手性中心
            ],
            'beta_amino': [
                r'\[NH3\]C(?!C)[^C]*\[C@@?H\].*C\(=O\)O',
                r'\[NH3\]C(?!C).*\[C\].*\(=O\)=O',
                r'C\(=O\)O.*\[C@@?H\].*C\[NH3\]',
                r'NCC.*C\(=O\)O',
                r'\[NH3\]C[^C\[].*C\(=O\)O',
            ],
            'gamma_amino': [
                r'\[NH3\]CC(?!C)[^C]*\[C@@?H\].*C\(=O\)O',
                r'\[NH3\]CC.*\[C\].*\(=O\)=O',
                r'\[NH3\]CC\[C@@?H\].*\[C\]\(=O\)=O',
                r'C\(=O\)O.*\[C@@?H\].*CC\[NH3\]',
                r'NCCC.*C\(=O\)O',
            ],
            'n_methyl': [
                r'NC\(=O\)',      # N-C=O键（酰胺）
                r'N\[C@@?H\]',    # N直接连手性碳
                r'CN\[C@@?H\]',   # 甲基-N-手性碳
                r'\[NH\]C',       # 仲胺
            ]
        }
    
    def classify_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """
        分类单个氨基酸
        
        Args:
            amino_acid: 氨基酸信息
        
        Returns:
            分类结果
        """
        smiles = amino_acid.smiles
        if not smiles:
            return {
                'amino_acid_id': amino_acid.id,
                'amino_acid_name': amino_acid.name,
                'classifications': ['unclassified'],
                'backbone_type': 'unknown',
                'details': {'error': 'No SMILES available'}
            }
        
        classifications = []
        details = {}
        
        # 检测各种类型
        for category, patterns in self.classification_rules.items():
            if self._matches_category(smiles, patterns):
                classifications.append(category)
                details[f'{category}_patterns'] = [p for p in patterns if re.search(p, smiles)]
        
        # 分析骨架类型
        backbone_type, backbone_info = self._analyze_backbone_structure(smiles)
        
        # 特殊处理：如果检测到beta或gamma，添加到分类中
        if backbone_type == 'beta' and 'beta_amino' not in classifications:
            classifications.append('beta_amino')
        elif backbone_type == 'gamma' and 'gamma_amino' not in classifications:
            classifications.append('gamma_amino')
        
        # 如果没有分类，标记为未分类
        if not classifications:
            classifications.append('unclassified')
        
        return {
            'amino_acid_id': amino_acid.id,
            'amino_acid_name': amino_acid.name,
            'smiles': smiles,
            'classifications': classifications,
            'backbone_type': backbone_type,
            'backbone_info': backbone_info,
            'details': details
        }
    
    def _matches_category(self, smiles: str, patterns: List[str]) -> bool:
        """检查SMILES是否匹配某个类别的模式"""
        return any(re.search(pattern, smiles) for pattern in patterns)
    
    def _analyze_backbone_structure(self, smiles: str) -> Tuple[str, str]:
        """
        分析氨基酸骨架结构
        
        Args:
            smiles: SMILES字符串
        
        Returns:
            (骨架类型, 详细信息)
        """
        if not smiles:
            return 'unknown', 'No SMILES available'
        
        # 特殊情况处理
        special_cases = {
            r'\[NH3\]CC\[C@@?H\]\(\[C\]\(=O\)=O\)\[NH3\]': ('gamma', 'DAB型γ氨基酸'),
            r'\[NH3\]CCC\[C@@?H\]\(\[C\]\(=O\)=O\)\[NH3\]': ('alpha', '鸟氨酸型（α氨基酸带侧链氨基）'),
            r'\[NH3\]CCCCCC\[C@@?H\]\(C\(=O\)O\)\[NH3\]': ('alpha', '赖氨酸衍生物（α氨基酸带侧链氨基）'),
        }
        
        for pattern, (backbone_type, info) in special_cases.items():
            if re.match(pattern, smiles):
                return backbone_type, info
        
        # γ氨基酸检测
        gamma_patterns = [
            r'\[NH3\]CC(?!C)[^C]*\[C@@?H\].*C\(=O\)O',
            r'\[NH3\]CC.*\[C\].*\(=O\)=O',
            r'\[NH3\]CC\[C@@?H\].*\[C\]\(=O\)=O',
            r'C\(=O\)O.*\[C@@?H\].*CC\[NH3\]',
            r'NCCC.*C\(=O\)O',
        ]
        
        for pattern in gamma_patterns:
            if re.search(pattern, smiles):
                return 'gamma', f'匹配模式: {pattern}'
        
        # β氨基酸检测
        beta_patterns = [
            r'\[NH3\]C(?!C)[^C]*\[C@@?H\].*C\(=O\)O',
            r'\[NH3\]C(?!C).*\[C\].*\(=O\)=O',
            r'C\(=O\)O.*\[C@@?H\].*C\[NH3\]',
            r'NCC.*C\(=O\)O',
            r'\[NH3\]C[^C\[].*C\(=O\)O',
        ]
        
        for pattern in beta_patterns:
            if re.search(pattern, smiles):
                # 确保不是γ氨基酸被误判
                if not any(re.search(p, smiles) for p in gamma_patterns):
                    return 'beta', f'匹配模式: {pattern}'
        
        # α氨基酸检测
        alpha_patterns = [
            r'\[NH3\]\[C@@?H\].*C\(=O\)O',
            r'\[NH3\]\[C@@?H\].*\[C\]\(=O\)=O',
            r'C\(=O\)O.*\[C@@?H\].*\[NH3\]',
            r'NCC\(=O\)O',  # 甘氨酸
        ]
        
        for pattern in alpha_patterns:
            if re.search(pattern, smiles):
                return 'alpha', f'匹配模式: {pattern}'
        
        return 'alpha', '默认分类'
    
    def batch_classify(self, amino_acids: List[AminoAcidInfo]) -> Dict[str, Any]:
        """
        批量分类氨基酸
        
        Args:
            amino_acids: 氨基酸信息列表
        
        Returns:
            批量分类结果
        """
        results = []
        classification_counts = defaultdict(int)
        detailed_classifications = defaultdict(list)
        
        print(f"🔬 开始批量分类 {len(amino_acids)} 个氨基酸...")
        
        for amino_acid in amino_acids:
            result = self.classify_amino_acid(amino_acid)
            results.append(result)
            
            # 统计分类
            for classification in result['classifications']:
                classification_counts[classification] += 1
                detailed_classifications[classification].append(amino_acid.id)
        
        # 生成报告
        report = {
            'total_amino_acids': len(amino_acids),
            'classification_counts': dict(classification_counts),
            'detailed_classifications': dict(detailed_classifications),
            'examples': {},
            'backbone_distribution': self._analyze_backbone_distribution(results)
        }
        
        # 添加示例
        for classification in classification_counts:
            for result in results:
                if classification in result['classifications']:
                    report['examples'][classification] = {
                        'amino_acid_id': result['amino_acid_id'],
                        'amino_acid_name': result['amino_acid_name'],
                        'smiles': result['smiles'],
                        'backbone_type': result['backbone_type']
                    }
                    break
        
        return {
            'results': results,
            'report': report
        }
    
    def _analyze_backbone_distribution(self, results: List[Dict[str, Any]]) -> Dict[str, int]:
        """分析骨架类型分布"""
        backbone_counts = defaultdict(int)
        for result in results:
            backbone_type = result.get('backbone_type', 'unknown')
            backbone_counts[backbone_type] += 1
        return dict(backbone_counts)
    
    def get_classification_summary(self, classification_results: Dict[str, Any]) -> str:
        """
        生成分类摘要报告
        
        Args:
            classification_results: 批量分类结果
        
        Returns:
            摘要报告字符串
        """
        report = classification_results['report']
        
        lines = [
            "氨基酸分类摘要报告",
            "=" * 50,
            "",
            f"总计氨基酸数量: {report['total_amino_acids']}",
            "",
            "分类统计:",
        ]
        
        # 分类统计
        for classification, count in sorted(report['classification_counts'].items()):
            if classification == 'unclassified':
                continue
            
            percentage = (count / report['total_amino_acids']) * 100
            lines.append(f"  {classification.replace('_', ' ').title()}: {count} ({percentage:.1f}%)")
            
            # 添加示例
            if classification in report['examples']:
                example = report['examples'][classification]
                lines.append(f"    示例: {example['amino_acid_id']} - {example['amino_acid_name']}")
        
        # 未分类
        if 'unclassified' in report['classification_counts']:
            count = report['classification_counts']['unclassified']
            percentage = (count / report['total_amino_acids']) * 100
            lines.append(f"  未分类: {count} ({percentage:.1f}%)")
        
        # 骨架类型分布
        lines.extend([
            "",
            "骨架类型分布:",
        ])
        
        for backbone_type, count in sorted(report['backbone_distribution'].items()):
            percentage = (count / report['total_amino_acids']) * 100
            lines.append(f"  {backbone_type.title()}: {count} ({percentage:.1f}%)")
        
        return "\n".join(lines)
    
    def find_special_amino_acids(self, amino_acids: List[AminoAcidInfo]) -> Dict[str, List[str]]:
        """
        查找特殊类型的氨基酸
        
        Args:
            amino_acids: 氨基酸信息列表
        
        Returns:
            特殊氨基酸字典
        """
        special_types = {
            'beta_amino_acids': [],
            'gamma_amino_acids': [],
            'd_amino_acids': [],
            'n_methyl_amino_acids': [],
            'aromatic_amino_acids': [],
            'cyclic_amino_acids': []
        }
        
        for amino_acid in amino_acids:
            result = self.classify_amino_acid(amino_acid)
            classifications = result['classifications']
            
            if 'beta_amino' in classifications:
                special_types['beta_amino_acids'].append(amino_acid.id)
            if 'gamma_amino' in classifications:
                special_types['gamma_amino_acids'].append(amino_acid.id)
            if 'd_amino' in classifications:
                special_types['d_amino_acids'].append(amino_acid.id)
            if 'n_methyl' in classifications:
                special_types['n_methyl_amino_acids'].append(amino_acid.id)
            if 'aromatic' in classifications:
                special_types['aromatic_amino_acids'].append(amino_acid.id)
            if 'cyclic' in classifications:
                special_types['cyclic_amino_acids'].append(amino_acid.id)
        
        return special_types
    
    def validate_classification_rules(self, test_cases: Dict[str, List[str]]) -> Dict[str, Any]:
        """
        验证分类规则的准确性
        
        Args:
            test_cases: 测试用例 {分类: [SMILES列表]}
        
        Returns:
            验证结果
        """
        validation_results = {}
        
        for category, smiles_list in test_cases.items():
            if category not in self.classification_rules:
                continue
            
            patterns = self.classification_rules[category]
            correct_matches = 0
            total_tests = len(smiles_list)
            
            for smiles in smiles_list:
                if self._matches_category(smiles, patterns):
                    correct_matches += 1
            
            accuracy = correct_matches / total_tests if total_tests > 0 else 0.0
            
            validation_results[category] = {
                'accuracy': accuracy,
                'correct_matches': correct_matches,
                'total_tests': total_tests,
                'patterns_used': patterns
            }
        
        return validation_results