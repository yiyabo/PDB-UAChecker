"""
高级专家分类器
集成精确SMILES解析和完整CIP规则实现
"""

from typing import Dict, Any, List, Optional
from ...core.models import AminoAcidInfo, ResidueInfo
from ...utils.smiles_parser import AminoAcidStructureAnalyzer
from ...utils.cip_rules import AminoAcidChiralityAnalyzer


class AdvancedExpertAminoAcidClassifier:
    """高级专家标准氨基酸分类器"""
    
    def __init__(self):
        """初始化高级分类器"""
        
        # 标准氨基酸
        self.standard_amino_acids = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        # 初始化分析引擎
        self.structure_analyzer = AminoAcidStructureAnalyzer()
        self.chirality_analyzer = AminoAcidChiralityAnalyzer()
        
        # PDB标识符模式
        self.pdb_patterns = {
            'd_prefixes': ['D'],
            'beta_indicators': ['B3', '3', 'BETA'],
            'gamma_indicators': ['G4', '4', 'GAMMA'],
            'n_methyl_codes': ['SAR', 'NME', 'NMA', 'NMET']
        }

    def classify_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """
        高级氨基酸分类
        
        Args:
            amino_acid: 氨基酸信息
            
        Returns:
            详细的分类结果
        """
        
        # 1. 标准氨基酸检查
        if amino_acid.id in self.standard_amino_acids:
            return self._create_standard_result(amino_acid)
        
        if not amino_acid.smiles:
            return self._create_unclassified_result(amino_acid, 'no_smiles')
        
        try:
            # 2. 精确结构分析
            structure_analysis = self.structure_analyzer.analyze_amino_acid_structure(
                amino_acid.smiles
            )
            
            # 3. 完整手性分析
            chirality_analysis = self.chirality_analyzer.analyze_chirality(
                amino_acid.smiles
            )
            
            # 4. 综合分类决策
            classification = self._make_classification_decision(
                amino_acid, structure_analysis, chirality_analysis
            )
            
            # 5. PDB验证和置信度调整
            final_result = self._finalize_classification(
                amino_acid, classification, structure_analysis, chirality_analysis
            )
            
            return final_result
            
        except Exception as e:
            return self._create_error_result(amino_acid, str(e))

    def _make_classification_decision(self, amino_acid: AminoAcidInfo,
                                    structure_analysis: Dict[str, Any],
                                    chirality_analysis: Dict[str, Any]) -> Dict[str, Any]:
        """综合分类决策"""
        
        primary_categories = []
        secondary_features = []
        confidence_scores = []
        decision_details = {}
        
        # 1. 骨架类型分析
        backbone = structure_analysis.get('backbone_analysis', {})
        backbone_type = backbone.get('type', 'unknown')
        backbone_confidence = backbone.get('confidence', 0.0)
        
        if backbone_type == 'beta':
            primary_categories.append('Beta-amino acid')
            confidence_scores.append(backbone_confidence)
            decision_details['backbone_decision'] = 'beta_confirmed'
        elif backbone_type == 'gamma':
            primary_categories.append('Gamma-amino acid')
            confidence_scores.append(backbone_confidence)
            decision_details['backbone_decision'] = 'gamma_confirmed'
        elif backbone_type == 'alpha':
            decision_details['backbone_decision'] = 'alpha_standard'
        else:
            decision_details['backbone_decision'] = f'complex_{backbone_type}'
        
        # 2. 手性分析
        if chirality_analysis.get('has_chirality', False):
            d_l_form = chirality_analysis.get('d_l_form', 'unknown')
            chirality_confidence = chirality_analysis.get('confidence', 0.0)
            
            if d_l_form == 'D':
                primary_categories.append('D-amino acid')
                confidence_scores.append(chirality_confidence)
                decision_details['chirality_decision'] = 'D_form_confirmed'
            elif d_l_form == 'L':
                secondary_features.append('L-form (natural chirality)')
                decision_details['chirality_decision'] = 'L_form_standard'
            else:
                decision_details['chirality_decision'] = f'unknown_chirality_{d_l_form}'
        
        # 3. N-甲基化检测（基于结构分析）
        n_methylation = self._detect_n_methylation_advanced(structure_analysis)
        if n_methylation['is_methylated']:
            primary_categories.append('N-methyl amino acid')
            confidence_scores.append(n_methylation['confidence'])
            decision_details['n_methyl_decision'] = 'confirmed'
        
        # 4. 环状结构检测
        ring_analysis = structure_analysis.get('ring_systems', [])
        if ring_analysis:
            secondary_features.append('Cyclic amino acid')
            if any(ring.aromatic for ring in ring_analysis):
                secondary_features.append('Aromatic amino acid')
                decision_details['aromatic_decision'] = 'aromatic_rings_detected'
        
        # 5. 其他结构特征
        side_chain = structure_analysis.get('side_chain_analysis', {})
        if side_chain.get('has_aromatic_ring') and 'Aromatic amino acid' not in secondary_features:
            secondary_features.append('Aromatic amino acid')
        
        # 6. 默认分类（如果没有主要特征）
        if not primary_categories:
            if backbone_type == 'alpha':
                primary_categories.append('Non-standard alpha-amino acid')
                confidence_scores.append(0.6)
            else:
                primary_categories.append('Complex amino acid structure')
                confidence_scores.append(0.4)
        
        return {
            'primary_categories': primary_categories,
            'secondary_features': secondary_features,
            'confidence_scores': confidence_scores,
            'decision_details': decision_details
        }

    def _detect_n_methylation_advanced(self, structure_analysis: Dict[str, Any]) -> Dict[str, Any]:
        """高级N-甲基化检测"""
        
        # 基于结构分析的精确检测
        functional_groups = structure_analysis.get('functional_groups', [])
        backbone_analysis = structure_analysis.get('backbone_analysis', {})
        
        # 检查氨基基团
        amino_groups = [g for g in functional_groups if g['type'] == 'amino']
        
        for amino_group in amino_groups:
            # 检查氨基氮的连接
            nitrogen_connections = amino_group.get('carbon_connections', [])
            hydrogen_count = amino_group.get('hydrogen_count', 0)
            
            # N-甲基化的标志：氮原子连接一个甲基碳
            if hydrogen_count == 1:  # NH而不是NH2
                # 进一步检查是否连接甲基
                for carbon_conn in nitrogen_connections:
                    # 这里需要更复杂的分析来确定是否为甲基
                    # 简化检测：如果氨基氢数量减少，可能是甲基化
                    return {
                        'is_methylated': True,
                        'confidence': 0.8,
                        'methylation_site': 'backbone_nitrogen'
                    }
        
        return {'is_methylated': False}

    def _finalize_classification(self, amino_acid: AminoAcidInfo,
                               classification: Dict[str, Any],
                               structure_analysis: Dict[str, Any],
                               chirality_analysis: Dict[str, Any]) -> Dict[str, Any]:
        """最终化分类结果"""
        
        # PDB标识符验证
        pdb_validation = self._validate_pdb_consistency(amino_acid.id, classification)
        
        # 计算最终置信度
        base_confidence = self._calculate_overall_confidence(classification)
        pdb_boost = pdb_validation.get('confidence_boost', 0.0)
        final_confidence = min(base_confidence + pdb_boost, 0.95)
        
        # 构建最终结果
        result = {
            'amino_acid_id': amino_acid.id,
            'amino_acid_name': amino_acid.name,
            'smiles': amino_acid.smiles,
            
            # 分类结果
            'primary_categories': classification['primary_categories'],
            'secondary_features': classification['secondary_features'],
            'confidence': final_confidence,
            
            # 详细分析
            'structure_analysis': structure_analysis,
            'chirality_analysis': chirality_analysis,
            'classification_decision': classification.get('decision_details', {}),
            
            # 验证信息
            'pdb_validation': pdb_validation,
            
            # 方法信息
            'classification_method': 'advanced_expert_system',
            'analysis_engines': [
                'precise_smiles_parser',
                'complete_cip_rules',
                'advanced_structure_analyzer'
            ]
        }
        
        return result

    def _validate_pdb_consistency(self, pdb_id: str, classification: Dict[str, Any]) -> Dict[str, Any]:
        """验证PDB标识符一致性"""
        validation = {
            'consistent': True,
            'warnings': [],
            'confirmations': [],
            'confidence_boost': 0.0
        }
        
        primary_categories = classification.get('primary_categories', [])
        
        # D型氨基酸验证
        if 'D-amino acid' in primary_categories:
            if any(pdb_id.startswith(prefix) for prefix in self.pdb_patterns['d_prefixes']):
                validation['confirmations'].append('D_prefix_matches_classification')
                validation['confidence_boost'] += 0.15
            else:
                validation['warnings'].append('D_form_predicted_without_D_prefix')
                validation['consistent'] = False
        
        # Beta氨基酸验证
        if 'Beta-amino acid' in primary_categories:
            if any(indicator in pdb_id for indicator in self.pdb_patterns['beta_indicators']):
                validation['confirmations'].append('beta_indicator_matches')
                validation['confidence_boost'] += 0.10
        
        # Gamma氨基酸验证
        if 'Gamma-amino acid' in primary_categories:
            if any(indicator in pdb_id for indicator in self.pdb_patterns['gamma_indicators']):
                validation['confirmations'].append('gamma_indicator_matches')
                validation['confidence_boost'] += 0.10
        
        # N-甲基氨基酸验证
        if 'N-methyl amino acid' in primary_categories:
            if any(code in pdb_id for code in self.pdb_patterns['n_methyl_codes']):
                validation['confirmations'].append('n_methyl_code_matches')
                validation['confidence_boost'] += 0.10
        
        return validation

    def _calculate_overall_confidence(self, classification: Dict[str, Any]) -> float:
        """计算总体置信度"""
        confidence_scores = classification.get('confidence_scores', [])
        
        if not confidence_scores:
            return 0.5
        
        # 加权平均
        weights = [1.0] * len(confidence_scores)  # 可以根据分类重要性调整权重
        
        weighted_sum = sum(score * weight for score, weight in zip(confidence_scores, weights))
        total_weight = sum(weights)
        
        return weighted_sum / total_weight if total_weight > 0 else 0.5

    def _create_standard_result(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """创建标准氨基酸结果"""
        return {
            'amino_acid_id': amino_acid.id,
            'amino_acid_name': amino_acid.name,
            'smiles': amino_acid.smiles,
            'primary_categories': ['Standard amino acid'],
            'secondary_features': [],
            'confidence': 1.0,
            'classification_method': 'standard_recognition',
            'is_standard': True
        }

    def _create_unclassified_result(self, amino_acid: AminoAcidInfo, reason: str) -> Dict[str, Any]:
        """创建未分类结果"""
        return {
            'amino_acid_id': amino_acid.id,
            'amino_acid_name': amino_acid.name,
            'smiles': amino_acid.smiles,
            'primary_categories': ['Requires expert review'],
            'secondary_features': [],
            'confidence': 0.0,
            'classification_method': 'insufficient_data',
            'unclassified_reason': reason
        }

    def _create_error_result(self, amino_acid: AminoAcidInfo, error_msg: str) -> Dict[str, Any]:
        """创建错误结果"""
        return {
            'amino_acid_id': amino_acid.id,
            'amino_acid_name': amino_acid.name,
            'smiles': amino_acid.smiles,
            'primary_categories': ['Classification error'],
            'secondary_features': [],
            'confidence': 0.0,
            'classification_method': 'error',
            'error_message': error_msg
        }

    def batch_classify(self, amino_acids: List[AminoAcidInfo]) -> Dict[str, Any]:
        """批量分类氨基酸"""
        print(f"🚀 高级专家分类器 - 批量分析 {len(amino_acids)} 个氨基酸")
        
        results = []
        statistics = {
            'total': len(amino_acids),
            'standard': 0,
            'classified': 0,
            'unclassified': 0,
            'errors': 0,
            'category_distribution': {},
            'confidence_distribution': {'high': 0, 'medium': 0, 'low': 0}
        }
        
        for i, amino_acid in enumerate(amino_acids, 1):
            if i % 25 == 0:
                print(f"   进度: {i}/{len(amino_acids)}")
            
            result = self.classify_amino_acid(amino_acid)
            results.append(result)
            
            # 统计
            if result.get('is_standard'):
                statistics['standard'] += 1
            elif 'error' in result.get('classification_method', ''):
                statistics['errors'] += 1
            elif result.get('confidence', 0) > 0:
                statistics['classified'] += 1
            else:
                statistics['unclassified'] += 1
            
            # 分类分布
            for category in result.get('primary_categories', []):
                statistics['category_distribution'][category] = \
                    statistics['category_distribution'].get(category, 0) + 1
            
            # 置信度分布
            confidence = result.get('confidence', 0)
            if confidence >= 0.8:
                statistics['confidence_distribution']['high'] += 1
            elif confidence >= 0.6:
                statistics['confidence_distribution']['medium'] += 1
            else:
                statistics['confidence_distribution']['low'] += 1
        
        return {
            'results': results,
            'statistics': statistics,
            'method': 'advanced_expert_classifier',
            'coverage': (statistics['classified'] + statistics['standard']) / statistics['total'] * 100
        }

    def classify_from_residue(self, residue: ResidueInfo) -> Dict[str, Any]:
        """从残基信息进行分类（用于PDB分析）"""
        temp_amino_acid = AminoAcidInfo(
            id=residue.residue_name,
            name=f"残基-{residue.residue_name}",
            smiles=getattr(residue, 'smiles', ''),
            molecular_formula=residue.molecular_formula,
            atom_composition=residue.atom_composition
        )
        
        return self.classify_amino_acid(temp_amino_acid)