"""
强化手性分析器
解决不同化学信息学工具间的手性判断不一致问题
整合多种方法提供最可靠的手性分析
"""

from typing import Dict, List, Optional, Tuple, Any
from enum import Enum

from .cip_rule_analyzer import CIPRuleAnalyzer


class ChiralityConsistencyLevel(Enum):
    """手性一致性等级"""
    HIGH = "high"           # 所有方法一致
    MODERATE = "moderate"   # 多数方法一致  
    LOW = "low"            # 方法间有分歧
    CONFLICTED = "conflicted"  # 严重分歧


class RobustChiralityAnalyzer:
    """
    强化手性分析器
    
    核心策略：
    1. 多方法验证：CIP规则 + RDKit + 名称解析 + 结构特征
    2. 一致性检测：检测不同方法间的分歧
    3. 智能权重：根据分子复杂度调整各方法权重
    4. 冲突解决：提供明确的冲突解决策略
    """
    
    def __init__(self):
        """初始化强化手性分析器"""
        # 核心分析器
        self.cip_analyzer = CIPRuleAnalyzer()
        
        # 各方法的基础权重
        self.method_weights = {
            'cip_rule_strict': 0.4,      # 我们的严格CIP实现
            'rdkit_standard': 0.3,       # RDKit标准实现
            'name_parsing': 0.2,         # 基于名称的判断
            'structural_features': 0.1   # 基于结构特征
        }
        
        # 分子复杂度调整因子
        self.complexity_adjustments = {
            'simple': {'cip_rule_strict': 1.0, 'rdkit_standard': 1.0},
            'moderate': {'cip_rule_strict': 1.2, 'rdkit_standard': 0.9},
            'complex': {'cip_rule_strict': 1.5, 'rdkit_standard': 0.7}
        }
    
    def analyze_chirality_robust(self, smiles: str, amino_acid_name: str = "",
                                amino_acid_code: str = "") -> Dict[str, Any]:
        """
        强化手性分析
        
        Args:
            smiles: SMILES字符串
            amino_acid_name: 氨基酸名称
            amino_acid_code: 氨基酸代码
            
        Returns:
            强化手性分析结果
        """
        if not smiles:
            return self._create_unknown_result("No SMILES provided")
        
        # Step 1: 使用多种方法分析
        method_results = self._run_multiple_methods(smiles, amino_acid_name, amino_acid_code)
        
        # Step 2: 评估分子复杂度
        complexity = self._assess_molecule_complexity(smiles)
        
        # Step 3: 检测一致性
        consistency_analysis = self._analyze_consistency(method_results)
        
        # Step 4: 计算加权结果
        weighted_result = self._calculate_weighted_result(
            method_results, complexity, consistency_analysis
        )
        
        # Step 5: 生成最终结果
        return self._generate_final_result(
            weighted_result, method_results, consistency_analysis, complexity
        )
    
    def _run_multiple_methods(self, smiles: str, amino_acid_name: str, 
                             amino_acid_code: str) -> Dict[str, Any]:
        """运行多种手性分析方法"""
        results = {}
        
        # 方法1: 严格CIP规则（我们的实现）
        try:
            cip_result = self.cip_analyzer.analyze_stereochemistry(smiles, amino_acid_name)
            results['cip_rule_strict'] = {
                'stereochemistry': cip_result.get('stereochemistry'),
                'confidence': cip_result.get('confidence', 0),
                'rs_configuration': cip_result.get('rs_configuration'),
                'evidence': cip_result.get('evidence', []),
                'success': True
            }
        except Exception as e:
            results['cip_rule_strict'] = {'success': False, 'error': str(e)}
        
        # 方法2: RDKit标准实现
        try:
            rdkit_result = self._analyze_with_rdkit(smiles)
            results['rdkit_standard'] = rdkit_result
        except Exception as e:
            results['rdkit_standard'] = {'success': False, 'error': str(e)}
        
        # 方法3: 基于名称解析
        try:
            name_result = self._analyze_by_name(amino_acid_name, amino_acid_code)
            results['name_parsing'] = name_result
        except Exception as e:
            results['name_parsing'] = {'success': False, 'error': str(e)}
        
        # 方法4: 结构特征分析
        try:
            structure_result = self._analyze_by_structure_features(smiles)
            results['structural_features'] = structure_result
        except Exception as e:
            results['structural_features'] = {'success': False, 'error': str(e)}
        
        return results
    
    def _analyze_with_rdkit(self, smiles: str) -> Dict[str, Any]:
        """使用RDKit分析手性"""
        try:
            from rdkit import Chem
            
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return {'success': False, 'error': 'Invalid SMILES'}
            
            # 标准化SMILES（减少表示差异的影响）
            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            canonical_mol = Chem.MolFromSmiles(canonical_smiles)
            
            # 获取手性中心
            chiral_centers = Chem.FindMolChiralCenters(canonical_mol, includeUnassigned=True)
            
            if chiral_centers:
                atom_idx, chirality = chiral_centers[0]  # 取第一个手性中心
                
                # 映射到D/L类型
                if chirality == 'R':
                    stereochemistry = 'D_form'
                elif chirality == 'S':
                    stereochemistry = 'L_form'
                else:
                    stereochemistry = 'unknown_chirality'
                
                return {
                    'success': True,
                    'stereochemistry': stereochemistry,
                    'rs_configuration': chirality,
                    'confidence': 0.9,  # RDKit通常很可靠
                    'canonical_smiles': canonical_smiles,
                    'chiral_atom_index': atom_idx
                }
            else:
                return {
                    'success': True,
                    'stereochemistry': 'achiral',
                    'confidence': 0.95
                }
                
        except ImportError:
            return {'success': False, 'error': 'RDKit not available'}
        except Exception as e:
            return {'success': False, 'error': str(e)}
    
    def _analyze_by_name(self, amino_acid_name: str, amino_acid_code: str) -> Dict[str, Any]:
        """基于名称分析手性"""
        evidence = []
        confidence = 0.0
        stereochemistry = 'unknown_chirality'
        
        # 检查氨基酸代码
        if amino_acid_code:
            if amino_acid_code.startswith('D') or amino_acid_code.startswith('d'):
                stereochemistry = 'D_form'
                confidence = 0.8
                evidence.append(f"代码前缀表明D型: {amino_acid_code}")
            elif amino_acid_code in ['ALA', 'PHE', 'TYR', 'TRP', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'SER', 'THR', 'CYS', 'ASN', 'GLN', 'ASP', 'GLU', 'LYS', 'ARG', 'HIS']:
                stereochemistry = 'L_form'
                confidence = 0.7
                evidence.append(f"标准氨基酸代码，通常为L型: {amino_acid_code}")
        
        # 检查名称
        if amino_acid_name:
            name_lower = amino_acid_name.lower()
            if name_lower.startswith('d-') or ' d-' in name_lower:
                stereochemistry = 'D_form'
                confidence = max(confidence, 0.8)
                evidence.append(f"名称前缀表明D型: {amino_acid_name}")
            elif name_lower.startswith('l-') or ' l-' in name_lower:
                stereochemistry = 'L_form'
                confidence = max(confidence, 0.8)
                evidence.append(f"名称前缀表明L型: {amino_acid_name}")
        
        return {
            'success': True,
            'stereochemistry': stereochemistry,
            'confidence': confidence,
            'evidence': evidence
        }
    
    def _analyze_by_structure_features(self, smiles: str) -> Dict[str, Any]:
        """基于结构特征分析手性"""
        # 简化的结构特征分析
        # 检查手性标记的密度和复杂度
        
        chiral_markers = smiles.count('@')
        total_atoms = len([c for c in smiles if c.isupper()])
        
        if chiral_markers == 0:
            return {
                'success': True,
                'stereochemistry': 'achiral',
                'confidence': 0.9,
                'evidence': ['未检测到手性标记']
            }
        
        # 简单的复杂度评估
        complexity_score = chiral_markers / max(total_atoms, 1)
        
        return {
            'success': True,
            'stereochemistry': 'has_chirality',
            'confidence': 0.6,  # 结构特征方法置信度较低
            'complexity_score': complexity_score,
            'evidence': [f'检测到{chiral_markers}个手性标记']
        }
    
    def _assess_molecule_complexity(self, smiles: str) -> str:
        """评估分子复杂度"""
        # 基于SMILES的复杂度指标
        chiral_centers = smiles.count('@')
        rings = smiles.count('1') + smiles.count('2') + smiles.count('3')
        branches = smiles.count('(')
        aromatic_atoms = len([c for c in smiles if c.islower() and c.isalpha()])
        
        complexity_score = (
            chiral_centers * 2 +
            rings * 1.5 +
            branches * 1 +
            aromatic_atoms * 0.5
        )
        
        if complexity_score <= 3:
            return 'simple'
        elif complexity_score <= 8:
            return 'moderate'
        else:
            return 'complex'
    
    def _analyze_consistency(self, method_results: Dict[str, Any]) -> Dict[str, Any]:
        """分析各方法间的一致性"""
        successful_results = {
            method: result for method, result in method_results.items()
            if result.get('success', False) and result.get('stereochemistry') != 'unknown_chirality'
        }
        
        if len(successful_results) < 2:
            return {
                'level': ChiralityConsistencyLevel.LOW,
                'agreement_count': len(successful_results),
                'total_methods': len(method_results),
                'conflicts': []
            }
        
        # 统计各种立体化学结果
        stereo_votes = {}
        for method, result in successful_results.items():
            stereo = result['stereochemistry']
            if stereo not in stereo_votes:
                stereo_votes[stereo] = []
            stereo_votes[stereo].append(method)
        
        # 分析一致性
        conflicts = []
        if len(stereo_votes) == 1:
            # 完全一致
            consistency_level = ChiralityConsistencyLevel.HIGH
        elif len(stereo_votes) == 2:
            # 有分歧
            vote_counts = [len(votes) for votes in stereo_votes.values()]
            if max(vote_counts) >= len(successful_results) * 0.7:
                consistency_level = ChiralityConsistencyLevel.MODERATE
            else:
                consistency_level = ChiralityConsistencyLevel.LOW
            
            # 记录冲突
            for stereo, methods in stereo_votes.items():
                for other_stereo, other_methods in stereo_votes.items():
                    if stereo != other_stereo:
                        conflicts.append(f"{stereo} ({', '.join(methods)}) vs {other_stereo} ({', '.join(other_methods)})")
        else:
            consistency_level = ChiralityConsistencyLevel.CONFLICTED
            conflicts = [f"多种不一致结果: {list(stereo_votes.keys())}"]
        
        return {
            'level': consistency_level,
            'stereo_votes': stereo_votes,
            'agreement_count': len(successful_results),
            'total_methods': len(method_results),
            'conflicts': conflicts
        }
    
    def _calculate_weighted_result(self, method_results: Dict[str, Any], 
                                  complexity: str, consistency_analysis: Dict[str, Any]) -> Dict[str, Any]:
        """计算加权结果"""
        # 根据分子复杂度调整权重
        adjusted_weights = self.method_weights.copy()
        complexity_adjustments = self.complexity_adjustments.get(complexity, {})
        
        for method, adjustment in complexity_adjustments.items():
            if method in adjusted_weights:
                adjusted_weights[method] *= adjustment
        
        # 标准化权重
        total_weight = sum(adjusted_weights.values())
        if total_weight > 0:
            adjusted_weights = {k: v / total_weight for k, v in adjusted_weights.items()}
        
        # 计算加权投票
        stereo_scores = {}
        total_confidence = 0
        
        for method, result in method_results.items():
            if not result.get('success', False):
                continue
            
            method_weight = adjusted_weights.get(method, 0)
            method_confidence = result.get('confidence', 0)
            stereo = result.get('stereochemistry')
            
            if stereo and stereo != 'unknown_chirality':
                weighted_score = method_weight * method_confidence
                stereo_scores[stereo] = stereo_scores.get(stereo, 0) + weighted_score
                total_confidence += weighted_score
        
        # 选择得分最高的结果
        if stereo_scores:
            best_stereo = max(stereo_scores, key=stereo_scores.get)
            best_score = stereo_scores[best_stereo]
            
            return {
                'stereochemistry': best_stereo,
                'confidence': best_score,
                'all_scores': stereo_scores,
                'total_confidence': total_confidence,
                'adjusted_weights': adjusted_weights
            }
        
        return {
            'stereochemistry': 'unknown_chirality',
            'confidence': 0.0
        }
    
    def _generate_final_result(self, weighted_result: Dict[str, Any], 
                              method_results: Dict[str, Any],
                              consistency_analysis: Dict[str, Any], 
                              complexity: str) -> Dict[str, Any]:
        """生成最终结果"""
        # 根据一致性调整置信度
        consistency_multiplier = {
            ChiralityConsistencyLevel.HIGH: 1.0,
            ChiralityConsistencyLevel.MODERATE: 0.8,
            ChiralityConsistencyLevel.LOW: 0.6,
            ChiralityConsistencyLevel.CONFLICTED: 0.4
        }
        
        consistency_level = consistency_analysis['level']
        final_confidence = weighted_result['confidence'] * consistency_multiplier[consistency_level]
        
        # 编译证据
        evidence = []
        for method, result in method_results.items():
            if result.get('success') and result.get('evidence'):
                evidence.extend([f"[{method}] {e}" for e in result['evidence'][:2]])
        
        # 添加一致性信息
        if consistency_analysis['conflicts']:
            evidence.append(f"检测到冲突: {'; '.join(consistency_analysis['conflicts'][:2])}")
        
        return {
            'stereochemistry': weighted_result['stereochemistry'],
            'confidence': final_confidence,
            'consistency_level': consistency_level.value,
            'complexity': complexity,
            'evidence': evidence,
            'method_results_summary': {
                method: {
                    'stereochemistry': result.get('stereochemistry'),
                    'confidence': result.get('confidence', 0),
                    'success': result.get('success', False)
                }
                for method, result in method_results.items()
            },
            'consistency_analysis': consistency_analysis,
            'analysis_method': 'robust_multi_method'
        }
    
    def _create_unknown_result(self, reason: str) -> Dict[str, Any]:
        """创建未知结果"""
        return {
            'stereochemistry': 'unknown_chirality',
            'confidence': 0.0,
            'consistency_level': 'unknown',
            'complexity': 'unknown',
            'evidence': [reason],
            'analysis_method': 'robust_multi_method',
            'failure_reason': reason
        }