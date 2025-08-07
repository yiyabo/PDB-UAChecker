"""
立体化学分析器
专门负责立体化学和对映异构体的分析
"""

from typing import Dict, List, Optional, Tuple, Any


class StereochemistryAnalyzer:
    """
    立体化学分析器
    
    分析功能：
    1. D/L型判断
    2. R/S构型分析
    3. 对映异构体识别
    4. 手性中心检测
    """
    
    def __init__(self):
        """初始化立体化学分析器"""
        # 立体化学指示符
        self.stereochemistry_indicators = {
            'd_form': {
                'smiles_patterns': [r'C@H(?!@)', r'\[C@H\]'],
                'name_patterns': [r'^d-', r'^D-', r'\bd\s', r'\bD\s'],
                'confidence': 0.85
            },
            'l_form': {
                'smiles_patterns': [r'C@@H', r'\[C@@H\]'],
                'name_patterns': [r'^l-', r'^L-', r'\bl\s', r'\bL\s'],
                'confidence': 0.85
            }
        }
    
    def analyze(self, smiles: str) -> Dict[str, Any]:
        """
        分析立体化学
        
        Args:
            smiles: SMILES字符串
            
        Returns:
            立体化学分析结果
        """
        if not smiles:
            return self._create_unknown_result("No SMILES provided")
        
        # 检测手性标记
        chiral_info = self._detect_chirality_markers(smiles)
        
        if chiral_info['has_chirality']:
            return chiral_info
        
        # 如果没有明确的手性标记，返回默认结果
        return {
            'stereochemistry': 'achiral_or_unknown',
            'confidence': 0.7,
            'chiral_centers': [],
            'evidence': ['未检测到手性标记'],
            'details': {
                'analysis_method': 'basic_detection',
                'stereochemical_type': 'unknown'
            }
        }
    
    def analyze_by_name(self, name: str) -> Dict[str, Any]:
        """
        基于名称分析立体化学
        
        Args:
            name: 氨基酸名称
            
        Returns:
            立体化学分析结果
        """
        import re
        
        name_lower = name.lower()
        
        # D型检测
        if any(re.search(pattern, name_lower) for pattern in self.stereochemistry_indicators['d_form']['name_patterns']):
            return {
                'stereochemistry': 'D_form',
                'confidence': 0.80,
                'evidence': [f"名称表明D型: {name}"],
                'details': {
                    'analysis_method': 'name_based',
                    'stereochemical_type': 'D'
                }
            }
        
        # L型检测
        if any(re.search(pattern, name_lower) for pattern in self.stereochemistry_indicators['l_form']['name_patterns']):
            return {
                'stereochemistry': 'L_form',
                'confidence': 0.80,
                'evidence': [f"名称表明L型: {name}"],
                'details': {
                    'analysis_method': 'name_based',
                    'stereochemical_type': 'L'
                }
            }
        
        return self._create_unknown_result("名称中未发现立体化学指示符")
    
    def _detect_chirality_markers(self, smiles: str) -> Dict[str, Any]:
        """检测手性标记"""
        import re
        
        evidence = []
        stereochemistry = 'achiral_or_unknown'
        confidence = 0.0
        chiral_centers = []
        
        # 检测@标记
        if '@' in smiles:
            # 检测L型（@@）
            if '@@' in smiles:
                l_matches = re.findall(r'C@@H?', smiles)
                if l_matches:
                    stereochemistry = 'L_suspected'
                    confidence = 0.75
                    evidence.append(f"检测到L型手性标记: {len(l_matches)}个")
                    chiral_centers = [{'type': 'L', 'count': len(l_matches)}]
            
            # 检测D型（@但不是@@）
            d_matches = re.findall(r'C@H(?!@)', smiles)
            if d_matches:
                if stereochemistry == 'achiral_or_unknown':
                    stereochemistry = 'D_suspected'
                    confidence = 0.75
                    evidence.append(f"检测到D型手性标记: {len(d_matches)}个")
                    chiral_centers = [{'type': 'D', 'count': len(d_matches)}]
                else:
                    # 混合手性
                    stereochemistry = 'mixed_chirality'
                    confidence = 0.70
                    evidence.append("检测到混合手性标记")
                    chiral_centers.append({'type': 'D', 'count': len(d_matches)})
        
        if confidence > 0:
            return {
                'has_chirality': True,
                'stereochemistry': stereochemistry,
                'confidence': confidence,
                'chiral_centers': chiral_centers,
                'evidence': evidence,
                'details': {
                    'analysis_method': 'smiles_chirality_markers',
                    'stereochemical_type': stereochemistry.replace('_suspected', '')
                }
            }
        
        return {'has_chirality': False}
    
    def _create_unknown_result(self, reason: str) -> Dict[str, Any]:
        """创建未知结果"""
        return {
            'stereochemistry': 'unknown',
            'confidence': 0.0,
            'chiral_centers': [],
            'evidence': [reason],
            'details': {
                'analysis_method': 'none',
                'reason': reason
            }
        }
    
    def compare_stereoisomers(self, smiles1: str, smiles2: str) -> Dict[str, Any]:
        """
        比较两个立体异构体
        
        Args:
            smiles1, smiles2: 要比较的SMILES
            
        Returns:
            比较结果
        """
        analysis1 = self.analyze(smiles1)
        analysis2 = self.analyze(smiles2)
        
        stereo1 = analysis1['stereochemistry']
        stereo2 = analysis2['stereochemistry']
        
        relationship = 'unknown'
        
        if stereo1 == stereo2:
            relationship = 'same_stereochemistry'
        elif ((stereo1 == 'D_suspected' and stereo2 == 'L_suspected') or
              (stereo1 == 'L_suspected' and stereo2 == 'D_suspected')):
            relationship = 'enantiomers'
        elif 'mixed' in stereo1 or 'mixed' in stereo2:
            relationship = 'complex_stereoisomers'
        
        return {
            'relationship': relationship,
            'structure1_stereochemistry': stereo1,
            'structure2_stereochemistry': stereo2,
            'confidence': min(analysis1['confidence'], analysis2['confidence']),
            'details': {
                'analysis1': analysis1,
                'analysis2': analysis2
            }
        }
    
    def identify_chiral_centers(self, smiles: str) -> List[Dict[str, Any]]:
        """
        识别手性中心
        
        Args:
            smiles: SMILES字符串
            
        Returns:
            手性中心列表
        """
        import re
        
        chiral_centers = []
        
        # 查找所有手性标记
        l_centers = re.finditer(r'C@@H?', smiles)
        for match in l_centers:
            chiral_centers.append({
                'position': match.start(),
                'type': 'L',
                'marker': match.group(),
                'confidence': 0.80
            })
        
        d_centers = re.finditer(r'C@H(?!@)', smiles)
        for match in d_centers:
            chiral_centers.append({
                'position': match.start(),
                'type': 'D', 
                'marker': match.group(),
                'confidence': 0.80
            })
        
        return chiral_centers
    
    def validate_stereochemistry(self, smiles: str, expected_type: str) -> Dict[str, Any]:
        """
        验证立体化学分配
        
        Args:
            smiles: SMILES字符串
            expected_type: 预期的立体化学类型
            
        Returns:
            验证结果
        """
        analysis = self.analyze(smiles)
        detected_type = analysis['stereochemistry']
        
        # 简化的匹配逻辑
        is_consistent = (
            expected_type.lower() in detected_type.lower() or
            detected_type == 'unknown' or
            detected_type == 'achiral_or_unknown'
        )
        
        return {
            'is_consistent': is_consistent,
            'detected_type': detected_type,
            'expected_type': expected_type,
            'confidence': analysis['confidence'],
            'evidence': analysis['evidence'],
            'recommendation': 'accept' if is_consistent else 'review_stereochemistry'
        }
    
    def get_stereochemical_statistics(self, amino_acids: List[Tuple[str, str]]) -> Dict[str, Any]:
        """
        获取立体化学统计
        
        Args:
            amino_acids: [(smiles, name), ...]
            
        Returns:
            统计结果
        """
        stats = {
            'total': len(amino_acids),
            'D_form': 0,
            'L_form': 0,
            'mixed': 0,
            'achiral': 0,
            'unknown': 0
        }
        
        for smiles, name in amino_acids:
            analysis = self.analyze(smiles)
            stereo_type = analysis['stereochemistry']
            
            if 'D' in stereo_type:
                stats['D_form'] += 1
            elif 'L' in stereo_type:
                stats['L_form'] += 1
            elif 'mixed' in stereo_type:
                stats['mixed'] += 1
            elif 'achiral' in stereo_type:
                stats['achiral'] += 1
            else:
                stats['unknown'] += 1
        
        # 计算百分比
        if stats['total'] > 0:
            for key in ['D_form', 'L_form', 'mixed', 'achiral', 'unknown']:
                stats[f'{key}_percentage'] = (stats[key] / stats['total']) * 100
        
        return stats
