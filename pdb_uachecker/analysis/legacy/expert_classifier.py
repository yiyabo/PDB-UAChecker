"""
符合专家修正标准的氨基酸分类器
基于修订版分类原则重新设计
"""

import re
from typing import Dict, Any, List, Tuple, Optional
from collections import defaultdict

from ...core.models import AminoAcidInfo, ResidueInfo


class ExpertAminoAcidClassifier:
    """专家标准氨基酸分类器"""
    
    def __init__(self):
        """初始化专家分类器"""
        
        # 标准氨基酸
        self.standard_amino_acids = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        }
        
        # PDB格式特殊标识
        self.pdb_identifiers = {
            'd_form_prefixes': ['D'],
            'beta_indicators': ['B3', '3'],
            'gamma_indicators': ['G4', '4'],
            'n_methyl_codes': ['SAR', 'NME', 'NMA']
        }

    def classify_amino_acid(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """
        根据专家修正标准分类氨基酸
        
        Args:
            amino_acid: 氨基酸信息
            
        Returns:
            详细分类结果
        """
        
        # 1. 检查标准氨基酸
        if amino_acid.id in self.standard_amino_acids:
            return self._create_standard_result(amino_acid)
        
        if not amino_acid.smiles:
            return self._create_unclassified_result(amino_acid, 'no_smiles')
        
        # 2. 确定骨架类型（最重要的分类）
        backbone_type = self._identify_backbone_type(amino_acid)
        
        # 3. 分析分子修饰
        modifications = self._analyze_modifications(amino_acid)
        
        # 4. 综合分类
        classification = self._build_classification(backbone_type, modifications)
        
        # 5. 验证PDB标识符一致性
        pdb_validation = self._validate_pdb_identifiers(amino_acid, classification)
        
        return {
            'amino_acid_id': amino_acid.id,
            'amino_acid_name': amino_acid.name,
            'smiles': amino_acid.smiles,
            'backbone_type': backbone_type,
            'modifications': modifications,
            'primary_categories': classification['primary'],
            'secondary_features': classification['secondary'],
            'confidence': classification['confidence'],
            'expert_validation': pdb_validation,
            'classification_details': {
                'method': 'expert_standard',
                'cip_analysis': modifications.get('chirality_details', {}),
                'backbone_analysis': backbone_type
            }
        }

    def _identify_backbone_type(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """识别氨基酸骨架类型（改进版）"""
        smiles = amino_acid.smiles
        
        # 标准α-氨基酸模式检测（最常见）
        alpha_patterns = [
            r'N\[C@@?H?\]\([^)]*\)C\(=O\)O',  # N[C@H](R)C(=O)O
            r'N\[C@@?H?\]\([^)]*C\(=O\)O\)',  # N[C@H](RC(=O)O)
            r'\[NH2\]\[C@@?H?\].*C\(=O\)O',   # [NH2][C@H]...C(=O)O
            r'N\[C@@?H?\].*C\(=O\)O',        # 通用α模式
        ]
        
        for pattern in alpha_patterns:
            if re.search(pattern, smiles, re.IGNORECASE):
                return {'type': 'alpha', 'confidence': 0.95, 'amino_position': 'alpha'}
        
        # Beta氨基酸模式：NH2-CH2-CH(R)-COOH
        beta_patterns = [
            r'N.*C.*\[C@@?H?\].*C\(=O\)O',  # N-C-[C@H]-COOH
            r'\[NH2\]C.*C.*C\(=O\)O',       # [NH2]C-C-COOH
        ]
        
        for pattern in beta_patterns:
            if re.search(pattern, smiles, re.IGNORECASE):
                return {'type': 'beta', 'confidence': 0.90, 'amino_position': 'beta'}
        
        # Gamma氨基酸模式：NH2-CH2-CH2-CH(R)-COOH  
        gamma_patterns = [
            r'N.*C.*C.*C.*C\(=O\)O',        # N-C-C-C-COOH
            r'\[NH2\]CC.*C.*C\(=O\)O',      # [NH2]CC-C-COOH
        ]
        
        for pattern in gamma_patterns:
            if re.search(pattern, smiles, re.IGNORECASE):
                return {'type': 'gamma', 'confidence': 0.90, 'amino_position': 'gamma'}
        
        # 如果有氨基和羧基但不匹配标准模式
        has_amino = bool(re.search(r'N(?!\d)|NH2|NH3', smiles, re.IGNORECASE))
        has_carboxyl = bool(re.search(r'C\(=O\)O', smiles, re.IGNORECASE))
        
        if has_amino and has_carboxyl:
            return {'type': 'alpha', 'confidence': 0.7, 'note': 'default_alpha_assumption'}
        
        return {'type': 'unknown', 'confidence': 0.0}

    def _analyze_carbon_chain(self, smiles: str, carboxyl_match) -> Dict[str, Any]:
        """分析从羧基开始的碳链结构"""
        # 简化版本 - 实际需要更复杂的SMILES解析
        carboxyl_pos = carboxyl_match.start()
        
        # 向前查找连接的碳原子
        alpha_carbon_pattern = r'(?:\[C@@?H?\]|C)(?=\(=O\)O|.*\(=O\)O)'
        
        return {
            'carboxyl_position': carboxyl_pos,
            'alpha_carbon': self._find_adjacent_carbon(smiles, carboxyl_pos),
            'chain_length': self._estimate_chain_length(smiles)
        }

    def _find_amino_positions(self, smiles: str, backbone_analysis: Dict) -> Dict[str, bool]:
        """找到氨基在主链中的位置"""
        
        # 识别氨基基团模式
        amino_patterns = [
            r'N(?!\[)',  # 简单氨基
            r'\[NH2\]',  # 明确氨基
            r'\[NH3\+\]',  # 质子化氨基
            r'N\[',  # 其他氨基形式
        ]
        
        amino_positions = []
        for pattern in amino_patterns:
            matches = list(re.finditer(pattern, smiles))
            amino_positions.extend([m.start() for m in matches])
        
        # 分析每个氨基相对于羧基的位置
        # 这里是简化版本，实际需要更精确的分子结构分析
        alpha_amino = self._check_alpha_amino(smiles, amino_positions)
        beta_amino = self._check_beta_amino(smiles, amino_positions)
        gamma_amino = self._check_gamma_amino(smiles, amino_positions)
        
        return {
            'alpha': alpha_amino,
            'beta': beta_amino,
            'gamma': gamma_amino,
            'positions': amino_positions
        }

    def _check_alpha_amino(self, smiles: str, amino_positions: List[int]) -> bool:
        """检查是否存在α位氨基"""
        # 查找典型的α-氨基酸模式
        alpha_patterns = [
            r'N\[C@@?H?\].*C\(=O\)O',  # N[C@H]...C(=O)O
            r'\[NH2\]\[C@@?H?\].*C\(=O\)O',  # [NH2][C@H]...C(=O)O
            r'N.*\[C@@?H?\].*C\(=O\)O',  # N...[C@H]...C(=O)O
        ]
        
        for pattern in alpha_patterns:
            if re.search(pattern, smiles):
                return True
        return False

    def _check_beta_amino(self, smiles: str, amino_positions: List[int]) -> bool:
        """检查是否存在β位氨基"""
        # β-氨基酸模式：NH2-CH2-CH(R)-COOH
        beta_patterns = [
            r'N.*C.*\[C@@?H?\].*C\(=O\)O',  # N-C-[C@H]-COOH
            r'\[NH2\]C.*\[C@@?H?\].*C\(=O\)O',  # [NH2]C-[C@H]-COOH
        ]
        
        for pattern in beta_patterns:
            if re.search(pattern, smiles):
                return True
        return False

    def _check_gamma_amino(self, smiles: str, amino_positions: List[int]) -> bool:
        """检查是否存在γ位氨基"""
        # γ-氨基酸模式：NH2-CH2-CH2-CH(R)-COOH
        gamma_patterns = [
            r'N.*C.*C.*\[C@@?H?\].*C\(=O\)O',  # N-C-C-[C@H]-COOH
        ]
        
        for pattern in gamma_patterns:
            if re.search(pattern, smiles):
                return True
        return False

    def _analyze_modifications(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """分析分子修饰"""
        modifications = {}
        
        # 1. D/L型分析（使用CIP规则）
        modifications['chirality'] = self._analyze_chirality_cip(amino_acid)
        
        # 2. N-甲基化检测（主链骨架）
        modifications['n_methylation'] = self._check_backbone_n_methyl(amino_acid)
        
        # 3. 环化结构检测
        modifications['cyclization'] = self._check_cyclic_structure(amino_acid)
        
        # 4. 芳香性检测
        modifications['aromaticity'] = self._check_aromatic_sidechain(amino_acid)
        
        return modifications

    def _analyze_chirality_cip(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """使用CIP规则分析手性"""
        smiles = amino_acid.smiles
        
        # 查找手性中心
        chiral_patterns = [
            r'\[C@H\]',   # R构型标记
            r'\[C@@H\]',  # S构型标记
        ]
        
        chirality_info = {
            'has_chirality': False,
            'stereo_centers': [],
            'd_l_form': 'unknown',
            'confidence': 0.0,
            'cip_analysis': {}
        }
        
        for pattern in chiral_patterns:
            matches = list(re.finditer(pattern, smiles))
            for match in matches:
                center_info = self._analyze_chiral_center_cip(smiles, match)
                chirality_info['stereo_centers'].append(center_info)
        
        if chirality_info['stereo_centers']:
            chirality_info['has_chirality'] = True
            # 基于CIP分析确定D/L型
            primary_center = chirality_info['stereo_centers'][0]  # 取主要手性中心
            
            if primary_center['cip_configuration'] == 'R':
                chirality_info['d_l_form'] = 'D'
                chirality_info['confidence'] = 0.85  # 需要验证CIP分析
            elif primary_center['cip_configuration'] == 'S':
                chirality_info['d_l_form'] = 'L'
                chirality_info['confidence'] = 0.85
        
        return chirality_info

    def _analyze_chiral_center_cip(self, smiles: str, match) -> Dict[str, Any]:
        """分析单个手性中心的CIP配置"""
        # 这是简化版本，实际需要完整的CIP优先级计算
        stereo_symbol = match.group()
        
        # 基本CIP分析（需要完整实现）
        cip_config = 'R' if '@' in stereo_symbol and '@@' not in stereo_symbol else 'S'
        
        return {
            'position': match.start(),
            'stereo_symbol': stereo_symbol,
            'cip_configuration': cip_config,
            'confidence': 0.7,  # 简化版本置信度较低
            'note': 'simplified_cip_analysis'
        }

    def _check_backbone_n_methyl(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """检查主链骨架N-甲基化"""
        smiles = amino_acid.smiles.upper()
        
        # 主链N-甲基化模式
        n_methyl_patterns = [
            r'CN\[C@@?H?\].*C\(=O\)O',  # CN[C@H]...COOH (α位)
            r'CNC.*\[C@@?H?\].*C\(=O\)O',  # CNC-[C@H]-COOH (β位)
            r'CNCC.*\[C@@?H?\].*C\(=O\)O',  # CNCC-[C@H]-COOH (γ位)
        ]
        
        for pattern in n_methyl_patterns:
            if re.search(pattern, smiles):
                return {
                    'is_n_methylated': True,
                    'confidence': 0.9,
                    'type': 'backbone_methylation'
                }
        
        return {'is_n_methylated': False}

    def _check_cyclic_structure(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """检查环化结构"""
        smiles = amino_acid.smiles
        
        # 检查环状标记
        ring_numbers = re.findall(r'\d', smiles)
        has_rings = len(set(ring_numbers)) > 0
        
        return {
            'is_cyclic': has_rings,
            'ring_count': len(set(ring_numbers)) // 2 if has_rings else 0,
            'confidence': 0.95 if has_rings else 1.0
        }

    def _check_aromatic_sidechain(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """检查芳香性侧链"""
        smiles = amino_acid.smiles.lower()
        
        # 精确芳香性检测
        aromatic_patterns = [
            r'c1ccccc1',        # 苯环
            r'c1cccc[nos]1',    # 五元芳香杂环
            r'c1ccc[nos]c1',    # 六元芳香杂环
        ]
        
        for pattern in aromatic_patterns:
            if re.search(pattern, smiles):
                return {
                    'is_aromatic': True,
                    'confidence': 0.95,
                    'aromatic_system': pattern
                }
        
        return {'is_aromatic': False}

    def _build_classification(self, backbone_type: Dict, modifications: Dict) -> Dict[str, Any]:
        """构建综合分类结果"""
        primary_categories = []
        secondary_features = []
        confidence_scores = []
        
        # 1. 骨架类型分类
        if backbone_type['type'] == 'beta':
            primary_categories.append('Beta-amino acid')
            confidence_scores.append(backbone_type['confidence'])
        elif backbone_type['type'] == 'gamma':
            primary_categories.append('Gamma-amino acid')
            confidence_scores.append(backbone_type['confidence'])
        
        # 2. 手性分类
        if modifications['chirality']['d_l_form'] == 'D':
            primary_categories.append('D-amino acid')
            confidence_scores.append(modifications['chirality']['confidence'])
        elif modifications['chirality']['d_l_form'] == 'L':
            # L型是标准构型，通常不特别标注，但在分析中记录
            secondary_features.append('L-form (standard chirality)')
        
        # 3. N-甲基化
        if modifications['n_methylation'].get('is_n_methylated'):
            primary_categories.append('N-methyl amino acid')
            confidence_scores.append(modifications['n_methylation']['confidence'])
        
        # 4. 次要特征
        if modifications['cyclization'].get('is_cyclic'):
            secondary_features.append('Cyclic amino acid')
        if modifications['aromaticity'].get('is_aromatic'):
            secondary_features.append('Aromatic amino acid')
        
        # 5. 如果没有主要分类，根据骨架类型提供默认分类
        if not primary_categories:
            if backbone_type['type'] == 'alpha':
                primary_categories.append('Non-standard alpha-amino acid')
                confidence_scores.append(0.6)
            elif backbone_type['type'] in ['complex', 'unknown']:
                primary_categories.append('Requires detailed analysis')
                confidence_scores.append(0.3)
            else:
                primary_categories.append('Unclassified amino acid')
                confidence_scores.append(0.2)
        
        # 计算总体置信度
        overall_confidence = sum(confidence_scores) / len(confidence_scores) if confidence_scores else 0.5
        
        return {
            'primary': primary_categories,
            'secondary': secondary_features,
            'confidence': min(overall_confidence, 0.95)  # 专家分类最高95%
        }

    def _validate_pdb_identifiers(self, amino_acid: AminoAcidInfo, classification: Dict) -> Dict[str, Any]:
        """验证PDB标识符一致性"""
        pdb_id = amino_acid.id
        validation = {'consistent': True, 'warnings': [], 'confidence_boost': 0.0}
        
        # D型验证
        if 'D-amino acid' in classification['primary']:
            if pdb_id.startswith('D'):
                validation['confidence_boost'] += 0.1
            else:
                validation['warnings'].append('D-form_predicted_but_no_D_prefix')
        
        # Beta/Gamma验证  
        if 'Beta-amino acid' in classification['primary']:
            if any(indicator in pdb_id for indicator in ['B3', '3']):
                validation['confidence_boost'] += 0.05
        
        if 'Gamma-amino acid' in classification['primary']:
            if any(indicator in pdb_id for indicator in ['G4', '4']):
                validation['confidence_boost'] += 0.05
        
        return validation

    def _create_standard_result(self, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """创建标准氨基酸结果"""
        return {
            'amino_acid_id': amino_acid.id,
            'amino_acid_name': amino_acid.name,
            'primary_categories': ['Standard amino acid'],
            'secondary_features': [],
            'confidence': 1.0,
            'method': 'standard_recognition'
        }

    def _create_unclassified_result(self, amino_acid: AminoAcidInfo, reason: str) -> Dict[str, Any]:
        """创建未分类结果"""
        return {
            'amino_acid_id': amino_acid.id,
            'amino_acid_name': amino_acid.name,
            'primary_categories': ['Requires expert review'],
            'confidence': 0.0,
            'unclassified_reason': reason
        }

    # 辅助方法（简化版本）
    def _find_adjacent_carbon(self, smiles: str, pos: int) -> Optional[int]:
        """找到相邻碳原子"""
        # 简化实现
        return pos - 1 if pos > 0 else None

    def _estimate_chain_length(self, smiles: str) -> int:
        """估算碳链长度"""
        # 简化实现
        return smiles.count('C')