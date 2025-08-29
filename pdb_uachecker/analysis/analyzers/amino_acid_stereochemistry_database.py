"""
标准氨基酸立体化学参考数据库
用于修正立体化学分析中的映射错误，特别是R/S到D/L的转换
"""

import re
from typing import Dict, List, Optional, Set, Any

class AminoAcidStereochemistryDatabase:
    """标准氨基酸立体化学参考数据库"""
    
    def __init__(self):
        """初始化数据库"""
        # 标准20种氨基酸的立体化学信息
        self.standard_amino_acids = {
            # 标准L型氨基酸的SMILES模式和预期立体化学
            'PRO': {
                'standard_form': 'L_form',
                'smiles_patterns': [
                    r'\[NH2\]1\[C@H\]\(C\(=O\)O\)CCC1',     # PRO的标准模式
                    r'N1\[C@@H\]\(C\(=O\)O\)CCC1',           # 替代表示
                    r'\[NH2\]1[C@H]\(.*C\(=O\)O.*\)CCC1'     # 通用模式
                ],
                'expected_config': 'S',  # PRO的L型对应S构型（特例）
                'notes': 'PRO是特殊的氮杂环氨基酸，L构型对应S构型'
            },
            
            'PHE': {
                'standard_form': 'L_form',
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(Cc1ccccc1\)C\(=O\)O',
                    r'N\[C@@H\]\(Cc1ccccc1\)C\(=O\)O'
                ],
                'expected_config': 'R',  # 标准L型对应R构型
                'notes': '苯丙氨酸，标准L型氨基酸'
            },
            
            'HIS': {
                'standard_form': 'L_form', 
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(Cc1c\[nH\]cn1\)C\(=O\)O',
                    r'N\[C@@H\]\(Cc1c\[nH\]cn1\)C\(=O\)O'
                ],
                'expected_config': 'R',
                'notes': '组氨酸，标准L型氨基酸'
            },
            
            'TRP': {
                'standard_form': 'L_form',
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(Cc1c\[nH\]c2c1cccc2\)C\(=O\)O',
                    r'N\[C@@H\]\(Cc1c\[nH\]c2c1cccc2\)C\(=O\)O'
                ],
                'expected_config': 'R',
                'notes': '色氨酸，标准L型氨基酸'
            },
            
            'CYS': {
                'standard_form': 'L_form',
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(CS\)C\(=O\)O',
                    r'N\[C@@H\]\(CS\)C\(=O\)O'
                ],
                'expected_config': 'R',
                'notes': '半胱氨酸，标准L型氨基酸'
            },
            
            'VAL': {
                'standard_form': 'L_form',
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(C\(C\)C\)C\(=O\)O',
                    r'N\[C@@H\]\(C\(C\)C\)C\(=O\)O'
                ],
                'expected_config': 'R',
                'notes': '缬氨酸，标准L型氨基酸'
            },
            
            'LEU': {
                'standard_form': 'L_form',
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(CC\(C\)C\)C\(=O\)O',
                    r'N\[C@@H\]\(CC\(C\)C\)C\(=O\)O'
                ],
                'expected_config': 'R',
                'notes': '亮氨酸，标准L型氨基酸'
            },
            
            'ASP': {
                'standard_form': 'L_form',
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(CC\(=O\)O\)C\(=O\)O',
                    r'N\[C@@H\]\(CC\(=O\)O\)C\(=O\)O'
                ],
                'expected_config': 'R',
                'notes': '天冬氨酸，标准L型氨基酸'
            },
            
            'GLU': {
                'standard_form': 'L_form',
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(CCC\(=O\)O\)C\(=O\)O',
                    r'N\[C@@H\]\(CCC\(=O\)O\)C\(=O\)O'
                ],
                'expected_config': 'R',
                'notes': '谷氨酸，标准L型氨基酸'
            },
            
            'ARG': {
                'standard_form': 'L_form',
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(CCCNC\(=\[NH2\]\)N\)C\(=O\)O',
                    r'N\[C@@H\]\(CCCNC.*N.*\)C\(=O\)O'
                ],
                'expected_config': 'R',
                'notes': '精氨酸，标准L型氨基酸'
            },
            
            'MET': {
                'standard_form': 'L_form', 
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(CCSC\)C\(=O\)O',
                    r'N\[C@@H\]\(CCSC\)C\(=O\)O'
                ],
                'expected_config': 'R',
                'notes': '蛋氨酸，标准L型氨基酸'
            },
            
            'GLN': {
                'standard_form': 'L_form',
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(CCC\(=O\)N\)C\(=O\)O',
                    r'N\[C@@H\]\(CCC\(=O\)N\)C\(=O\)O'
                ],
                'expected_config': 'R',
                'notes': '谷氨酰胺，标准L型氨基酸'
            },
            
            'ASN': {
                'standard_form': 'L_form',
                'smiles_patterns': [
                    r'\[NH3\]\[C@@H\]\(CC\(=O\)N\)C\(=O\)O',
                    r'N\[C@@H\]\(CC\(=O\)N\)C\(=O\)O'
                ],
                'expected_config': 'R',
                'notes': '天冬酰胺，标准L型氨基酸'
            }
        }
        
        # D型氨基酸对应的SMILES模式
        self.d_type_patterns = {
            'general_d_markers': [
                r'\[C@H\]',        # S构型标记（一般对应D型）
                r'C@H',            # 简化S构型标记
            ],
            'specific_d_patterns': [
                r'\[NH3\]\[C@H\].*C\(=O\)O',     # 标准D型氨基酸模式
                r'N\[C@H\].*C\(=O\)O',           # 清理后D型模式
                r'OC\(=O\)\[C@H\]'               # 反向D型模式
            ]
        }
    
    def identify_stereochemistry_by_database(self, smiles: str, amino_acid_code: str = "") -> Dict[str, Any]:
        """
        基于数据库识别立体化学构型
        
        Args:
            smiles: SMILES字符串
            amino_acid_code: 氨基酸代码
            
        Returns:
            立体化学分析结果
        """
        if not smiles:
            return self._create_unknown_result("空SMILES")
        
        # 标准化SMILES用于匹配
        normalized_smiles = self._normalize_smiles_for_matching(smiles)
        
        # 优先检查标准氨基酸
        if amino_acid_code and amino_acid_code.upper() in self.standard_amino_acids:
            return self._check_standard_amino_acid(normalized_smiles, amino_acid_code.upper())
        
        # 通用立体化学识别
        return self._identify_general_stereochemistry(normalized_smiles, smiles)
    
    def _check_standard_amino_acid(self, normalized_smiles: str, code: str) -> Dict[str, Any]:
        """检查标准氨基酸的立体化学"""
        amino_acid_info = self.standard_amino_acids[code]
        
        # 检查SMILES是否匹配标准氨基酸模式
        pattern_matched = False
        for pattern in amino_acid_info['smiles_patterns']:
            if re.search(pattern, normalized_smiles, re.IGNORECASE):
                pattern_matched = True
                break
        
        if pattern_matched:
            # 检查实际构型标记
            config_from_smiles = self._extract_configuration_from_smiles(normalized_smiles)
            expected_config = amino_acid_info['expected_config']
            
            # 对于特殊情况（如PRO），使用数据库预期值
            if code == 'PRO':
                # PRO的特殊处理：L型PRO对应S构型但仍然是L型氨基酸
                if '[C@H]' in normalized_smiles or 'C@H' in normalized_smiles:
                    return {
                        'success': True,
                        'stereochemistry': 'L_form',  # 强制L型（数据库校正）
                        'rs_configuration': 'S',
                        'confidence': 0.95,
                        'source': 'database_override',
                        'evidence': [f'{code}是标准L型氨基酸（数据库校正）', amino_acid_info['notes']],
                        'database_matched': True
                    }
            
            # 标准逻辑
            if config_from_smiles == expected_config:
                return {
                    'success': True,
                    'stereochemistry': amino_acid_info['standard_form'],
                    'rs_configuration': config_from_smiles,
                    'confidence': 0.95,
                    'source': 'database_match',
                    'evidence': [f'{code}匹配标准{amino_acid_info["standard_form"]}模式'],
                    'database_matched': True
                }
            else:
                # 构型不匹配，可能是D型
                opposite_form = 'D_form' if amino_acid_info['standard_form'] == 'L_form' else 'L_form'
                return {
                    'success': True,
                    'stereochemistry': opposite_form,
                    'rs_configuration': config_from_smiles,
                    'confidence': 0.90,
                    'source': 'database_inferred',
                    'evidence': [f'{code}的{opposite_form}构型（基于构型标记推断）'],
                    'database_matched': True
                }
        
        return self._create_unknown_result(f"未匹配到{code}的标准模式")
    
    def _identify_general_stereochemistry(self, normalized_smiles: str, original_smiles: str) -> Dict[str, Any]:
        """通用立体化学识别"""
        config = self._extract_configuration_from_smiles(normalized_smiles)
        
        if config == 'R':
            return {
                'success': True,
                'stereochemistry': 'L_form',
                'rs_configuration': 'R',
                'confidence': 0.80,
                'source': 'general_mapping',
                'evidence': ['R构型通常对应L型氨基酸'],
                'database_matched': False
            }
        elif config == 'S':
            return {
                'success': True,
                'stereochemistry': 'D_form',
                'rs_configuration': 'S', 
                'confidence': 0.80,
                'source': 'general_mapping',
                'evidence': ['S构型通常对应D型氨基酸'],
                'database_matched': False
            }
        else:
            return self._create_unknown_result("无法确定立体化学构型")
    
    def _extract_configuration_from_smiles(self, smiles: str) -> str:
        """从SMILES中提取R/S构型"""
        if '[C@@H]' in smiles or 'C@@H' in smiles:
            return 'R'
        elif '[C@H]' in smiles or 'C@H' in smiles:
            return 'S'
        else:
            return 'unknown'
    
    def _normalize_smiles_for_matching(self, smiles: str) -> str:
        """标准化SMILES以便模式匹配"""
        # 移除多余的空格和特殊字符
        normalized = smiles.strip()
        
        # 统一离子表示
        normalized = re.sub(r'\[NH3\+\]', '[NH3]', normalized)
        normalized = re.sub(r'\[O-\]', 'O', normalized)
        
        return normalized
    
    def _create_unknown_result(self, reason: str) -> Dict[str, Any]:
        """创建未知结果"""
        return {
            'success': True,
            'stereochemistry': 'unknown_chirality',
            'rs_configuration': 'unknown',
            'confidence': 0.0,
            'source': 'database_unknown',
            'evidence': [reason],
            'database_matched': False
        }
    
    def get_supported_amino_acids(self) -> List[str]:
        """获取支持的标准氨基酸列表"""
        return list(self.standard_amino_acids.keys())
    
    def add_custom_amino_acid(self, code: str, info: Dict[str, Any]):
        """添加自定义氨基酸到数据库"""
        self.standard_amino_acids[code] = info

# 创建全局实例
amino_acid_stereo_db = AminoAcidStereochemistryDatabase()