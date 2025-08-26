"""
增强特征识别器
专门检测芳香族、卤素、官能团等结构特征
"""

import re
from typing import Dict, List, Any, Set, Optional
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors


class EnhancedFeatureDetector:
    """
    增强特征识别器
    
    主要功能：
    1. 芳香族特征识别
    2. 卤素取代基检测
    3. 官能团识别
    4. 结构特征分类
    """
    
    def __init__(self):
        """初始化增强特征识别器"""
        # 芳香族模式
        self.aromatic_patterns = [
            r'c1c+c+c+c+c+1',     # 苯环
            r'c1c+c+c+n+1',       # 吡啶环
            r'c1c+c+n+n+1',       # 咪唑环 
            r'c1c+n+c+c+1',       # 吡嗪环
            r'c1n+c+c+c+1',       # 吡咯环
            r'c1[nH]c+c+c+1',     # 吲哚环
            r'c1c+c+c+s+1',       # 噻吩环
            r'c1c+c+o+c+1',       # 呋喃环
        ]
        
        # 卤素原子
        self.halogens = ['F', 'Cl', 'Br', 'I']
        
        # 官能团SMARTS模式
        self.functional_groups = {
            'hydroxyl': ['[OH]', '[CX4][OH]'],
            'carbonyl': ['[CX3]=[OX1]'],
            'carboxyl': ['C(=O)[OH]', '[CX3](=[OX1])[OH]'],
            'amino': ['[NX3][H2]', '[NX3][H][H]', '[NH3+]'],
            'amide': ['[NX3][CX3]=[OX1]'],
            'ether': ['[OX2][CX4][CX4]'],
            'thiol': ['[SH]', '[CX4][SH]'],
            'sulfide': ['[SX2][CX4]'],
            'nitrile': ['[CX2]#[NX1]', 'C#N'],
            'nitro': ['[NX3+]([OX1-])[OX1-]', '[N+](=O)[O-]'],
            'aldehyde': ['[CX3H1](=O)'],
            'ketone': ['[CX3](=[OX1])[CX4]'],
            'ester': ['[CX3](=[OX1])[OX2][CX4]'],
            'phosphate': ['[PX4](=[OX1])([OX2])([OX2])'],
            'halogen': ['[F,Cl,Br,I]']
        }
        
        # 结构分类
        self.structural_classifications = {
            'aromatic': ['含有芳香环结构'],
            'halogen_substituted': ['含有卤素取代基'],
            'hydroxyl_substituted': ['含有羟基取代基'],
            'sulfur_containing': ['含有硫原子'],
            'nitrogen_heterocycle': ['含有氮杂环'],
            'phosphorus_containing': ['含有磷原子'],
            'cyclic': ['含有环状结构'],
            'branched': ['含有支链结构']
        }
    
    def detect_features(self, smiles: str, amino_acid_id: Optional[str] = None) -> Dict[str, Any]:
        """
        检测分子的结构特征
        
        Args:
            smiles: SMILES字符串
            amino_acid_id: 氨基酸ID（可选）
            
        Returns:
            特征检测结果
        """
        if not smiles:
            return self._create_empty_result("空SMILES字符串")
        
        try:
            # 使用RDKit进行详细分析
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                rdkit_features = self._analyze_with_rdkit(mol, smiles)
                if rdkit_features['feature_count'] > 0:
                    return rdkit_features
            
            # RDKit失败时，使用字符串模式匹配
            pattern_features = self._analyze_with_patterns(smiles)
            return pattern_features
            
        except Exception as e:
            # 如果所有方法都失败，至少尝试基本的字符串检查
            return self._basic_string_analysis(smiles)
    
    def _analyze_with_rdkit(self, mol, smiles: str) -> Dict[str, Any]:
        """使用RDKit进行详细特征分析"""
        detected_features = []
        feature_details = {}
        evidence = []
        
        # 1. 芳香族特征检测
        aromatic_info = self._detect_aromatic_features(mol)
        if aromatic_info['is_aromatic']:
            detected_features.append('aromatic')
            feature_details['aromatic'] = aromatic_info
            evidence.extend(aromatic_info['evidence'])
        
        # 2. 卤素检测
        halogen_info = self._detect_halogens(mol)
        if halogen_info['has_halogens']:
            detected_features.append('halogen_substituted')
            feature_details['halogens'] = halogen_info
            evidence.extend(halogen_info['evidence'])
        
        # 3. 官能团检测
        functional_groups = self._detect_functional_groups(mol, smiles)
        detected_features.extend(functional_groups['groups'])
        feature_details['functional_groups'] = functional_groups
        evidence.extend(functional_groups['evidence'])
        
        # 4. 环状结构检测
        ring_info = self._detect_ring_systems(mol)
        if ring_info['has_rings']:
            detected_features.append('cyclic')
            feature_details['rings'] = ring_info
            evidence.extend(ring_info['evidence'])
        
        # 5. 杂原子检测
        heteroatom_info = self._detect_heteroatoms(mol)
        detected_features.extend(heteroatom_info['features'])
        feature_details['heteroatoms'] = heteroatom_info
        evidence.extend(heteroatom_info['evidence'])
        
        # 6. 计算置信度
        confidence = self._calculate_feature_confidence(detected_features, evidence)
        
        return {
            'detected_features': list(set(detected_features)),  # 去重
            'feature_count': len(set(detected_features)),
            'feature_details': feature_details,
            'confidence': confidence,
            'method': 'rdkit_analysis',
            'evidence': evidence
        }
    
    def _detect_aromatic_features(self, mol) -> Dict[str, Any]:
        """检测芳香族特征"""
        aromatic_atoms = []
        aromatic_rings = []
        
        # 查找芳香原子
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                aromatic_atoms.append(atom.GetIdx())
        
        # 分析环系统
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                aromatic_rings.append(ring)
        
        is_aromatic = len(aromatic_atoms) > 0
        evidence = []
        
        if is_aromatic:
            evidence.append(f"检测到{len(aromatic_atoms)}个芳香原子")
            evidence.append(f"检测到{len(aromatic_rings)}个芳香环")
        
        return {
            'is_aromatic': is_aromatic,
            'aromatic_atom_count': len(aromatic_atoms),
            'aromatic_ring_count': len(aromatic_rings),
            'evidence': evidence
        }
    
    def _detect_halogens(self, mol) -> Dict[str, Any]:
        """检测卤素原子"""
        halogen_atoms = {}
        total_halogens = 0
        
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol in self.halogens:
                halogen_atoms[symbol] = halogen_atoms.get(symbol, 0) + 1
                total_halogens += 1
        
        evidence = []
        if total_halogens > 0:
            for halogen, count in halogen_atoms.items():
                evidence.append(f"检测到{count}个{halogen}原子")
        
        return {
            'has_halogens': total_halogens > 0,
            'halogen_types': list(halogen_atoms.keys()),
            'halogen_counts': halogen_atoms,
            'total_count': total_halogens,
            'evidence': evidence
        }
    
    def _detect_functional_groups(self, mol, smiles: str) -> Dict[str, Any]:
        """检测官能团"""
        detected_groups = []
        group_details = {}
        evidence = []
        
        for group_name, patterns in self.functional_groups.items():
            for pattern in patterns:
                try:
                    if group_name == 'halogen':
                        continue  # 卤素已经单独处理
                    
                    # 使用SMARTS模式匹配
                    if pattern.startswith('['):
                        # SMARTS模式
                        smarts = Chem.MolFromSmarts(pattern)
                        if smarts and mol.HasSubstructMatch(smarts):
                            detected_groups.append(group_name)
                            matches = mol.GetSubstructMatches(smarts)
                            group_details[group_name] = {
                                'pattern': pattern,
                                'match_count': len(matches)
                            }
                            evidence.append(f"检测到{group_name}: {len(matches)}个匹配")
                            break
                    else:
                        # 字符串模式
                        if pattern in smiles:
                            detected_groups.append(group_name)
                            group_details[group_name] = {
                                'pattern': pattern,
                                'match_count': smiles.count(pattern)
                            }
                            evidence.append(f"检测到{group_name}: {pattern}")
                            break
                            
                except Exception:
                    continue
        
        return {
            'groups': detected_groups,
            'group_details': group_details,
            'evidence': evidence
        }
    
    def _detect_ring_systems(self, mol) -> Dict[str, Any]:
        """检测环系统"""
        ring_info = mol.GetRingInfo()
        ring_count = ring_info.NumRings()
        ring_sizes = []
        
        for ring in ring_info.AtomRings():
            ring_sizes.append(len(ring))
        
        evidence = []
        if ring_count > 0:
            evidence.append(f"检测到{ring_count}个环")
            if ring_sizes:
                evidence.append(f"环大小: {', '.join(map(str, sorted(ring_sizes)))}")
        
        return {
            'has_rings': ring_count > 0,
            'ring_count': ring_count,
            'ring_sizes': ring_sizes,
            'evidence': evidence
        }
    
    def _detect_heteroatoms(self, mol) -> Dict[str, Any]:
        """检测杂原子"""
        heteroatom_counts = {}
        features = []
        evidence = []
        
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            if symbol not in ['C', 'H']:  # 非C、H原子
                heteroatom_counts[symbol] = heteroatom_counts.get(symbol, 0) + 1
        
        # 基于杂原子添加特征
        if 'S' in heteroatom_counts:
            features.append('sulfur_containing')
            evidence.append(f"含有{heteroatom_counts['S']}个硫原子")
        
        if 'P' in heteroatom_counts:
            features.append('phosphorus_containing')
            evidence.append(f"含有{heteroatom_counts['P']}个磷原子")
        
        if 'N' in heteroatom_counts and heteroatom_counts['N'] > 1:  # 除了氨基氮外的氮
            # 检查是否在环中
            ring_info = mol.GetRingInfo()
            nitrogen_in_ring = False
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == 'N' and atom.IsInRing():
                    nitrogen_in_ring = True
                    break
            
            if nitrogen_in_ring:
                features.append('nitrogen_heterocycle')
                evidence.append("含有氮杂环")
        
        return {
            'heteroatom_counts': heteroatom_counts,
            'features': features,
            'evidence': evidence
        }
    
    def _analyze_with_patterns(self, smiles: str) -> Dict[str, Any]:
        """使用字符串模式进行特征分析"""
        detected_features = []
        evidence = []
        
        # 检查芳香族模式
        aromatic_found = False
        for pattern in self.aromatic_patterns:
            if re.search(pattern, smiles):
                aromatic_found = True
                break
        
        if aromatic_found or any(c.islower() for c in smiles):  # 小写字母表示芳香
            detected_features.append('aromatic')
            evidence.append("检测到芳香族模式")
        
        # 检查卤素
        halogen_found = []
        for halogen in self.halogens:
            if halogen in smiles:
                halogen_found.append(halogen)
        
        if halogen_found:
            detected_features.append('halogen_substituted')
            evidence.append(f"检测到卤素: {', '.join(halogen_found)}")
        
        # 检查基本官能团
        if 'OH' in smiles or ')O' in smiles:
            detected_features.append('hydroxyl_substituted')
            evidence.append("检测到羟基")
        
        if 'S' in smiles:
            detected_features.append('sulfur_containing')
            evidence.append("检测到硫原子")
        
        if 'C#N' in smiles or 'N#C' in smiles:
            detected_features.append('nitrile_containing')
            evidence.append("检测到氰基")
        
        if 'P' in smiles:
            detected_features.append('phosphorus_containing')
            evidence.append("检测到磷原子")
        
        # 检查环状结构（通过数字配对）
        if re.search(r'\d', smiles):
            detected_features.append('cyclic')
            evidence.append("检测到环状结构标记")
        
        confidence = min(0.8, len(detected_features) * 0.2 + 0.4)
        
        return {
            'detected_features': detected_features,
            'feature_count': len(detected_features),
            'confidence': confidence,
            'method': 'pattern_analysis',
            'evidence': evidence
        }
    
    def _basic_string_analysis(self, smiles: str) -> Dict[str, Any]:
        """基本的字符串分析（最后的回退方法）"""
        detected_features = []
        evidence = []
        
        # 最基本的特征检查
        if any(c.islower() for c in smiles):
            detected_features.append('aromatic')
            evidence.append("包含小写字母（芳香标记）")
        
        if any(hal in smiles for hal in ['F', 'Cl', 'Br', 'I']):
            detected_features.append('halogen_substituted')
            evidence.append("包含卤素原子")
        
        if 'O' in smiles:
            detected_features.append('oxygen_containing')
            evidence.append("包含氧原子")
        
        return {
            'detected_features': detected_features,
            'feature_count': len(detected_features),
            'confidence': 0.3,  # 低置信度
            'method': 'basic_string_analysis',
            'evidence': evidence
        }
    
    def _calculate_feature_confidence(self, features: List[str], evidence: List[str]) -> float:
        """计算特征检测的置信度"""
        base_confidence = 0.5
        
        # 基于特征数量的置信度
        feature_bonus = min(0.3, len(features) * 0.1)
        
        # 基于证据数量的置信度
        evidence_bonus = min(0.2, len(evidence) * 0.05)
        
        return min(1.0, base_confidence + feature_bonus + evidence_bonus)
    
    def _create_empty_result(self, reason: str) -> Dict[str, Any]:
        """创建空的特征检测结果"""
        return {
            'detected_features': [],
            'feature_count': 0,
            'confidence': 0.0,
            'method': 'failed_analysis',
            'evidence': [f"特征检测失败: {reason}"]
        }


def create_enhanced_feature_detector():
    """创建增强特征识别器实例"""
    return EnhancedFeatureDetector()