"""
分子结构分析器
使用RDKit进行精确的分子结构分析，替代正则表达式方法
"""

from typing import Dict, List, Optional, Tuple, Any
import re

from ..core.models import MolecularAnalysis, AminoAcidInfo


class MolecularStructureAnalyzer:
    """分子结构分析器 - 基于RDKit的化学分析"""
    
    def __init__(self):
        """初始化分析器"""
        self.rdkit_available = self._check_rdkit_availability()
        
        if self.rdkit_available:
            print("🧪 RDKit可用 - 启用高精度分子分析")
        else:
            print("⚠️ RDKit不可用 - 使用基础分析模式")
    
    def _check_rdkit_availability(self) -> bool:
        """检查RDKit是否可用"""
        try:
            import rdkit
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors, Descriptors
            return True
        except ImportError:
            return False
    
    def analyze(self, smiles: str) -> MolecularAnalysis:
        """
        分析分子结构
        
        Args:
            smiles: SMILES字符串
            
        Returns:
            分子分析结果
        """
        if not smiles:
            return MolecularAnalysis(smiles="", is_valid=False)
        
        if self.rdkit_available:
            return self._analyze_with_rdkit(smiles)
        else:
            return self._analyze_basic(smiles)
    
    def _analyze_with_rdkit(self, smiles: str) -> MolecularAnalysis:
        """使用RDKit进行完整分子分析"""
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors, Descriptors
            
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return MolecularAnalysis(smiles=smiles, is_valid=False)
            
            # 基本验证
            analysis = MolecularAnalysis(smiles=smiles, is_valid=True)
            
            # 1. 芳香性分析
            aromatic_atoms = []
            for i, atom in enumerate(mol.GetAtoms()):
                if atom.GetIsAromatic():
                    aromatic_atoms.append(i)
            analysis.aromatic_atoms = aromatic_atoms
            
            # 2. 环系统分析
            ring_info = mol.GetRingInfo()
            analysis.ring_systems = list(ring_info.AtomRings())
            
            # 3. 手性中心检测
            chiral_centers = []
            for center in Chem.FindMolChiralCenters(mol, includeUnassigned=True):
                atom_idx, chirality = center
                chiral_centers.append((atom_idx, chirality))
            analysis.chiral_centers = chiral_centers
            
            # 4. 官能团识别
            analysis.functional_groups = self._identify_functional_groups_rdkit(mol)
            
            # 5. 主链分析
            analysis.backbone_analysis = self._analyze_backbone_rdkit(mol)
            
            # 6. 分子描述符
            analysis.molecular_descriptors = {
                'molecular_weight': Descriptors.MolWt(mol),
                'logp': Descriptors.MolLogP(mol),
                'tpsa': rdMolDescriptors.CalcTPSA(mol),
                'num_rotatable_bonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
                'num_hbd': rdMolDescriptors.CalcNumHBD(mol),
                'num_hba': rdMolDescriptors.CalcNumHBA(mol),
            }
            
            return analysis
            
        except Exception as e:
            print(f"⚠️ RDKit分析失败: {e}")
            return MolecularAnalysis(smiles=smiles, is_valid=False)
    
    def _analyze_basic(self, smiles: str) -> MolecularAnalysis:
        """基础分析模式（不使用RDKit）"""
        analysis = MolecularAnalysis(smiles=smiles, is_valid=True)
        
        # 基本的芳香性检测（启发式）
        aromatic_patterns = ['c1', 'c2', 'n1', 's1', 'o1']
        if any(pattern in smiles.lower() for pattern in aromatic_patterns):
            # 粗略估计芳香原子（仅用于fallback）
            analysis.aromatic_atoms = [i for i, c in enumerate(smiles) if c.islower()]
        
        # 基本的环检测
        ring_numbers = re.findall(r'\d+', smiles)
        if ring_numbers:
            analysis.ring_systems = [list(range(len(ring_numbers)))]  # 简化表示
        
        # 基本的官能团识别
        analysis.functional_groups = self._identify_functional_groups_basic(smiles)
        
        # 基本的主链分析
        analysis.backbone_analysis = self._analyze_backbone_basic(smiles)
        
        return analysis
    
    def _identify_functional_groups_rdkit(self, mol) -> Dict[str, List[int]]:
        """使用RDKit识别官能团"""
        from rdkit import Chem
        
        functional_groups = {}
        
        # 预定义的官能团模式
        patterns = {
            'amino': '[NH2]',
            'carboxyl': 'C(=O)O',
            'hydroxyl': '[OH]',
            'thiol': '[SH]',
            'amide': 'C(=O)N',
            'phenol': 'c[OH]',
            'imidazole': 'c1c[nH]cn1',
            'indole': 'c1c[nH]c2ccccc12',
            'guanidinium': 'NC(=N)N'
        }
        
        for group_name, pattern in patterns.items():
            try:
                pattern_mol = Chem.MolFromSmarts(pattern)
                if pattern_mol:
                    matches = mol.GetSubstructMatches(pattern_mol)
                    if matches:
                        functional_groups[group_name] = [list(match) for match in matches]
            except:
                continue
        
        return functional_groups
    
    def _identify_functional_groups_basic(self, smiles: str) -> Dict[str, List[int]]:
        """基础官能团识别"""
        functional_groups = {}
        
        # 简单的模式匹配
        patterns = {
            'amino': ['N', '[NH2]', '[NH3]'],
            'carboxyl': ['C(=O)O', 'COOH'],
            'hydroxyl': ['OH'],
            'thiol': ['SH'],
            'amide': ['C(=O)N'],
        }
        
        for group_name, group_patterns in patterns.items():
            matches = []
            for pattern in group_patterns:
                for match in re.finditer(re.escape(pattern), smiles):
                    matches.append([match.start()])
            if matches:
                functional_groups[group_name] = matches
        
        return functional_groups
    
    def _analyze_backbone_rdkit(self, mol) -> Dict[str, Any]:
        """使用RDKit分析氨基酸主链结构"""
        from rdkit import Chem
        
        try:
            # 寻找氨基和羧基
            amino_pattern = Chem.MolFromSmarts('[NH2,NH3]')
            carboxyl_pattern = Chem.MolFromSmarts('C(=O)[OH]')
            
            amino_matches = mol.GetSubstructMatches(amino_pattern)
            carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
            
            if not amino_matches or not carboxyl_matches:
                return {'backbone_type': 'unknown', 'reason': 'no_amino_carboxyl_groups'}
            
            # 计算最短路径
            min_distance = float('inf')
            best_path = None
            
            for amino_match in amino_matches:
                amino_atom = amino_match[0]  # 氮原子
                for carboxyl_match in carboxyl_matches:
                    carboxyl_carbon = carboxyl_match[0]  # 羧基碳原子
                    
                    try:
                        path = Chem.GetShortestPath(mol, amino_atom, carboxyl_carbon)
                        if path and len(path) < min_distance:
                            min_distance = len(path)
                            best_path = path
                    except:
                        continue
            
            if best_path is None:
                return {'backbone_type': 'unknown', 'reason': 'no_path_found'}
            
            # 根据路径长度判断主链类型
            path_length = len(best_path) - 1  # 减去起点
            
            if path_length <= 2:
                backbone_type = 'alpha'
            elif path_length == 3:
                backbone_type = 'beta'
            elif path_length == 4:
                backbone_type = 'gamma'
            else:
                backbone_type = 'extended'
            
            return {
                'backbone_type': backbone_type,
                'path_length': path_length,
                'amino_carboxyl_path': best_path,
                'confidence': 0.95
            }
            
        except Exception as e:
            return {'backbone_type': 'unknown', 'reason': f'analysis_error: {e}'}
    
    def _analyze_backbone_basic(self, smiles: str) -> Dict[str, Any]:
        """基础主链分析"""
        # 简化的启发式方法
        if 'NCC(=O)O' in smiles:
            return {'backbone_type': 'alpha', 'confidence': 0.7}
        elif 'NCCC(=O)O' in smiles:
            return {'backbone_type': 'beta', 'confidence': 0.7}
        elif 'NCCCC(=O)O' in smiles:
            return {'backbone_type': 'gamma', 'confidence': 0.7}
        else:
            return {'backbone_type': 'unknown', 'confidence': 0.0}
    
    def detect_stereochemistry(self, smiles: str) -> Dict[str, Any]:
        """检测立体化学信息"""
        if not self.rdkit_available:
            return self._detect_stereochemistry_basic(smiles)
        
        try:
            from rdkit import Chem
            
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return {'stereochemistry': 'unknown', 'confidence': 0.0}
            
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            
            if not chiral_centers:
                return {'stereochemistry': 'achiral', 'confidence': 1.0}
            
            # 对于氨基酸，通常只关心α碳的立体化学
            # 这里需要更复杂的逻辑来确定D/L构型
            stereochemistry_info = {
                'chiral_centers_count': len(chiral_centers),
                'chiral_centers': chiral_centers,
                'stereochemistry': 'chiral_detected',
                'confidence': 0.8,
                'note': 'D/L assignment requires additional analysis'
            }
            
            return stereochemistry_info
            
        except Exception as e:
            return {'stereochemistry': 'unknown', 'confidence': 0.0, 'error': str(e)}
    
    def _detect_stereochemistry_basic(self, smiles: str) -> Dict[str, Any]:
        """基础立体化学检测"""
        if '@' in smiles:
            if 'C@@H' in smiles:
                return {'stereochemistry': 'L_suspected', 'confidence': 0.6}
            elif 'C@H' in smiles:
                return {'stereochemistry': 'D_suspected', 'confidence': 0.6}
            else:
                return {'stereochemistry': 'chiral_detected', 'confidence': 0.5}
        else:
            return {'stereochemistry': 'achiral_or_unknown', 'confidence': 0.7}
    
    def is_aromatic(self, smiles: str) -> Tuple[bool, float]:
        """检测是否为芳香化合物"""
        if not self.rdkit_available:
            return self._is_aromatic_basic(smiles)
        
        try:
            from rdkit import Chem
            
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False, 0.0
            
            # 检查是否有芳香原子
            aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
            
            if aromatic_atoms:
                confidence = min(1.0, len(aromatic_atoms) / mol.GetNumAtoms() + 0.5)
                return True, confidence
            else:
                return False, 1.0
                
        except Exception:
            return self._is_aromatic_basic(smiles)
    
    def _is_aromatic_basic(self, smiles: str) -> Tuple[bool, float]:
        """基础芳香性检测"""
        aromatic_indicators = ['c1', 'c2', 'cccc', 'benzene', 'phenyl', 'indole', 'imidazole']
        
        for indicator in aromatic_indicators:
            if indicator in smiles.lower():
                return True, 0.7
        
        return False, 0.8
    
    def has_rings(self, smiles: str) -> Tuple[bool, int, float]:
        """检测环状结构"""
        if not self.rdkit_available:
            return self._has_rings_basic(smiles)
        
        try:
            from rdkit import Chem
            
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return False, 0, 0.0
            
            ring_info = mol.GetRingInfo()
            num_rings = ring_info.NumRings()
            
            return num_rings > 0, num_rings, 1.0
            
        except Exception:
            return self._has_rings_basic(smiles)
    
    def _has_rings_basic(self, smiles: str) -> Tuple[bool, int, float]:
        """基础环检测"""
        ring_numbers = re.findall(r'\d+', smiles)
        unique_numbers = set(ring_numbers)
        ring_count = len([n for n in unique_numbers if ring_numbers.count(n) >= 2])
        
        return ring_count > 0, ring_count, 0.8
