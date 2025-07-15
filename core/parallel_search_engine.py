#!/usr/bin/env python3
"""
多重验证并行搜索引擎
实现科学严谨的氨基酸识别策略
"""

import sys
import json
import sqlite3
import numpy as np
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass

sys.path.insert(0, '.')
sys.path.insert(0, 'core')
sys.path.insert(0, 'legacy')

try:
    from unified_pdb_analyzer import UnifiedSearchEngine, MatchResult, ResidueInfo
except ImportError:
    from core.unified_pdb_analyzer import UnifiedSearchEngine, MatchResult, ResidueInfo

@dataclass
class ParallelMatchResult(MatchResult):
    """并行验证匹配结果"""
    formula_score: Optional[float] = None
    atom_score: Optional[float] = None
    fingerprint_score: Optional[float] = None
    structure_3d_score: Optional[float] = None
    verification_details: Optional[Dict] = None

class ParallelSearchEngine(UnifiedSearchEngine):
    """多重验证并行搜索引擎"""
    
    def __init__(self):
        super().__init__()
        self.thresholds = {
            'molecular_formula': 1.0,    # 完全匹配
            'atom_composition': 0.9,     # 允许小幅偏差
            'fingerprint_similarity': 0.7,  # 化学信息学标准
            'structure_3d': 0.5          # 提高到更协调的水平
        }
        print("🚀 多重验证并行搜索引擎已初始化")
        print(f"📊 验证阈值: {self.thresholds}")
    
    def parallel_search(self, residue: ResidueInfo, 
                       enable_fingerprint: bool = True,
                       enable_3d: bool = True) -> List[ParallelMatchResult]:
        """多重验证并行搜索"""
        
        print(f"\n🔍 开始并行验证搜索: {residue.residue_name}")
        
        # 第一步：残基名精确匹配（快速通道）
        exact_matches = self.search_by_residue_name(residue.residue_name)
        if exact_matches:
            print(f"   ✅ 残基名精确匹配成功: {exact_matches[0].amino_acid_name}")
            # 转换为ParallelMatchResult
            result = ParallelMatchResult(
                amino_acid_id=exact_matches[0].amino_acid_id,
                amino_acid_name=exact_matches[0].amino_acid_name,
                confidence_score=1.0,
                match_method="residue_name_exact",
                molecular_formula=exact_matches[0].molecular_formula,
                smiles=exact_matches[0].smiles,
                residue_info=residue,
                verification_details={'method': 'exact_match', 'bypass_parallel': True}
            )
            return [result]
        
        print(f"   ⚠️ 残基名精确匹配失败，启动并行验证...")
        
        # 第二步：并行验证策略
        return self._parallel_verification(residue, enable_fingerprint, enable_3d)
    
    def _parallel_verification(self, residue: ResidueInfo, 
                              enable_fingerprint: bool,
                              enable_3d: bool) -> List[ParallelMatchResult]:
        """执行并行验证"""
        
        # 获取所有候选氨基酸
        candidates = self._get_candidates(residue)
        
        if not candidates:
            print(f"   ❌ 未找到候选氨基酸")
            return []
        
        print(f"   📋 找到 {len(candidates)} 个候选氨基酸")
        
        # 对每个候选进行四重验证
        verified_results = []
        
        for candidate_id in candidates:
            verification_result = self._verify_candidate(
                residue, candidate_id, enable_fingerprint, enable_3d
            )
            
            if verification_result:
                verified_results.append(verification_result)
        
        # 按综合置信度排序
        verified_results.sort(key=lambda x: x.confidence_score, reverse=True)
        
        print(f"   ✅ 并行验证完成，通过验证: {len(verified_results)} 个")
        
        return verified_results
    
    def _get_candidates(self, residue: ResidueInfo) -> List[str]:
        """获取候选氨基酸列表"""
        candidates = set()
        
        # 基于分子式筛选
        if residue.molecular_formula:
            formula_matches = self.search_by_molecular_formula(residue.molecular_formula)
            candidates.update([m.amino_acid_id for m in formula_matches])
        
        # 基于原子组成筛选
        if residue.atom_composition:
            atom_matches = self.search_by_atom_composition(residue.atom_composition)
            candidates.update([m.amino_acid_id for m in atom_matches])
        
        # 如果没有基础筛选结果，考虑所有氨基酸（用于指纹相似性搜索）
        if not candidates:
            candidates = set(self.amino_acids.keys())
        
        return list(candidates)
    
    def _verify_candidate(self, residue: ResidueInfo, candidate_id: str,
                         enable_fingerprint: bool, enable_3d: bool) -> Optional[ParallelMatchResult]:
        """对候选氨基酸进行四重验证"""
        
        candidate_data = self.amino_acids[candidate_id]
        verification_scores = {}
        verification_details = {}
        
        # 验证1：分子式匹配
        formula_score = self._verify_molecular_formula(
            residue.molecular_formula, candidate_data['molecular_formula']
        )
        verification_scores['formula'] = formula_score
        verification_details['formula_match'] = formula_score >= self.thresholds['molecular_formula']
        
        # 验证2：原子组成匹配
        atom_score = self._verify_atom_composition(
            residue.atom_composition, candidate_data['atom_composition']
        )
        verification_scores['atom'] = atom_score
        verification_details['atom_match'] = atom_score >= self.thresholds['atom_composition']
        
        # 验证3：指纹相似性匹配
        fingerprint_score = 0.0
        if enable_fingerprint and candidate_data['smiles']:
            fingerprint_score = self._verify_fingerprint_similarity(
                residue, candidate_data['smiles']
            )
        verification_scores['fingerprint'] = fingerprint_score
        verification_details['fingerprint_match'] = fingerprint_score >= self.thresholds['fingerprint_similarity']
        
        # 验证4：3D结构验证
        structure_score = 0.0
        if enable_3d:
            structure_score = self._verify_3d_structure(residue, candidate_id)
        verification_scores['structure'] = structure_score
        verification_details['structure_match'] = structure_score >= self.thresholds['structure_3d']
        
        # 检查是否所有验证都通过
        required_verifications = ['formula', 'atom']
        if enable_fingerprint:
            required_verifications.append('fingerprint')
        if enable_3d:
            required_verifications.append('structure')
        
        all_passed = all(
            verification_details[f'{v}_match'] for v in required_verifications
        )
        
        if not all_passed:
            return None
        
        # 计算综合置信度（乘积模型）
        confidence_score = self._calculate_combined_confidence(verification_scores)
        
        # 创建结果
        result = ParallelMatchResult(
            amino_acid_id=candidate_id,
            amino_acid_name=candidate_data['name'],
            confidence_score=confidence_score,
            match_method="parallel_verification",
            molecular_formula=candidate_data['molecular_formula'],
            smiles=candidate_data['smiles'],
            residue_info=residue,
            formula_score=formula_score,
            atom_score=atom_score,
            fingerprint_score=fingerprint_score,
            structure_3d_score=structure_score,
            verification_details=verification_details
        )
        
        return result
    
    def _verify_molecular_formula(self, residue_formula: Optional[str],
                                 candidate_formula: str) -> float:
        """
        验证分子式匹配（智能氢原子处理）

        策略：
        1. 优先进行完整分子式比较（包含氢原子）
        2. 如果不匹配，尝试重原子分子式比较（忽略氢原子）
        3. 这样既保持化学完整性，又兼容数据不一致情况
        """
        if not residue_formula:
            return 0.0

        # 策略1: 完整分子式比较
        if residue_formula == candidate_formula:
            return 1.0

        # 策略2: 重原子分子式比较（fallback）
        residue_no_h = self._remove_hydrogen_from_formula(residue_formula)
        candidate_no_h = self._remove_hydrogen_from_formula(candidate_formula)

        if residue_no_h == candidate_no_h:
            return 0.9  # 稍微降低分数，因为不是完整匹配

        return 0.0

    def _remove_hydrogen_from_formula(self, formula: str) -> str:
        """从分子式中移除氢原子部分"""
        try:
            import re
            # 移除H和其后的数字（如H7, H13等）
            no_hydrogen = re.sub(r'H\d*', '', formula)
            return no_hydrogen
        except Exception as e:
            print(f"氢原子移除失败: {e}")
            return formula
    
    def _verify_atom_composition(self, residue_atoms: Optional[Dict],
                                candidate_atoms: Dict) -> float:
        """
        验证原子组成匹配（智能氢原子处理）

        策略：
        1. 检测数据库是否包含氢原子
        2. 如果都包含或都不包含氢原子，直接比较
        3. 如果一个包含一个不包含，忽略氢原子进行比较
        """
        if not residue_atoms:
            return 0.0

        # 检测氢原子情况
        residue_has_h = 'H' in residue_atoms
        candidate_has_h = 'H' in candidate_atoms

        # 选择比较策略
        if residue_has_h == candidate_has_h:
            # 氢原子情况一致，直接比较
            comp1, comp2 = residue_atoms, candidate_atoms
        else:
            # 氢原子情况不一致，忽略氢原子比较
            comp1 = {k: v for k, v in residue_atoms.items() if k != 'H'}
            comp2 = {k: v for k, v in candidate_atoms.items() if k != 'H'}

        # 完全匹配检查
        if comp1 == comp2:
            return 1.0

        # 计算Jaccard相似性
        all_elements = set(comp1.keys()) | set(comp2.keys())
        if not all_elements:
            return 0.0

        intersection = 0
        union = 0

        for element in all_elements:
            res_count = comp1.get(element, 0)
            cand_count = comp2.get(element, 0)
            intersection += min(res_count, cand_count)
            union += max(res_count, cand_count)

        return intersection / union if union > 0 else 0.0
    
    def _verify_fingerprint_similarity(self, residue: ResidueInfo,
                                     candidate_smiles: str) -> float:
        """验证指纹相似性 - 改进实现"""
        try:
            # 方法1: 尝试真正的分子指纹计算
            try:
                from legacy.isomer_identifier import IsomerIdentifier
                if not hasattr(self, '_isomer_identifier_cache'):
                    self._isomer_identifier_cache = IsomerIdentifier()

                # 如果residue有SMILES信息，使用真正的指纹相似性
                if hasattr(residue, 'smiles') and residue.smiles:
                    relationship = self._isomer_identifier_cache.identify_relationship(
                        residue.smiles, candidate_smiles
                    )
                    return relationship.get('structural_similarity', 0.0)

            except ImportError:
                pass  # 继续使用fallback方法

            # 方法2: 改进的加权组成相似性
            if residue.molecular_formula and residue.atom_composition:
                weighted_similarity = self._calculate_weighted_composition_similarity(
                    residue, candidate_smiles
                )
                return weighted_similarity

            # 方法3: 基础组成相似性（最后的fallback）
            if residue.atom_composition:
                basic_similarity = self._calculate_basic_composition_similarity(
                    residue, candidate_smiles
                )
                return basic_similarity * 0.8  # 降低置信度

            # 如果没有任何信息，返回低分
            return 0.1

        except Exception as e:
            print(f"指纹相似性计算失败: {e}")
            return 0.0

    def _calculate_weighted_composition_similarity(self, residue: ResidueInfo, candidate_smiles: str) -> float:
        """
        基于加权原子组成计算相似性（重原子骨架为主）

        化学信息学依据：
        - 分子指纹主要基于重原子拓扑结构
        - 氢原子对分子骨架特征影响较小
        - 这是ECFP、MACCS等主流指纹算法的标准做法
        """
        try:
            # 获取候选氨基酸的原子组成（使用缓存索引）
            candidate_atoms = self._get_candidate_atom_composition(candidate_smiles)
            if not candidate_atoms or not residue.atom_composition:
                return 0.0

            # 提取重原子组成（忽略氢原子，符合分子指纹标准）
            residue_heavy = {k: v for k, v in residue.atom_composition.items() if k != 'H'}
            candidate_heavy = {k: v for k, v in candidate_atoms.items() if k != 'H'}

            # 重原子权重（基于化学重要性）
            element_weights = {
                'C': 1.0,   # 碳骨架最重要
                'N': 0.9,   # 氮原子很重要（氨基酸特征）
                'O': 0.9,   # 氧原子很重要（羧基等）
                'S': 0.8,   # 硫原子重要（半胱氨酸等）
                'P': 0.8,   # 磷原子重要（磷酸化等）
                'F': 0.8,   # 氟原子重要（强电负性）
                'Cl': 0.7,  # 氯原子
                'Br': 0.6,  # 溴原子
                'I': 0.5,   # 碘原子
            }

            # 计算加权Jaccard相似性（仅重原子）
            all_elements = set(residue_heavy.keys()) | set(candidate_heavy.keys())
            weighted_intersection = 0.0
            weighted_union = 0.0

            for element in all_elements:
                weight = element_weights.get(element, 0.5)  # 默认权重0.5
                res_count = residue_heavy.get(element, 0)
                cand_count = candidate_heavy.get(element, 0)

                weighted_intersection += weight * min(res_count, cand_count)
                weighted_union += weight * max(res_count, cand_count)

            return weighted_intersection / weighted_union if weighted_union > 0 else 0.0

        except Exception as e:
            print(f"加权组成相似性计算失败: {e}")
            return 0.0

    def _calculate_basic_composition_similarity(self, residue: ResidueInfo, candidate_smiles: str) -> float:
        """基础原子组成相似性（fallback方法）"""
        try:
            candidate_atoms = self._get_candidate_atom_composition(candidate_smiles)
            if not candidate_atoms or not residue.atom_composition:
                return 0.0

            residue_atoms = residue.atom_composition

            # 标准Jaccard相似性
            all_elements = set(residue_atoms.keys()) | set(candidate_atoms.keys())
            intersection = sum(min(residue_atoms.get(e, 0), candidate_atoms.get(e, 0)) for e in all_elements)
            union = sum(max(residue_atoms.get(e, 0), candidate_atoms.get(e, 0)) for e in all_elements)

            return intersection / union if union > 0 else 0.0

        except Exception as e:
            print(f"基础组成相似性计算失败: {e}")
            return 0.0

    def _get_candidate_atom_composition(self, candidate_smiles: str) -> Dict:
        """获取候选氨基酸的原子组成（使用缓存索引）"""
        try:
            # 建立SMILES到原子组成的索引缓存
            if not hasattr(self, '_smiles_to_composition_cache'):
                self._smiles_to_composition_cache = {}
                for amino_id, aa_data in self.amino_acids.items():
                    smiles = aa_data.get('smiles', '')
                    if smiles:
                        self._smiles_to_composition_cache[smiles] = aa_data.get('atom_composition', {})

            return self._smiles_to_composition_cache.get(candidate_smiles, {})

        except Exception as e:
            print(f"获取候选原子组成失败: {e}")
            return {}

    def _verify_3d_structure(self, residue: ResidueInfo, candidate_id: str) -> float:
        """验证3D结构匹配 - 改进实现"""
        try:
            # 加载3D结构数据（使用缓存）
            if not hasattr(self, '_structure_data_cache'):
                self._load_3d_structure_data()

            if candidate_id not in self._structure_data_cache:
                return 0.0  # 无3D结构数据

            # 提取PDB坐标
            pdb_coords = self._extract_pdb_coordinates(residue)
            if not pdb_coords:
                return 0.0  # 坐标提取失败

            # 获取标准3D结构坐标
            standard_structure = self._structure_data_cache[candidate_id]
            standard_coords = self._extract_standard_coordinates(standard_structure)

            if not standard_coords:
                return 0.0  # 标准坐标不可用

            # 原子匹配和对齐
            aligned_coords = self._align_atom_coordinates(pdb_coords, standard_coords)
            if not aligned_coords:
                return 0.0  # 原子匹配失败

            # 使用Kabsch算法计算最优RMSD
            rmsd = self._calculate_kabsch_rmsd(aligned_coords[0], aligned_coords[1])

            # 基于分子大小的动态阈值转换为分数
            score = self._rmsd_to_score_dynamic(rmsd, len(aligned_coords[0]))

            return score

        except Exception as e:
            print(f"3D结构验证失败: {e}")
            return 0.0

    def _load_3d_structure_data(self):
        """加载3D结构数据（改进缓存版本）"""
        try:
            import sqlite3
            import json

            conn = sqlite3.connect('core/amino_acids.db')
            cursor = conn.cursor()

            cursor.execute('''
                SELECT id, structure_data
                FROM amino_acids
                WHERE structure_data IS NOT NULL
            ''')

            self._structure_data_cache = {}
            for amino_id, structure_json in cursor.fetchall():
                try:
                    structure = json.loads(structure_json)
                    self._structure_data_cache[amino_id] = structure
                except:
                    continue

            conn.close()
            print(f"✅ 加载了 {len(self._structure_data_cache)} 种氨基酸的3D结构到缓存")

        except Exception as e:
            print(f"⚠️ 3D结构数据加载失败: {e}")
            self._structure_data_cache = {}

    def _extract_pdb_coordinates(self, residue: ResidueInfo) -> List[List[float]]:
        """
        提取PDB残基的坐标（标准化原子顺序，忽略氢原子）

        注意：忽略氢原子是为了与数据库中的3D结构数据保持一致。
        数据库中的3D结构通常不包含氢原子坐标，这样可以确保
        PDB残基与标准结构的原子数匹配，提高3D结构验证的准确性。
        """
        try:
            if not residue.atoms or len(residue.atoms) < 2:
                return []

            # 过滤掉氢原子，保持与数据库3D结构一致性
            non_hydrogen_atoms = [atom for atom in residue.atoms if atom.get('element') != 'H']

            if len(non_hydrogen_atoms) < 2:
                return []

            # 标准化原子顺序：按元素类型排序，然后按坐标排序
            sorted_atoms = self._standardize_atom_order(non_hydrogen_atoms)

            coordinates = []
            for atom in sorted_atoms:
                if all(key in atom for key in ['x', 'y', 'z']):
                    coordinates.append([float(atom['x']), float(atom['y']), float(atom['z'])])

            return coordinates

        except Exception as e:
            print(f"PDB坐标提取失败: {e}")
            return []

    def _standardize_atom_order(self, atoms: List[Dict]) -> List[Dict]:
        """标准化原子顺序，确保结果一致性"""
        try:
            # 定义元素优先级（常见氨基酸元素）
            element_priority = {
                'N': 1, 'C': 2, 'O': 3, 'S': 4, 'P': 5,
                'F': 6, 'Cl': 7, 'Br': 8, 'I': 9, 'H': 10
            }

            def sort_key(atom):
                element = atom.get('element', 'X')
                priority = element_priority.get(element, 99)
                x = float(atom.get('x', 0))
                y = float(atom.get('y', 0))
                z = float(atom.get('z', 0))

                # 排序键：元素优先级 -> x坐标 -> y坐标 -> z坐标
                return (priority, round(x, 3), round(y, 3), round(z, 3))

            return sorted(atoms, key=sort_key)

        except Exception as e:
            print(f"原子顺序标准化失败: {e}")
            return atoms

    def _extract_standard_coordinates(self, standard_structure: Dict) -> List[List[float]]:
        """提取标准结构的坐标"""
        try:
            if 'atoms_info' not in standard_structure:
                return []

            coordinates = []
            for atom in standard_structure['atoms_info']:
                if 'coordinates' in atom and len(atom['coordinates']) >= 3:
                    coord = atom['coordinates']
                    coordinates.append([float(coord[0]), float(coord[1]), float(coord[2])])

            return coordinates

        except Exception as e:
            print(f"标准坐标提取失败: {e}")
            return []

    def _align_atom_coordinates(self, pdb_coords: List[List[float]],
                               standard_coords: List[List[float]]) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """原子坐标对齐和匹配"""
        try:
            if len(pdb_coords) != len(standard_coords):
                # 尝试基于距离的最佳匹配
                if len(pdb_coords) < len(standard_coords):
                    # PDB原子较少，选择最接近的标准原子
                    matched_standard = self._find_closest_atoms(pdb_coords, standard_coords)
                    return np.array(pdb_coords), np.array(matched_standard)
                else:
                    # 标准原子较少，选择最接近的PDB原子
                    matched_pdb = self._find_closest_atoms(standard_coords, pdb_coords)
                    return np.array(matched_pdb), np.array(standard_coords)
            else:
                # 原子数相同，直接使用
                return np.array(pdb_coords), np.array(standard_coords)

        except Exception as e:
            print(f"原子对齐失败: {e}")
            return None

    def _find_closest_atoms(self, reference_coords: List[List[float]],
                           target_coords: List[List[float]]) -> List[List[float]]:
        """找到最接近的原子匹配"""
        try:
            ref_array = np.array(reference_coords)
            target_array = np.array(target_coords)

            matched_coords = []
            used_indices = set()

            for ref_coord in ref_array:
                # 计算到所有目标原子的距离
                distances = np.linalg.norm(target_array - ref_coord, axis=1)

                # 找到最近的未使用原子
                sorted_indices = np.argsort(distances)
                for idx in sorted_indices:
                    if idx not in used_indices:
                        matched_coords.append(target_coords[idx])
                        used_indices.add(idx)
                        break

            return matched_coords

        except Exception as e:
            print(f"原子匹配失败: {e}")
            return target_coords[:len(reference_coords)]

    def _calculate_kabsch_rmsd(self, coords1: np.ndarray, coords2: np.ndarray) -> float:
        """使用Kabsch算法计算最优叠合RMSD"""
        try:
            if coords1.shape != coords2.shape:
                return float('inf')

            # 中心化坐标
            centroid1 = np.mean(coords1, axis=0)
            centroid2 = np.mean(coords2, axis=0)

            coords1_centered = coords1 - centroid1
            coords2_centered = coords2 - centroid2

            # Kabsch算法：计算最优旋转矩阵
            H = coords1_centered.T @ coords2_centered
            U, S, Vt = np.linalg.svd(H)
            R = Vt.T @ U.T

            # 确保是右手坐标系（行列式为正）
            if np.linalg.det(R) < 0:
                Vt[-1, :] *= -1
                R = Vt.T @ U.T

            # 应用最优旋转
            coords1_rotated = coords1_centered @ R.T

            # 计算最优RMSD
            diff = coords1_rotated - coords2_centered
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

            return float(rmsd)

        except Exception as e:
            print(f"Kabsch RMSD计算失败: {e}")
            return float('inf')

    def _rmsd_to_score_dynamic(self, rmsd: float, atom_count: int) -> float:
        """基于分子大小的动态RMSD评分"""
        try:
            # 基于原子数的动态阈值
            if atom_count <= 10:
                max_rmsd = 1.0  # 小分子更严格
            elif atom_count <= 20:
                max_rmsd = 1.5  # 中等分子
            else:
                max_rmsd = 2.0  # 大分子稍微宽松

            # 使用指数衰减函数
            score = np.exp(-rmsd / max_rmsd)

            return float(np.clip(score, 0.0, 1.0))

        except Exception as e:
            print(f"RMSD评分转换失败: {e}")
            return 0.0

    def _extract_pdb_geometry(self, residue: ResidueInfo) -> Dict:
        """从PDB残基提取几何信息"""
        try:
            atoms = residue.atoms
            if len(atoms) < 2:
                return {}

            # 计算键长（简化：相邻原子距离）
            bond_lengths = []
            for i in range(len(atoms) - 1):
                atom1 = atoms[i]
                atom2 = atoms[i + 1]

                dist = ((atom1['x'] - atom2['x'])**2 +
                       (atom1['y'] - atom2['y'])**2 +
                       (atom1['z'] - atom2['z'])**2)**0.5
                bond_lengths.append(dist)

            # 计算键角（简化：三个连续原子的角度）
            bond_angles = []
            for i in range(len(atoms) - 2):
                atom1 = atoms[i]
                atom2 = atoms[i + 1]
                atom3 = atoms[i + 2]

                # 向量
                v1 = np.array([atom1['x'] - atom2['x'],
                              atom1['y'] - atom2['y'],
                              atom1['z'] - atom2['z']])
                v2 = np.array([atom3['x'] - atom2['x'],
                              atom3['y'] - atom2['y'],
                              atom3['z'] - atom2['z']])

                # 计算角度
                cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                cos_angle = np.clip(cos_angle, -1.0, 1.0)
                angle = np.arccos(cos_angle) * 180.0 / np.pi
                bond_angles.append(angle)

            return {
                'bond_lengths': bond_lengths,
                'bond_angles': bond_angles,
                'coordinates': [[atom['x'], atom['y'], atom['z']] for atom in atoms]
            }

        except Exception as e:
            print(f"PDB几何提取失败: {e}")
            return {}

    def _compare_bond_lengths(self, pdb_lengths: List[float], standard_lengths: List[Dict]) -> float:
        """比较键长"""
        try:
            if not pdb_lengths or not standard_lengths:
                return 0.0

            pdb_avg = np.mean(pdb_lengths) if pdb_lengths else 0
            standard_values = [bond['length'] for bond in standard_lengths if 'length' in bond]
            standard_avg = np.mean(standard_values) if standard_values else 0

            if standard_avg == 0:
                return 0.0

            # 计算相对误差
            relative_error = abs(pdb_avg - standard_avg) / standard_avg
            score = max(0, 1 - relative_error)

            return score

        except Exception as e:
            return 0.0

    def _compare_bond_angles(self, pdb_angles: List[float], standard_angles: List[Dict]) -> float:
        """比较键角"""
        try:
            if not pdb_angles or not standard_angles:
                return 0.0

            pdb_avg = np.mean(pdb_angles) if pdb_angles else 0
            standard_values = [angle['angle'] for angle in standard_angles if 'angle' in angle]
            standard_avg = np.mean(standard_values) if standard_values else 0

            if standard_avg == 0:
                return 0.0

            # 键角差异（度数）
            angle_diff = abs(pdb_avg - standard_avg)
            score = max(0, 1 - angle_diff / 180.0)

            return score

        except Exception as e:
            return 0.0

    def _calculate_rmsd_score(self, pdb_coords: List[List[float]],
                             standard_atoms: List[Dict]) -> float:
        """计算RMSD分数"""
        try:
            standard_coords = [[atom['coordinates'][0], atom['coordinates'][1], atom['coordinates'][2]]
                              for atom in standard_atoms]

            if len(pdb_coords) != len(standard_coords):
                return 0.0

            pdb_coords = np.array(pdb_coords)
            standard_coords = np.array(standard_coords)

            # 中心化坐标
            pdb_centered = pdb_coords - np.mean(pdb_coords, axis=0)
            standard_centered = standard_coords - np.mean(standard_coords, axis=0)

            # 计算RMSD
            diff = pdb_centered - standard_centered
            rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

            # RMSD转换为相似性分数（RMSD越小，分数越高）
            score = max(0, 1 - rmsd / 5.0)  # 假设5Å为最大可接受RMSD

            return score

        except Exception as e:
            return 0.0

    def _calculate_combined_confidence(self, scores: Dict[str, float]) -> float:
        """计算综合置信度（乘积模型）"""
        
        # 基础分数（分子式和原子组成）
        base_score = scores.get('formula', 0.0) * scores.get('atom', 0.0)
        
        # 高级分数（指纹和3D结构）
        advanced_scores = []
        if scores.get('fingerprint', 0.0) > 0:
            advanced_scores.append(scores['fingerprint'])
        if scores.get('structure', 0.0) > 0:
            advanced_scores.append(scores['structure'])
        
        if advanced_scores:
            advanced_score = np.mean(advanced_scores)
            # 乘积模型：基础分数 × 高级分数
            combined_score = base_score * advanced_score
        else:
            combined_score = base_score
        
        # 确保分数在合理范围内
        return min(combined_score, 1.0)
    
    def get_search_statistics(self) -> Dict:
        """获取搜索统计信息"""
        return {
            'total_amino_acids': len(self.amino_acids),
            'thresholds': self.thresholds,
            'verification_methods': [
                'molecular_formula',
                'atom_composition', 
                'fingerprint_similarity',
                '3d_structure'
            ]
        }

def test_parallel_search():
    """测试并行搜索引擎"""
    print("🧪 测试多重验证并行搜索引擎")
    print("=" * 50)
    
    engine = ParallelSearchEngine()
    
    # 创建测试残基
    test_residue = ResidueInfo(
        chain_id="A",
        residue_number=1,
        residue_name="ALA",
        atoms=[],
        molecular_formula="C3H7NO2",
        atom_composition={"C": 3, "H": 7, "N": 1, "O": 2}
    )
    
    # 测试并行搜索
    results = engine.parallel_search(test_residue)
    
    print(f"\n📊 搜索结果:")
    for result in results:
        print(f"  {result.amino_acid_id}: {result.amino_acid_name}")
        print(f"    置信度: {result.confidence_score:.3f}")
        print(f"    方法: {result.match_method}")
        if result.verification_details:
            print(f"    验证详情: {result.verification_details}")

if __name__ == "__main__":
    test_parallel_search()
