"""
指纹相似性验证器
基于分子指纹进行氨基酸验证，采用重原子骨架方法
"""

from typing import Dict, Any, Optional

from .base import BaseVerifier
from ..models import ResidueInfo, AminoAcidInfo, VerificationMethod
from ...utils.chemistry import ChemistryUtils, MolecularFingerprint


class FingerprintSimilarityVerifier(BaseVerifier):
    """指纹相似性验证器"""
    
    def __init__(self, threshold: float = 0.7):
        super().__init__(threshold)
        self.chemistry_utils = ChemistryUtils()
        self.fingerprint_utils = MolecularFingerprint()
    
    def get_method(self) -> VerificationMethod:
        return VerificationMethod.FINGERPRINT_SIMILARITY
    
    def calculate_score(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> float:
        """
        计算指纹相似性验证分数
        
        验证逻辑：
        1. 基于重原子组成计算加权相似性
        2. 如果有SMILES，计算Morgan指纹相似性
        3. 返回更高的分数
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            验证分数 (0.0-1.0)
        """
        # 方法1: 基于重原子组成的加权相似性
        composition_score = self._calculate_weighted_composition_similarity(residue, amino_acid)
        
        # 方法2: 基于Morgan指纹的相似性（如果可用）
        fingerprint_score = self._calculate_morgan_fingerprint_similarity(residue, amino_acid)
        
        # 返回更高的分数
        if fingerprint_score is not None:
            return max(composition_score, fingerprint_score)
        else:
            return composition_score
    
    def _calculate_weighted_composition_similarity(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> float:
        """
        基于重原子组成计算加权相似性
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            加权相似性分数 (0.0-1.0)
        """
        # 提取重原子组成（符合分子指纹标准）
        residue_heavy = residue.heavy_atom_composition
        candidate_heavy = amino_acid.heavy_atom_composition
        
        if not residue_heavy or not candidate_heavy:
            return 0.0
        
        # 计算加权Jaccard相似性（仅重原子）
        return self.chemistry_utils.calculate_weighted_jaccard_similarity(
            residue_heavy, candidate_heavy, self.chemistry_utils.ELEMENT_WEIGHTS
        )
    
    def _calculate_morgan_fingerprint_similarity(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> Optional[float]:
        """
        基于Morgan指纹计算相似性
        
        策略：
        1. 优先使用数据库中的预计算指纹
        2. 如果没有预计算指纹，尝试从SMILES计算
        3. 从残基生成SMILES（暂未实现，返回None）
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            指纹相似性分数，失败返回None
        """
        # 策略1: 使用数据库中的预计算指纹
        if amino_acid.fingerprints and 'ecfp2' in amino_acid.fingerprints:
            # 尝试通过残基名匹配获取指纹
            residue_fingerprint = self._get_residue_fingerprint_from_database(residue)
            if residue_fingerprint:
                candidate_fingerprint = amino_acid.fingerprints['ecfp2']
                return self.fingerprint_utils.calculate_tanimoto_similarity(
                    residue_fingerprint, candidate_fingerprint
                )
        
        # 策略2: 从SMILES计算指纹
        candidate_smiles = amino_acid.smiles
        if not candidate_smiles:
            return None
        
        # 从残基生成SMILES（暂未实现）
        residue_smiles = self._generate_residue_smiles(residue)
        if not residue_smiles:
            return None
        
        # 计算Morgan指纹
        residue_fp = self.fingerprint_utils.calculate_morgan_fingerprint(residue_smiles)
        candidate_fp = self.fingerprint_utils.calculate_morgan_fingerprint(candidate_smiles)
        
        if not residue_fp or not candidate_fp:
            return None
        
        # 计算Tanimoto相似性
        return self.fingerprint_utils.calculate_tanimoto_similarity(residue_fp, candidate_fp)
    
    def _get_residue_fingerprint_from_database(self, residue: ResidueInfo) -> Optional[str]:
        """
        从数据库获取残基的指纹
        
        Args:
            residue: 残基信息
        
        Returns:
            指纹字符串，失败返回None
        """
        try:
            # 导入数据库管理器
            from ..database import DatabaseManager
            from ...utils.config import default_config
            
            db_manager = DatabaseManager(default_config)
            amino_acid = db_manager.get_amino_acid_by_id(residue.residue_name)
            
            if amino_acid and amino_acid.fingerprints and 'ecfp2' in amino_acid.fingerprints:
                return amino_acid.fingerprints['ecfp2']
            
            return None
        except Exception:
            return None
    
    def _generate_residue_smiles(self, residue: ResidueInfo) -> Optional[str]:
        """
        从残基信息生成SMILES（简化实现）
        
        注意：这是一个简化的实现，实际应用中可能需要更复杂的算法
        
        Args:
            residue: 残基信息
        
        Returns:
            SMILES字符串，失败返回None
        """
        # 这里可以实现从3D坐标生成SMILES的算法
        # 目前返回None，表示不支持
        return None
    
    def _get_verification_details(self, residue: ResidueInfo, amino_acid: AminoAcidInfo, score: float) -> Dict[str, Any]:
        """获取指纹相似性验证详情"""
        details = super()._get_verification_details(residue, amino_acid, score)
        
        # 计算各种相似性分数
        composition_score = self._calculate_weighted_composition_similarity(residue, amino_acid)
        fingerprint_score = self._calculate_morgan_fingerprint_similarity(residue, amino_acid)
        
        details.update({
            'residue_heavy_composition': residue.heavy_atom_composition,
            'candidate_heavy_composition': amino_acid.heavy_atom_composition,
            'weighted_composition_similarity': composition_score,
            'morgan_fingerprint_similarity': fingerprint_score,
            'candidate_smiles': amino_acid.smiles,
            'method_used': self._get_method_used(composition_score, fingerprint_score)
        })
        
        return details
    
    def _get_method_used(self, composition_score: float, fingerprint_score: Optional[float]) -> str:
        """获取使用的方法"""
        if fingerprint_score is not None:
            if fingerprint_score >= composition_score:
                return "morgan_fingerprint"
            else:
                return "weighted_composition"
        else:
            return "weighted_composition_only"
    
    def calculate_detailed_similarity(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """
        计算详细的相似性分析
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            详细相似性分析结果
        """
        result = {
            'residue_heavy_composition': residue.heavy_atom_composition,
            'candidate_heavy_composition': amino_acid.heavy_atom_composition,
            'candidate_smiles': amino_acid.smiles
        }
        
        # 加权组成相似性
        composition_score = self._calculate_weighted_composition_similarity(residue, amino_acid)
        result['weighted_composition_similarity'] = composition_score
        
        # Morgan指纹相似性
        fingerprint_score = self._calculate_morgan_fingerprint_similarity(residue, amino_acid)
        result['morgan_fingerprint_similarity'] = fingerprint_score
        
        # 简单Jaccard相似性（作为对比）
        simple_jaccard = self.chemistry_utils.calculate_jaccard_similarity(
            residue.heavy_atom_composition, amino_acid.heavy_atom_composition
        )
        result['simple_jaccard_similarity'] = simple_jaccard
        
        # 最终分数
        if fingerprint_score is not None:
            result['final_score'] = max(composition_score, fingerprint_score)
            result['best_method'] = 'morgan_fingerprint' if fingerprint_score >= composition_score else 'weighted_composition'
        else:
            result['final_score'] = composition_score
            result['best_method'] = 'weighted_composition'
        
        # 元素权重分析
        result['element_weights_used'] = self.chemistry_utils.ELEMENT_WEIGHTS
        
        return result
    
    def verify_with_database_fingerprint(self, residue: ResidueInfo, amino_acid: AminoAcidInfo, 
                                       fingerprint_type: str = 'ecfp2') -> Dict[str, Any]:
        """
        使用数据库中的预计算指纹进行验证
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
            fingerprint_type: 指纹类型
        
        Returns:
            验证结果
        """
        result = {
            'fingerprint_type': fingerprint_type,
            'database_fingerprint_available': False,
            'similarity_score': 0.0
        }
        
        # 检查数据库中是否有预计算的指纹
        if fingerprint_type in amino_acid.fingerprints:
            db_fingerprint = amino_acid.fingerprints[fingerprint_type]
            result['database_fingerprint_available'] = True
            result['database_fingerprint'] = db_fingerprint
            
            # 生成残基指纹（如果可能）
            residue_smiles = self._generate_residue_smiles(residue)
            if residue_smiles:
                residue_fingerprint = self.fingerprint_utils.calculate_morgan_fingerprint(residue_smiles)
                if residue_fingerprint:
                    similarity = self.fingerprint_utils.calculate_tanimoto_similarity(
                        residue_fingerprint, db_fingerprint
                    )
                    result['similarity_score'] = similarity
                    result['residue_fingerprint'] = residue_fingerprint
        
        return result