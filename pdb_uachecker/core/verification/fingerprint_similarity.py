"""
指纹相似性验证器
基于分子指纹进行氨基酸验证，采用重原子骨架方法
增强版：支持从PDB坐标生成SMILES和分子指纹
"""

from typing import Dict, Any, Optional, List
import logging

from .base import BaseVerifier
from ..models import ResidueInfo, AminoAcidInfo, VerificationMethod
from ...utils.chemistry import ChemistryUtils, MolecularFingerprint
from .coordinate_to_smiles import CoordinateToSmilesConverter, EnhancedFingerprintGenerator
from .stereochemistry_detector import StereochemistryDetector


class FingerprintSimilarityVerifier(BaseVerifier):
    """指纹相似性验证器"""
    
    def __init__(self, threshold: float = 0.7):
        super().__init__(threshold)
        self.chemistry_utils = ChemistryUtils()
        self.fingerprint_utils = MolecularFingerprint()
        
        # 新增的增强组件
        self.coord_converter = CoordinateToSmilesConverter()
        self.enhanced_fingerprint = EnhancedFingerprintGenerator(self.coord_converter)
        self.stereo_detector = StereochemistryDetector()
    
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
        从残基信息生成SMILES（增强实现）
        
        使用新的坐标到SMILES转换器，支持立体化学
        
        Args:
            residue: 残基信息
        
        Returns:
            SMILES字符串，失败返回None
        """
        try:
            # 使用增强的坐标转换器
            smiles = self.coord_converter.convert_residue_to_smiles(
                residue, include_stereochemistry=True
            )
            
            if smiles and self.coord_converter.validate_generated_smiles(smiles):
                logging.info(f"成功从残基 {residue.residue_name} 生成SMILES: {smiles}")
                return smiles
            else:
                logging.warning(f"残基 {residue.residue_name} SMILES生成失败或无效")
                return None
                
        except Exception as e:
            logging.error(f"残基SMILES生成异常: {e}")
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
    
    def calculate_enhanced_fingerprint_similarity(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> Dict[str, Any]:
        """
        计算增强的指纹相似性（包含立体化学）
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
            
        Returns:
            增强的相似性分析结果
        """
        result = {
            'residue_smiles': None,
            'candidate_smiles': amino_acid.smiles,
            'morgan_fingerprint_similarity': None,
            'stereochemistry_match': None,
            'enhanced_score': 0.0,
            'method_used': 'enhanced_fingerprint',
            'confidence': 0.0,
            'details': {}
        }
        
        try:
            # 1. 生成残基SMILES
            residue_smiles = self._generate_residue_smiles(residue)
            result['residue_smiles'] = residue_smiles
            
            if residue_smiles and amino_acid.smiles:
                # 2. 计算Morgan指纹相似性
                residue_fp = self.enhanced_fingerprint.generate_fingerprint_from_smiles(residue_smiles)
                candidate_fp = self.enhanced_fingerprint.generate_fingerprint_from_smiles(amino_acid.smiles)
                
                if residue_fp and candidate_fp:
                    similarity = self.fingerprint_utils.calculate_tanimoto_similarity(residue_fp, candidate_fp)
                    result['morgan_fingerprint_similarity'] = similarity
                
                # 3. 立体化学比较
                residue_stereo = self.stereo_detector.detect_stereochemistry_from_residue(residue)
                candidate_stereo = self.stereo_detector.detect_stereochemistry_from_smiles(amino_acid.smiles)
                
                stereo_comparison = self.stereo_detector.compare_stereochemistry(residue_stereo, candidate_stereo)
                result['stereochemistry_match'] = stereo_comparison
                
                # 4. 综合评分
                fingerprint_score = result['morgan_fingerprint_similarity'] or 0.0
                stereo_score = stereo_comparison.get('similarity_score', 0.5)
                
                # 加权综合分数
                result['enhanced_score'] = (fingerprint_score * 0.7 + stereo_score * 0.3)
                
                # 5. 置信度评估
                result['confidence'] = self._calculate_enhanced_confidence(result)
                
                result['details'] = {
                    'residue_stereochemistry': self.stereo_detector.generate_stereochemistry_report(residue_stereo),
                    'candidate_stereochemistry': self.stereo_detector.generate_stereochemistry_report(candidate_stereo),
                    'fingerprint_available': bool(residue_fp and candidate_fp),
                    'stereochemistry_analyzed': bool(stereo_comparison.get('match_type') != 'unknown')
                }
                
            else:
                # 降级到组成指纹
                composition_fp = self.enhanced_fingerprint.generate_composition_fingerprint(residue)
                result['enhanced_score'] = self._calculate_weighted_composition_similarity(residue, amino_acid)
                result['method_used'] = 'composition_fallback'
                result['confidence'] = 0.3
                result['details']['fallback_reason'] = 'SMILES生成失败'
                
        except Exception as e:
            logging.error(f"增强指纹相似性计算失败: {e}")
            result['enhanced_score'] = 0.0
            result['details']['error'] = str(e)
        
        return result
    
    def _calculate_enhanced_confidence(self, result: Dict[str, Any]) -> float:
        """
        计算增强指纹相似性的置信度
        
        Args:
            result: 相似性分析结果
            
        Returns:
            置信度分数
        """
        confidence_factors = []
        
        # 因子1: SMILES生成成功
        if result['residue_smiles']:
            confidence_factors.append(0.8)
        else:
            confidence_factors.append(0.2)
        
        # 因子2: Morgan指纹计算成功
        if result['morgan_fingerprint_similarity'] is not None:
            confidence_factors.append(0.9)
        else:
            confidence_factors.append(0.3)
        
        # 因子3: 立体化学分析质量
        stereo_match = result.get('stereochemistry_match', {})
        if stereo_match.get('match_type') in ['same_configuration', 'opposite_configuration']:
            confidence_factors.append(0.9)
        elif stereo_match.get('match_type') in ['both_achiral', 'mixed_configuration']:
            confidence_factors.append(0.7)
        else:
            confidence_factors.append(0.4)
        
        # 因子4: 分子复杂度
        if result['residue_smiles']:
            smiles_length = len(result['residue_smiles'])
            if smiles_length > 20:  # 复杂分子
                confidence_factors.append(0.8)
            elif smiles_length > 10:  # 中等复杂度
                confidence_factors.append(0.9)
            else:  # 简单分子
                confidence_factors.append(0.7)
        else:
            confidence_factors.append(0.5)
        
        # 计算综合置信度
        return sum(confidence_factors) / len(confidence_factors)
    
    def verify_isomer_distinction(self, residue: ResidueInfo, isomer_candidates: List[AminoAcidInfo]) -> Dict[str, Any]:
        """
        验证异构体区分能力
        
        专门用于测试系统对异构体的区分效果
        
        Args:
            residue: 残基信息
            isomer_candidates: 异构体候选列表
            
        Returns:
            异构体区分分析结果
        """
        analysis = {
            'total_candidates': len(isomer_candidates),
            'fingerprint_scores': [],
            'stereochemistry_analysis': [],
            'best_match': None,
            'distinction_quality': 'unknown',
            'confidence': 0.0
        }
        
        try:
            scores_with_details = []
            
            for candidate in isomer_candidates:
                enhanced_result = self.calculate_enhanced_fingerprint_similarity(residue, candidate)
                
                score_info = {
                    'candidate_id': candidate.id,
                    'candidate_name': candidate.name,
                    'enhanced_score': enhanced_result['enhanced_score'],
                    'fingerprint_similarity': enhanced_result.get('morgan_fingerprint_similarity'),
                    'stereochemistry_similarity': enhanced_result.get('stereochemistry_match', {}).get('similarity_score'),
                    'confidence': enhanced_result['confidence'],
                    'method_used': enhanced_result['method_used']
                }
                
                scores_with_details.append(score_info)
                analysis['fingerprint_scores'].append(score_info)
            
            # 按分数排序
            scores_with_details.sort(key=lambda x: x['enhanced_score'], reverse=True)
            
            if scores_with_details:
                analysis['best_match'] = scores_with_details[0]
                
                # 评估区分质量
                if len(scores_with_details) > 1:
                    best_score = scores_with_details[0]['enhanced_score']
                    second_score = scores_with_details[1]['enhanced_score']
                    score_gap = best_score - second_score
                    
                    if score_gap > 0.3:
                        analysis['distinction_quality'] = 'excellent'
                    elif score_gap > 0.15:
                        analysis['distinction_quality'] = 'good'
                    elif score_gap > 0.05:
                        analysis['distinction_quality'] = 'fair'
                    else:
                        analysis['distinction_quality'] = 'poor'
                else:
                    analysis['distinction_quality'] = 'single_candidate'
                
                # 计算整体置信度
                avg_confidence = sum(score['confidence'] for score in scores_with_details) / len(scores_with_details)
                analysis['confidence'] = avg_confidence
            
        except Exception as e:
            logging.error(f"异构体区分验证失败: {e}")
            analysis['error'] = str(e)
        
        return analysis
    
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