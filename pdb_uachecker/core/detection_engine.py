"""
检测引擎
整合异构体检测、结构匹配和数据库查询的主引擎
"""

import time
from typing import List, Dict, Optional
import logging

from .models_simple import (
    ResidueInfo, AminoAcidInfo, VerificationResult, VerificationScore, 
    VerificationMethod, IsomerDetectionResult, MatchResult, AnalysisResult
)
from .isomer_detector import IsomerDetector, IsomerMatch
from .structure_matcher import StructureMatcher, StereochemistryAnalyzer
from .database_simple import get_database


class DetectionEngine:
    """
    检测引擎
    核心功能：给定PDB残基，识别数据库中的非天然氨基酸
    """
    
    def __init__(self, db_path: Optional[str] = None):
        """
        初始化检测引擎
        
        Args:
            db_path: 数据库路径
        """
        self.database = get_database(db_path)
        self.isomer_detector = IsomerDetector()
        self.structure_matcher = StructureMatcher()
        self.stereochemistry_analyzer = StereochemistryAnalyzer()
        
        # 检测阈值（降低用于测试）
        self.min_confidence_threshold = 0.3
        self.high_confidence_threshold = 0.7
    
    def detect_amino_acid(self, residue: ResidueInfo) -> IsomerDetectionResult:
        """
        检测单个残基对应的非天然氨基酸
        
        Args:
            residue: PDB残基信息
            
        Returns:
            检测结果
        """
        start_time = time.time()
        
        try:
            # 1. 从数据库获取候选氨基酸
            candidates = self._get_candidates(residue)
            
            if not candidates:
                logging.info(f"残基 {residue.residue_name} 没有找到候选氨基酸")
                return IsomerDetectionResult(
                    residue_info=residue,
                    detection_confidence=0.0,
                    analysis_time=time.time() - start_time
                )
            
            # 2. 对每个候选进行异构体检测
            matches = []
            for candidate in candidates:
                isomer_match = self.isomer_detector.detect_isomer(residue, candidate)
                
                # 创建详细的匹配结果
                match_result = self._create_match_result(isomer_match, candidate, residue)
                matches.append(match_result)
            
            # 3. 按置信度排序
            matches.sort(key=lambda x: x.confidence_score, reverse=True)
            
            # 4. 选择最佳匹配
            best_match = matches[0] if matches and matches[0].confidence_score >= self.min_confidence_threshold else None
            
            detection_confidence = best_match.confidence_score if best_match else 0.0
            
            return IsomerDetectionResult(
                residue_info=residue,
                best_match=best_match,
                all_candidates=matches,
                detection_confidence=detection_confidence,
                analysis_time=time.time() - start_time
            )
            
        except Exception as e:
            logging.error(f"残基检测失败 {residue.residue_name}: {e}")
            return IsomerDetectionResult(
                residue_info=residue,
                detection_confidence=0.0,
                analysis_time=time.time() - start_time
            )
    
    def analyze_pdb_file(self, pdb_residues: List[ResidueInfo], 
                        pdb_filename: str = "unknown") -> AnalysisResult:
        """
        分析整个PDB文件中的非天然氨基酸
        
        Args:
            pdb_residues: PDB中的残基列表
            pdb_filename: PDB文件名
            
        Returns:
            分析结果
        """
        start_time = time.time()
        detected_nna = []
        
        logging.info(f"开始分析PDB文件: {pdb_filename}, 残基数量: {len(pdb_residues)}")
        
        for i, residue in enumerate(pdb_residues):
            logging.info(f"正在处理残基 {i+1}/{len(pdb_residues)}: {residue.residue_name}")
            
            detection_result = self.detect_amino_acid(residue)
            
            # 只保存有检测结果的残基
            if detection_result.best_match or detection_result.detection_confidence > 0:
                detected_nna.append(detection_result)
        
        analysis_time = time.time() - start_time
        
        return AnalysisResult(
            pdb_file=pdb_filename,
            total_residues=len(pdb_residues),
            detected_nna=detected_nna,
            analysis_time=analysis_time
        )
    
    def _get_candidates(self, residue: ResidueInfo) -> List[AminoAcidInfo]:
        """获取候选氨基酸"""
        # 基于分子式查找候选
        candidates = self.database.get_candidates_for_residue(
            molecular_formula=residue.molecular_formula,
            molecular_weight=self._estimate_molecular_weight(residue),
            weight_tolerance=2.0
        )
        
        logging.debug(f"为残基 {residue.residue_name} 找到 {len(candidates)} 个候选氨基酸")
        
        return candidates
    
    def _estimate_molecular_weight(self, residue: ResidueInfo) -> float:
        """估算分子量"""
        # 简化的原子量映射
        atomic_weights = {
            'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999,
            'S': 32.06, 'P': 30.974, 'F': 18.998, 'Cl': 35.45,
            'Br': 79.904, 'I': 126.904
        }
        
        total_weight = 0.0
        for element, count in residue.atom_composition.items():
            weight = atomic_weights.get(element, 12.011)  # 默认使用碳的原子量
            total_weight += weight * count
        
        return total_weight
    
    def _create_match_result(self, isomer_match: IsomerMatch, 
                           candidate: AminoAcidInfo, residue: ResidueInfo) -> MatchResult:
        """创建详细的匹配结果"""
        
        # 创建验证分数
        scores = [
            VerificationScore(
                method=VerificationMethod.FINGERPRINT_SIMILARITY,
                score=isomer_match.fingerprint_similarity,
                passed=isomer_match.fingerprint_similarity >= 0.8
            ),
            VerificationScore(
                method=VerificationMethod.STRUCTURE_3D,
                score=max(0, 1 - isomer_match.structure_3d_rmsd / 2.0),
                passed=isomer_match.structure_3d_rmsd < 1.0
            ),
            VerificationScore(
                method=VerificationMethod.MOLECULAR_FORMULA,
                score=1.0 if residue.molecular_formula == candidate.molecular_formula else 0.0,
                passed=residue.molecular_formula == candidate.molecular_formula
            )
        ]
        
        # 立体化学分析
        stereo_analysis = self.stereochemistry_analyzer.analyze_stereochemistry(residue, candidate)
        
        # 结构分析
        structural_analysis = {
            'isomer_type': isomer_match.isomer_type,
            'smiles_similarity': isomer_match.smiles_similarity,
            'structure_rmsd': isomer_match.structure_3d_rmsd,
            'stereochemistry': stereo_analysis,
            'molecular_formula_match': residue.molecular_formula == candidate.molecular_formula,
            'heavy_atom_count_match': residue.heavy_atom_count == len(candidate.heavy_atom_composition)
        }
        
        verification_result = VerificationResult(
            amino_acid_id=candidate.id,
            amino_acid_name=candidate.name,
            scores=scores,
            overall_confidence=isomer_match.overall_confidence,
            isomer_type=isomer_match.isomer_type
        )
        
        return MatchResult(
            amino_acid_info=candidate,
            verification_result=verification_result,
            structural_analysis=structural_analysis
        )
    
    def get_detection_statistics(self, results: List[IsomerDetectionResult]) -> Dict[str, any]:
        """获取检测统计信息"""
        total_residues = len(results)
        successful_detections = len([r for r in results if r.best_match])
        high_confidence_detections = len([r for r in results 
                                        if r.best_match and r.best_match.is_high_confidence])
        
        # 异构体类型统计
        isomer_types = {}
        for result in results:
            if result.best_match:
                isomer_type = result.best_match.verification_result.isomer_type
                isomer_types[isomer_type] = isomer_types.get(isomer_type, 0) + 1
        
        # 平均检测时间
        avg_time = sum(r.analysis_time for r in results) / len(results) if results else 0
        
        return {
            'total_residues': total_residues,
            'successful_detections': successful_detections,
            'detection_rate': successful_detections / total_residues * 100 if total_residues > 0 else 0,
            'high_confidence_detections': high_confidence_detections,
            'high_confidence_rate': high_confidence_detections / total_residues * 100 if total_residues > 0 else 0,
            'isomer_type_distribution': isomer_types,
            'average_detection_time': avg_time,
            'total_analysis_time': sum(r.analysis_time for r in results)
        }
    
    def update_thresholds(self, min_confidence: float = None, 
                         high_confidence: float = None,
                         rmsd_threshold: float = None,
                         fingerprint_threshold: float = None):
        """更新检测阈值"""
        if min_confidence is not None:
            self.min_confidence_threshold = min_confidence
            
        if high_confidence is not None:
            self.high_confidence_threshold = high_confidence
            
        if rmsd_threshold is not None:
            self.structure_matcher.rmsd_threshold = rmsd_threshold
            
        if fingerprint_threshold is not None:
            self.isomer_detector.fingerprint_threshold = fingerprint_threshold
            
        logging.info(f"阈值已更新: min_conf={self.min_confidence_threshold}, "
                    f"high_conf={self.high_confidence_threshold}")
    
    def validate_system(self) -> Dict[str, bool]:
        """验证系统可用性"""
        validation = {
            'database_available': False,
            'rdkit_available': False,
            'all_components_ready': False
        }
        
        try:
            # 测试数据库
            stats = self.database.get_statistics()
            validation['database_available'] = stats.get('total_amino_acids', 0) > 0
            
            # 测试RDKit
            from rdkit import Chem
            test_mol = Chem.MolFromSmiles('CCO')
            validation['rdkit_available'] = test_mol is not None
            
            validation['all_components_ready'] = (validation['database_available'] and 
                                               validation['rdkit_available'])
            
        except Exception as e:
            logging.error(f"系统验证失败: {e}")
        
        return validation