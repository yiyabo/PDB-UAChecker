"""
简化的数据模型
专注于科研需求，移除工业级的复杂性
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any
from enum import Enum


class VerificationMethod(Enum):
    """验证方法枚举"""
    MOLECULAR_FORMULA = "molecular_formula"
    FINGERPRINT_SIMILARITY = "fingerprint_similarity"
    STRUCTURE_3D = "structure_3d"


@dataclass
class AtomInfo:
    """原子信息 - 简化版"""
    element: str
    x: float
    y: float
    z: float
    atom_name: Optional[str] = None
    formal_charge: int = 0


@dataclass
class ResidueInfo:
    """残基信息 - 简化版"""
    residue_name: str
    residue_number: int
    chain_id: str
    atoms: List[AtomInfo]
    
    @property
    def molecular_formula(self) -> str:
        """计算分子式"""
        composition = {}
        for atom in self.atoms:
            element = atom.element.strip()
            if element:
                composition[element] = composition.get(element, 0) + 1
        
        # 按标准顺序生成分子式 (Hill系统：C, H, 然后按字母顺序)
        formula_parts = []
        
        # 首先C和H
        for element in ['C', 'H']:
            if element in composition:
                count = composition[element]
                if count == 1:
                    formula_parts.append(element)
                else:
                    formula_parts.append(f"{element}{count}")
        
        # 其他元素按字母顺序
        other_elements = sorted([e for e in composition.keys() if e not in ['C', 'H']])
        for element in other_elements:
            count = composition[element]
            if count == 1:
                formula_parts.append(element)
            else:
                formula_parts.append(f"{element}{count}")
        
        return ''.join(formula_parts)
    
    @property
    def atom_composition(self) -> Dict[str, int]:
        """原子组成"""
        composition = {}
        for atom in self.atoms:
            element = atom.element.strip()
            if element:
                composition[element] = composition.get(element, 0) + 1
        return composition
    
    @property
    def heavy_atom_count(self) -> int:
        """重原子数量"""
        return len([atom for atom in self.atoms if atom.element != 'H'])


@dataclass
class AminoAcidInfo:
    """氨基酸信息 - 简化版"""
    id: str
    name: str
    molecular_formula: str
    molecular_weight: float
    smiles: str
    atom_composition: Dict[str, int]
    
    @property
    def heavy_atom_composition(self) -> Dict[str, int]:
        """重原子组成"""
        return {k: v for k, v in self.atom_composition.items() if k != 'H'}


@dataclass
class VerificationScore:
    """验证分数 - 简化版"""
    method: VerificationMethod
    score: float
    passed: bool
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class VerificationResult:
    """验证结果 - 简化版"""
    amino_acid_id: str
    amino_acid_name: str
    scores: List[VerificationScore]
    overall_confidence: float
    isomer_type: str = "unknown"  # "identical", "structural", "stereoisomer", "different"
    
    @property
    def is_match(self) -> bool:
        """是否匹配成功"""
        return self.overall_confidence >= 0.6


@dataclass  
class IsomerDetectionResult:
    """异构体检测结果"""
    residue_info: ResidueInfo
    best_match: Optional['MatchResult'] = None
    all_candidates: List['MatchResult'] = field(default_factory=list)
    detection_confidence: float = 0.0
    analysis_time: float = 0.0


@dataclass
class MatchResult:
    """匹配结果 - 简化版"""
    amino_acid_info: AminoAcidInfo
    verification_result: VerificationResult
    structural_analysis: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def confidence_score(self) -> float:
        """置信度分数"""
        return self.verification_result.overall_confidence
    
    @property
    def is_high_confidence(self) -> bool:
        """是否为高置信度匹配"""
        return self.confidence_score >= 0.8


@dataclass
class AnalysisResult:
    """分析结果 - 简化版"""
    pdb_file: str
    total_residues: int
    detected_nna: List[IsomerDetectionResult]
    analysis_time: float
    
    @property
    def detection_rate(self) -> float:
        """检测率"""
        detected_count = len([r for r in self.detected_nna if r.best_match])
        return (detected_count / self.total_residues * 100) if self.total_residues > 0 else 0.0


# 异常类 - 简化版
class PDBUACheckerError(Exception):
    """基础异常类"""
    pass


class DetectionError(PDBUACheckerError):
    """检测相关异常"""
    pass