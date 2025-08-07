"""
数据模型定义
定义系统中使用的核心数据结构
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any, Tuple
from enum import Enum

# 在模块顶层导入，避免运行时导入的性能开销
try:
    from ..utils.chemistry import ChemistryUtils
except ImportError:
    # 兼容性处理，如果导入失败则延迟导入
    ChemistryUtils = None


class VerificationMethod(Enum):
    """验证方法枚举"""
    RESIDUE_NAME = "residue_name"
    MOLECULAR_FORMULA = "molecular_formula"
    ATOM_COMPOSITION = "atom_composition"
    FINGERPRINT_SIMILARITY = "fingerprint_similarity"
    STRUCTURE_3D = "structure_3d"


@dataclass
class AtomInfo:
    """原子信息"""
    atom_name: str
    element: str
    x: float
    y: float
    z: float
    residue_name: str
    residue_number: int
    chain_id: str


@dataclass
class ResidueInfo:
    """残基信息"""
    residue_name: str
    residue_number: int
    chain_id: str
    atoms: List[AtomInfo]
    molecular_formula: str = field(default="", init=False)
    atom_composition: Dict[str, int] = field(default_factory=dict, init=False)
    
    def __post_init__(self):
        """自动计算分子式和原子组成"""
        self._calculate_molecular_properties()
    
    def _calculate_molecular_properties(self):
        """计算分子属性"""
        # 使用模块级别的ChemistryUtils，避免运行时导入
        if ChemistryUtils is None:
            # 如果模块导入时失败，这里再次尝试导入
            try:
                from ..utils.chemistry import ChemistryUtils as CU
                chemistry_utils = CU
            except ImportError as e:
                print(f"⚠️ ChemistryUtils导入失败: {e}")
                self.molecular_formula = ""
                self.atom_composition = {}
                return
        else:
            chemistry_utils = ChemistryUtils
        
        try:
            # 确保atoms列表存在且不为空
            if not self.atoms:  # 简化检查，不使用hasattr
                self.molecular_formula = ""
                self.atom_composition = {}
                return
            
            # 转换原子信息为字典格式
            atom_dicts = []
            for atom in self.atoms:
                if hasattr(atom, 'element') and atom.element and atom.element.strip():
                    atom_dicts.append({
                        'element': atom.element.strip(),
                        'atom_name': getattr(atom, 'atom_name', ''),
                        'x': getattr(atom, 'x', 0.0),
                        'y': getattr(atom, 'y', 0.0),
                        'z': getattr(atom, 'z', 0.0)
                    })
            
            if atom_dicts:
                self.atom_composition = chemistry_utils.calculate_atom_composition(atom_dicts)
                self.molecular_formula = chemistry_utils.calculate_molecular_formula(self.atom_composition)
            else:
                self.molecular_formula = ""
                self.atom_composition = {}
                
        except Exception as e:
            print(f"⚠️ 分子式计算异常: {e}")
            self.molecular_formula = ""
            self.atom_composition = {}
    
    @property
    def residue_key(self) -> str:
        """残基唯一标识"""
        return f"{self.chain_id}:{self.residue_name}{self.residue_number}"
    
    @property
    def heavy_atom_composition(self) -> Dict[str, int]:
        """重原子组成（不含氢）"""
        if ChemistryUtils is None:
            # 如果ChemistryUtils不可用，手动过滤氢原子
            return {k: v for k, v in (self.atom_composition or {}).items() if k != 'H'}
        return ChemistryUtils.extract_heavy_atoms(self.atom_composition or {})
    
    @property
    def coordinates(self) -> List[List[float]]:
        """原子坐标列表"""
        return [[atom.x, atom.y, atom.z] for atom in self.atoms]


@dataclass
class AminoAcidInfo:
    """氨基酸信息"""
    id: str
    name: str
    molecular_formula: str
    molecular_weight: float
    smiles: str
    atom_composition: Dict[str, int]
    key_features: List[str] = field(default_factory=list)
    fingerprints: Dict[str, str] = field(default_factory=dict)
    
    @property
    def heavy_atom_composition(self) -> Dict[str, int]:
        """重原子组成（不含氢）"""
        if ChemistryUtils is None:
            # 如果ChemistryUtils不可用，手动过滤氢原子
            return {k: v for k, v in self.atom_composition.items() if k != 'H'}
        return ChemistryUtils.extract_heavy_atoms(self.atom_composition)


@dataclass
class VerificationScore:
    """单项验证分数"""
    method: VerificationMethod
    score: float
    passed: bool
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class VerificationResult:
    """验证结果"""
    amino_acid_id: str
    amino_acid_name: str
    scores: List[VerificationScore]
    overall_confidence: float
    passed_verifications: int
    total_verifications: int
    
    @property
    def is_match(self) -> bool:
        """是否匹配成功"""
        return self.overall_confidence > 0.0
    
    @property
    def verification_summary(self) -> Dict[str, bool]:
        """验证摘要"""
        return {score.method.value: score.passed for score in self.scores}


@dataclass
class MatchResult:
    """匹配结果"""
    amino_acid_info: AminoAcidInfo
    residue_info: ResidueInfo
    verification_result: VerificationResult
    match_method: str
    
    @property
    def confidence_score(self) -> float:
        """置信度分数"""
        return self.verification_result.overall_confidence
    
    @property
    def amino_acid_id(self) -> str:
        """氨基酸ID"""
        return self.amino_acid_info.id
    
    @property
    def amino_acid_name(self) -> str:
        """氨基酸名称"""
        return self.amino_acid_info.name


@dataclass
class AnalysisResult:
    """分析结果"""
    pdb_file: str
    total_residues: int
    identified_residues: int
    matches: List[MatchResult]
    analysis_time: float
    capabilities: Dict[str, bool] = field(default_factory=dict)
    
    @property
    def identification_rate(self) -> float:
        """识别率"""
        return (self.identified_residues / self.total_residues * 100) if self.total_residues > 0 else 0.0
    
    @property
    def summary(self) -> Dict[str, Any]:
        """分析摘要"""
        return {
            'pdb_file': self.pdb_file,
            'total_residues': self.total_residues,
            'identified_residues': self.identified_residues,
            'identification_rate': self.identification_rate,
            'analysis_time': self.analysis_time,
            'capabilities': self.capabilities
        }


@dataclass
class Structure3DInfo:
    """3D结构信息"""
    coordinates: List[List[float]]
    elements: List[str]
    rmsd: Optional[float] = None
    kabsch_score: Optional[float] = None
    
    @property
    def atom_count(self) -> int:
        """原子数量"""
        return len(self.coordinates)


@dataclass
class ClassificationResult:
    """高精度分类结果"""
    amino_acid_id: str
    amino_acid_name: str
    categories: List[str]
    confidence: float
    classification_method: str
    evidence: List[str]
    details: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def is_high_confidence(self) -> bool:
        """是否为高置信度分类"""
        return self.confidence >= 0.9
    
    @property
    def requires_review(self) -> bool:
        """是否需要人工审查"""
        return self.confidence < 0.8 or 'requires_manual_review' in self.categories


@dataclass
class ChemicalKnowledge:
    """化学知识库条目"""
    amino_acid_id: str
    canonical_smiles: str
    stereochemistry: Optional[str] = None  # 'D', 'L', 或 None
    backbone_type: Optional[str] = None    # 'alpha', 'beta', 'gamma'
    functional_groups: List[str] = field(default_factory=list)
    structural_features: List[str] = field(default_factory=list)
    alternative_names: List[str] = field(default_factory=list)
    confidence: float = 1.0
    source: str = "curated_database"


@dataclass
class MolecularAnalysis:
    """分子结构分析结果"""
    smiles: str
    is_valid: bool
    aromatic_atoms: List[int] = field(default_factory=list)
    ring_systems: List[List[int]] = field(default_factory=list)
    chiral_centers: List[Tuple[int, str]] = field(default_factory=list)
    functional_groups: Dict[str, List[int]] = field(default_factory=dict)
    backbone_analysis: Optional[Dict[str, Any]] = None
    molecular_descriptors: Dict[str, float] = field(default_factory=dict)


# 异常类定义
class PDBUACheckerError(Exception):
    """基础异常类"""
    pass


class DatabaseError(PDBUACheckerError):
    """数据库相关异常"""
    pass


class VerificationError(PDBUACheckerError):
    """验证相关异常"""
    pass


class ParsingError(PDBUACheckerError):
    """解析相关异常"""
    pass