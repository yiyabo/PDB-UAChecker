"""
验证器基类
定义验证器的通用接口和行为
"""

from abc import ABC, abstractmethod
from typing import Dict, Any

from ..models import ResidueInfo, AminoAcidInfo, VerificationScore, VerificationMethod


class BaseVerifier(ABC):
    """验证器基类"""
    
    def __init__(self, threshold: float):
        self.threshold = threshold
    
    @abstractmethod
    def get_method(self) -> VerificationMethod:
        """获取验证方法类型"""
        pass
    
    @abstractmethod
    def calculate_score(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> float:
        """
        计算验证分数
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            验证分数 (0.0-1.0)
        """
        pass
    
    def verify(self, residue: ResidueInfo, amino_acid: AminoAcidInfo) -> VerificationScore:
        """
        执行验证
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
        
        Returns:
            验证结果
        """
        try:
            score = self.calculate_score(residue, amino_acid)
            passed = score >= self.threshold
            details = self._get_verification_details(residue, amino_acid, score)
            
            return VerificationScore(
                method=self.get_method(),
                score=score,
                passed=passed,
                details=details
            )
        
        except Exception as e:
            # 验证失败时返回0分
            return VerificationScore(
                method=self.get_method(),
                score=0.0,
                passed=False,
                details={'error': str(e)}
            )
    
    def _get_verification_details(self, residue: ResidueInfo, amino_acid: AminoAcidInfo, score: float) -> Dict[str, Any]:
        """
        获取验证详情（子类可重写）
        
        Args:
            residue: 残基信息
            amino_acid: 氨基酸信息
            score: 计算得到的分数
        
        Returns:
            验证详情字典
        """
        return {
            'threshold': self.threshold,
            'score': score,
            'passed': score >= self.threshold
        }
    
    def set_threshold(self, threshold: float):
        """设置阈值"""
        if not 0.0 <= threshold <= 1.0:
            raise ValueError("阈值必须在0.0-1.0之间")
        self.threshold = threshold
    
    def get_threshold(self) -> float:
        """获取阈值"""
        return self.threshold