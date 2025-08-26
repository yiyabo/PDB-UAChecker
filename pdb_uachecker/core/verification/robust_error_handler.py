"""
鲁棒错误处理机制
提供智能的降级策略、错误恢复和详细诊断
"""

import logging
import traceback
from typing import Dict, List, Optional, Any, Callable, Tuple
from dataclasses import dataclass, field
from enum import Enum
import time
from functools import wraps


class ErrorSeverity(Enum):
    """错误严重程度"""
    CRITICAL = "critical"  # 系统无法继续运行
    HIGH = "high"         # 功能无法正常工作
    MEDIUM = "medium"     # 功能部分受限
    LOW = "low"          # 轻微问题，不影响主要功能
    INFO = "info"        # 信息性提示


class RecoveryStrategy(Enum):
    """恢复策略"""
    RETRY = "retry"              # 重试操作
    FALLBACK = "fallback"        # 使用降级方法
    SKIP = "skip"               # 跳过当前操作
    DEFAULT_VALUE = "default"    # 返回默认值
    USER_INTERVENTION = "user"   # 需要用户干预


@dataclass
class ErrorContext:
    """错误上下文信息"""
    operation_name: str
    input_data_type: str
    input_data_size: Optional[int] = None
    timestamp: float = field(default_factory=time.time)
    component: str = "unknown"
    user_context: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ErrorRecord:
    """错误记录"""
    error_type: str
    error_message: str
    severity: ErrorSeverity
    context: ErrorContext
    stack_trace: str
    recovery_strategy: RecoveryStrategy
    recovery_success: bool = False
    recovery_details: str = ""
    occurrence_count: int = 1
    
    def __post_init__(self):
        self.timestamp = self.context.timestamp


class RobustErrorHandler:
    """鲁棒错误处理器"""
    
    def __init__(self, max_retry_attempts: int = 3, retry_delay: float = 0.1):
        """
        初始化错误处理器
        
        Args:
            max_retry_attempts: 最大重试次数
            retry_delay: 重试延迟（秒）
        """
        self.max_retry_attempts = max_retry_attempts
        self.retry_delay = retry_delay
        self.error_history = []
        self.recovery_handlers = {}
        self.fallback_strategies = {}
        
        # 注册默认恢复处理器
        self._register_default_handlers()
    
    def _register_default_handlers(self):
        """注册默认的错误恢复处理器"""
        
        # RDKit相关错误
        self.register_recovery_handler(
            "rdkit.*", 
            self._handle_rdkit_errors,
            ErrorSeverity.HIGH
        )
        
        # 价态验证错误
        self.register_recovery_handler(
            "valence.*",
            self._handle_valence_errors,
            ErrorSeverity.MEDIUM
        )
        
        # 数据格式错误
        self.register_recovery_handler(
            "format.*",
            self._handle_format_errors,
            ErrorSeverity.MEDIUM
        )
        
        # 化学键推断错误
        self.register_recovery_handler(
            "bond_inference.*",
            self._handle_bond_inference_errors,
            ErrorSeverity.MEDIUM
        )
        
        # 内存错误
        self.register_recovery_handler(
            "memory.*",
            self._handle_memory_errors,
            ErrorSeverity.HIGH
        )
    
    def register_recovery_handler(self, error_pattern: str, handler: Callable, severity: ErrorSeverity):
        """
        注册错误恢复处理器
        
        Args:
            error_pattern: 错误模式（支持通配符）
            handler: 处理函数
            severity: 错误严重程度
        """
        self.recovery_handlers[error_pattern] = {
            'handler': handler,
            'severity': severity
        }
    
    def register_fallback_strategy(self, operation_name: str, fallback_func: Callable):
        """
        注册降级策略
        
        Args:
            operation_name: 操作名称
            fallback_func: 降级函数
        """
        self.fallback_strategies[operation_name] = fallback_func
    
    def safe_execute(self, func: Callable, *args, context: ErrorContext = None, **kwargs) -> Tuple[Any, Optional[ErrorRecord]]:
        """
        安全执行函数，带错误处理和恢复
        
        Args:
            func: 要执行的函数
            *args: 位置参数
            context: 错误上下文
            **kwargs: 关键字参数
            
        Returns:
            (result, error_record) 元组
        """
        if context is None:
            context = ErrorContext(
                operation_name=func.__name__,
                input_data_type=str(type(args[0]) if args else "None"),
                component="unknown"
            )
        
        last_error_record = None
        
        # 尝试多次执行
        for attempt in range(self.max_retry_attempts + 1):
            try:
                result = func(*args, **kwargs)
                
                # 执行成功，更新恢复状态
                if last_error_record:
                    last_error_record.recovery_success = True
                    last_error_record.recovery_details = f"重试{attempt}次后成功"
                
                return result, last_error_record
                
            except Exception as e:
                error_record = self._create_error_record(e, context)
                last_error_record = error_record
                
                # 记录错误
                self._record_error(error_record)
                
                # 尝试恢复
                recovery_result = self._attempt_recovery(error_record, func, args, kwargs)
                
                if recovery_result is not None:
                    error_record.recovery_success = True
                    return recovery_result, error_record
                
                # 如果不是最后一次尝试，继续重试
                if attempt < self.max_retry_attempts:
                    logging.warning(f"第{attempt + 1}次尝试失败，{self.retry_delay}秒后重试: {e}")
                    time.sleep(self.retry_delay)
                    continue
                else:
                    logging.error(f"所有重试尝试失败: {e}")
                    break
        
        # 所有尝试都失败，尝试降级策略
        fallback_result = self._try_fallback_strategy(context, args, kwargs)
        if fallback_result is not None:
            last_error_record.recovery_success = True
            last_error_record.recovery_strategy = RecoveryStrategy.FALLBACK
            last_error_record.recovery_details = "使用降级策略成功"
            return fallback_result, last_error_record
        
        # 完全失败
        return None, last_error_record
    
    def _create_error_record(self, error: Exception, context: ErrorContext) -> ErrorRecord:
        """创建错误记录"""
        error_type = type(error).__name__
        error_message = str(error)
        stack_trace = traceback.format_exc()
        
        # 确定错误严重程度和恢复策略
        severity, recovery_strategy = self._analyze_error(error, context)
        
        return ErrorRecord(
            error_type=error_type,
            error_message=error_message,
            severity=severity,
            context=context,
            stack_trace=stack_trace,
            recovery_strategy=recovery_strategy
        )
    
    def _analyze_error(self, error: Exception, context: ErrorContext) -> Tuple[ErrorSeverity, RecoveryStrategy]:
        """分析错误类型并确定严重程度和恢复策略"""
        error_type = type(error).__name__.lower()
        error_message = str(error).lower()
        
        # RDKit错误
        if "rdkit" in error_message or "sanitization" in error_message:
            return ErrorSeverity.HIGH, RecoveryStrategy.FALLBACK
        
        # 价态错误
        if "valence" in error_message or "bond" in error_message:
            return ErrorSeverity.MEDIUM, RecoveryStrategy.RETRY
        
        # 数据格式错误
        if error_type in ['keyerror', 'attributeerror', 'typeerror']:
            return ErrorSeverity.MEDIUM, RecoveryStrategy.DEFAULT_VALUE
        
        # 内存错误
        if error_type in ['memoryerror', 'outofmemoryerror']:
            return ErrorSeverity.HIGH, RecoveryStrategy.FALLBACK
        
        # 值错误
        if error_type in ['valueerror', 'indexerror']:
            return ErrorSeverity.MEDIUM, RecoveryStrategy.DEFAULT_VALUE
        
        # 默认情况
        return ErrorSeverity.MEDIUM, RecoveryStrategy.RETRY
    
    def _attempt_recovery(self, error_record: ErrorRecord, func: Callable, args: tuple, kwargs: dict) -> Any:
        """尝试错误恢复"""
        error_pattern = f"{error_record.context.component}.{error_record.error_type.lower()}"
        
        # 寻找匹配的恢复处理器
        for pattern, handler_info in self.recovery_handlers.items():
            if self._pattern_matches(pattern, error_pattern):
                try:
                    logging.info(f"尝试使用恢复处理器: {pattern}")
                    result = handler_info['handler'](error_record, func, args, kwargs)
                    if result is not None:
                        error_record.recovery_details = f"使用处理器{pattern}恢复成功"
                        return result
                except Exception as recovery_error:
                    logging.warning(f"恢复处理器{pattern}执行失败: {recovery_error}")
        
        return None
    
    def _pattern_matches(self, pattern: str, text: str) -> bool:
        """简单的模式匹配（支持*通配符）"""
        import re
        pattern_regex = pattern.replace('*', '.*')
        return bool(re.match(pattern_regex, text))
    
    def _try_fallback_strategy(self, context: ErrorContext, args: tuple, kwargs: dict) -> Any:
        """尝试降级策略"""
        operation_name = context.operation_name
        
        if operation_name in self.fallback_strategies:
            try:
                logging.info(f"尝试降级策略: {operation_name}")
                fallback_func = self.fallback_strategies[operation_name]
                return fallback_func(*args, **kwargs)
            except Exception as fallback_error:
                logging.error(f"降级策略执行失败: {fallback_error}")
        
        return None
    
    def _record_error(self, error_record: ErrorRecord):
        """记录错误到历史中"""
        # 检查是否是重复错误
        for existing_record in self.error_history:
            if (existing_record.error_type == error_record.error_type and
                existing_record.context.operation_name == error_record.context.operation_name):
                existing_record.occurrence_count += 1
                return
        
        # 新错误，添加到历史中
        self.error_history.append(error_record)
        
        # 限制历史记录长度
        if len(self.error_history) > 1000:
            self.error_history = self.error_history[-500:]  # 保留最近500条
    
    def _handle_rdkit_errors(self, error_record: ErrorRecord, func: Callable, args: tuple, kwargs: dict) -> Any:
        """处理RDKit相关错误"""
        error_message = error_record.error_message.lower()
        
        # 价态超限错误
        if "valence" in error_message:
            logging.info("检测到RDKit价态错误，尝试放宽分子清理选项")
            # 如果是分子清理错误，尝试部分清理
            try:
                # 假设args[0]是分子对象
                if args and hasattr(args[0], 'GetNumAtoms'):
                    from rdkit import Chem
                    mol = args[0]
                    Chem.SanitizeMol(mol, sanitizeOps=(
                        Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                        Chem.SanitizeFlags.SANITIZE_KEKULIZE
                    ))
                    return mol
            except:
                pass
        
        # 键重复错误
        if "bond already exists" in error_message:
            logging.info("检测到键重复错误，使用保守键推断")
            # 返回信号，让调用者使用保守策略
            return "USE_CONSERVATIVE_BONDING"
        
        return None
    
    def _handle_valence_errors(self, error_record: ErrorRecord, func: Callable, args: tuple, kwargs: dict) -> Any:
        """处理价态验证错误"""
        logging.info("尝试自动修正价态问题")
        
        # 如果有价态验证器可用，尝试自动修正
        try:
            if args and hasattr(args[0], 'GetNumAtoms'):
                from .valence_validator import SmartValenceValidator
                validator = SmartValenceValidator(tolerance_level="loose")
                mol, corrections = validator.apply_automatic_corrections(args[0], [])
                if corrections:
                    logging.info(f"应用了价态修正: {corrections}")
                    return mol
        except:
            pass
        
        return None
    
    def _handle_format_errors(self, error_record: ErrorRecord, func: Callable, args: tuple, kwargs: dict) -> Any:
        """处理数据格式错误"""
        error_type = error_record.error_type
        
        # KeyError处理
        if error_type == "KeyError":
            logging.info("检测到数据字段缺失，提供默认值")
            # 根据上下文提供合理的默认值
            operation_name = error_record.context.operation_name
            if "coordinate" in operation_name.lower():
                return {"x": 0.0, "y": 0.0, "z": 0.0, "element": "C"}
            elif "atom" in operation_name.lower():
                return None  # 跳过这个原子
        
        # AttributeError处理
        if error_type == "AttributeError":
            logging.info("检测到属性访问错误，使用兼容性包装器")
            # 使用统一数据处理器
            try:
                from .unified_data_processor import UnifiedDataProcessor
                processor = UnifiedDataProcessor()
                if args:
                    return processor.standardize_atom(args[0])
            except:
                pass
        
        return None
    
    def _handle_bond_inference_errors(self, error_record: ErrorRecord, func: Callable, args: tuple, kwargs: dict) -> Any:
        """处理化学键推断错误"""
        logging.info("化学键推断失败，尝试保守策略")
        
        # 使用保守的化学键推断
        try:
            if args and len(args) >= 2:
                mol, atoms_info = args[0], args[1]
                from .advanced_bond_inference import AdvancedBondInference
                bond_inference = AdvancedBondInference()
                result = bond_inference._conservative_inference(mol, atoms_info)
                return result
        except:
            pass
        
        return None
    
    def _handle_memory_errors(self, error_record: ErrorRecord, func: Callable, args: tuple, kwargs: dict) -> Any:
        """处理内存错误"""
        logging.warning("检测到内存错误，尝试减少计算复杂度")
        
        # 减少计算量的策略
        if error_record.context.input_data_size and error_record.context.input_data_size > 100:
            logging.info("输入数据较大，尝试简化处理")
            # 返回信号让调用者使用简化算法
            return "USE_SIMPLIFIED_ALGORITHM"
        
        return None
    
    def get_error_statistics(self) -> Dict[str, Any]:
        """获取错误统计信息"""
        if not self.error_history:
            return {"total_errors": 0}
        
        stats = {
            "total_errors": len(self.error_history),
            "by_severity": {},
            "by_type": {},
            "by_component": {},
            "recovery_success_rate": 0.0,
            "most_common_errors": []
        }
        
        # 按严重程度统计
        for record in self.error_history:
            severity = record.severity.value
            stats["by_severity"][severity] = stats["by_severity"].get(severity, 0) + record.occurrence_count
        
        # 按错误类型统计
        for record in self.error_history:
            error_type = record.error_type
            stats["by_type"][error_type] = stats["by_type"].get(error_type, 0) + record.occurrence_count
        
        # 按组件统计
        for record in self.error_history:
            component = record.context.component
            stats["by_component"][component] = stats["by_component"].get(component, 0) + record.occurrence_count
        
        # 计算恢复成功率
        total_recoveries = sum(1 for record in self.error_history if record.recovery_success)
        stats["recovery_success_rate"] = total_recoveries / len(self.error_history) if self.error_history else 0.0
        
        # 最常见错误
        error_counts = [(record.error_type, record.occurrence_count) for record in self.error_history]
        stats["most_common_errors"] = sorted(error_counts, key=lambda x: x[1], reverse=True)[:5]
        
        return stats
    
    def generate_diagnostic_report(self) -> str:
        """生成诊断报告"""
        stats = self.get_error_statistics()
        
        report = f"""
🔧 错误处理诊断报告
{'=' * 50}

📊 总体统计:
  • 总错误数: {stats['total_errors']}
  • 恢复成功率: {stats['recovery_success_rate']:.1%}

🚨 按严重程度分布:
"""
        
        for severity, count in stats.get("by_severity", {}).items():
            report += f"  • {severity.upper()}: {count}\n"
        
        report += "\n🔍 最常见错误:\n"
        for error_type, count in stats.get("most_common_errors", []):
            report += f"  • {error_type}: {count}次\n"
        
        report += f"\n📈 组件错误分布:\n"
        for component, count in stats.get("by_component", {}).items():
            report += f"  • {component}: {count}\n"
        
        # 添加建议
        report += self._generate_recommendations(stats)
        
        return report
    
    def _generate_recommendations(self, stats: Dict) -> str:
        """基于统计信息生成改进建议"""
        recommendations = "\n💡 改进建议:\n"
        
        # 基于错误类型的建议
        common_errors = dict(stats.get("most_common_errors", []))
        
        if "ValueError" in common_errors and common_errors["ValueError"] > 5:
            recommendations += "  • 考虑增强输入数据验证\n"
        
        if "KeyError" in common_errors and common_errors["KeyError"] > 3:
            recommendations += "  • 改进数据格式标准化\n"
        
        if "AttributeError" in common_errors:
            recommendations += "  • 使用统一数据接口\n"
        
        # 基于严重程度的建议
        severity_stats = stats.get("by_severity", {})
        if severity_stats.get("high", 0) > 2:
            recommendations += "  • 优先修复高严重程度错误\n"
        
        if stats.get("recovery_success_rate", 0) < 0.7:
            recommendations += "  • 改进错误恢复机制\n"
        
        return recommendations


# 装饰器：自动错误处理
def robust_operation(component: str = "unknown", max_retries: int = 3):
    """
    装饰器：为函数添加鲁棒错误处理
    
    Args:
        component: 组件名称
        max_retries: 最大重试次数
    """
    def decorator(func: Callable):
        @wraps(func)
        def wrapper(*args, **kwargs):
            error_handler = RobustErrorHandler(max_retry_attempts=max_retries)
            
            context = ErrorContext(
                operation_name=func.__name__,
                input_data_type=str(type(args[0]) if args else "None"),
                input_data_size=len(args[0]) if args and hasattr(args[0], '__len__') else None,
                component=component
            )
            
            result, error_record = error_handler.safe_execute(func, *args, context=context, **kwargs)
            
            if error_record and not error_record.recovery_success:
                logging.error(f"函数 {func.__name__} 执行失败: {error_record.error_message}")
            
            return result
        
        return wrapper
    return decorator