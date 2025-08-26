"""
智能价态验证器
处理原子价态验证、异常检测和自动修正
"""

import logging
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass
from collections import defaultdict

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    # 创建空的Chem类用于类型提示
    class Chem:
        class Mol:
            pass


@dataclass
class ValenceRule:
    """价态规则"""
    element: str
    common_valences: List[int]
    max_valence: int
    formal_charge_range: Tuple[int, int] = (-3, 3)
    special_cases: Dict[str, int] = None  # 特殊情况下的价态


@dataclass
class ValenceIssue:
    """价态问题"""
    atom_idx: int
    element: str
    current_valence: float
    expected_valences: List[int]
    severity: str  # "error", "warning", "info"
    suggestion: str = ""


class SmartValenceValidator:
    """智能价态验证器"""
    
    # 扩展的价态规则表
    VALENCE_RULES = {
        'H': ValenceRule('H', [1], 1, (-1, 1)),
        'C': ValenceRule('C', [4], 4, (-4, 4), {'aromatic': 4, 'radical': 3}),
        'N': ValenceRule('N', [3, 5], 5, (-3, 5), {'aromatic': 3, 'quaternary': 4}),
        'O': ValenceRule('O', [2], 6, (-2, 2), {'peroxide': 1, 'radical': 1}),
        'F': ValenceRule('F', [1], 1, (-1, 0)),
        'P': ValenceRule('P', [3, 5], 5, (-3, 5), {'phosphonate': 4}),
        'S': ValenceRule('S', [2, 4, 6], 6, (-2, 6), {'sulfoxide': 4, 'sulfone': 6}),
        'Cl': ValenceRule('Cl', [1, 3, 5, 7], 7, (-1, 7)),
        'Br': ValenceRule('Br', [1, 3, 5, 7], 7, (-1, 7)),
        'I': ValenceRule('I', [1, 3, 5, 7], 7, (-1, 7)),
        'Se': ValenceRule('Se', [2, 4, 6], 6, (-2, 6)),
    }
    
    # 常见分子片段的价态模式
    FRAGMENT_PATTERNS = {
        'carboxyl': {'C': 4, 'O': 2, 'O': 2},  # -COOH
        'amino': {'N': 3},  # -NH2
        'hydroxyl': {'O': 2},  # -OH
        'carbonyl': {'C': 4, 'O': 2},  # C=O
        'phosphate': {'P': 5, 'O': 2, 'O': 2, 'O': 2, 'O': 2},  # PO4
        'sulfate': {'S': 6, 'O': 2, 'O': 2, 'O': 2, 'O': 2},  # SO4
    }
    
    def __init__(self, tolerance_level: str = "moderate"):
        """
        初始化价态验证器
        
        Args:
            tolerance_level: 容忍级别 ("strict", "moderate", "loose")
        """
        self.tolerance_level = tolerance_level
        self.tolerance_factors = {
            "strict": 0.1,
            "moderate": 0.2,
            "loose": 0.3
        }
        self.current_tolerance = self.tolerance_factors.get(tolerance_level, 0.2)
    
    def validate_molecule_valences(self, mol) -> Tuple[bool, List[ValenceIssue], Dict]:
        """
        全面验证分子中的原子价态
        
        Args:
            mol: RDKit分子对象
            
        Returns:
            (is_valid, issues_list, analysis_report)
        """
        if not RDKIT_AVAILABLE or mol is None:
            return False, [], {"error": "RDKit不可用或分子为空"}
        
        issues = []
        analysis_report = {
            "total_atoms": mol.GetNumAtoms(),
            "validated_atoms": 0,
            "warnings": 0,
            "errors": 0,
            "valence_distribution": defaultdict(int),
            "unusual_patterns": []
        }
        
        try:
            for atom_idx in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(atom_idx)
                atom_issues = self._validate_single_atom(mol, atom_idx)
                issues.extend(atom_issues)
                
                # 更新统计信息
                analysis_report["validated_atoms"] += 1
                for issue in atom_issues:
                    if issue.severity == "error":
                        analysis_report["errors"] += 1
                    elif issue.severity == "warning":
                        analysis_report["warnings"] += 1
                
                # 记录价态分布
                actual_valence = self._calculate_atom_valence(atom)
                element = atom.GetSymbol()
                analysis_report["valence_distribution"][f"{element}({actual_valence})"] += 1
            
            # 检测不寻常的分子模式
            unusual_patterns = self._detect_unusual_patterns(mol)
            analysis_report["unusual_patterns"] = unusual_patterns
            
            # 判断整体有效性
            error_count = sum(1 for issue in issues if issue.severity == "error")
            is_valid = error_count == 0
            
            return is_valid, issues, analysis_report
            
        except Exception as e:
            logging.error(f"价态验证失败: {e}")
            return False, [], {"error": str(e)}
    
    def _validate_single_atom(self, mol, atom_idx: int) -> List[ValenceIssue]:
        """
        验证单个原子的价态
        
        Args:
            mol: 分子对象
            atom_idx: 原子索引
            
        Returns:
            问题列表
        """
        issues = []
        atom = mol.GetAtomWithIdx(atom_idx)
        element = atom.GetSymbol()
        
        # 获取价态规则
        valence_rule = self.VALENCE_RULES.get(element)
        if not valence_rule:
            return issues  # 未知元素，跳过验证
        
        # 计算实际价态
        actual_valence = self._calculate_atom_valence(atom)
        formal_charge = atom.GetFormalCharge()
        
        # 调整后的有效价态
        effective_valence = actual_valence - formal_charge
        
        # 检查是否在允许范围内
        is_common = effective_valence in valence_rule.common_valences
        is_within_max = effective_valence <= valence_rule.max_valence
        is_formal_charge_ok = (valence_rule.formal_charge_range[0] <= 
                              formal_charge <= valence_rule.formal_charge_range[1])
        
        if not is_within_max:
            # 严重错误：超过最大价态
            issues.append(ValenceIssue(
                atom_idx=atom_idx,
                element=element,
                current_valence=actual_valence,
                expected_valences=valence_rule.common_valences,
                severity="error",
                suggestion=f"价态{actual_valence}超过{element}的最大价态{valence_rule.max_valence}"
            ))
        elif not is_common and not self._is_acceptable_special_case(atom, mol):
            # 警告：不常见的价态但在最大范围内
            issues.append(ValenceIssue(
                atom_idx=atom_idx,
                element=element,
                current_valence=actual_valence,
                expected_valences=valence_rule.common_valences,
                severity="warning",
                suggestion=f"不常见的价态{actual_valence}，常见价态为{valence_rule.common_valences}"
            ))
        
        if not is_formal_charge_ok:
            issues.append(ValenceIssue(
                atom_idx=atom_idx,
                element=element,
                current_valence=actual_valence,
                expected_valences=valence_rule.common_valences,
                severity="warning",
                suggestion=f"形式电荷{formal_charge}超出正常范围{valence_rule.formal_charge_range}"
            ))
        
        return issues
    
    def _calculate_atom_valence(self, atom) -> float:
        """
        计算原子的实际价态
        
        Args:
            atom: 原子对象
            
        Returns:
            价态数值
        """
        try:
            valence = 0.0
            for bond in atom.GetBonds():
                valence += bond.GetBondTypeAsDouble()
            
            # 加上显式氢原子
            valence += atom.GetTotalNumHs()
            
            return valence
        except Exception as e:
            logging.warning(f"价态计算失败: {e}")
            return 0.0
    
    def _is_acceptable_special_case(self, atom, mol) -> bool:
        """
        检查是否是可接受的特殊情况
        
        Args:
            atom: 原子对象
            mol: 分子对象
            
        Returns:
            是否可接受
        """
        element = atom.GetSymbol()
        
        try:
            # 检查芳香性
            if atom.GetIsAromatic():
                return True
            
            # 检查是否在已知的分子片段中
            if self._is_in_known_fragment(atom, mol):
                return True
            
            # 检查是否是自由基
            if self._is_radical_center(atom):
                return True
            
            # 检查是否是配位化合物
            if self._is_coordination_complex(atom, mol):
                return True
            
            return False
            
        except Exception as e:
            logging.warning(f"特殊情况检查失败: {e}")
            return False
    
    def _is_in_known_fragment(self, atom, mol) -> bool:
        """检查原子是否在已知的分子片段中"""
        # 简化实现：检查周围原子模式
        atom_idx = atom.GetIdx()
        neighbors = [mol.GetAtomWithIdx(neighbor.GetIdx()) 
                    for neighbor in atom.GetNeighbors()]
        
        neighbor_elements = [n.GetSymbol() for n in neighbors]
        
        # 检查常见模式
        if atom.GetSymbol() == 'C':
            # 羧基碳
            if 'O' in neighbor_elements and neighbor_elements.count('O') >= 2:
                return True
            # 酰胺碳
            if 'N' in neighbor_elements and 'O' in neighbor_elements:
                return True
        
        elif atom.GetSymbol() == 'N':
            # 季铵氮
            if len(neighbors) == 4:
                return True
            # 硝基氮
            if neighbor_elements.count('O') >= 2:
                return True
        
        elif atom.GetSymbol() == 'S':
            # 亚砜或砜
            if neighbor_elements.count('O') >= 1:
                return True
        
        return False
    
    def _is_radical_center(self, atom) -> bool:
        """检查是否是自由基中心"""
        try:
            # 检查未成对电子
            return atom.GetNumRadicalElectrons() > 0
        except:
            return False
    
    def _is_coordination_complex(self, atom, mol) -> bool:
        """检查是否是配位化合物中心"""
        # 简化实现：检查是否有异常多的配位原子
        if atom.GetSymbol() in ['Fe', 'Cu', 'Zn', 'Mn', 'Co', 'Ni', 'Cr']:
            neighbor_count = len(list(atom.GetNeighbors()))
            if neighbor_count > 4:  # 可能是八面体或其他配位几何
                return True
        
        return False
    
    def _detect_unusual_patterns(self, mol) -> List[str]:
        """检测不寻常的分子模式"""
        patterns = []
        
        try:
            # 检测高价态原子
            for atom in mol.GetAtoms():
                valence = self._calculate_atom_valence(atom)
                element = atom.GetSymbol()
                rule = self.VALENCE_RULES.get(element)
                
                if rule and valence > max(rule.common_valences):
                    patterns.append(f"高价态{element}({valence})")
            
            # 检测长键
            for bond in mol.GetBonds():
                try:
                    bond_length = self._estimate_bond_length(bond, mol)
                    if bond_length and bond_length > 2.5:  # 异常长的键
                        atom1 = bond.GetBeginAtom().GetSymbol()
                        atom2 = bond.GetEndAtom().GetSymbol()
                        patterns.append(f"长键{atom1}-{atom2}({bond_length:.2f}Å)")
                except:
                    continue
            
            # 检测环张力
            rings = mol.GetRingInfo()
            for ring in rings.AtomRings():
                if len(ring) < 5:  # 小环可能有张力
                    patterns.append(f"{len(ring)}元环(可能有张力)")
        
        except Exception as e:
            logging.warning(f"异常模式检测失败: {e}")
        
        return patterns
    
    def _estimate_bond_length(self, bond, mol) -> Optional[float]:
        """估算键长"""
        try:
            if mol.GetNumConformers() == 0:
                return None
            
            conf = mol.GetConformer()
            begin_pos = conf.GetAtomPosition(bond.GetBeginAtomIdx())
            end_pos = conf.GetAtomPosition(bond.GetEndAtomIdx())
            
            return begin_pos.Distance(end_pos)
        except:
            return None
    
    def suggest_corrections(self, mol, issues: List[ValenceIssue]) -> List[Dict]:
        """
        建议修正方案
        
        Args:
            mol: 分子对象
            issues: 问题列表
            
        Returns:
            修正建议列表
        """
        corrections = []
        
        for issue in issues:
            if issue.severity == "error":
                correction = self._generate_correction_for_error(mol, issue)
                if correction:
                    corrections.append(correction)
        
        return corrections
    
    def _generate_correction_for_error(self, mol, issue: ValenceIssue) -> Optional[Dict]:
        """为价态错误生成修正建议"""
        try:
            atom = mol.GetAtomWithIdx(issue.atom_idx)
            excess_valence = issue.current_valence - max(issue.expected_valences)
            
            correction = {
                "atom_idx": issue.atom_idx,
                "element": issue.element,
                "problem": f"价态超限{excess_valence}",
                "suggestions": []
            }
            
            # 建议1: 调整形式电荷
            if excess_valence <= 2:
                correction["suggestions"].append({
                    "type": "formal_charge",
                    "action": f"设置形式电荷为+{int(excess_valence)}",
                    "feasibility": "high"
                })
            
            # 建议2: 移除部分化学键
            if excess_valence >= 1:
                correction["suggestions"].append({
                    "type": "remove_bonds",
                    "action": f"移除{int(excess_valence)}个化学键",
                    "feasibility": "medium"
                })
            
            # 建议3: 改变键级
            correction["suggestions"].append({
                "type": "bond_order",
                "action": "将部分单键改为双键或三键",
                "feasibility": "low"
            })
            
            return correction
            
        except Exception as e:
            logging.error(f"修正建议生成失败: {e}")
            return None
    
    def apply_automatic_corrections(self, mol, issues: List[ValenceIssue]) -> Tuple[Chem.Mol, List[str]]:
        """
        自动应用修正
        
        Args:
            mol: 分子对象
            issues: 问题列表
            
        Returns:
            (修正后的分子, 应用的修正列表)
        """
        if not RDKIT_AVAILABLE:
            return mol, ["RDKit不可用"]
        
        mol_copy = Chem.RWMol(mol)
        applied_corrections = []
        
        try:
            for issue in issues:
                if issue.severity == "error":
                    success = self._apply_correction_to_atom(mol_copy, issue)
                    if success:
                        applied_corrections.append(
                            f"原子{issue.atom_idx}({issue.element}): {success}"
                        )
            
            # 尝试清理分子
            try:
                Chem.SanitizeMol(mol_copy, sanitizeOps=(
                    Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                    Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                ))
            except:
                pass  # 清理失败不影响返回
            
            return mol_copy, applied_corrections
            
        except Exception as e:
            logging.error(f"自动修正失败: {e}")
            return mol, [f"修正失败: {e}"]
    
    def _apply_correction_to_atom(self, mol, issue: ValenceIssue) -> Optional[str]:
        """对单个原子应用修正"""
        try:
            atom = mol.GetAtomWithIdx(issue.atom_idx)
            excess_valence = issue.current_valence - max(issue.expected_valences)
            
            if excess_valence <= 0:
                return None
            
            # 方法1: 设置形式电荷
            if excess_valence <= 2 and abs(atom.GetFormalCharge()) < 2:
                new_charge = atom.GetFormalCharge() + int(excess_valence)
                atom.SetFormalCharge(new_charge)
                return f"设置形式电荷为{new_charge}"
            
            # 方法2: 移除最长的键
            bonds = list(atom.GetBonds())
            if bonds and excess_valence >= 1:
                # 按键长排序（如果有坐标信息）
                longest_bond = max(bonds, key=lambda b: self._get_bond_priority_for_removal(b, mol))
                if longest_bond:
                    mol.RemoveBond(longest_bond.GetBeginAtomIdx(), longest_bond.GetEndAtomIdx())
                    return f"移除键到原子{longest_bond.GetOtherAtomIdx(issue.atom_idx)}"
            
            return None
            
        except Exception as e:
            logging.warning(f"原子修正失败: {e}")
            return None
    
    def _get_bond_priority_for_removal(self, bond, mol) -> float:
        """获取键移除的优先级（返回值越大越优先移除）"""
        try:
            # 优先移除较长的键
            if mol.GetNumConformers() > 0:
                bond_length = self._estimate_bond_length(bond, mol)
                if bond_length:
                    return bond_length
            
            # 其次优先移除单键而不是双键、三键
            return -bond.GetBondTypeAsDouble()
            
        except:
            return 0.0
    
    def generate_valence_report(self, mol) -> Dict:
        """
        生成详细的价态分析报告
        
        Args:
            mol: 分子对象
            
        Returns:
            分析报告字典
        """
        is_valid, issues, analysis = self.validate_molecule_valences(mol)
        
        report = {
            "is_valid": is_valid,
            "summary": {
                "total_atoms": analysis.get("total_atoms", 0),
                "error_count": analysis.get("errors", 0),
                "warning_count": analysis.get("warnings", 0),
                "unusual_patterns": len(analysis.get("unusual_patterns", []))
            },
            "details": {
                "issues": [
                    {
                        "atom_idx": issue.atom_idx,
                        "element": issue.element,
                        "current_valence": issue.current_valence,
                        "expected_valences": issue.expected_valences,
                        "severity": issue.severity,
                        "suggestion": issue.suggestion
                    }
                    for issue in issues
                ],
                "valence_distribution": dict(analysis.get("valence_distribution", {})),
                "unusual_patterns": analysis.get("unusual_patterns", [])
            }
        }
        
        # 添加修正建议
        if not is_valid:
            corrections = self.suggest_corrections(mol, issues)
            report["correction_suggestions"] = corrections
        
        return report