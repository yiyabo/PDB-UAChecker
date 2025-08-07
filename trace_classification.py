#!/usr/bin/env python3
"""
追踪具体的分类决策过程
"""

import sys
from pathlib import Path

# 添加项目路径
sys.path.insert(0, str(Path(__file__).parent))

from pdb_uachecker.analysis.intelligent_classifier import IntelligentAminoAcidClassifier
from pdb_uachecker.core.database.enhanced_manager import EnhancedDatabaseManager
from pdb_uachecker.utils.config import Config

def trace_classification_decision():
    """追踪分类决策过程"""
    print("🔍 追踪分类决策过程")
    print("=" * 60)
    
    # 初始化
    config = Config()
    db = EnhancedDatabaseManager(config)
    classifier = IntelligentAminoAcidClassifier()
    
    # 获取问题案例
    problem_cases = ['NLE', 'ABA', 'CSA']
    
    for aa_id in problem_cases:
        amino_acids = [aa for aa in db.get_all_amino_acids() if aa.id == aa_id]
        if amino_acids:
            aa = amino_acids[0]
            print(f"\n🧪 追踪 {aa.id} 的分类过程:")
            print(f"   SMILES: {aa.smiles}")
            
            # 手动执行分类逻辑
            result = {'standard_category': 'requires_expert_review', 'confidence': 0.0}
            
            # Tier 1 检查
            print("   Tier 1 检查:")
            if aa.id in classifier.standard_amino_acids:
                print(f"     ✓ 标准氨基酸: {aa.id}")
                result = {'standard_category': 'Standard', 'confidence': 1.0}
            elif classifier._is_d_amino_by_naming(aa):
                print("     ✓ D型氨基酸 (命名)")
                result = {'standard_category': 'D-amino acid', 'confidence': 1.0}
            elif classifier._is_d_amino_by_chirality(aa):
                print("     ✓ D型氨基酸 (手性)")
                result = {'standard_category': 'D-amino acid', 'confidence': 1.0}
            else:
                print("     ✗ 无Tier 1匹配")
            
            # 如果没有Tier 1匹配，检查Tier 2
            if result['confidence'] < 1.0:
                print("   Tier 2 检查:")
                aromatic_result = classifier._is_aromatic_amino_acid(aa)
                print(f"     芳香性检测: {aromatic_result}")
                
                if aromatic_result['is_match']:
                    print("     ✓ 芳香性氨基酸匹配")
                    result = {'standard_category': 'Aromatic', 'confidence': 0.95}
                else:
                    print("     ✗ 无芳香性匹配")
            
            print(f"   最终分类: {result}")
            
            # 实际运行分类器对比
            actual_result = classifier.classify_amino_acid(aa)
            print(f"   分类器结果: {actual_result['standard_category']}, 置信度: {actual_result['confidence']}")

if __name__ == "__main__":
    trace_classification_decision()