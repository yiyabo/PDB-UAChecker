"""
分析器测试
"""

import unittest
import tempfile
import os
from pathlib import Path

from ..analysis import PDBAnalyzer
from ..core.models import ResidueInfo, AtomInfo, AminoAcidInfo
from ..utils.config import Config


class TestPDBAnalyzer(unittest.TestCase):
    """PDB分析器测试类"""
    
    def setUp(self):
        """测试前准备"""
        self.config = Config()
        self.analyzer = PDBAnalyzer(self.config)
    
    def test_analyzer_initialization(self):
        """测试分析器初始化"""
        self.assertIsNotNone(self.analyzer.parser)
        self.assertIsNotNone(self.analyzer.database)
        self.assertIsNotNone(self.analyzer.verification_engine)
    
    def test_create_sample_residue(self):
        """测试创建样本残基"""
        # 创建丙氨酸残基
        atoms = [
            AtomInfo('N', 'N', 0.0, 0.0, 0.0, 'ALA', 1, 'A'),
            AtomInfo('CA', 'C', 1.5, 0.0, 0.0, 'ALA', 1, 'A'),
            AtomInfo('C', 'C', 1.5, 1.5, 0.0, 'ALA', 1, 'A'),
            AtomInfo('O', 'O', 1.5, 2.5, 0.0, 'ALA', 1, 'A'),
            AtomInfo('CB', 'C', 2.5, 0.0, 0.0, 'ALA', 1, 'A'),
        ]
        
        residue = ResidueInfo('ALA', 1, 'A', atoms)
        
        self.assertEqual(residue.residue_name, 'ALA')
        self.assertEqual(len(residue.atoms), 5)
        self.assertIsNotNone(residue.molecular_formula)
        self.assertIsNotNone(residue.atom_composition)
    
    def test_database_stats(self):
        """测试数据库统计"""
        stats = self.analyzer.get_database_stats()
        self.assertIsInstance(stats, dict)
        self.assertIn('total_amino_acids', stats)
    
    def test_system_validation(self):
        """测试系统验证"""
        validation = self.analyzer.validate_system()
        self.assertIsInstance(validation, dict)
        self.assertIn('overall_status', validation)
        self.assertIn('components', validation)
    
    def create_test_pdb_file(self) -> str:
        """创建测试PDB文件"""
        pdb_content = """ATOM      1  N   ALA A   1      20.154  16.967  14.365  1.00 20.00           N  
ATOM      2  CA  ALA A   1      19.030  16.101  14.618  1.00 20.00           C  
ATOM      3  C   ALA A   1      17.664  16.849  14.897  1.00 20.00           C  
ATOM      4  O   ALA A   1      17.764  18.067  15.086  1.00 20.00           O  
ATOM      5  CB  ALA A   1      18.756  15.178  13.425  1.00 20.00           C  
END
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            return f.name
    
    def test_pdb_analysis(self):
        """测试PDB分析（需要有效的数据库）"""
        # 创建测试PDB文件
        pdb_file = self.create_test_pdb_file()
        
        try:
            # 执行分析
            result = self.analyzer.analyze_pdb(pdb_file)
            
            # 验证结果
            self.assertIsNotNone(result)
            self.assertEqual(result.pdb_file, pdb_file)
            self.assertGreaterEqual(result.total_residues, 1)
            self.assertGreaterEqual(result.analysis_time, 0)
        
        finally:
            # 清理测试文件
            os.unlink(pdb_file)


if __name__ == '__main__':
    unittest.main()