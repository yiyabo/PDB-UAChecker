#!/usr/bin/env python3
import sys
import os

# 添加项目路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 逐步测试导入
print("Testing imports...")

try:
    print("1. Testing core models import...")
    from src.core.models import AminoAcidRecord
    print("   ✓ AminoAcidRecord imported successfully")
except Exception as e:
    print(f"   ✗ AminoAcidRecord import failed: {e}")

try:
    print("2. Testing core database import...")
    from src.core.database import AminoAcidDatabase
    print("   ✓ AminoAcidDatabase imported successfully")
except Exception as e:
    print(f"   ✗ AminoAcidDatabase import failed: {e}")

try:
    print("3. Testing search indexing import...")
    from src.search.indexing import IndexManager
    print("   ✓ IndexManager imported successfully")
except Exception as e:
    print(f"   ✗ IndexManager import failed: {e}")

try:
    print("4. Testing search strategies import...")
    from src.search.strategies import ResidueNameMatcher
    print("   ✓ ResidueNameMatcher imported successfully")
except Exception as e:
    print(f"   ✗ ResidueNameMatcher import failed: {e}")

try:
    print("5. Testing search engine import...")
    from src.search.engine import ScalableSearchEngine
    print("   ✓ ScalableSearchEngine imported successfully")
except Exception as e:
    print(f"   ✗ ScalableSearchEngine import failed: {e}")

print("\nImport testing complete!")
