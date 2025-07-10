#!/usr/bin/env python3
"""
PDB-UAChecker 安装配置
"""

import os
from setuptools import setup, find_packages

# 读取README文件
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

# 读取requirements文件
def read_requirements():
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith('#')]

setup(
    name="pdb-uachecker",
    version="1.0.0",
    author="PDB-UAChecker Team",
    author_email="your.email@example.com",
    description="Advanced Non-Natural Amino Acid PDB Search Engine",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/yiyabo/PDB-UAChecker",
    project_urls={
        "Bug Tracker": "https://github.com/yiyabo/PDB-UAChecker/issues",
        "Documentation": "https://github.com/yiyabo/PDB-UAChecker/blob/main/README.md",
        "Source Code": "https://github.com/yiyabo/PDB-UAChecker",
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=2.12.0", 
            "black>=21.0.0",
            "flake8>=3.9.0",
            "mypy>=0.900",
        ],
        "docs": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=0.5.0",
        ],
        "web": [
            "fastapi>=0.68.0",
            "uvicorn>=0.15.0",
        ]
    },
    entry_points={
        "console_scripts": [
            "pdb-uachecker=src.api.cli:main",
            "pdb-uachecker-web=src.api.web_server:main",
        ],
    },
    include_package_data=True,
    package_data={
        "src": ["data/*", "docs/*"],
    },
    keywords=[
        "bioinformatics",
        "chemistry", 
        "pdb",
        "amino-acids",
        "isomers",
        "molecular-fingerprints",
        "search-engine",
        "machine-learning"
    ],
    zip_safe=False,
)