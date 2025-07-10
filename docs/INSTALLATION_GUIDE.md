# 可扩展非天然氨基酸PDB搜索引擎 - 安装和使用指南

## 📋 **系统要求**

### 基础要求
- **Python**: 3.8+ (推荐 3.9+)
- **操作系统**: Windows 10+, macOS 10.15+, Linux (Ubuntu 18.04+)
- **内存**: 最少 4GB RAM (推荐 8GB+)
- **存储**: 最少 1GB 可用空间
- **网络**: 用于下载依赖包

### 推荐配置
- **CPU**: 4核心+ (支持并行处理)
- **内存**: 16GB+ (大规模数据处理)
- **存储**: SSD硬盘 (提升I/O性能)

## 🔧 **安装步骤**

### 1. 克隆或下载项目
```bash
# 如果使用Git
git clone <repository-url>
cd amino_acids_search_engine

# 或直接下载并解压项目文件
```

### 2. 创建虚拟环境 (推荐)
```bash
# 使用venv
python -m venv amino_acids_env
source amino_acids_env/bin/activate  # Linux/macOS
# 或
amino_acids_env\Scripts\activate     # Windows

# 使用conda
conda create -n amino_acids_env python=3.9
conda activate amino_acids_env
```

### 3. 安装依赖包

#### 基础依赖 (第一阶段)
```bash
pip install numpy pandas sqlite3
```

#### 性能优化依赖 (第二阶段)
```bash
pip install scikit-learn scipy psutil
```

#### 高级功能依赖 (第三阶段)
```bash
# Web API和界面
pip install fastapi uvicorn pydantic

# 科学计算和可视化
pip install matplotlib seaborn plotly

# 可选：RDKit (分子处理)
conda install -c conda-forge rdkit  # 推荐使用conda安装
# 或
pip install rdkit-pypi
```

#### 一键安装所有依赖
```bash
pip install -r requirements.txt
```

### 4. 验证安装
```bash
python -c "
import numpy, pandas, sklearn, fastapi
print('所有依赖包安装成功!')
"
```

## 🚀 **快速开始**

### 1. 基础搜索 (第一阶段)
```bash
# 运行基础演示
python scalable_search_demo.py

# 使用命令行接口
python scalable_search_cli.py search --input test --report
```

### 2. 性能优化搜索 (第二阶段)
```bash
# 运行性能演示
python phase2_demo.py

# 性能对比测试
python performance_comparison.py
```

### 3. 高级功能 (第三阶段)
```bash
# 运行高级功能演示
python phase3_demo.py

# 启动Web服务器
python web_api_server.py
# 访问 http://localhost:8000
```

## 📚 **使用示例**

### 基础搜索
```python
from scalable_search_engine import ScalableSearchEngine

# 初始化搜索引擎
engine = ScalableSearchEngine()

# 搜索氨基酸
results = engine.search({
    'residue_name': '0A1',
    'molecular_formula': 'C10H13NO3'
})

for result in results:
    print(f"找到: {result.amino_acid_id} - {result.amino_acid_record.name}")
```

### 性能优化搜索
```python
from performance_optimized_engine import PerformanceOptimizedSearchEngine, PerformanceConfig

# 配置性能参数
config = PerformanceConfig(
    lsh_num_tables=10,
    memory_cache_size=10000,
    max_workers=8
)

# 初始化优化引擎
engine = PerformanceOptimizedSearchEngine(config=config)

# 执行优化搜索
results = engine.optimized_search({
    'residue_name': '0A1'
}, methods=['residue_name', 'ecfp_similarity'])
```

### 3D结构分析
```python
from advanced_features_engine import Structure3DAnalyzer

# 初始化分析器
analyzer = Structure3DAnalyzer()

# 解析PDB结构
with open('structure.pdb', 'r') as f:
    pdb_content = f.read()

structure = analyzer.parse_pdb_structure(pdb_content, "0A1")

# 提取几何特征
features = analyzer.extract_geometric_features(structure)
print(f"键长数量: {len(features.bond_lengths)}")
print(f"手性中心: {len(features.chiral_centers)}")
```

### 机器学习分类
```python
from advanced_features_engine import MLClassifier

# 初始化分类器
ml_classifier = MLClassifier()

# 训练模型 (需要训练数据)
# X, y = ml_classifier.prepare_training_data(amino_acid_records)
# results = ml_classifier.train_models(X, y)

# 预测 (需要已训练的模型)
# prediction = ml_classifier.predict(amino_acid_record)
```

### Web API使用
```bash
# 启动服务器
python web_api_server.py

# 在另一个终端测试API
curl -X GET "http://localhost:8000/api/health"

curl -X POST "http://localhost:8000/api/search" \
     -H "Content-Type: application/json" \
     -d '{"residue_name": "0A1"}'
```

## 🔧 **配置选项**

### 性能配置
```python
from performance_optimized_engine import PerformanceConfig

config = PerformanceConfig(
    # LSH配置
    lsh_num_tables=10,          # LSH哈希表数量
    lsh_hash_size=16,           # 哈希大小
    
    # 缓存配置
    memory_cache_size=10000,    # 内存缓存大小
    disk_cache_size=100000,     # 磁盘缓存大小
    
    # 并行配置
    max_workers=8,              # 最大工作线程
    batch_size=100,             # 批处理大小
    
    # 内存配置
    max_memory_usage=200*1024*1024,  # 最大内存使用 (200MB)
    compression_enabled=True,    # 启用压缩
)
```

### Web服务配置
```python
# 在web_api_server.py中修改
api_server.run_server(
    host="0.0.0.0",    # 监听地址
    port=8000,         # 端口号
    debug=False        # 调试模式
)
```

## 📁 **项目结构**

```
amino_acids_search_engine/
├── 第一阶段 (核心架构)
│   ├── scalable_search_engine.py
│   ├── scalable_search_cli.py
│   └── scalable_search_demo.py
│
├── 第二阶段 (性能优化)
│   ├── performance_optimized_engine.py
│   ├── performance_comparison.py
│   └── phase2_demo.py
│
├── 第三阶段 (高级功能)
│   ├── advanced_features_engine.py
│   ├── web_api_server.py
│   ├── advanced_analytics.py
│   └── phase3_demo.py
│
├── 文档
│   ├── README_scalable.md
│   ├── README_phase2.md
│   ├── README_phase3.md
│   └── INSTALLATION_GUIDE.md
│
├── 数据
│   ├── amino_acids.db          # SQLite数据库
│   ├── test/                   # 测试PDB文件
│   └── cache/                  # 缓存目录
│
└── 配置
    ├── requirements.txt        # Python依赖
    └── config.json            # 配置文件
```

## 🐛 **常见问题**

### 1. 导入错误
```
ImportError: No module named 'sklearn'
```
**解决方案**: 安装scikit-learn
```bash
pip install scikit-learn
```

### 2. 内存不足
```
MemoryError: Unable to allocate array
```
**解决方案**: 调整内存配置
```python
config = PerformanceConfig(
    max_memory_usage=100*1024*1024,  # 减少到100MB
    compression_enabled=True
)
```

### 3. Web服务启动失败
```
ImportError: No module named 'fastapi'
```
**解决方案**: 安装Web依赖
```bash
pip install fastapi uvicorn
```

### 4. PDB文件解析失败
```
ValueError: 无法解析PDB结构
```
**解决方案**: 检查PDB文件格式
- 确保文件包含ATOM或HETATM记录
- 检查残基名称是否正确
- 验证坐标格式

### 5. 数据库连接错误
```
sqlite3.OperationalError: database is locked
```
**解决方案**: 
- 关闭其他使用数据库的进程
- 删除.db-journal文件
- 重新初始化数据库

## 📊 **性能调优**

### 1. 内存优化
```python
# 启用压缩
config.compression_enabled = True

# 调整缓存大小
config.memory_cache_size = 5000  # 减少内存使用

# 使用内存映射
config.use_memory_mapping = True
```

### 2. 并行优化
```python
# 调整工作线程数
import multiprocessing
config.max_workers = min(8, multiprocessing.cpu_count())

# 调整批处理大小
config.batch_size = 50  # 减少批处理大小
```

### 3. 缓存优化
```python
# 预加载热点查询
config.preload_popular_queries = True

# 调整缓存TTL
config.cache_ttl = 1800  # 30分钟
```

## 🔄 **更新和维护**

### 更新依赖包
```bash
pip install --upgrade -r requirements.txt
```

### 清理缓存
```bash
# 删除缓存目录
rm -rf cache/

# 重建索引
python -c "
from performance_optimized_engine import PerformanceOptimizedSearchEngine
engine = PerformanceOptimizedSearchEngine()
print('索引重建完成')
"
```

### 备份数据
```bash
# 备份数据库
cp amino_acids.db amino_acids_backup.db

# 备份配置
cp config.json config_backup.json
```

## 📞 **技术支持**

如果遇到问题，请：

1. **检查日志**: 查看控制台输出的错误信息
2. **验证环境**: 确认Python版本和依赖包
3. **查看文档**: 参考相应阶段的README文档
4. **运行测试**: 使用demo脚本验证功能
5. **重新安装**: 清理环境后重新安装

## 🎉 **开始使用**

现在您已经完成了安装，可以开始使用可扩展非天然氨基酸PDB搜索引擎了！

建议的学习路径：
1. 运行 `scalable_search_demo.py` 了解基础功能
2. 运行 `phase2_demo.py` 体验性能优化
3. 运行 `phase3_demo.py` 探索高级功能
4. 启动 `web_api_server.py` 使用Web界面

祝您使用愉快！🚀
