#!/usr/bin/env python3
"""
可扩展非天然氨基酸PDB搜索引擎 - Web API服务
基于FastAPI的RESTful API和Web界面
"""

import os
import json
import time
import asyncio
from typing import Dict, List, Optional, Any, Union
from datetime import datetime
import uuid

# 科学计算库
try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False
    print("警告: NumPy未安装，某些统计功能将不可用")

# 科学计算库
try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False
    print("警告: NumPy未安装，某些统计功能将不可用")

# Web框架
try:
    from fastapi import FastAPI, HTTPException, UploadFile, File, Form, BackgroundTasks
    from fastapi.responses import HTMLResponse, JSONResponse, FileResponse
    from fastapi.staticfiles import StaticFiles
    from fastapi.middleware.cors import CORSMiddleware
    from pydantic import BaseModel, Field
    import uvicorn
    FASTAPI_AVAILABLE = True
except ImportError:
    FASTAPI_AVAILABLE = False
    print("警告: FastAPI未安装，Web API功能将不可用")

# 导入搜索引擎
from performance_optimized_engine import (
    PerformanceOptimizedSearchEngine,
    PerformanceConfig,
    OptimizedCompatibilityWrapper
)
from advanced_features_engine import (
    Structure3DAnalyzer,
    MLClassifier,
    StereochemistryAnalyzer,
    Structure3D
)

# 数据模型
class SearchQuery(BaseModel):
    """搜索查询模型"""
    residue_name: Optional[str] = None
    molecular_formula: Optional[str] = None
    atom_composition: Optional[Dict[str, int]] = None
    molecular_weight: Optional[float] = None
    weight_tolerance: Optional[float] = 5.0
    smiles: Optional[str] = None
    methods: List[str] = ['residue_name', 'molecular_formula', 'ecfp_similarity']
    max_results: int = 10

class BatchSearchQuery(BaseModel):
    """批量搜索查询模型"""
    queries: List[SearchQuery]
    parallel: bool = True

class MLPredictionQuery(BaseModel):
    """机器学习预测查询模型"""
    amino_acid_id: str
    molecular_formula: Optional[str] = None
    atom_composition: Optional[Dict[str, int]] = None
    molecular_weight: Optional[float] = None
    smiles: Optional[str] = None
    key_features: List[str] = []
    model_name: str = 'random_forest'

class StructureAnalysisQuery(BaseModel):
    """结构分析查询模型"""
    pdb_content: str
    residue_name: str
    analysis_type: str = 'full'  # 'geometry', 'chirality', 'stereoisomers', 'full'

class TaskStatus(BaseModel):
    """任务状态模型"""
    task_id: str
    status: str  # 'pending', 'running', 'completed', 'failed'
    progress: float = 0.0
    result: Optional[Any] = None
    error: Optional[str] = None
    created_at: datetime
    updated_at: datetime

class AdvancedSearchAPI:
    """高级搜索API服务"""
    
    def __init__(self):
        # 初始化搜索引擎
        config = PerformanceConfig(
            lsh_num_tables=10,
            memory_cache_size=10000,
            max_workers=8,
            compression_enabled=True
        )
        
        self.search_engine = PerformanceOptimizedSearchEngine(config=config)
        self.structure_analyzer = Structure3DAnalyzer()
        self.ml_classifier = MLClassifier()
        self.stereochemistry_analyzer = StereochemistryAnalyzer()
        
        # 任务管理
        self.tasks = {}
        
        # 初始化FastAPI应用
        if FASTAPI_AVAILABLE:
            self.app = FastAPI(
                title="非天然氨基酸PDB搜索引擎API",
                description="可扩展的非天然氨基酸PDB搜索引擎，支持3D结构分析、机器学习分类和立体化学分析",
                version="3.0.0"
            )
            
            # 配置CORS
            self.app.add_middleware(
                CORSMiddleware,
                allow_origins=["*"],
                allow_credentials=True,
                allow_methods=["*"],
                allow_headers=["*"],
            )
            
            # 注册路由
            self._register_routes()
        else:
            self.app = None
    
    def _register_routes(self):
        """注册API路由"""
        
        @self.app.get("/", response_class=HTMLResponse)
        async def root():
            """主页"""
            return self._get_web_interface()
        
        @self.app.get("/api/health")
        async def health_check():
            """健康检查"""
            return {
                "status": "healthy",
                "timestamp": datetime.now().isoformat(),
                "version": "3.0.0",
                "components": {
                    "search_engine": "active",
                    "ml_classifier": "active" if self.ml_classifier.is_trained else "inactive",
                    "structure_analyzer": "active",
                    "stereochemistry_analyzer": "active"
                }
            }
        
        @self.app.get("/api/stats")
        async def get_statistics():
            """获取系统统计信息"""
            try:
                stats = self.search_engine.get_comprehensive_statistics()
                return {
                    "success": True,
                    "data": stats,
                    "timestamp": datetime.now().isoformat()
                }
            except Exception as e:
                raise HTTPException(status_code=500, detail=str(e))
        
        @self.app.post("/api/search")
        async def search_amino_acids(query: SearchQuery):
            """搜索氨基酸"""
            try:
                # 构建查询数据
                query_data = {}
                if query.residue_name:
                    query_data['residue_name'] = query.residue_name
                if query.molecular_formula:
                    query_data['molecular_formula'] = query.molecular_formula
                if query.atom_composition:
                    query_data['atom_composition'] = query.atom_composition
                if query.molecular_weight:
                    query_data['molecular_weight'] = query.molecular_weight
                    query_data['weight_tolerance'] = query.weight_tolerance
                if query.smiles:
                    query_data['smiles'] = query.smiles
                
                # 执行搜索
                start_time = time.time()
                results = self.search_engine.optimized_search(
                    query_data, methods=query.methods, max_results=query.max_results
                )
                search_time = time.time() - start_time
                
                # 转换结果格式
                formatted_results = []
                for result in results:
                    formatted_result = {
                        "amino_acid_id": result.amino_acid_id,
                        "match_method": result.match_method,
                        "confidence_score": result.confidence_score,
                        "amino_acid_record": result.amino_acid_record.to_dict(),
                        "additional_info": result.additional_info or {}
                    }
                    formatted_results.append(formatted_result)
                
                return {
                    "success": True,
                    "data": {
                        "results": formatted_results,
                        "search_time_ms": search_time * 1000,
                        "total_results": len(formatted_results),
                        "query": query_data,
                        "methods_used": query.methods
                    },
                    "timestamp": datetime.now().isoformat()
                }
                
            except Exception as e:
                raise HTTPException(status_code=500, detail=str(e))
        
        @self.app.post("/api/search/batch")
        async def batch_search(query: BatchSearchQuery, background_tasks: BackgroundTasks):
            """批量搜索"""
            try:
                task_id = str(uuid.uuid4())
                
                # 创建任务
                task = TaskStatus(
                    task_id=task_id,
                    status="pending",
                    created_at=datetime.now(),
                    updated_at=datetime.now()
                )
                self.tasks[task_id] = task
                
                # 添加后台任务
                background_tasks.add_task(
                    self._process_batch_search, task_id, query.queries, query.parallel
                )
                
                return {
                    "success": True,
                    "data": {
                        "task_id": task_id,
                        "status": "pending",
                        "total_queries": len(query.queries)
                    },
                    "timestamp": datetime.now().isoformat()
                }
                
            except Exception as e:
                raise HTTPException(status_code=500, detail=str(e))
        
        @self.app.get("/api/task/{task_id}")
        async def get_task_status(task_id: str):
            """获取任务状态"""
            if task_id not in self.tasks:
                raise HTTPException(status_code=404, detail="任务不存在")
            
            task = self.tasks[task_id]
            return {
                "success": True,
                "data": {
                    "task_id": task.task_id,
                    "status": task.status,
                    "progress": task.progress,
                    "result": task.result,
                    "error": task.error,
                    "created_at": task.created_at.isoformat(),
                    "updated_at": task.updated_at.isoformat()
                },
                "timestamp": datetime.now().isoformat()
            }
        
        @self.app.post("/api/upload/pdb")
        async def upload_pdb_file(file: UploadFile = File(...), 
                                 background_tasks: BackgroundTasks = None):
            """上传PDB文件进行分析"""
            try:
                if not file.filename.endswith('.pdb'):
                    raise HTTPException(status_code=400, detail="只支持PDB文件")
                
                # 读取文件内容
                content = await file.read()
                pdb_content = content.decode('utf-8')
                
                task_id = str(uuid.uuid4())
                
                # 创建任务
                task = TaskStatus(
                    task_id=task_id,
                    status="pending",
                    created_at=datetime.now(),
                    updated_at=datetime.now()
                )
                self.tasks[task_id] = task
                
                # 添加后台任务
                if background_tasks:
                    background_tasks.add_task(
                        self._process_pdb_file, task_id, pdb_content, file.filename
                    )
                
                return {
                    "success": True,
                    "data": {
                        "task_id": task_id,
                        "filename": file.filename,
                        "status": "pending"
                    },
                    "timestamp": datetime.now().isoformat()
                }
                
            except Exception as e:
                raise HTTPException(status_code=500, detail=str(e))
        
        @self.app.post("/api/ml/predict")
        async def ml_predict(query: MLPredictionQuery):
            """机器学习预测"""
            try:
                if not self.ml_classifier.is_trained:
                    raise HTTPException(status_code=400, detail="机器学习模型未训练")
                
                # 创建临时氨基酸记录
                from scalable_search_engine import AminoAcidRecord
                
                record = AminoAcidRecord(
                    id=query.amino_acid_id,
                    name="Query Record",
                    molecular_formula=query.molecular_formula or "",
                    molecular_weight=query.molecular_weight or 0.0,
                    smiles=query.smiles or "",
                    atom_composition=query.atom_composition or {},
                    key_features=query.key_features
                )
                
                # 执行预测
                prediction = self.ml_classifier.predict(record, model_name=query.model_name)
                
                return {
                    "success": True,
                    "data": prediction,
                    "timestamp": datetime.now().isoformat()
                }
                
            except Exception as e:
                raise HTTPException(status_code=500, detail=str(e))
        
        @self.app.post("/api/structure/analyze")
        async def analyze_structure(query: StructureAnalysisQuery):
            """结构分析"""
            try:
                # 解析PDB结构
                structure = self.structure_analyzer.parse_pdb_structure(
                    query.pdb_content, query.residue_name
                )
                
                if not structure:
                    raise HTTPException(status_code=400, detail="无法解析PDB结构")
                
                analysis_result = {}
                
                if query.analysis_type in ['geometry', 'full']:
                    # 几何分析
                    geometric_features = self.structure_analyzer.extract_geometric_features(structure)
                    validation = self.structure_analyzer.validate_geometry(structure)
                    
                    analysis_result['geometry'] = {
                        'bond_lengths': geometric_features.bond_lengths,
                        'bond_angles': geometric_features.bond_angles,
                        'dihedral_angles': geometric_features.dihedral_angles,
                        'ring_systems': geometric_features.ring_systems,
                        'validation': validation
                    }
                
                if query.analysis_type in ['chirality', 'full']:
                    # 手性分析
                    chirality = self.stereochemistry_analyzer.analyze_chirality(structure)
                    analysis_result['chirality'] = chirality
                
                if query.analysis_type in ['stereoisomers', 'full']:
                    # E/Z异构体分析
                    ez_isomers = self.stereochemistry_analyzer.detect_ez_isomers(structure)
                    analysis_result['ez_isomers'] = ez_isomers
                
                # 基本结构信息
                analysis_result['structure_info'] = {
                    'amino_acid_id': structure.amino_acid_id,
                    'num_atoms': len(structure.atoms),
                    'num_bonds': len(structure.bonds),
                    'molecular_volume': structure.molecular_volume,
                    'surface_area': structure.surface_area,
                    'center_of_mass': structure.center_of_mass
                }
                
                return {
                    "success": True,
                    "data": analysis_result,
                    "timestamp": datetime.now().isoformat()
                }
                
            except Exception as e:
                raise HTTPException(status_code=500, detail=str(e))

    async def _process_batch_search(self, task_id: str, queries: List[SearchQuery], parallel: bool):
        """处理批量搜索任务"""
        try:
            task = self.tasks[task_id]
            task.status = "running"
            task.updated_at = datetime.now()

            results = []
            total_queries = len(queries)

            for i, query in enumerate(queries):
                # 构建查询数据
                query_data = {}
                if query.residue_name:
                    query_data['residue_name'] = query.residue_name
                if query.molecular_formula:
                    query_data['molecular_formula'] = query.molecular_formula
                if query.atom_composition:
                    query_data['atom_composition'] = query.atom_composition
                if query.molecular_weight:
                    query_data['molecular_weight'] = query.molecular_weight
                    query_data['weight_tolerance'] = query.weight_tolerance
                if query.smiles:
                    query_data['smiles'] = query.smiles

                # 执行搜索
                search_results = self.search_engine.optimized_search(
                    query_data, methods=query.methods, max_results=query.max_results
                )

                # 转换结果格式
                formatted_results = []
                for result in search_results:
                    formatted_result = {
                        "amino_acid_id": result.amino_acid_id,
                        "match_method": result.match_method,
                        "confidence_score": result.confidence_score,
                        "amino_acid_record": result.amino_acid_record.to_dict(),
                        "additional_info": result.additional_info or {}
                    }
                    formatted_results.append(formatted_result)

                results.append({
                    "query_index": i,
                    "query": query_data,
                    "results": formatted_results
                })

                # 更新进度
                task.progress = (i + 1) / total_queries
                task.updated_at = datetime.now()

            # 完成任务
            task.status = "completed"
            task.result = {
                "batch_results": results,
                "total_queries": total_queries,
                "completed_queries": len(results)
            }
            task.progress = 1.0
            task.updated_at = datetime.now()

        except Exception as e:
            task.status = "failed"
            task.error = str(e)
            task.updated_at = datetime.now()

    async def _process_pdb_file(self, task_id: str, pdb_content: str, filename: str):
        """处理PDB文件分析任务"""
        try:
            task = self.tasks[task_id]
            task.status = "running"
            task.updated_at = datetime.now()

            # 使用兼容包装器处理PDB文件
            wrapper = OptimizedCompatibilityWrapper()

            # 临时保存文件
            temp_filename = f"temp_{task_id}.pdb"
            with open(temp_filename, 'w') as f:
                f.write(pdb_content)

            try:
                # 搜索PDB文件
                search_results = wrapper.search_pdb_files([temp_filename])

                # 结构分析
                structure_analyses = []

                # 解析所有残基
                lines = pdb_content.strip().split('\n')
                residue_names = set()

                for line in lines:
                    if line.startswith(('ATOM', 'HETATM')):
                        residue_name = line[17:20].strip()
                        residue_names.add(residue_name)

                # 分析每个残基
                for residue_name in residue_names:
                    structure = self.structure_analyzer.parse_pdb_structure(
                        pdb_content, residue_name
                    )

                    if structure:
                        # 几何分析
                        geometric_features = self.structure_analyzer.extract_geometric_features(structure)
                        validation = self.structure_analyzer.validate_geometry(structure)

                        # 手性分析
                        chirality = self.stereochemistry_analyzer.analyze_chirality(structure)

                        structure_analysis = {
                            'residue_name': residue_name,
                            'structure_info': {
                                'num_atoms': len(structure.atoms),
                                'num_bonds': len(structure.bonds),
                                'molecular_volume': structure.molecular_volume,
                                'surface_area': structure.surface_area,
                                'center_of_mass': structure.center_of_mass
                            },
                            'geometry': {
                                'bond_lengths_stats': {
                                    'mean': float(np.mean(geometric_features.bond_lengths)) if NUMPY_AVAILABLE and geometric_features.bond_lengths else 0,
                                    'std': float(np.std(geometric_features.bond_lengths)) if NUMPY_AVAILABLE and geometric_features.bond_lengths else 0,
                                    'count': len(geometric_features.bond_lengths) if geometric_features.bond_lengths else 0
                                },
                                'bond_angles_stats': {
                                    'mean': float(np.mean(geometric_features.bond_angles)) if NUMPY_AVAILABLE and geometric_features.bond_angles else 0,
                                    'std': float(np.std(geometric_features.bond_angles)) if NUMPY_AVAILABLE and geometric_features.bond_angles else 0,
                                    'count': len(geometric_features.bond_angles) if geometric_features.bond_angles else 0
                                },
                                'ring_systems': geometric_features.ring_systems,
                                'validation': validation
                            },
                            'chirality': chirality
                        }
                        structure_analyses.append(structure_analysis)

                # 完成任务
                task.status = "completed"
                task.result = {
                    "filename": filename,
                    "search_results": search_results,
                    "structure_analyses": structure_analyses,
                    "total_residues": len(structure_analyses)
                }
                task.progress = 1.0
                task.updated_at = datetime.now()

            finally:
                # 清理临时文件
                if os.path.exists(temp_filename):
                    os.remove(temp_filename)

        except Exception as e:
            task.status = "failed"
            task.error = str(e)
            task.updated_at = datetime.now()

    def _get_web_interface(self) -> str:
        """获取Web界面HTML"""
        return """
        <!DOCTYPE html>
        <html lang="zh-CN">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>非天然氨基酸PDB搜索引擎</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }
                .container { max-width: 1000px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; }
                h1 { color: #2c3e50; text-align: center; }
                .section { margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }
                .form-group { margin: 10px 0; }
                label { display: block; margin-bottom: 5px; font-weight: bold; }
                input, select, textarea { width: 100%; padding: 8px; border: 1px solid #ddd; border-radius: 4px; box-sizing: border-box; }
                button { background: #3498db; color: white; padding: 10px 20px; border: none; border-radius: 4px; cursor: pointer; }
                button:hover { background: #2980b9; }
                .results { margin-top: 15px; padding: 10px; background: #f8f9fa; border-radius: 5px; }
                .api-docs { background: #e9ecef; padding: 15px; border-radius: 5px; margin-top: 20px; }
            </style>
        </head>
        <body>
            <div class="container">
                <h1>🧬 非天然氨基酸PDB搜索引擎</h1>

                <div class="section">
                    <h2>🔍 氨基酸搜索</h2>
                    <div class="form-group">
                        <label>残基名称:</label>
                        <input type="text" id="residue_name" placeholder="例如: 0A1">
                    </div>
                    <div class="form-group">
                        <label>分子式:</label>
                        <input type="text" id="molecular_formula" placeholder="例如: C10H13NO3">
                    </div>
                    <div class="form-group">
                        <label>分子量:</label>
                        <input type="number" id="molecular_weight" placeholder="例如: 195.22">
                    </div>
                    <button onclick="searchAminoAcids()">搜索</button>
                    <div id="search_results" class="results" style="display:none;"></div>
                </div>

                <div class="section">
                    <h2>📁 PDB文件上传</h2>
                    <div class="form-group">
                        <label>选择PDB文件:</label>
                        <input type="file" id="pdb_file" accept=".pdb">
                    </div>
                    <button onclick="uploadPDBFile()">上传并分析</button>
                    <div id="upload_results" class="results" style="display:none;"></div>
                </div>

                <div class="api-docs">
                    <h2>📚 API文档</h2>
                    <p>可用的API端点：</p>
                    <div style="font-family: monospace;">
                        <div>GET /api/health - 健康检查</div>
                        <div>GET /api/stats - 系统统计</div>
                        <div>POST /api/search - 搜索氨基酸</div>
                        <div>POST /api/upload/pdb - 上传PDB文件</div>
                    </div>
                </div>
            </div>

            <script>
                async function searchAminoAcids() {
                    const query = {
                        residue_name: document.getElementById('residue_name').value || null,
                        molecular_formula: document.getElementById('molecular_formula').value || null,
                        molecular_weight: parseFloat(document.getElementById('molecular_weight').value) || null,
                        methods: ['residue_name', 'molecular_formula'],
                        max_results: 10
                    };

                    try {
                        const response = await fetch('/api/search', {
                            method: 'POST',
                            headers: {'Content-Type': 'application/json'},
                            body: JSON.stringify(query)
                        });

                        const result = await response.json();
                        displaySearchResults(result);
                    } catch (error) {
                        alert('搜索失败: ' + error.message);
                    }
                }

                function displaySearchResults(result) {
                    const resultsDiv = document.getElementById('search_results');

                    if (result.success && result.data.results.length > 0) {
                        let html = `<h3>找到 ${result.data.total_results} 个结果</h3>`;

                        result.data.results.forEach((item) => {
                            html += `
                                <div style="margin: 10px 0; padding: 10px; background: white; border-radius: 3px;">
                                    <strong>${item.amino_acid_id}</strong> - ${item.amino_acid_record.name}<br>
                                    分子式: ${item.amino_acid_record.molecular_formula}<br>
                                    分子量: ${item.amino_acid_record.molecular_weight}<br>
                                    置信度: ${item.confidence_score.toFixed(3)}
                                </div>
                            `;
                        });

                        resultsDiv.innerHTML = html;
                    } else {
                        resultsDiv.innerHTML = '<h3>未找到匹配结果</h3>';
                    }

                    resultsDiv.style.display = 'block';
                }

                async function uploadPDBFile() {
                    const fileInput = document.getElementById('pdb_file');
                    const file = fileInput.files[0];

                    if (!file) {
                        alert('请选择PDB文件');
                        return;
                    }

                    const formData = new FormData();
                    formData.append('file', file);

                    try {
                        const response = await fetch('/api/upload/pdb', {
                            method: 'POST',
                            body: formData
                        });

                        const result = await response.json();

                        if (result.success) {
                            document.getElementById('upload_results').innerHTML =
                                `<h3>文件上传成功</h3><p>任务ID: ${result.data.task_id}</p><p>状态: 处理中...</p>`;
                            document.getElementById('upload_results').style.display = 'block';

                            // 开始轮询任务状态
                            pollTaskStatus(result.data.task_id, 'upload_results');
                        }
                    } catch (error) {
                        alert('上传失败: ' + error.message);
                    }
                }

                async function pollTaskStatus(taskId, resultDivId) {
                    const maxAttempts = 30; // 最多轮询30次
                    let attempts = 0;

                    const poll = async () => {
                        try {
                            const response = await fetch(`/api/task/${taskId}`);
                            const result = await response.json();

                            if (result.success) {
                                const task = result.data;

                                if (task.status === 'completed') {
                                    displayTaskResult(task, resultDivId);
                                } else if (task.status === 'failed') {
                                    document.getElementById(resultDivId).innerHTML =
                                        `<h3>任务失败</h3><p>错误: ${task.error}</p>`;
                                } else if (attempts < maxAttempts) {
                                    // 更新进度
                                    document.getElementById(resultDivId).innerHTML =
                                        `<h3>处理中...</h3><p>进度: ${(task.progress * 100).toFixed(1)}%</p><p>状态: ${task.status}</p>`;

                                    // 2秒后继续轮询
                                    setTimeout(poll, 2000);
                                } else {
                                    document.getElementById(resultDivId).innerHTML =
                                        `<h3>处理超时</h3><p>请手动查询任务状态: ${taskId}</p>`;
                                }
                            }

                            attempts++;
                        } catch (error) {
                            console.error('轮询任务状态失败:', error);
                            document.getElementById(resultDivId).innerHTML =
                                `<h3>查询失败</h3><p>请手动查询任务: ${taskId}</p>`;
                        }
                    };

                    // 开始轮询
                    poll();
                }

                function displayTaskResult(task, resultDivId) {
                    let html = '<h3>分析完成！</h3>';

                    if (task.result.search_results && task.result.search_results.length > 0) {
                        html += `<h4>🔍 搜索结果 (${task.result.search_results.length} 个匹配)</h4>`;
                        task.result.search_results.forEach(result => {
                            html += `
                                <div style="margin: 10px 0; padding: 10px; background: #e8f4fd; border-radius: 5px;">
                                    <strong>PDB: ${result.pdb_id}</strong><br>
                                    氨基酸: ${result.amino_acid_id}<br>
                                    置信度: ${result.confidence_score.toFixed(3)}<br>
                                    匹配方法: ${result.match_method}
                                </div>
                            `;
                        });
                    }

                    if (task.result.structure_analyses && task.result.structure_analyses.length > 0) {
                        html += `<h4>🧬 结构分析 (${task.result.structure_analyses.length} 个残基)</h4>`;
                        task.result.structure_analyses.forEach(analysis => {
                            html += `
                                <div style="margin: 10px 0; padding: 10px; background: #f0f8f0; border-radius: 5px;">
                                    <strong>残基: ${analysis.residue_name}</strong><br>
                                    原子数: ${analysis.structure_info.num_atoms}<br>
                                    键数: ${analysis.structure_info.num_bonds}<br>
                                    分子体积: ${analysis.structure_info.molecular_volume.toFixed(2)} Ų<br>
                                    手性中心: ${analysis.chirality.num_chiral_centers}<br>
                                    几何验证: ${analysis.geometry.validation.overall_valid ? '✅ 通过' : '❌ 未通过'}
                                </div>
                            `;
                        });
                    }

                    // 添加手动查询链接
                    html += `
                        <div style="margin-top: 15px; padding: 10px; background: #fff3cd; border-radius: 5px;">
                            <strong>📋 详细结果查询:</strong><br>
                            <a href="/api/task/${task.task_id}" target="_blank">查看完整JSON结果</a>
                        </div>
                    `;

                    document.getElementById(resultDivId).innerHTML = html;
                }
            </script>
        </body>
        </html>
        """

    def run_server(self, host: str = "127.0.0.1", port: int = 8000, debug: bool = False):
        """运行Web服务器"""
        if not FASTAPI_AVAILABLE:
            print("错误: FastAPI未安装，无法启动Web服务器")
            print("请安装: pip install fastapi uvicorn")
            return

        print(f"启动Web服务器...")

        # 显示正确的访问地址
        display_host = "localhost" if host in ["127.0.0.1", "0.0.0.0"] else host
        print(f"访问地址: http://{display_host}:{port}")
        print(f"API文档: http://{display_host}:{port}/docs")

        # 启动服务器
        try:
            uvicorn.run(self.app, host=host, port=port, reload=debug)
        except Exception as e:
            print(f"启动服务器失败: {e}")
            print("尝试使用以下命令直接启动:")
            print(f"uvicorn web_api_server:api_server.app --host {host} --port {port}")
            if debug:
                print(" --reload")

def main():
    """主程序"""
    import argparse

    parser = argparse.ArgumentParser(description="可扩展非天然氨基酸PDB搜索引擎 - Web API服务")
    parser.add_argument("--host", default="127.0.0.1", help="服务器主机地址 (默认: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=8000, help="服务器端口 (默认: 8000)")
    parser.add_argument("--debug", action="store_true", help="启用调试模式")

    args = parser.parse_args()

    print("可扩展非天然氨基酸PDB搜索引擎 - Web API服务")

    api_server = AdvancedSearchAPI()

    if api_server.app:
        api_server.run_server(host=args.host, port=args.port, debug=args.debug)
    else:
        print("Web API服务不可用，请安装FastAPI: pip install fastapi uvicorn")

if __name__ == "__main__":
    main()
