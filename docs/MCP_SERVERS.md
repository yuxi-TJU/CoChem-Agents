# MCP Servers for ChemAgent

## 概述

ChemAgent 支持通过 MCP (Model Context Protocol) 集成各种化学工具。我们优先使用官方的 MCP 实现，为不存在官方实现的工具提供自定义集成。

## 已支持的 MCP 服务器

### 1. RDKit MCP Server (官方)

**包名**: `mcp-rdkit`  
**状态**: ✅ 官方支持  
**PyPI**: https://pypi.org/project/mcp-rdkit/

#### 功能
- 分子解析（SMILES, MOL文件）
- 2D/3D结构生成
- 分子描述符计算（MW, LogP, TPSA等）
- 分子指纹（Morgan, MACCS, RDKit等）
- 相似性计算（Tanimoto, Dice, Cosine）
- 分子可视化（2D图像, 3D结构, SVG输出）

#### 安装
```bash
# 通过 ChemAgent 安装
python chemagent_install.py mcp

# 或单独安装
pip install mcp-rdkit

# 或使用安装脚本
./install_rdkit_mcp.sh
```

#### 使用示例
```python
# MCP 请求示例
{
  "method": "calculate_properties",
  "params": {
    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "properties": ["MW", "LogP", "TPSA", "HBD", "HBA"]
  }
}
```

### 2. PubChem MCP Server

**状态**: ⚠️ 使用内部集成（官方MCP待发布）  
**后备方案**: PubChemPy

#### 功能
- 化合物搜索
- 生物活性数据
- 专利信息查询
- 同义词查找

### 3. ChEMBL MCP Server

**状态**: ⚠️ 使用内部集成（官方MCP待发布）  
**后备方案**: ChEMBL Web Services

#### 功能
- 生物活性数据
- 靶点搜索
- 药物信息
- 临床试验数据

### 4. USPTO Patents MCP

**状态**: ⚠️ 使用内部集成  
**后备方案**: USPTO API

#### 功能
- 专利搜索
- 先有技术查询
- 专利状态检查

### 5. Chemical Suppliers MCP

**状态**: ⚠️ 使用内部集成  
**后备方案**: 供应商API

#### 功能
- 价格查询
- 可用性检查
- 批量定价

## MCP 服务器架构

```
ChemAgent
    │
    ├── MCP Orchestrator
    │   ├── 检查官方MCP可用性
    │   ├── 路由到官方MCP（如果存在）
    │   └── 使用内部工具（作为后备）
    │
    ├── 官方 MCP 服务器
    │   ├── mcp-rdkit ✅
    │   ├── mcp-pubchem (待发布)
    │   └── mcp-chembl (待发布)
    │
    └── 内部工具（后备方案）
        ├── RDKit Python
        ├── PubChemPy
        ├── ChEMBL API
        └── 自定义实现
```

## 配置文件

MCP 服务器配置存储在 `~/.chemagent/mcp_config.json`:

```json
{
  "mcp_servers": {
    "rdkit": {
      "type": "official",
      "name": "RDKit MCP Server",
      "endpoint": "http://localhost:8766/mcp",
      "capabilities": [
        "molecular_analysis",
        "descriptor_calculation",
        "fingerprint_generation"
      ],
      "status": "installed"
    },
    "pubchem": {
      "type": "fallback",
      "name": "PubChem (Internal)",
      "module": "chemagent.tools.pubchem",
      "status": "available"
    }
  }
}
```

## 开发新的 MCP 服务器

如果您想为某个化学工具开发官方 MCP 服务器：

1. **遵循 MCP 协议规范**
2. **实现标准方法**：
   - `initialize`
   - `execute`
   - `list_capabilities`
3. **发布到 PyPI**
4. **提交 PR 到 ChemAgent**

## 故障排除

### RDKit MCP 无法启动
```bash
# 检查安装
python -c "import mcp_rdkit"

# 手动启动服务器
python -m mcp_rdkit

# 查看日志
tail -f ~/.chemagent/logs/rdkit_mcp.log
```

### 端口冲突
```bash
# 修改配置文件中的端口
vim ~/.chemagent/mcp_config.json
# 将 8766 改为其他端口
```

### 权限问题
```bash
# 确保脚本可执行
chmod +x ~/.chemagent/mcp_servers/start_rdkit_mcp.sh
```

## 未来计划

我们正在积极推动更多官方 MCP 服务器的开发：

- [ ] OpenBabel MCP
- [ ] Gaussian/ORCA MCP
- [ ] AutoDock MCP
- [ ] PyMOL MCP
- [ ] ChemDraw MCP

## 贡献

欢迎贡献新的 MCP 服务器集成！请查看 [CONTRIBUTING.md](../CONTRIBUTING.md) 了解详情。
