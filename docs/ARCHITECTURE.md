# ChemAgent 新架构设计

## 核心理念：编排而非重新包装

### 问题
当前我们重新包装了 RDKit、PubChem 等工具，这带来了：
- 🔴 **维护负担**：需要跟踪上游更新
- 🔴 **功能限制**：只暴露了部分功能
- 🔴 **重复工作**：如果 RDKit 发布官方 MCP，我们的工作就浪费了

### 解决方案：使用官方 MCP + 智能编排

## 新架构

```yaml
# ChemAgent 配置
mcp_servers:
  # 使用官方 MCP 服务器
  rdkit:
    source: "official"  # 官方 RDKit MCP
    url: "github.com/rdkit/rdkit-mcp"
    
  pubchem:
    source: "official"  # 官方 PubChem MCP
    url: "github.com/pubchem/pubchem-mcp"
    
  # 我们只提供没有官方 MCP 的工具
  chemagent_extras:
    source: "custom"
    tools:
      - literature_search  # 如果没有官方的
      - workflow_orchestrator
```

## ChemAgent 的价值定位

### 1. **智能工作流编排**

```python
# chemagent/workflows/drug_discovery.yaml
name: Drug Discovery Workflow
steps:
  - name: "Structure Validation"
    tool: "rdkit.validate_smiles"  # 使用官方 RDKit MCP
    
  - name: "Property Calculation"
    tool: "rdkit.calculate_descriptors"
    params:
      descriptors: ["MW", "LogP", "TPSA"]
      
  - name: "ADMET Prediction"
    tool: "admetlab.predict"  # 使用官方 ADMETlab MCP
    
  - name: "Patent Check"
    tool: "surechem.search"  # 使用官方专利 MCP
    
  - name: "Optimization Suggestions"
    orchestrator: true  # 我们的智能编排
    combine_results: true
    generate_report: true
```

### 2. **提示词工程和最佳实践**

```python
# chemagent/prompts/chemistry_guide.md
## 分子分析最佳实践

当用户要求分析分子时：

1. **首先验证结构**
   使用 `rdkit.validate_smiles` 确保输入有效

2. **基础属性**
   调用 `rdkit.calculate_descriptors` 获取:
   - 分子量 (MW)
   - 亲脂性 (LogP)
   - 极性表面积 (TPSA)

3. **药物相似性**
   使用 `rdkit.lipinski_filter` 检查五规则

4. **ADMET 预测**
   如果是药物分子，调用 `admetlab.predict`

5. **结构优化建议**
   基于上述结果，提供改进建议
```

### 3. **MCP 工具发现和管理**

```python
# chemagent/mcp_registry.py
class MCPRegistry:
    """注册和管理可用的 MCP 工具"""
    
    def discover_tools(self):
        """自动发现可用的 MCP 工具"""
        tools = {
            "rdkit": self.check_rdkit_mcp(),
            "pubchem": self.check_pubchem_mcp(),
            "chembl": self.check_chembl_mcp(),
            # ...
        }
        return tools
    
    def install_missing(self):
        """安装缺失的 MCP 服务器"""
        missing = self.find_missing_tools()
        for tool in missing:
            if tool.has_official_mcp:
                self.install_official(tool)
            else:
                self.suggest_alternative(tool)
```

### 4. **只提供真正需要的工具**

```python
# chemagent/custom_tools/
├── workflow_orchestrator.py  # 编排多个 MCP 工具
├── report_generator.py       # 生成综合报告
└── missing_tools/           # 只实现没有官方 MCP 的工具
    ├── literature_aggregator.py  # 如果没有官方的
    └── safety_checker.py         # 如果没有官方的
```

## 实施方案

### 第一步：调研现有 MCP 工具

```markdown
## 已有官方 MCP 的工具
- [ ] RDKit - 检查是否有官方 MCP
- [ ] PubChem - 检查 API 是否支持 MCP
- [ ] ChEMBL - 检查官方支持
- [ ] UniProt - 蛋白质数据库

## 需要我们实现的
- [ ] 文献聚合（跨 PubMed、ChemRxiv）
- [ ] 工作流编排
- [ ] 综合报告生成
```

### 第二步：创建工具适配层

```python
# chemagent/adapters/tool_adapter.py
class ToolAdapter:
    """适配不同来源的 MCP 工具"""
    
    def adapt_rdkit(self, official_mcp):
        """为官方 RDKit MCP 添加化学领域特定的包装"""
        return ChemistryAwareWrapper(official_mcp)
    
    def create_workflow(self, tools):
        """组合多个工具创建工作流"""
        return Workflow(tools)
```

### 第三步：提供集成脚本

```bash
#!/bin/bash
# install_mcp_tools.sh

echo "Installing official MCP tools..."

# RDKit MCP
if [ ! -d "$HOME/.mcp/rdkit" ]; then
    echo "Installing RDKit MCP..."
    git clone https://github.com/rdkit/rdkit-mcp.git
    # ... 安装步骤
fi

# PubChem MCP
if [ ! -d "$HOME/.mcp/pubchem" ]; then
    echo "Installing PubChem MCP..."
    # ... 安装步骤
fi

# ChemAgent 编排器
echo "Installing ChemAgent Orchestrator..."
pip install chemagent-orchestrator
```

## 优势

### 对比当前方案

| 方面 | 当前（重新包装） | 新方案（编排） |
|------|-----------------|---------------|
| 维护成本 | 高（需要跟踪更新） | 低（使用官方） |
| 功能完整性 | 受限（只暴露部分） | 完整（官方全功能） |
| 扩展性 | 差（需要修改代码） | 好（配置即可） |
| 社区贡献 | 重复劳动 | 真正的价值 |

### 真正的价值

1. **工作流编排** - 组合工具完成复杂任务
2. **领域知识** - 化学专业的提示词和指导
3. **最佳实践** - 告诉 AI 如何正确使用工具
4. **缺失工具** - 只实现真正缺失的功能

## 示例：新的使用方式

### Claude Code 中

```markdown
User: 分析这个分子 CC(=O)OC1=CC=CC=C1C(=O)O

Claude (with ChemAgent):
我将使用化学工具分析这个分子（阿司匹林）。

[通过 ChemAgent 编排器调用官方 RDKit MCP]
1. 验证 SMILES ✓
2. 计算属性：
   - MW: 180.16
   - LogP: 1.19
   
[调用官方 ADMETlab MCP]
3. ADMET 预测：
   - 口服吸收: 良好
   - BBB 透过: 中等

[ChemAgent 编排器生成综合报告]
综合分析：这是阿司匹林，具有良好的药物特性...
```

## 行动计划

1. **调研阶段**
   - 列出所有需要的化学工具
   - 检查哪些有官方 MCP
   - 确定真正需要实现的工具

2. **重构阶段**
   - 移除重复包装的代码
   - 实现编排器
   - 创建工作流系统

3. **优化阶段**
   - 完善提示词工程
   - 添加更多工作流模板
   - 创建最佳实践文档

这样 ChemAgent 就成为一个真正有价值的项目：
- ✅ 不重复造轮子
- ✅ 专注于编排和智能
- ✅ 易于维护和扩展
- ✅ 为社区提供真正的价值
