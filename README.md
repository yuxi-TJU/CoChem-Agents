# CoChem Agents

A collaborative framework for building chemistry & materials agents with general agent framework, like Gemini CLI, Claude code, or Codex, +**MCP**.

**CoChem Agents:** an open framework for chemistry & materials AI. Use CodeX as the general agent shell and add capabilities via the Model Context Protocol (MCP). Skip one-off agents—publish MCP servers (RDKit, Materials Project, sims, lab APIs) and compose them. Integrate once, reuse everywhere.

**CoChem Agents** turns the “one-agent-per-domain” pattern on its head. Instead of crafting bespoke chemistry or materials agents, we use **general agent framework, like Codex** as the general-purpose agent shell and plug in domain tools via the **Model Context Protocol (MCP)**. Anyone can contribute an MCP server—RDKit, Materials Project, your internal pipeline—and it becomes instantly usable by the same agent. This creates an open, extensible ecosystem rather than a zoo of siloed agents. 


## Why chemistry & materials agents matter

AI is rapidly accelerating discovery across chemistry and materials—from structure/property prediction to polymer and crystal modeling—pushing research beyond static prediction toward **agentic** workflows that plan, act, and iterate. Surveys and community reports document both the momentum and the need for robust tooling to make these systems practical in the lab and in silico. 

## What hasn’t worked

Most prior efforts ship a **standalone agent per subfield** (drug design, catalysis, crystals…), each with custom glue code, brittle integrations, and duplicated effort. Evaluations often emphasize reasoning but struggle with **reproducibility and tool-generalization**, so systems don’t travel well between tasks or labs. Meanwhile, tool access (APIs, DBs, codes) is fragmented and hard to standardize across agents. 


## The core technical obstacles

 - **Heterogeneous tools & schemas:** cheminformatics libs, materials databases, simulation engines—all different call patterns and data models. 

 - **Agent–tool wiring & maintenance:** each agent re-implements connectors and auth, leading to drift and duplication. 

 - **Security & governance:** opening tools to agents raises questions around auth, data access, and isolation. 

 - **Evaluation & provenance:** agent benchmarks underweight reproducibility and end-to-end paper-to-protocol faithfulness. 

## Our approach (what’s different)

**1.One agent framework to rule them all**
Use Codex as the generic agent runtime (chat + tools + prompts). No more domain-specific shells. 

**2.Tools as MCP servers**
Expose chemistry/materials capabilities as MCP tools (standardized names, schemas, metadata). Any MCP-compatible client (like Codex) can discover and call them—zero bespoke glue in the agent. 


**3.Open ecosystem, not one-off agents**

 - **Cheminformatics:** community MCP servers for **RDKit** provide descriptor calc, substructure search, rendering, and more. Plug and use. 


 - **Materials data:** connect to **Materials Project** via its public API (or an MCP wrapper) for structures, formation energies, and band gaps. 


 - **Custom science:** **FastMCP** + Gemini CLI make it straightforward to publish your lab’s pipeline as a reusable tool, not a bespoke agent. 

**Result:** we standardize the interface (MCP) and reuse the agent shell (Gemini CLI). New capabilities arrive as new MCP servers, instantly composable with existing ones. This reduces integration tax, improves reproducibility, and channels community energy into shared tooling rather than parallel agent rewrites. 


## TL;DR (project intent)

 - **Mission:** build an open, multi-tool ecosystem for chemistry & materials agents by unifying on Codex + MCP.

 - **Why it matters:** agentic science needs interoperable tools, not more siloed agents. 

 - **How you can help:** contribute or refine an MCP server (RDKit, Materials Project, simulations, ELN/SDMS, robo-lab APIs). The agent comes for free. 

## 快速开始
### 安装
### 🚀 快速安装（推荐）
#### 一键安装（最简单）
```bash
curl -fsSL https://raw.githubusercontent.com/dazhaolang/ai-chemkit/main/install.sh | bash
```
#### 或者手动安装
# 克隆仓库
```bash
git clone https://github.com/dazhaolang/ai-chemkit.git
cd ai-chemkit

# 运行交互式安装器
python chemagent_install.py
```

#### 安装RDKit MCP服务器（官方支持）
```bash
# 安装官方的mcp-rdkit包
python chemagent_install.py mcp

# 或单独安装
./install_rdkit_mcp.sh
```

### 📦 安装模式
ChemAgent 提供多种安装模式，满足不同需求：

#### 1. **快速安装** - 自动检测并安装所有功能
```bash
python chemagent_install.py
# 或安装后使用: chemagent install
```

#### 2. **交互式安装** - 选择要安装的组件
```bash
python chemagent_install.py --interactive
# 选择平台、功能、工具等
```

#### 3. **最小安装** - 仅安装核心功能
```bash
python chemagent_install.py --minimal
# 轻量级安装，适合资源受限环境
```

#### 4. **开发者模式** - 包含开发工具
```bash
python chemagent_install.py --profile developer
# 包含测试、代码格式化、类型检查等工具
```

### 🎯 其他选项
```bash
# 查看安装状态
python chemagent_install.py status

# 安装示例文件
python chemagent_install.py examples

# 更新到最新版本
python chemagent_install.py update

# 静默安装（自动化）
python chemagent_install.py --yes --quiet

# 仅安装特定平台
python chemagent_install.py --platform claude-code
python chemagent_install.py --platform gemini-cli

# Gemini CLI 专用安装（基于提示词）
chmod +x install_gemini_simple.sh
./install_gemini_simple.sh

# 查看所有选项
python chemagent_install.py --help
```

## 命令系统

ChemAgent 采用类似 SuperClaude_Framework 的 Markdown 命令定义方式，命令定义简单灵活：

### 命令定义位置
 - **项目级:** `.claude/commands/` - 项目专用命令
 - **用户级:** `~/.claude/commands/` - 个人全局命令
 - **系统级:** ChemAgent 自带的默认命令

### 创建自定义命令
```bash
# 在项目中创建命令
mkdir -p .claude/commands
cat > .claude/commands/my-analysis.md << EOF
---
description: 我的分子分析流程
tools: [read_file, web_search]
---

请执行以下分析步骤：
1. 验证分子结构
2. 计算基本性质
3. 预测ADMET
4. 生成报告
EOF
```

### Claude Code (Cursor) 中使用
安装后，在 Cursor 中直接使用化学命令：

```markdown
# 分析分子
cc-analyze aspirin

# 设计合成路线  
cc-synthesize "CC(=O)OC1=CC=CC=C1C(=O)O"

# 调用化学专家
@chemist 请帮我设计一个环保的阿司匹林合成路线

# 执行工作流
cc-workflow drug-discovery target.pdb
```

### Codex 中使用
```bash
# 分析分子
gemini chem:analyze "CC(=O)OC1=CC=CC=C1C(=O)O" --properties

# 批处理
gemini chem:batch molecules.csv --operation analyze

# 图像识别（Gemini 特色）
gemini chem:analyze molecule.png --image-to-structure

# 使用别名（source ~/.gemini/chemagent_aliases.sh 后）
chem aspirin  # 快速分析
synthesize ibuprofen  # 合成规划
```
```
“Integrate once, compose everywhere.”
