# CoChem Agents



<img width="2816" height="1536" alt="Gemini_Generated_Image_9vntz29vntz29vnt" src="https://github.com/user-attachments/assets/b20a37a8-f1a6-43b2-a34f-db14cf4b87e0" />


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



## TL;DR (project intent)

 - **Mission:** build an open, multi-tool ecosystem for chemistry & materials agents by unifying on Codex + MCP.

 - **Why it matters:** agentic science needs interoperable tools, not more siloed agents. 

 - **How you can help:** contribute or refine an MCP server (RDKit, Materials Project, simulations, ELN/SDMS, robo-lab APIs). The agent comes for free. 

## 快速开始
### 安装
#### 克隆仓库
```bash
git clone https://github.com/yuxi-TJU/CoChem-Agents.git
cd ai-chemkit
```

#### 安装依赖
```bash
npm install
npm run build #生成 dist/.
```
#### 配置各 API 密钥
CHEMSPIDER_API_KEY、MATERIALS_PROJECT_API_KEY 等，通过环境变量或 .env.

#### CLI 安装到 Codex
##### 自动安装（推荐）
```bash
npx chemagent
```
##### 手动安装
```bash
npx chemagent-cli install --platform codex --home <路径> --configure-mcp
```

#### 配置 & 注册 MCP 服务
~/.chemagent/mcp_config.json，里面列着所有内置和社区 MCP
```bash
npx chemagent-<name> #注册你需要的服务器
```

#### 安装RDKit MCP服务器（官方支持）
```bash
# 安装官方的mcp-rdkit包
python chemagent_install.py mcp

# 或单独安装
./install_rdkit_mcp.sh
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

### Codex 中使用
```bash
# 分析分子
codex exec --skip-git-repo-check -- "请调用 chemagent_admet_predict 工具, 使用 SMILES CC(=O)OC1=CC=CC=C1C(=O)O, includeRiskDetails true"

# 批处理
batch molecules.csv --operation analyze

# 图像识别
analyze molecule.png --image-to-structure

# 使用别名
analyze aspirin  # 快速分析
synthesize ibuprofen  # 合成规划
```

## Case Study
### Case1
```bash
# 分析分子
codex exec --skip-git-repo-check -- "请调用 chemagent_admet_predict 工具, 使用 SMILES CC(=O)OC1=CC=CC=C1C(=O)O, includeRiskDetails true"
```
```bash
tool chemagent-admet.chemagent_admet_predict({"smiles":"CC(=O)OC1=CC=CC=C1C(=O)O","includeRiskDetails":true})
codex
Ran chemagent_admet_predict for `CC(=O)OC1=CC=CC=C1C(=O)O`. Key outputs:

- Basic props: formula `C9H8O4`, MW 180.04, logP 1.13, logS −1.93, TPSA 63.6 Å², 3 rotatable bonds, 1 HBD/4 HBA.
- Lipinski: no violations; passes.
- Scores: druglikeness −0.478 (moderate), drugScore 0.143.
- Toxicity: high-risk flags for mutagenicity, tumorigenicity, reproductive effects (all cite known issue `dklB@@QmR[fUxUZBBF@@`); no irritation risk noted.
- Warning: overall high-risk toxicity signals present; manual review recommended.
tokens used
1,360
Ran chemagent_admet_predict for `CC(=O)OC1=CC=CC=C1C(=O)O`. Key outputs:

- Basic props: formula `C9H8O4`, MW 180.04, logP 1.13, logS −1.93, TPSA 63.6 Å², 3 rotatable bonds, 1 HBD/4 HBA.
- Lipinski: no violations; passes.
- Scores: druglikeness −0.478 (moderate), drugScore 0.143.
- Toxicity: high-risk flags for mutagenicity, tumorigenicity, reproductive effects (all cite known issue `dklB@@QmR[fUxUZBBF@@`); no irritation risk noted.
- Warning: overall high-risk toxicity signals present; manual review recommended.
```

```
“Integrate once, compose everywhere.”
