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

### Codex 中使用
```bash
# 分析分子
analyze "CC(=O)OC1=CC=CC=C1C(=O)O" --properties

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
### Case2
```bash
# 使用别名
/sc:chem-analyze aspirin --source name --focus drug
```
```bash
Aspirin (acetylsalicylic acid, PubChem CID 2244) standardizes to SMILES CC(=O)OC1=CC=CC=C1C(=O)O and InChIKey
  BSYNRYMUTXBXSQ-UHFFFAOYSA-N. Canonical properties pulled from PubChem/PUG:

  - Formula C9H8O4, MW 180.16 g/mol.
  - TPSA 63.6 Å²; LogP (XLogP3-AA) 1.2.
  - H-bond donors 1, acceptors 4; rotatable bonds 3.

  Lipinski screen is fully satisfied (MW<500, HBD≤5, HBA≤10, LogP<5, TPSA<140, <10 rotatable bonds), aligning with
  aspirin’s established oral bioavailability. Chemically the molecule combines a salicylate core with an acetyl ester;
  hydrolysis (esterases or alkaline media) releases salicylic acid, the primary COX inhibitor. The aromatic/acetyl motif
  yields modest lipophilicity while the carboxylate gives aqueous solubility when ionized (pKa ≈ 3.5).

  Safety considerations are dominated by systemic COX-1 inhibition (GI irritation/bleeding), salicylate toxicity at
  high dose (tinnitus, metabolic acidosis), Reye’s syndrome in children with viral illness, and bronchospasm in aspirin-
  exacerbated respiratory disease; platelet inhibition prolongs bleeding. Storage needs to minimize hydrolysis (dry,
  cool, pH-neutral conditions).

  Suggested follow-ups if you need more depth:

  1. Run chemagent-admet MCP on the SMILES for fresh descriptor/tox panels tailored to your models.
  2. Use /sc:chem-intel "aspirin" --scope safety if you want summarized regulatory/adverse-event notes tied to the
     latest labels.
```
### Case3
```bash
# 合成规划
/sc:chem-intel "Ibuprofen synthesis" --scope synthesis
```
```bash
Route Overview

  - Standard Boots/Hoechst process begins from isobutylbenzene: para-selective Friedel–Crafts acylation with acetyl
    chloride/AlCl₃ furnishes 4-isobutylacetophenone in ~90 % p-selectivity after distillation.
  - Darzens-like α-halogenation (Br₂/HBr) of the acetyl side chain followed by base-induced hydrolysis rearranges the α-
    bromo ketone to the corresponding α-hydroxy acid; dehydration yields 4-isobutyl-α-methylstyrene oxide intermediates.
  - Carbonylation (HF/BF₃ or HCN/H₂SO₄ route) inserts the extra carbon: formation of the nitrile (via HCN addition) and
    subsequent hydrolysis under acidic or basic conditions leads to ibuprofen’s propionic acid motif.
  - Final steps: Racemate resolution optional (common product is racemic); neutralize/isolate as the free acid, then
    crystallize to pharmacopeial spec (drying ≤1 % moisture).

  Practical Notes

  - Maintain strict temperature control during Friedel–Crafts (≤10 °C charge, raise to 40 °C) to minimize ortho-isomer
    formation.
  - Handle HF/BF₃ or HCN sequences with full alkylation-level PPE and scrubbers; consider the “green” BHC catalytic
    carbonylation (Pd-catalyzed carbonylation of 4-isobutylstyrene with CO and HCl) if hazardous reagents are
    problematic.
  - Typical yields over the 3-step Boots sequence are 60–65 %; BHC process improves overall yield to >75 % and reduces
    waste.

  Next Steps

  1. Confirm available feedstocks (isobutylbenzene, acetyl chloride, HF/BF₃ or Pd/CO setup) and decide whether to run
     Boots or BHC variant.
  2. Draft detailed batch instructions (stoichiometry, solvent charge, quench) for each step and run a 100 g pilot to
     lock selectivity and impurity profile.
  3. Plan workup/purification (crystallization from heptane or IPA/water) and analytical release (HPLC assay, chiral
     excess if needed).
```

```
“Integrate once, compose everywhere.”
