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

## Case Study
### Case1
```codex
请帮我设计一个环保的阿司匹林合成路线（请调用mcp工具并告诉我你调用了哪些mcp工具）
```
```codex
• Green Aspirin Route

  - Used MCP tool chemagent-extras: chem_synthesize to frame the eco-design constraints and ensure route planning stayed
    aligned with renewable feedstocks/recyclable catalysts.
  - Start from bio-derived phenol converted to sodium phenolate; perform a Kolbe–Schmitt fixation with recycled/biogenic
    CO₂ at 130 °C and 8–10 bar in water or glycerol/water to give sodium salicylate, then mildly acidify with recyclable
    dilute acetic acid to obtain salicylic acid.
  - Dry salicylic acid and acetylate with bio-ethanol–derived acetic anhydride (0.95–1.0 equiv) in a recyclable solvent
    such as 2-MeTHF or EtOAc at 70–80 °C, using a reusable solid acid catalyst (Amberlyst-15 or supported succinic acid/
    Fe³⁺); catalyst is filtered, washed, and reused.
  - Quench with a small amount of reclaimed water to hydrolyze excess anhydride; phase-separate, recover acetic acid for
    reconversion to acetic anhydride, and recycle organic solvent by distillation.
  - Crude aspirin is recrystallized from ethanol/water (both bio-sourced) with mother liquors redistilled; solid is
    dried under reduced pressure powered by low-carbon electricity.
  - Waste minimization: no mineral acids beyond trace acetic acid, closed-loop solvents, solid acid catalyst
    regeneration, and CO₂ capture integration keep the E-factor low and the process compliant with green-chemistry
    principles.

  Next steps: 1) quantify mass balance/E-factor for target scale; 2) run small Kolbe–Schmitt trials to confirm CO₂
  recycle efficiency.
```

### Case2
```codex
请帮我分析作为有机中间合成体，邻苯二甲酸单乙酸酯可以制备什么材料，并给出制备方法（请调用mcp工具并告诉我你调用了哪些mcp工具）
```
```codex
• Downstream Materials

  - chemagent-extras: chem_synthesize (targets: cellulose acetate phthalate; DEHP plasticizer) flagged phthalic
    monoacetate as a versatile half-ester that can be driven either toward anhydrides for acylation chemistry or toward
    fully esterified products with long-chain alcohols. Building on those prompts, here are practical materials you
    can access:
  - Cellulose acetate phthalate (enteric coating polymer) — Convert phthalic monoacetate to the mixed anhydride by
    treating with acetic anhydride (1.2 eq) and catalytic pyridine at 80 °C; strip excess anhydride to get a reactive
    acylation mixture still containing the monoacetate. Add pre-activated cellulose acetate (DS≈1.7) in acetone and
    raise to 50 °C while adding triethylamine as HCl scavenger; the phthalyl group installs at remaining hydroxyls,
    giving cellulose acetate phthalate after neutralization and precipitation into water. Wash, dry, then mill for
    tablet coating use.
  - Di-2-ethylhexyl phthalate (DEHP) or other PVC plasticizer-grade diesters — Hydrolyze the monoacetate (50 % aq. NaOH,
    80 °C) to revert the free phthalic acid half, re-esterify via acid chloride formation (SOCl₂, catalytic DMF, 70 °C)
    to produce phthalyl chloride; react with 2-ethylhexanol (2.2 eq) and pyridine at 0–25 °C, then heat to 100 °C for
    completion. Alternatively, run a direct Fischer esterification on the monoacetate with excess 2-ethylhexanol under
    Dean–Stark to displace the terminal acetate and generate the diester; vacuum-strip and polish through a neutral
    alumina bed to reach plasticizer purity.
  - Alkyd/unsaturated polyester resins (coil coatings, wire enamels) — Use the monoacetate as a “half-ester” feed in an
    alkyd cook: charge it with glycerol (or pentaerythritol) and fatty-acid-modified monoglycerides at 220 °C under N₂.
    The monoacetate first cleaves to phthalic acid in situ, then esterifies to form the phthalate-rich segments that
    lend hardness and gloss. Apply vacuum toward the end to pull off acetic acid and drive polycondensation to acid
    value <15 mg KOH/g; the resulting resin dissolves in high-boiling esters for coating formulations.
  - Specialty imide monomers (e.g., phthalic-imide curing agents) — Convert monoacetate to phthalic anhydride via
    vacuum pyrolysis (180 °C, 20 mmHg) or acetic anhydride reflux with catalytic ZnCl₂; react with ammonia or amines
    (e.g., 4,4′-diaminodiphenyl ether) in glacial acetic acid at 120 °C to form amic acids, then dehydrate with acetic
    anhydride to imide monomers used in high-Tg polyimides.

  Each route leverages the monoacetate either as a latent anhydride (for polymer acylations) or as a controllable half-
  ester that can be selectively converted to acid chlorides/diesters, aligning with the MCP guidance on using it as a
  modular intermediate.
```
