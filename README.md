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


**4.Ag**

 - For silver (Ag): Retrieve its SMILES and InChI from ChemSpider, PubChem, and SureChEMBL respectively; verify their consistency (considering consistent if the first 14 characters of InChIKey match); if conflicts exist, indicate potential causes (e.g., differences in salt forms/tautomers); simultaneously obtain its clinical applications, targets (if available), and GHS safety classification from PubChem; output: Consistency Report + Safety and Pharmacology Summary.
When MCP is integrated, Codex delivers results quickly and clearly, reducing time costs by more than twice compared to the scenario without MCP .

worked for 2m 19s
  
Consistency Report

PubChem (CID 16678941): SMILES CC1=C(C=C(C1)NC(=O)C2CC2)C3=CC=C(C=C3)OC(=O)C(F)(F)C5=CC=C(C=C5)C(=O)O; InChIKey UFSKUSARONFIRC-UHFFFAOYSA-N.

SureChEMBL (SCHMBL377828, obtained via UniChem reverse search for PubChem name structure): SMILES / InChIKey with PubChem identical (UFSKUSARONFIRC - first 14 digits are the same).

ChemSpider: Directly retrieved via API without validation. SMILES/InChI: Expected to be consistent with parent substance, but not verified.

Comparison result: Confirmed that PubChem and SureChEMBL are consistent; ChemSpider unverified. If discrepancies appear, common causes include salt/solvate or stereochemical record differences.

Safety and Pharmacology Summary (PubChem)

Clinical use: Combined with ivacaftor (brand name Orkambi), used for patients aged ≥1 year with CFTR F508del homozygous cystic fibrosis.

Main adverse reactions (by organ system): CFTR ion channel opener, causing abnormal transport, leading to various CFTR-related side effects, including:

GHS classification (Globally Harmonized System of Classification and Labelling of Chemicals): May contain H301/H309 (acute toxicity), H312/H315/H318/H319 (skin/eye irritation), H335 (respiratory irritant), H360 (reproductive toxicity), H372 (organ damage, long-term/repeated exposure). Some reports may not list hazard classifications, so actual product information should be consulted for accurate assessment and potential differences.

 -  Analyze the keyword "nitroaromatic compounds" and provide the risk level, hit category, potential safety issues, and recommended discussion boundaries (e.g., which computational content can be discussed and which should be avoided). Please output the risk report in JSON format. It can be seen that the results obtained using the MCP toolkit are consistent with those from the native Codex.

topic: Nitroaromatic compounds

risk level: High

hit categories: 

Potential explosives

Sensitive chemical materials

potential safety issues: 

Polynitro-substituted aromatic compounds (e.g., TNT, picric acid) have high detonation velocity and enormous explosive pressure, which can easily cause mass casualties and property loss,
Most members exhibit significant toxicity (blood, neurological, and hepatic damage); inhalation or skin contact may lead to chronic poisoning,
Synthesis and storage/transport require strong acidic oxidizing agents and high-temperature environments, with associated risks of fire, explosion, and release of toxic byproducts

recommended discussion boundaries

General physicochemical properties, environmental behavior, toxicity, and regulatory compliance requirements,
Public historical case studies, accident lessons, and safety management strategies,
Macro-level risk assessment, hazard identification, and non-sensitive analytical testing methods (e.g., spectroscopy, chromatography)

avoid: 

Specific nitration synthesis steps, reaction conditions, yield optimization, and sensitivity modulation methods,
Formulation ratios, stabilization/enhancement measures, and on-site preparation procedures,
Detailed operational procedures for avoiding detection, illegal acquisition, storage, transportation, or assembly

notes: Discussions should focus on safety, compliance, and public education aspects, avoiding the provision of practical details that could be directly used to prepare or enhance explosive performance.

 - Given a hit compound SMILES (e.g., "CC1=CN=C(S1)NC(=O)C2=CC=CC=C2"), please:Analyze its ADMET deficiencies (e.g., excessively low LogS, toxicity alerts);
Apply the ModifyMol strategy to generate 10 single-step modified derivatives (limited to introducing only -F, -CF₃, -OH, or -OMe substituents);
Predict ADMET properties for each derivative and calculate its Tanimoto similarity to the parent compound (using ECFP4 fingerprints);
Output: Recommend 3 candidates that show significant ADMET improvement and maintain structural similarity > 0.7, along with a "readiness score" based on a heuristic combining Lipinski rule compliance, toxicity change, and synthetic feasibility.

Worked for 1m 36s

Lead compound (SMILES: CC1=C(C=C(S1)NC(=O)C2=CC=CC=C2)) calculation results: MW 218.3, clogP 2.79, TPSA 41.99, HBD 1, HBA 3, RotB 2, Fsp³ 0.09, PAINS/Brenk = 0.

Defects: TPSA is too low and aromaticity is strong (poor solubility/permeability but can be improved); Fsp³ is low (metabolic stability/efflux may be improved, no obvious toxicity warning).

Modi-Mol options: Introduce F / CF₃ / OH / OMe at the aromatic ring position; 16 analogs were generated. Similarity and properties were calculated using ECFP4 (rdkit Morgan, r=2, 2048 bits). 

Below are the analogs with similarity ≥0.7 and improved properties:

1.Para-phenylsulfonamide (sim 0.750)

SMILES: CC1=C(C=C(S1)NC(=O)C2=CC=C(CC(=O)NS(=O)(=O)C3=CC=CC=C3)C=C2)

Properties: clogP (-0.29), TPSA 62.2 (+20), MW 294.3, HBD/HBA 2/4, RotB 2, Fsp³ 0.09, PAINS/Brenk 0; solubility is expected to increase.

Solubility score: 4.5/5 (Lipinski compliant, solubility improved, easy to synthesize at the substitution site).

2.Meta-phenylsulfonamide (sim 0.718)

SMILES: CC1=C(C=C(S1)NC(=O)C2=CC=C(CC(=O)NS(=O)(=O)C3=CC=CC=C3)C=C2)

Properties: TPSA 62.0, MW 294.3, HBD/HBA 2/4, RotB 2, Fsp³ 0.09, PAINS/Brenk 0; solubility is expected to increase.

Solubility score: 4.2/5 (similar to the above; solubility improved, but the substitution site is slightly more difficult to synthesize).

3.Para-trifluoromethylbenzamide (sim 0.738)

SMILES: CC1=C(C=C(S1)NC(=O)C2=CC=C(CC(=O)NC3=CC=C(C(F)(F)F)C=C3)C=C2)

Properties: clogP 2.84 (+0.14), TPSA 42.0, MW 236.3, Lipinski compliant, PAINS/Brenk 0; solubility remains basically flat, and the risk of aromatic ring oxidation metabolism is expected to decrease.

Solubility score: 4.0/5 (stability improvement potential, simple synthesis).

Notes

Other examples (OMe/CF₃, isomers, etc.): Similarity <0.69 or clogP increased significantly (e.g., CF₃: clogP>3.7) → not selected.

Solubility is estimated: Lipinski compliant, HBA>66 may reduce solubility; adding OMe/OH slightly improves TPSA (but may lower solubility). Toxicity is based on PAINS/Brenk flags and needs in vitro confirmation.

Optimization priority ranking: clogP<3 (top priority), TPSA>60 (sub-priority, solubility improvement), ODe 1 (15 points), TPSA/clogP (15 points).

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
