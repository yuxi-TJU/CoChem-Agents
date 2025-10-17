# ChemAgent - 化学增强包 for Claude Code & Gemini CLI

## 概述

ChemAgent 是一个专为 **Claude Code (Cursor)** 和 **Gemini CLI** 设计的化学领域增强包，通过预设命令、子代理（sub-agents）、专业工具和优化提示词，让这些强大的 AI 编程助手在化学研究、药物设计和材料科学等领域发挥更大作用。

灵感来源于 [SuperClaude_Framework](https://github.com/SuperClaude-Org/SuperClaude_Framework)，ChemAgent 采用类似的安装包模式，为 AI 编程助手提供即插即用的化学能力增强。

## 核心特性

### Claude Code (Cursor) 增强
- 📝 **cc-系列命令**：`cc-analyze`、`cc-synthesize`、`cc-predict` 等化学专用命令
- 🎭 **化学子代理**：`@organic-chemist`、`@drug-designer` 等专业角色
- 🔧 **自动上下文**：`.cursorrules` 文件自动加载化学上下文
- 🧪 **集成工具**：RDKit、PubChem、ChEMBL 等工具无缝集成

### Gemini CLI 增强  
- 🌟 **chem:系列命令**：`gemini chem:analyze`、`gemini chem:synthesize` 等
- 📋 **基于提示词**：通过规则文件和提示词增强，轻量级集成
- 🖼️ **多模态支持**：利用 Gemini 的图像识别能力处理分子结构图
- 📦 **批处理**：`gemini chem:batch` 批量处理分子数据
- 🎭 **化学角色**：`@drug-designer`、`@synthesist` 等专业角色

### 通用特性
- 🚀 **一键安装**：自动检测并安装到 Claude Code 或 Gemini CLI
- 🧪 **专业工具**：RDKit、PubChem、ChemSpace、OpenBabel 集成
- 📚 **预设提示词**：药物设计、材料科学、有机合成等领域优化
- 🔬 **MCP 支持**：通过 MCP 协议扩展更多工具

## 快速开始

### 安装

### 🚀 快速安装（推荐）

#### 一键安装（最简单）
```bash
curl -fsSL https://raw.githubusercontent.com/dazhaolang/ai-chemkit/main/install.sh | bash
```

#### 或者手动安装
```bash
# 克隆仓库
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
- **项目级**: `.claude/commands/` - 项目专用命令
- **用户级**: `~/.claude/commands/` - 个人全局命令  
- **系统级**: ChemAgent 自带的默认命令

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

### Gemini CLI 中使用

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

## 项目结构

```
chemagent/
├── commands/            # Markdown 命令定义（如 SuperClaude）
│   ├── cc-analyze.md   # 分析命令
│   ├── cc-synthesize.md # 合成命令
│   ├── cc-check.md     # 安全/专利检查
│   └── ...            # 更多命令
├── chemagent/
│   ├── enhancers/      # 平台增强器
│   │   ├── claude_enhancer.py
│   │   └── gemini_enhancer.py
│   ├── commands/       # 命令系统
│   │   └── loader.py   # Markdown 命令加载器
│   ├── roles/          # 化学专家（简化为5个）
│   │   ├── chemist.py
│   │   ├── drug_designer.py
│   │   ├── synthesist.py
│   │   ├── safety_expert.py
│   │   └── data_analyst.py
│   ├── tools/          # 增强工具（ChemCrow功能）
│   │   ├── patent_search.py    # 专利检查
│   │   ├── literature_enhanced.py # 文献搜索
│   │   ├── safety_assessment.py # 安全评估
│   │   └── price_lookup.py     # 价格查询
│   └── mcp_tools/      # MCP 编排器
│       └── orchestrator.py  # 优先调用官方MCP
├── .cursorrules        # Claude Code 规则文件
└── install.py          # 智能安装脚本
```

## 化学命令列表

### Claude Code (cc-系列)
| 命令 | 描述 | 示例 |
|------|------|------|
| `cc-analyze` | 分子分析 | `cc-analyze aspirin --properties` |
| `cc-synthesize` | 合成设计 | `cc-synthesize target.smi --retro` |
| `cc-predict` | 反应预测 | `cc-predict "A + B" --mechanism` |
| `cc-optimize` | 性质优化 | `cc-optimize mol.sdf --target logP=2.5` |
| `cc-drug-design` | 药物设计流程 | `cc-drug-design lead.smi --admet` |
| `cc-screen` | 虚拟筛选 | `cc-screen library.sdf --target protein.pdb` |

### Gemini CLI (chem:系列)
| 命令 | 描述 | 示例 |
|------|------|------|
| `chem:analyze` | 分子分析 | `gemini chem:analyze "CCO"` |
| `chem:synthesize` | 合成规划 | `gemini chem:synthesize aspirin` |
| `chem:batch` | 批量处理 | `gemini chem:batch mols.csv` |
| `chem:visualize` | 结构可视化 | `gemini chem:visualize "c1ccccc1"` |

## 化学子代理（Sub-Agents）

| 子代理 | 专长 | 调用方式 |
|--------|------|----------|
| `@organic-chemist` | 有机合成、机理 | Claude Code 中直接 @ |
| `@drug-designer` | 药物设计、ADMET | Claude Code 中直接 @ |
| `@materials-scientist` | 材料、聚合物 | Claude Code 中直接 @ |
| `@comp-chemist` | 量子化学计算 | Claude Code 中直接 @ |
| `@analytical-chemist` | 光谱、分析方法 | Claude Code 中直接 @ |

## 与其他项目的对比

| 特性 | ChemAgent | ChemCrow | SuperClaude |
|------|-----------|----------|-------------|
| 定位 | AI 助手增强包 | 独立 Agent | Claude 增强框架 |
| 平台支持 | Claude Code + Gemini CLI | LangChain | Claude Code |
| 化学工具 | ✅ 完整集成 | ✅ 完整集成 | ❌ 通用工具 |
| 安装方式 | 一键安装脚本 | pip 安装 | 一键安装脚本 |
| 命令系统 | cc-/chem: 系列 | Python API | sc- 系列 |
| 子代理 | 5 个化学专家 | 无 | 15 个通用代理 |

## 开发路线图

- [ ] 支持更多 AI 平台（GitHub Copilot、Codeium）
- [ ] 集成更多化学数据库（ChemSpider、ZINC）
- [ ] 添加量子化学计算接口
- [ ] 支持分子对接和 MD 模拟
- [ ] 开发 VS Code 扩展
- [ ] 创建 Web UI 界面

## 贡献

欢迎贡献代码、报告问题或提出建议！特别欢迎：
- 新的化学命令实现
- 更多平台的适配器
- 化学工具集成
- 文档和示例

## 致谢

- [SuperClaude_Framework](https://github.com/SuperClaude-Org/SuperClaude_Framework) - 架构灵感
- [ChemCrow](https://github.com/ur-whitelab/chemcrow-public) - 化学工具参考
- RDKit、PubChem 等开源化学工具

## 许可

MIT License
