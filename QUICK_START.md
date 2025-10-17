# ChemAgent 快速开始指南

## 什么是 ChemAgent？

ChemAgent 是一个**化学增强包**，为 Claude Code (Cursor) 和 Gemini CLI 等 AI 编程助手添加专业的化学能力。

**不是**：独立的应用程序或框架
**而是**：AI 助手的化学工具包

## 架构概览

```
用户 → AI助手(Claude/Gemini) → MCP协议 → ChemAgent工具服务器 → 化学工具(RDKit等)
```

## 快速安装

```bash
# 1. 克隆项目
git clone https://github.com/yourusername/chemagent.git
cd chemagent

# 2. 运行安装脚本
python install.py  # 自动检测并安装

# 3. 启动 MCP 服务器
python -m chemagent.mcp_server
```

## 使用方式

### 在 Claude Code (Cursor) 中

```markdown
# 直接使用化学命令
cc-analyze aspirin --properties

# 调用化学专家
@organic-chemist 设计布洛芬的合成路线
```

### 在 Gemini CLI 中

```bash
# 使用 chem: 命令
gemini chem:analyze "CCO"

# 批处理
gemini chem:batch molecules.csv
```

## 核心功能

1. **分子分析** - 属性计算、药物相似性评估
2. **反应预测** - 产物预测、机理生成
3. **合成规划** - 逆合成分析、路线设计
4. **文献搜索** - PubMed、ChemRxiv 检索
5. **批量处理** - 大规模分子筛选

## 工作原理

1. **AI 理解意图** - Claude/Gemini 理解你的化学问题
2. **MCP 调用** - 通过标准协议调用化学工具
3. **工具执行** - RDKit、PubChem 等执行计算
4. **AI 解释结果** - 生成易懂的解释和建议

## 为什么选择这种设计？

- ✅ **利用 AI 的理解能力** - 不需要学习复杂的命令
- ✅ **标准化接口** - 通过 MCP 协议统一工具调用
- ✅ **易于扩展** - 添加新工具只需定义 MCP 接口
- ✅ **跨平台** - 同一套工具可用于多个 AI 平台

## 文件结构

```
chemagent/
├── mcp_tools/          # MCP 工具定义和实现
│   ├── tool_definitions.py  # 工具接口定义
│   └── chemistry_tools.py   # 工具实现
├── enhancers/          # 平台增强器
│   ├── claude_enhancer.py   # Claude Code 集成
│   └── gemini_enhancer.py   # Gemini CLI 集成
├── commands/           # 化学命令实现
├── tools/              # 底层化学工具
│   ├── rdkit_tools.py      # RDKit 封装
│   ├── pubchem.py          # PubChem API
│   └── literature.py       # 文献搜索
└── roles/              # 化学专家角色定义
```

## 常见问题

**Q: 需要安装哪些依赖？**
A: 主要是 RDKit 和一些 Python 库，安装脚本会自动处理。

**Q: 可以离线使用吗？**
A: 基础功能（RDKit）可以离线，数据库查询需要网络。

**Q: 如何添加新工具？**
A: 在 `mcp_tools/tool_definitions.py` 添加定义，在 `chemistry_tools.py` 实现。

**Q: 支持哪些 AI 平台？**
A: 目前支持 Claude Code 和 Gemini CLI，计划支持更多。

## 获取帮助

- 查看 `examples/` 目录的使用示例
- 阅读 `docs/` 目录的详细文档
- 提交 Issue 到 GitHub

## 下一步

1. 尝试基础命令
2. 探索高级功能
3. 根据需求定制工具
4. 贡献新功能

---

记住：ChemAgent 是你的 AI 助手的化学增强，不是替代品。充分利用 AI 的理解能力 + 专业的化学工具 = 强大的化学研究助手！
