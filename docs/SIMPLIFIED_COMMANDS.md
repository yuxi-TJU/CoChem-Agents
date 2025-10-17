# ChemAgent 精简命令和子代理系统

参考 SuperClaude_Framework 的设计理念：少而精，灵活组合。

## 🎭 核心子代理（5个足够）

### 1. @chemist
```yaml
角色: 通用化学专家
覆盖: 有机化学、无机化学、分析化学
适用: 大部分化学问题的首选
```

### 2. @drug-designer  
```yaml
角色: 药物设计专家
覆盖: ADMET、先导优化、靶点分析
适用: 药物研发相关
```

### 3. @synthesist
```yaml
角色: 合成专家
覆盖: 逆合成、反应设计、工艺优化
适用: 合成路线规划
```

### 4. @safety-expert
```yaml
角色: 安全与合规专家
覆盖: 毒性、法规、专利、环保
适用: 安全评估和合规检查
```

### 5. @data-analyst
```yaml
角色: 化学数据分析师
覆盖: SAR分析、化学信息学、机器学习
适用: 数据处理和模式识别
```

## 📋 核心命令（15个精选）

### 分析类（3个）

#### cc-analyze
```bash
# 万能分析命令
cc-analyze <input> [--deep] [--quick] [--focus area]

# 自动识别输入类型（SMILES、名称、文件）
# --deep: 深度分析（所有可能的属性）
# --quick: 快速分析（只算关键属性）
# --focus: 聚焦特定领域（drug/material/synthesis）
```

#### cc-compare
```bash
# 比较命令（支持多个分子）
cc-compare <mol1> <mol2> [...] [--aspect property/structure/activity]
```

#### cc-search
```bash
# 统一搜索命令
cc-search <query> [--in database/literature/patents] [--similar] [--substructure]
```

### 设计类（3个）

#### cc-design
```bash
# 分子设计（根据上下文自动选择策略）
cc-design <objective> [--from scaffold] [--optimize properties] [--novel]

# 示例：
cc-design "improve solubility" --from aspirin
cc-design "new kinase inhibitor" --novel
```

#### cc-optimize
```bash
# 优化命令（自动识别优化目标）
cc-optimize <molecule> [--for admet/synthesis/cost] [--constraints]
```

#### cc-synthesize
```bash
# 合成规划（整合逆合成和路线设计）
cc-synthesize <target> [--practical] [--green] [--steps N]
```

### 预测类（2个）

#### cc-predict
```bash
# 通用预测命令
cc-predict <what> <for molecule/reaction> 

# 示例：
cc-predict properties <molecule>
cc-predict products <reaction>
cc-predict activity <molecule> --target kinase
```

#### cc-simulate
```bash
# 模拟命令（对接、动力学等）
cc-simulate <type> <input> [--conditions]

# 示例：
cc-simulate docking <ligand> <protein>
cc-simulate reaction <reactants> --temperature 100C
```

### 工作流类（3个）

#### cc-workflow
```bash
# 执行预定义工作流
cc-workflow <drug-discovery|lead-opt|safety-check> <input>
```

#### cc-batch
```bash
# 批处理命令
cc-batch <operation> <input-file> [--parallel] [--filter]
```

#### cc-report
```bash
# 生成报告
cc-report <type> <data> [--format pdf/html] [--sections]
```

### 辅助类（4个）

#### cc-explain
```bash
# 解释概念或结果
cc-explain <concept/result> [--level beginner/expert] [--visual]
```

#### cc-suggest
```bash
# 智能建议
cc-suggest [--next-steps] [--alternatives] [--improvements]
```

#### cc-check
```bash
# 检查命令（安全、专利、质量等）
cc-check <safety|patent|quality|compliance> <input>
```

#### cc-help
```bash
# 上下文相关的帮助
cc-help [command] [--examples] [--tips]
```

## 🔗 命令组合哲学

### 简单任务 = 单个命令
```bash
cc-analyze aspirin
cc-search "kinase inhibitors"
cc-synthesize ibuprofen
```

### 复杂任务 = 命令组合
```bash
# 药物发现流程
cc-search "EGFR inhibitors" | 
cc-analyze --focus drug |
cc-optimize --for admet |
cc-synthesize --practical

# 安全评估
cc-analyze <molecule> |
cc-check safety |
cc-check patent |
cc-report safety-assessment
```

### 让 AI 智能选择
用户只需描述目标，AI 自动组合合适的命令：

```
用户: "帮我设计一个更好的阿司匹林"

AI 自动执行:
1. cc-analyze aspirin
2. cc-design "improve aspirin" --optimize "reduce side effects"
3. cc-check safety <new-molecule>
4. cc-synthesize <new-molecule> --practical
```

## 💡 设计原则

### 1. 命令通用化
- 一个命令覆盖多种相关功能
- 通过参数控制具体行为
- AI 根据上下文智能选择

### 2. 自动识别
- 自动识别输入类型（SMILES、名称、文件）
- 自动识别任务类型（药物、材料、合成）
- 自动选择合适的工具

### 3. 渐进式复杂度
- 基础用法简单：`cc-analyze aspirin`
- 高级用法灵活：`cc-analyze aspirin --deep --focus drug --export json`
- AI 帮助用户逐步深入

### 4. 组合优于配置
- 少量通用命令
- 通过组合实现复杂功能
- 管道和链式调用

## 🎯 与 SuperClaude_Framework 对比

| 方面 | SuperClaude | ChemAgent |
|------|-------------|-----------|
| 命令数量 | ~25个 | 15个 |
| 子代理数 | ~15个 | 5个 |
| 设计理念 | 通用开发 | 化学专用 |
| 命令前缀 | sc- | cc- |
| 核心价值 | 开发效率 | 化学研究 |

## 📝 实际使用示例

### 场景1：快速分子分析
```bash
# 最简单
cc-analyze caffeine

# 需要详细信息
cc-analyze caffeine --deep

# 关注特定方面
cc-analyze caffeine --focus drug
```

### 场景2：药物优化
```bash
# 一步到位
cc-optimize "lead-compound.smi" --for admet

# 或者分步控制
cc-analyze "lead-compound.smi" |
cc-design "reduce toxicity" |
cc-check safety
```

### 场景3：文献研究
```bash
# 简单搜索
cc-search "PROTAC degraders"

# 深度研究
cc-search "PROTAC degraders" --in literature |
cc-analyze --batch |
cc-report literature-review
```

## 🚀 优势

1. **易学易用** - 15个命令容易记忆
2. **功能完整** - 通过组合覆盖所有需求  
3. **AI 友好** - 命令语义清晰，AI容易理解
4. **灵活扩展** - 参数系统支持新功能
5. **保持简洁** - 不会命令爆炸

这样的设计既保持了 SuperClaude_Framework 的简洁优雅，又针对化学领域做了优化。
