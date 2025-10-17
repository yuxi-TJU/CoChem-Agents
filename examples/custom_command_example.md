# 创建自定义化学命令示例

ChemAgent 采用 SuperClaude_Framework 的命令定义方式，让您可以轻松创建自己的化学工作流。

## 示例 1：创建药物筛选命令

在您的项目中创建 `.claude/commands/drug-screen.md`：

```markdown
---
description: 对候选药物进行全面筛选
tools: [read_file, web_search]
parameters:
  target:
    description: 靶点名称或PDB文件
    required: true
  library:
    description: 化合物库文件（SDF/CSV）
    required: true
---

# 药物筛选工作流 - ${target}

对化合物库 ${library} 进行以下筛选：

## 1. 初步过滤
- 分子量：200-500 Da
- LogP：0-5
- 氢键供体：≤5
- 氢键受体：≤10
- 可旋转键：≤10

## 2. ADMET 预筛选
- 预测口服吸收
- 血脑屏障透过性
- CYP450 抑制
- hERG 毒性

## 3. 分子对接
- 准备受体结构
- 对接所有通过初筛的分子
- 按对接分数排序

## 4. 视觉检查
- 生成前10个分子的对接构象图
- 分析关键相互作用
- 标注氢键和疏水作用

## 5. 生成报告
创建包含以下内容的报告：
- 筛选统计
- 前20个候选物
- SAR 分析
- 建议的后续实验
```

使用命令：
```
cc-drug-screen PDB:1A2B compounds.sdf
```

## 示例 2：创建绿色化学评估命令

创建 `.claude/commands/green-check.md`：

```markdown
---
description: 评估反应的绿色化学指标
tools: [read_file]
---

# 绿色化学评估

评估以下反应的环保性：

## 1. 原子经济性
- 计算原子利用率
- 识别浪费的原子
- 建议改进方案

## 2. E-因子分析
- 计算总废物量
- 分类废物类型
- 与行业标准比较

## 3. 溶剂评估
- 检查溶剂毒性
- 建议绿色替代品
- 评估回收可能性

## 4. 能源效率
- 反应温度要求
- 反应时间
- 能源强度计算

## 5. 安全性评估
- 试剂危害性
- 反应条件风险
- 废物处理要求

## 总体评分
基于12条绿色化学原则给出评分和改进建议。
```

## 示例 3：创建文献综述命令

创建 `.claude/commands/lit-review.md`：

```markdown
---
description: 生成特定主题的文献综述
tools: [web_search]
parameters:
  topic:
    description: 研究主题
    required: true
  years:
    description: 时间范围（如 2020-2024）
    required: false
    default: "last 5 years"
---

# 文献综述：${topic}

时间范围：${years}

## 1. 文献搜索
搜索以下数据库：
- PubMed/MEDLINE
- Web of Science  
- Google Scholar
- ChemRxiv

## 2. 筛选标准
- 相关性评分 > 0.7
- 优先综述文章
- 高影响因子期刊
- 引用次数 > 10

## 3. 内容提取
对每篇文献提取：
- 主要发现
- 使用方法
- 关键数据
- 局限性

## 4. 趋势分析
- 研究热点演变
- 技术发展路线
- 主要研究团队
- 未解决的问题

## 5. 综述撰写
生成包含以下部分的综述：
1. 摘要（200字）
2. 引言
3. 主要进展
4. 方法比较
5. 挑战与机遇
6. 结论
7. 参考文献（ACS格式）
```

## 高级功能

### 参数化命令
命令支持参数，可以创建灵活的模板：

```markdown
---
parameters:
  molecule:
    description: 输入分子
    required: true
  property:
    description: 要优化的性质
    required: true
  constraint:
    description: 约束条件
    required: false
---

优化 ${molecule} 的 ${property} 性质
${constraint ? `同时保持约束：${constraint}` : ''}
```

### 条件逻辑
使用简单的条件渲染：

```markdown
${high_throughput ? '执行批量处理模式' : '执行详细分析'}
${temperature > 100 ? '警告：高温反应，注意安全' : ''}
```

### 工具集成
指定命令所需的工具：

```markdown
---
tools: [read_file, web_search, grep]
---
```

## 最佳实践

1. **命名规范**：使用描述性名称，如 `solubility-optimize.md`
2. **参数验证**：在命令开始时检查必需参数
3. **错误处理**：包含失败情况的处理步骤
4. **输出格式**：明确指定期望的输出格式
5. **文档化**：在 description 中清楚说明命令用途

## 共享命令

您可以通过以下方式共享命令：

1. **Git仓库**：将 `.claude/commands/` 提交到版本控制
2. **命令包**：创建命令集合并发布
3. **团队共享**：通过共享目录分发命令

这种方式让化学工作流的自动化变得简单而强大！
