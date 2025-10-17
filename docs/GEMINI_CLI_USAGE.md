# ChemAgent for Gemini CLI - 使用指南

## 概述

ChemAgent 通过提示词（prompts）和规则文件为 Gemini CLI 提供化学专业能力增强。与 Claude Code 的增强方式类似，这是一种轻量级的集成方式，无需复杂的代码集成。

## 安装方式

### 方法1：简单脚本安装（推荐）

```bash
# 使用提示词基础的安装脚本
chmod +x install_gemini_simple.sh
./install_gemini_simple.sh
```

这将安装：
- 化学规则文件到 `~/.gemini/rules/`
- 化学提示词到 `~/.gemini/prompts/`
- 命令别名和配置

### 方法2：手动安装

1. 复制规则文件：
```bash
mkdir -p ~/.gemini/rules
cp gemini_rules.md ~/.gemini/rules/chemistry.md
```

2. 在 Gemini CLI 配置中引用规则：
```json
{
  "system_prompts": {
    "chemistry": "~/.gemini/rules/chemistry.md"
  }
}
```

## 使用方式

### 1. 直接命令模式

在 Gemini CLI 中直接使用化学命令：

```bash
# 分析分子
gemini "chem:analyze aspirin"
gemini "分析阿司匹林的分子性质"

# 合成规划
gemini "chem:synthesize ibuprofen"
gemini "如何合成布洛芬？"

# 反应预测
gemini "chem:predict benzene + Br2 with FeBr3"
gemini "预测苯与溴在三溴化铁催化下的反应"
```

### 2. 角色模式

使用特定的化学专家角色：

```bash
# 药物设计师
gemini "@drug-designer 设计一个选择性COX-2抑制剂"

# 合成专家
gemini "@synthesist 设计紫杉醇的全合成路线"

# 安全专家
gemini "@safety-expert 评估二氯甲烷的安全风险"

# 数据分析师
gemini "@data-analyst 分析这组化合物的构效关系"
```

### 3. 自然语言模式

直接用自然语言提问，Gemini 会自动识别化学相关查询：

```bash
gemini "这个分子的LogP是多少：CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
gemini "比较阿司匹林和布洛芬的药理作用"
gemini "解释什么是Diels-Alder反应"
```

## 命令参考

### 核心命令（15个）

#### 分析类
- `chem:analyze <molecule>` - 分析分子结构和性质
- `chem:compare <mol1> <mol2>` - 比较多个分子
- `chem:search <query>` - 搜索数据库和文献

#### 设计类
- `chem:design <requirements>` - 基于需求设计分子
- `chem:optimize <molecule> <property>` - 优化分子性质
- `chem:synthesize <target>` - 设计合成路线

#### 预测类
- `chem:predict <reaction/property>` - 预测反应或性质
- `chem:simulate <molecule>` - 运行分子模拟

#### 工作流类
- `chem:workflow <name>` - 执行预定义工作流
- `chem:batch <file>` - 批量处理分子
- `chem:report <data>` - 生成综合报告

#### 辅助类
- `chem:explain <concept>` - 解释化学概念
- `chem:suggest <context>` - 提供建议
- `chem:check <type> <target>` - 检查安全/专利/合规
- `chem:help` - 获取帮助

## 示例用例

### 药物发现工作流

```bash
# 1. 搜索已知抑制剂
gemini "chem:search PDE5 inhibitors"

# 2. 分析先导化合物
gemini "chem:analyze sildenafil"

# 3. 设计类似物
gemini "@drug-designer 基于西地那非设计新的PDE5抑制剂"

# 4. 优化ADMET
gemini "chem:optimize [新分子SMILES] solubility"

# 5. 评估安全性
gemini "chem:check safety [新分子SMILES]"

# 6. 规划合成
gemini "chem:synthesize [新分子SMILES]"
```

### 批量分析

```bash
# 准备分子列表文件 molecules.txt
echo "aspirin
ibuprofen
paracetamol" > molecules.txt

# 批量分析
gemini "chem:batch molecules.txt --operation analyze"

# 或直接在命令中列出
gemini "分析这些分子的药物相似性：阿司匹林、布洛芬、对乙酰氨基酚"
```

### 反应机理研究

```bash
# 预测反应
gemini "chem:predict CH3CHO + CH3MgBr"

# 解释机理
gemini "chem:explain Grignard反应的机理"

# 优化条件
gemini "@synthesist 如何优化Grignard反应的产率"
```

## 输出格式

ChemAgent 会根据查询类型自动格式化输出：

### 分子分析输出
```
Structure: [SMILES/InChI]
Properties:
  - MW: 180.16 g/mol
  - LogP: 1.19
  - TPSA: 63.60 Ų
  - Drug-likeness: Pass
Visualization: [2D structure]
Notes: [observations]
```

### 合成路线输出
```
Target: [structure]
Retrosynthesis:
  Step 1: [disconnection]
  Step 2: [disconnection]
Forward Synthesis:
  1. [reaction + conditions]
  2. [reaction + conditions]
Considerations: [practical notes]
```

## 高级功能

### 1. 多模态输入

Gemini CLI 支持图像输入，可以识别分子结构图：

```bash
gemini "分析这个分子结构" --image molecule.png
```

### 2. 上下文保持

Gemini 会记住对话上下文：

```bash
gemini "chem:analyze caffeine"
gemini "它的溶解度如何？"  # 自动理解指咖啡因
gemini "设计一个类似物"     # 基于咖啡因设计
```

### 3. 自定义提示词

创建自定义化学提示词：

```bash
# 创建自定义提示词
cat > ~/.gemini/prompts/my_analysis.md << EOF
分析分子时额外包括：
- 专利状态
- 商业供应商
- 价格范围
- 相关文献
EOF

# 使用自定义提示词
gemini --prompt my_analysis "analyze aspirin"
```

## 配置选项

编辑 `~/.gemini/config/chemistry_config.json` 自定义行为：

```json
{
  "chemistry": {
    "default_role": "chemist",
    "auto_validate_structures": true,
    "include_3d": false,
    "citation_style": "ACS",
    "safety_warnings": true,
    "language": "zh-CN"  // 中文输出
  }
}
```

## 故障排除

### 问题：命令未识别
确保已加载化学增强：
```bash
source ~/.gemini/load_chemistry.sh
```

### 问题：结构无法解析
使用标准SMILES格式或化合物名称：
```bash
gemini "chem:analyze CC(=O)OC1=CC=CC=C1C(=O)O"  # SMILES
gemini "chem:analyze aspirin"                   # 名称
```

### 问题：角色未激活
确保使用正确的角色标记：
```bash
gemini "@drug-designer ..."  # 正确
gemini "drug-designer ..."   # 错误
```

## 最佳实践

1. **明确指定格式**：指明输入是SMILES、InChI还是名称
2. **使用标准命名**：优先使用IUPAC名称或常用名
3. **批量处理大数据**：对多个分子使用batch命令
4. **保存重要结果**：使用 `--output` 保存分析结果
5. **验证关键数据**：对关键结果进行交叉验证

## 与其他工具集成

ChemAgent 可以与其他工具配合使用：

```bash
# 结合 RDKit
gemini "chem:analyze" | rdkit-process

# 导出到 ChemDraw
gemini "chem:analyze aspirin --format mol" > aspirin.mol

# 生成报告
gemini "chem:report" --format pdf > report.pdf
```

## 更多资源

- 示例集合：`~/.gemini/chemistry_examples.md`
- 完整规则：`~/.gemini/rules/chemistry.md`
- 命令参考：`gemini chem:help`
- GitHub：https://github.com/yourusername/chemagent
