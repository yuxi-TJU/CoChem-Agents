# ChemAgent for Gemini CLI - 集成方案

## 设计理念

与Claude Code的集成方式类似，ChemAgent对Gemini CLI的增强采用**基于提示词（Prompt-based）**的轻量级方案，而非复杂的代码集成。这种方式的优势：

1. **轻量级**：无需修改Gemini CLI核心代码
2. **灵活性**：用户可以轻松自定义和扩展
3. **兼容性**：适配不同版本的Gemini CLI
4. **易维护**：提示词和规则文件易于更新

## 集成架构

```
Gemini CLI
    ├── 系统提示词 (System Prompts)
    │   └── chemistry.md (化学规则文件)
    ├── 命令提示词 (Command Prompts)
    │   ├── chem_analyze.md
    │   ├── chem_synthesize.md
    │   └── chem_*.md
    ├── 角色定义 (Role Definitions)
    │   ├── @chemist
    │   ├── @drug-designer
    │   └── @synthesist
    └── 命令别名 (Command Aliases)
        └── chem:* 系列命令
```

## 安装方式对比

### 1. 简单脚本安装（推荐）

**文件**: `install_gemini_simple.sh`

这是最轻量的安装方式，仅安装：
- 提示词文件
- 规则文件
- Shell别名

```bash
./install_gemini_simple.sh
```

**优点**：
- 纯提示词驱动，无代码依赖
- 安装快速，占用空间小
- 易于理解和自定义

### 2. Python完整安装

**文件**: `install_gemini.py`

包含更多功能：
- 提示词和规则
- Python处理脚本（可选）
- 扩展配置
- 示例文件

```bash
python install_gemini.py
```

**优点**：
- 功能更完整
- 支持高级特性
- 可与ChemAgent核心集成

## 核心文件说明

### 1. gemini_rules.md
主要的化学规则文件，定义了：
- 化学命令的解释方式
- 响应格式
- 专业知识库
- 安全和伦理准则

### 2. 提示词文件 (Prompts)
每个命令对应的详细提示词：
- `chem_analyze.md` - 分子分析提示
- `chem_synthesize.md` - 合成规划提示
- `chem_drug_design.md` - 药物设计提示
- `chem_safety.md` - 安全评估提示

### 3. 配置文件
`chemistry_config.json` - 自定义行为配置

### 4. Shell集成
`load_chemistry.sh` - 命令别名和便捷函数

## 工作原理

### 1. 命令识别
当用户输入化学相关命令时，Gemini CLI会：
1. 识别`chem:`前缀或化学关键词
2. 加载对应的提示词
3. 应用化学规则文件
4. 生成专业响应

### 2. 角色切换
使用`@role`语法激活特定专家角色：
```bash
gemini "@drug-designer 设计一个激酶抑制剂"
```

### 3. 上下文增强
规则文件提供持续的化学上下文：
- 自动验证分子结构
- 应用化学约束
- 考虑安全因素
- 引用专业知识

## 与Claude Code对比

| 特性 | Claude Code | Gemini CLI |
|------|------------|------------|
| 集成方式 | .cursorrules文件 | 规则+提示词文件 |
| 命令前缀 | cc- | chem: |
| 角色定义 | Markdown文件 | 提示词内嵌 |
| 多模态 | 文本为主 | 支持图像输入 |
| 安装方式 | 自动检测 | Shell脚本 |

## 自定义和扩展

### 添加新命令
创建新的提示词文件：
```bash
cat > ~/.gemini/prompts/chem_custom.md << EOF
# 自定义化学命令
执行特定的分析...
EOF
```

### 添加新角色
在规则文件中定义：
```markdown
### @custom-role
专门处理...
```

### 修改行为
编辑配置文件：
```json
{
  "chemistry": {
    "custom_setting": "value"
  }
}
```

## 最佳实践

1. **保持提示词简洁**：每个提示词专注单一功能
2. **使用标准格式**：统一输入输出格式
3. **版本控制**：跟踪提示词变更
4. **测试验证**：定期测试命令响应
5. **用户反馈**：根据使用情况优化

## 故障排除

### 常见问题

1. **命令未识别**
   - 检查是否加载了chemistry规则
   - 验证命令前缀是否正确

2. **响应不专业**
   - 确认规则文件路径正确
   - 检查Gemini配置

3. **性能问题**
   - 减少提示词复杂度
   - 使用缓存机制

## 未来发展

1. **MCP集成**：当Gemini支持MCP时集成官方工具
2. **插件系统**：开发Gemini原生插件
3. **云函数**：利用Google Cloud Functions
4. **API接口**：提供REST API访问

## 贡献指南

欢迎贡献新的：
- 提示词模板
- 化学规则
- 命令示例
- 使用案例

提交PR到：https://github.com/yourusername/chemagent
