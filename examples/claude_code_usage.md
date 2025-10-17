# ChemAgent for Claude Code (Cursor) 使用示例

## 安装后的使用

### 1. 化学命令使用

在 Cursor 中，安装 ChemAgent 后，你可以直接使用 cc- 系列命令：

```markdown
# 分析分子
cc-analyze aspirin --properties --druglike --visualize

# 输出示例：
## Molecular Analysis
**Structure**: CC(=O)OC1=CC=CC=C1C(=O)O
**Properties**:
- MW: 180.16
- LogP: 1.19
- TPSA: 63.6
- HBD: 1, HBA: 4
**Drug-likeness**: Passes Lipinski's Rule of Five ✓
```

### 2. 子代理调用

使用 @ 符号调用专业化学子代理：

```markdown
@organic-chemist 请设计一个环保的布洛芬合成路线

# 子代理会提供：
1. 逆合成分析
2. 前向合成步骤
3. 绿色化学考虑
4. 反应条件优化
```

### 3. 复杂工作流

#### 药物设计工作流
```markdown
cc-drug-design "CN1C=NC2=C1C(=O)N(C(=O)N2C)C" \
  --optimize-admet \
  --check-patents \
  --suggest-analogs

# 会执行：
1. ADMET 属性分析
2. 专利检索
3. 类似物生成
4. 优化建议
```

#### 反应预测
```markdown
cc-predict \
  --reactants "CC=O" "CCO" \
  --conditions "H+, 60°C" \
  --mechanism \
  --yield-estimate

# 输出：
- 产物：CH3CH(OCH2CH3)2 (缩醛)
- 机理：酸催化缩醛化
- 预期产率：85-90%
```

### 4. 批量处理

```markdown
cc-screen molecules.sdf \
  --filter "lipinski" \
  --filter "mw<500" \
  --filter "logp<5" \
  --output filtered.sdf
```

### 5. 与 Claude 的自然交互

```markdown
User: 我有一个分子 SMILES: CC(C)Cc1ccc(cc1)C(C)C(O)=O，
      请帮我分析它的药物性质并提出优化建议。

Claude (with ChemAgent): 
让我使用 ChemAgent 工具来分析这个分子...

[自动调用 cc-analyze]
这是布洛芬（Ibuprofen）的结构。

**分析结果**：
- 分子量：206.28 Da ✓
- LogP：3.97 (良好的亲脂性) ✓
- TPSA：37.3 Å² (良好的口服吸收) ✓
- 可旋转键：4 ✓

**优化建议**：
1. 降低 LogP：可以引入极性基团
2. 提高选择性：考虑加入氟原子
3. 改善半衰期：考虑氘代修饰

[自动调用 cc-optimize]
生成了 5 个优化的类似物...
```

## 高级功能

### 1. 自定义工作流

```python
# .cursor/workflows/lead_optimization.yaml
name: Lead Optimization Workflow
steps:
  - command: cc-analyze
    params: 
      - properties
      - druglike
  - command: cc-optimize
    params:
      - target: "qed>0.7"
      - maintain_scaffold: true
  - command: cc-synthesize
    params:
      - max_steps: 5
      - green_chemistry: true
```

### 2. 上下文感知

ChemAgent 会根据项目类型自动调整：

- **药物发现项目**：自动关注 ADMET、毒性
- **材料科学项目**：关注聚合物性质、热稳定性
- **有机合成项目**：关注反应条件、产率

### 3. 集成外部工具

```markdown
# 与 AutoDock Vina 集成
cc-dock ligand.mol2 \
  --protein 4LDE.pdb \
  --exhaustiveness 32 \
  --num_poses 10

# 与量子化学软件集成
cc-quantum optimize.xyz \
  --method "B3LYP/6-31G*" \
  --solvent water
```

## 最佳实践

1. **始终验证结构**：使用 `cc-analyze` 验证输入的 SMILES
2. **批量操作**：对大量分子使用批处理命令
3. **保存会话**：重要的分析结果会自动保存
4. **利用子代理**：复杂问题调用专业子代理
5. **组合命令**：使用管道组合多个命令

## 故障排除

### MCP 服务器未启动
```bash
# 手动启动 MCP 服务器
python -m chemagent.mcp_tools.server
```

### 工具未找到
```markdown
# 重新加载工具
cc-reload-tools
```

### 性能优化
```markdown
# 使用缓存
cc-config --enable-cache --cache-ttl 3600
```
