# ChemAgent for Gemini CLI 使用示例

## 安装后的使用

### 1. 基本化学命令

```bash
# 分析分子
gemini chem:analyze "CC(=O)OC1=CC=CC=C1C(=O)O"

# 输出：
Analyzing molecule...
Structure: Aspirin (Acetylsalicylic acid)
Properties:
  - MW: 180.16
  - LogP: 1.19
  - Drug-like: Yes ✓
```

### 2. 批量处理（Gemini 特色）

```bash
# 处理 CSV 文件中的多个分子
gemini chem:batch molecules.csv \
  --operation analyze \
  --output results.csv

# 处理 100 个分子仅需几秒
Processing 100 molecules...
[████████████████████] 100/100 complete
Results saved to results.csv
```

### 3. 多模态分析（Gemini 独特功能）

```bash
# 从图片识别分子结构
gemini chem:analyze molecule_structure.png \
  --image-to-smiles \
  --properties

# Gemini 的视觉能力可以：
# - 识别手绘结构
# - 解析文献中的分子图
# - 分析反应机理图
```

### 4. 合成路线规划

```bash
# 设计合成路线
gemini chem:synthesize "Ibuprofen" \
  --max-steps 5 \
  --green-chemistry \
  --cost-optimize

# 输出：
Retrosynthetic Analysis:
Step 1: Friedel-Crafts acylation
Step 2: Wolff-Kishner reduction
Step 3: Carboxylation
Estimated yield: 65%
Green score: 8/10
```

### 5. 使用别名（配置后）

```bash
# Source 别名文件
source ~/.gemini/chemagent_aliases.sh

# 使用简化命令
chem aspirin              # 快速分析
synthesize paracetamol    # 合成规划
chembatch drugs.sdf       # 批处理
```

## 高级功能

### 1. 云函数集成

```bash
# 部署为 Google Cloud Function
gemini deploy chem:analyze \
  --project my-project \
  --region us-central1

# 远程调用
gemini cloud:chem:analyze "CCO" \
  --async \
  --webhook https://myapp.com/callback
```

### 2. 流式处理

```bash
# 实时处理数据流
cat molecules.txt | gemini chem:stream \
  --operation filter \
  --criteria "mw<500 AND logp<5" \
  | tee filtered.txt
```

### 3. 与 Jupyter 集成

```python
# 在 Jupyter 中使用
!gemini chem:analyze "c1ccccc1" --format json > result.json

import json
with open('result.json') as f:
    data = json.load(f)
    print(f"LogP: {data['properties']['logp']}")
```

### 4. 并行处理

```bash
# 利用 Gemini 的并行能力
gemini chem:parallel \
  --input large_library.sdf \
  --operation "analyze,optimize,dock" \
  --workers 10 \
  --output results/

# 处理 10,000 个分子
Time: 2 minutes (vs 30 minutes sequential)
```

## 特色工作流

### 1. 文献挖掘 + 分子分析

```bash
# 从文献中提取分子并分析
gemini chem:literature-mining \
  --query "COVID-19 inhibitors" \
  --extract-molecules \
  --analyze-properties \
  --output covid_drugs.csv
```

### 2. 图像批处理

```bash
# 处理文件夹中的所有分子图片
gemini chem:image-batch \
  --input-dir ./structures/ \
  --recognize-structure \
  --calculate-properties \
  --output-format sdf
```

### 3. AI 辅助优化

```bash
# 使用 Gemini 的 AI 能力优化分子
gemini chem:ai-optimize \
  --input lead_compound.mol \
  --target "improve solubility" \
  --constraints "maintain activity" \
  --suggestions 10
```

## 配置文件示例

```yaml
# ~/.gemini/chemagent_config.yaml
chemistry:
  default_format: smiles
  auto_validate: true
  cache_results: true
  
tools:
  rdkit:
    enabled: true
    compute_3d: false
  
  pubchem:
    enabled: true
    cache_ttl: 3600
    
batch:
  chunk_size: 100
  parallel_workers: 4
  
visualization:
  default_format: png
  dpi: 300
```

## 性能优化技巧

### 1. 使用缓存
```bash
gemini chem:config --enable-cache --cache-dir ~/.gemini/cache
```

### 2. 批量操作
```bash
# 不要这样（慢）
for mol in $(cat molecules.txt); do
  gemini chem:analyze "$mol"
done

# 这样做（快）
gemini chem:batch molecules.txt --operation analyze
```

### 3. 并行处理
```bash
gemini chem:parallel --workers $(nproc) --input big_file.sdf
```

## 故障排除

### MCP 服务器连接
```bash
# 检查 MCP 服务器状态
gemini chem:status

# 手动指定服务器
gemini chem:analyze "CCO" --server ws://localhost:8765
```

### 依赖问题
```bash
# 检查依赖
gemini chem:check-deps

# 重新安装
gemini chem:reinstall
```

## 与其他工具集成

### 1. 与 RDKit 直接交互
```bash
gemini chem:rdkit --command "MolFromSmiles('CCO')"
```

### 2. 导出到其他格式
```bash
gemini chem:analyze "aspirin" --export pymol > aspirin.pml
gemini chem:analyze "aspirin" --export gaussian > aspirin.gjf
```

### 3. 管道组合
```bash
# 复杂管道
gemini chem:search "antibiotic" |
gemini chem:filter "mw<600" |
gemini chem:optimize "logp" |
gemini chem:synthesize --check-availability
```
