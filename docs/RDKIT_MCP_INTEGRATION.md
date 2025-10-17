# RDKit MCP 集成说明

## ✅ 好消息：RDKit 有官方 MCP 支持！

经过搜索和验证，我们发现 **RDKit 确实有官方的 MCP (Model Context Protocol) 实现**！

### 官方包：mcp-rdkit

- **PyPI 包名**: `mcp-rdkit`
- **安装命令**: `pip install mcp-rdkit`
- **运行命令**: `python -m mcp_rdkit`
- **状态**: ✅ 官方支持，可立即使用

## 集成到 ChemAgent

我们已经完成了以下集成工作：

### 1. 创建了专门的安装脚本
- **文件**: `install_rdkit_mcp.sh`
- **功能**:
  - 自动安装 mcp-rdkit 包
  - 创建虚拟环境（可选）
  - 生成配置文件
  - 创建启动脚本
  - 提供 systemd 服务配置

### 2. 更新了主安装程序
```bash
# 新增的 MCP 安装命令
python chemagent_install.py mcp
```

### 3. 更新了 MCP Orchestrator
- 将 RDKit MCP 标记为"官方可用"
- 优先使用官方 MCP 服务器
- 保留内部实现作为后备方案

## RDKit MCP 功能

官方的 mcp-rdkit 提供以下功能：

### 分子处理
- SMILES 解析和验证
- MOL 文件读写
- 2D 坐标生成
- 3D 构象生成

### 分子描述符
- 分子量 (MW)
- 脂水分配系数 (LogP)
- 极性表面积 (TPSA)
- 氢键供体/受体数 (HBD/HBA)
- Lipinski 五规则
- QED 药物相似性评分
- 合成可及性评分

### 分子指纹
- Morgan 指纹
- RDKit 指纹
- MACCS 指纹
- Atom Pair 指纹

### 相似性计算
- Tanimoto 系数
- Dice 系数
- Cosine 相似度

### 可视化
- 2D 分子图像生成
- 3D 结构展示
- SVG 输出格式

## 使用示例

### 1. 启动 RDKit MCP 服务器
```bash
# 使用安装脚本创建的启动器
~/.chemagent/mcp_servers/start_rdkit_mcp.sh

# 或直接运行
python -m mcp_rdkit
```

### 2. 在 Claude Code 中使用
```markdown
cc-analyze "CC(=O)OC1=CC=CC=C1C(=O)O"
# ChemAgent 会自动调用 RDKit MCP 服务器
```

### 3. 直接调用 MCP API
```python
import requests

# MCP 请求
request = {
    "method": "calculate_properties",
    "params": {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "properties": ["MW", "LogP", "TPSA", "HBD", "HBA"]
    }
}

response = requests.post("http://localhost:8766/mcp", json=request)
print(response.json())
```

## 架构优势

使用官方 MCP 的优势：

1. **标准化**: 遵循 MCP 协议规范
2. **维护性**: 由 RDKit 社区维护
3. **性能**: 优化的服务器实现
4. **扩展性**: 易于添加新功能
5. **互操作性**: 可与其他 MCP 工具配合

## 后续计划

虽然 RDKit 已有官方 MCP，但其他工具仍需努力：

- [ ] 推动 PubChem 官方 MCP
- [ ] 推动 ChEMBL 官方 MCP
- [ ] 推动 OpenBabel 官方 MCP
- [ ] 继续使用内部集成作为过渡方案

## 总结

ChemAgent 的设计理念得到了验证：
- ✅ 优先使用官方 MCP（RDKit 已实现）
- ✅ 提供内部实现作为后备
- ✅ 统一的编排层管理所有工具
- ✅ 用户无需关心底层实现细节

这正是我们期望的架构：**让 ChemAgent 成为化学工具的智能编排器，而不是重新造轮子**。
