# ChemAgent 后端工具和软件说明

## 🧪 分子分析功能的实现

### 完整属性计算
**实际使用的软件**：
- **RDKit** (开源，Python) - 主要工具
  - 分子量、LogP、TPSA、HBD/HBA
  - 芳香环数、可旋转键数
  - 摩尔折射率、拓扑极性表面积
- **Mordred** (开源，Python) - 扩展描述符
  - 1800+ 分子描述符
- **OpenBabel** (开源) - 格式转换和基础属性

**实现代码示例**：
```python
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

mol = Chem.MolFromSmiles(smiles)
properties = {
    "MW": Descriptors.MolWt(mol),
    "LogP": Descriptors.MolLogP(mol),
    "TPSA": Descriptors.TPSA(mol)
}
```

### 药物相似性评估
**实际使用的软件**：
- **RDKit** - Lipinski's Rule of Five
- **QED (Quantitative Estimate of Drug-likeness)** - Bickerton et al. 算法
- **PAINS (Pan-Assay Interference Compounds)** - 通过 SMARTS 模式匹配

**需要集成的商业/学术工具**：
- **SwissADME** (免费 web 服务) - 更准确的药物相似性
- **DataWarrior** (开源) - 药物相似性评分

### 毒性预测
**当前实现**：
- **结构警告** - 基于 SMARTS 的毒性片段检测
- **PAINS 过滤** - 问题化合物识别

**可以集成的专业工具**：
- **ProTox-II** (免费 web 服务) - 急性毒性、器官毒性
- **ToxTree** (开源，JRC) - 决策树毒性预测
- **VEGA** (开源) - QSAR 毒性模型
- **DeepTox** - 深度学习毒性预测

## 🔬 反应工具的实现

### 反应预测
**可用的开源工具**：
- **RDChiral** - 手性反应预测
- **ASKCOS** (MIT, 免费 API) - 基于 AI 的反应预测
- **RXNMapper** (IBM) - 反应原子映射
- **Molecular Transformer** - 基于 SMILES 的反应预测

**商业工具**（需要许可）：
- **Synthia/Chematica** (Merck) - 专业反应预测
- **Reaxys** (Elsevier) - 反应数据库和预测

### 合成规划
**开源方案**：
- **AiZynthFinder** (AstraZeneca, 开源) - 逆合成规划
- **ASKCOS** - 多步合成路线规划
- **RetroPath2.0** - 代谢和合成路径

**实现示例**：
```python
# 使用 AiZynthFinder
from aizynthfinder import AiZynthFinder

finder = AiZynthFinder()
finder.target_smiles = target_molecule
finder.expansion_policy = "uspto"
finder.run()
routes = finder.routes
```

### 逆合成分析
**主要工具**：
- **AiZynthFinder** (开源，AstraZeneca)
- **Retrosynthetic Planning** (开源模板)
- **Computer-Aided Synthesis Planning (CASP)** 工具集

### 反应机理生成
**当前方案**：
- **基于规则的机理生成** - 使用反应模板
- **Arrow Pushing** - 电子流动可视化

**可集成的工具**：
- **ReactionMechanismGenerator (RMG)** - 自动机理生成
- **ORCA/Gaussian** - 过渡态计算（需要量子化学）

## 💊 药物设计工具的实现

### ADMET 预测
**开源工具**：
- **ADMETlab 2.0** (免费 web 服务，有 Python API)
- **pkCSM** (免费 web 服务) - ADMET 预测
- **SwissADME** (免费) - ADME 参数

**实现方案**：
```python
# 可以通过 API 调用
import requests

def predict_admet(smiles):
    # ADMETlab API (示例)
    response = requests.post(
        "http://admetlab.scbdd.com/api/predict",
        json={"smiles": smiles}
    )
    return response.json()
```

### 分子对接
**开源工具**：
- **AutoDock Vina** (开源) - 最流行的对接软件
- **Smina** (AutoDock Vina 分支) - 更快速
- **rDock** (开源) - 高通量虚拟筛选
- **GNINA** - 基于 CNN 的对接打分

**Python 集成**：
```python
# 使用 Meeko + AutoDock Vina
from meeko import MoleculePreparation
from vina import Vina

v = Vina(sf_name='vina')
v.set_receptor('protein.pdbqt')
v.set_ligand_from_string(ligand_pdbqt)
v.compute_vina_maps(center=[x, y, z], box_size=[20, 20, 20])
v.dock(exhaustiveness=32)
```

### 先导优化
**开源方案**：
- **REINVENT** (AstraZeneca, 开源) - 基于 RL 的分子生成
- **GuacaMol** (BenevolentAI) - 分子优化基准
- **STONED** - 分子突变算法
- **Molecule.one** - 逆合成可行性评分

### 虚拟筛选
**工具链**：
- **RDKit** - 药效团筛选、子结构搜索
- **Open Drug Discovery Toolkit (ODDT)** - 完整筛选流程
- **PyRx** - AutoDock 的 GUI 包装
- **Pharmit** - 药效团搜索

## 📚 文献搜索能力集成方案

### 1. PubMed/PMC 集成
```python
from Bio import Entrez
import requests

class LiteratureSearch:
    def __init__(self):
        Entrez.email = "your-email@example.com"
    
    def search_pubmed(self, query, max_results=10):
        """搜索 PubMed 文献"""
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results
        )
        record = Entrez.read(handle)
        ids = record["IdList"]
        
        # 获取摘要
        handle = Entrez.efetch(
            db="pubmed",
            id=ids,
            rettype="abstract",
            retmode="text"
        )
        return handle.read()
```

### 2. ChemRxiv/bioRxiv 预印本
```python
def search_chemrxiv(query):
    """搜索 ChemRxiv 预印本"""
    api_url = "https://chemrxiv.org/engage/chemrxiv/public-api/v1/items"
    params = {"term": query}
    response = requests.get(api_url, params=params)
    return response.json()
```

### 3. Google Scholar (通过 scholarly)
```python
from scholarly import scholarly

def search_google_scholar(query):
    """搜索 Google Scholar"""
    search_query = scholarly.search_pubs(query)
    papers = []
    for i in range(10):  # 获取前10篇
        paper = next(search_query)
        papers.append({
            'title': paper['bib']['title'],
            'author': paper['bib']['author'],
            'year': paper['bib'].get('pub_year'),
            'citations': paper.get('num_citations')
        })
    return papers
```

### 4. 专利文献搜索
```python
class PatentSearch:
    def search_uspto(self, query):
        """搜索美国专利"""
        # 使用 USPTO API
        pass
    
    def search_espacenet(self, query):
        """搜索欧洲专利"""
        # 使用 EPO OPS API
        pass
```

### 5. 化学文献数据库
- **Reaxys** (需要订阅) - 化学反应和物质数据
- **SciFinder** (CAS, 需要订阅) - 化学文献和反应
- **Web of Science** (需要订阅) - 引文数据库

## 🔧 实际集成建议

### 核心依赖（必需）
```bash
pip install rdkit-pypi
pip install pubchempy
pip install biopython  # for PubMed
pip install requests
```

### 扩展依赖（推荐）
```bash
pip install aizynthfinder  # 逆合成
pip install meeko vina     # 分子对接
pip install mordred        # 描述符
pip install scholarly       # Google Scholar
```

### MCP 工具定义示例
```python
# 添加文献搜索工具到 MCP
{
    "name": "search_literature",
    "description": "Search scientific literature",
    "inputSchema": {
        "type": "object",
        "properties": {
            "query": {"type": "string"},
            "database": {
                "type": "string",
                "enum": ["pubmed", "chemrxiv", "scholar", "all"]
            },
            "limit": {"type": "integer", "default": 10},
            "filters": {
                "type": "object",
                "properties": {
                    "year_from": {"type": "integer"},
                    "year_to": {"type": "integer"},
                    "journal": {"type": "string"}
                }
            }
        }
    }
}
```

## 📊 工具成熟度评估

| 功能类别 | 开源方案成熟度 | 商业方案必要性 | 推荐实现 |
|---------|--------------|---------------|---------|
| 基础属性计算 | ⭐⭐⭐⭐⭐ 非常成熟 | 低 | RDKit |
| ADMET 预测 | ⭐⭐⭐⭐ 成熟 | 中 | ADMETlab API |
| 反应预测 | ⭐⭐⭐ 发展中 | 高 | ASKCOS + 规则 |
| 逆合成 | ⭐⭐⭐⭐ 成熟 | 中 | AiZynthFinder |
| 分子对接 | ⭐⭐⭐⭐⭐ 非常成熟 | 低 | AutoDock Vina |
| 文献搜索 | ⭐⭐⭐⭐ 成熟 | 中 | PubMed + Scholar |
| 毒性预测 | ⭐⭐⭐ 发展中 | 高 | ProTox-II API |

## 💡 实施优先级建议

### 第一阶段（立即可实现）
1. ✅ RDKit 全功能集成
2. ✅ PubChem API 集成
3. ✅ 基础 ADMET (通过 RDKit)
4. ➕ PubMed 文献搜索
5. ➕ Google Scholar 集成

### 第二阶段（需要额外开发）
1. ➕ ADMETlab 2.0 API 集成
2. ➕ AiZynthFinder 逆合成
3. ➕ AutoDock Vina 对接
4. ➕ ChemRxiv 预印本搜索
5. ➕ ProTox-II 毒性预测

### 第三阶段（需要商业许可或大量开发）
1. ⏳ ASKCOS 完整集成
2. ⏳ 专利数据库接入
3. ⏳ Reaxys/SciFinder（如有订阅）
4. ⏳ 量子化学计算（Gaussian/ORCA）
