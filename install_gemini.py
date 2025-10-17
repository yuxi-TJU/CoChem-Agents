#!/usr/bin/env python3
"""
ChemAgent Installer for Gemini CLI
专门用于将ChemAgent安装到Gemini CLI的安装脚本
"""

import os
import sys
import json
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Any, Optional, List
import argparse

# 尝试导入rich用于更好的输出
try:
    from rich import print
    from rich.console import Console
    from rich.table import Table
    from rich.progress import Progress, SpinnerColumn, TextColumn
    from rich.prompt import Prompt, Confirm
    RICH_AVAILABLE = True
    console = Console()
except ImportError:
    RICH_AVAILABLE = False
    console = None
    print("[警告] rich库未安装，使用基础输出模式")


class GeminiCLIInstaller:
    """Gemini CLI专用安装器"""
    
    def __init__(self, quiet: bool = False, yes: bool = False):
        self.quiet = quiet
        self.yes = yes
        self.gemini_home = Path.home() / ".gemini"
        self.extensions_dir = self.gemini_home / "extensions"
        self.chemagent_dir = self.extensions_dir / "chemagent"
        self.config_file = self.gemini_home / "config.json"
        
    def check_gemini_cli(self) -> bool:
        """检查Gemini CLI是否已安装"""
        # 检查gemini命令
        result = subprocess.run(
            ["which", "gemini"],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            # 检查常见安装位置
            common_paths = [
                Path.home() / ".local" / "bin" / "gemini",
                Path("/usr/local/bin/gemini"),
                Path("/opt/gemini/bin/gemini")
            ]
            
            for path in common_paths:
                if path.exists():
                    return True
            
            return False
        
        return True
    
    def detect_gemini_version(self) -> Optional[str]:
        """检测Gemini CLI版本"""
        try:
            result = subprocess.run(
                ["gemini", "--version"],
                capture_output=True,
                text=True
            )
            if result.returncode == 0:
                # 解析版本信息
                version_line = result.stdout.strip()
                if "gemini" in version_line.lower():
                    return version_line
        except:
            pass
        
        return None
    
    def create_extension_structure(self):
        """创建Gemini扩展目录结构"""
        if not self.quiet:
            print("\n📁 创建Gemini扩展目录结构...")
        
        # 创建主目录
        self.chemagent_dir.mkdir(parents=True, exist_ok=True)
        
        # 创建子目录
        subdirs = ["functions", "prompts", "models", "data", "cache", "logs"]
        for subdir in subdirs:
            (self.chemagent_dir / subdir).mkdir(exist_ok=True)
        
        if not self.quiet:
            print("✅ 扩展目录结构创建完成")
    
    def install_manifest(self):
        """安装扩展清单文件"""
        if not self.quiet:
            print("\n📝 创建扩展清单...")
        
        manifest = {
            "name": "ChemAgent",
            "version": "1.0.0",
            "description": "Chemistry Enhancement Package for Gemini CLI",
            "description_zh": "Gemini CLI化学增强包",
            "author": "ChemAgent Team",
            "license": "MIT",
            "homepage": "https://github.com/yourusername/chemagent",
            "entry_point": "main.py",
            "requirements": [
                "rdkit>=2023.0.0",
                "pubchempy>=1.0.4",
                "chembl_webresource_client>=0.10.8"
            ],
            "commands": {
                "chem:analyze": "分析分子结构和性质",
                "chem:synthesize": "设计合成路线",
                "chem:predict": "预测反应和性质",
                "chem:optimize": "优化分子结构",
                "chem:search": "搜索化学数据库",
                "chem:batch": "批量处理分子",
                "chem:visualize": "可视化分子结构",
                "chem:dock": "分子对接",
                "chem:safety": "安全性评估",
                "chem:patent": "专利检索"
            },
            "capabilities": [
                "molecular_analysis",
                "synthesis_planning",
                "reaction_prediction",
                "property_optimization",
                "database_search",
                "batch_processing",
                "structure_visualization",
                "image_to_structure",
                "multimodal_input",
                "cloud_computing"
            ],
            "settings": {
                "auto_load": True,
                "cache_enabled": True,
                "default_model": "gemini-pro",
                "vision_model": "gemini-pro-vision",
                "max_batch_size": 100,
                "timeout": 30
            }
        }
        
        manifest_file = self.chemagent_dir / "manifest.json"
        with open(manifest_file, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2, ensure_ascii=False)
        
        if not self.quiet:
            print("✅ 扩展清单已创建")
    
    def install_main_script(self):
        """安装主入口脚本"""
        if not self.quiet:
            print("\n🔧 安装主入口脚本...")
        
        main_content = '''#!/usr/bin/env python3
"""
ChemAgent Main Entry Point for Gemini CLI
Gemini CLI的ChemAgent主入口
"""

import sys
import json
import os
from pathlib import Path
from typing import Dict, Any, List, Optional

# 添加ChemAgent到Python路径
chemagent_root = Path(__file__).parent.parent.parent.parent
if chemagent_root.exists():
    sys.path.insert(0, str(chemagent_root))

try:
    from chemagent.mcp_tools.orchestrator import MCPOrchestrator
    from chemagent.commands.loader import get_command_prompt
    CHEMAGENT_AVAILABLE = True
except ImportError:
    CHEMAGENT_AVAILABLE = False
    print("警告: ChemAgent核心模块未找到，部分功能可能不可用", file=sys.stderr)


class GeminiChemHandler:
    """Gemini CLI化学命令处理器"""
    
    def __init__(self):
        self.orchestrator = MCPOrchestrator() if CHEMAGENT_AVAILABLE else None
        self.commands = self._load_commands()
    
    def _load_commands(self) -> Dict[str, Any]:
        """加载可用命令"""
        commands = {
            "analyze": self.analyze_molecule,
            "synthesize": self.plan_synthesis,
            "predict": self.predict_reaction,
            "optimize": self.optimize_molecule,
            "search": self.search_database,
            "batch": self.batch_process,
            "visualize": self.visualize_structure,
            "dock": self.molecular_docking,
            "safety": self.assess_safety,
            "patent": self.check_patent
        }
        return commands
    
    def handle_command(self, command: str, args: List[str]) -> Dict[str, Any]:
        """处理化学命令"""
        # 移除chem:前缀
        if command.startswith("chem:"):
            command = command[5:]
        
        if command in self.commands:
            return self.commands[command](args)
        else:
            return {
                "error": f"未知命令: {command}",
                "available_commands": list(self.commands.keys())
            }
    
    def analyze_molecule(self, args: List[str]) -> Dict[str, Any]:
        """分析分子"""
        if not args:
            return {"error": "请提供分子SMILES或名称"}
        
        molecule = args[0]
        options = args[1:] if len(args) > 1 else []
        
        if self.orchestrator:
            result = self.orchestrator.execute_command(
                "cc-analyze",
                {"molecule": molecule, "options": options}
            )
            return result
        else:
            return {
                "molecule": molecule,
                "message": "ChemAgent核心未加载，返回基础分析",
                "properties": {
                    "input": molecule,
                    "type": "SMILES" if "C" in molecule else "name"
                }
            }
    
    def plan_synthesis(self, args: List[str]) -> Dict[str, Any]:
        """规划合成路线"""
        if not args:
            return {"error": "请提供目标分子"}
        
        target = args[0]
        
        if self.orchestrator:
            result = self.orchestrator.execute_command(
                "cc-synthesize",
                {"target": target}
            )
            return result
        else:
            return {
                "target": target,
                "message": "合成路线规划需要ChemAgent核心模块"
            }
    
    def predict_reaction(self, args: List[str]) -> Dict[str, Any]:
        """预测反应"""
        if not args:
            return {"error": "请提供反应物"}
        
        if self.orchestrator:
            result = self.orchestrator.execute_command(
                "cc-predict",
                {"reactants": args}
            )
            return result
        else:
            return {
                "reactants": args,
                "message": "反应预测需要ChemAgent核心模块"
            }
    
    def optimize_molecule(self, args: List[str]) -> Dict[str, Any]:
        """优化分子"""
        if not args:
            return {"error": "请提供分子结构"}
        
        molecule = args[0]
        target_property = args[1] if len(args) > 1 else "druglike"
        
        if self.orchestrator:
            result = self.orchestrator.execute_command(
                "cc-optimize",
                {"molecule": molecule, "target": target_property}
            )
            return result
        else:
            return {
                "molecule": molecule,
                "target": target_property,
                "message": "分子优化需要ChemAgent核心模块"
            }
    
    def search_database(self, args: List[str]) -> Dict[str, Any]:
        """搜索数据库"""
        if not args:
            return {"error": "请提供搜索关键词"}
        
        query = " ".join(args)
        
        if self.orchestrator:
            result = self.orchestrator.execute_command(
                "cc-search",
                {"query": query}
            )
            return result
        else:
            return {
                "query": query,
                "message": "数据库搜索需要ChemAgent核心模块"
            }
    
    def batch_process(self, args: List[str]) -> Dict[str, Any]:
        """批量处理"""
        if not args:
            return {"error": "请提供输入文件"}
        
        input_file = args[0]
        operation = args[1] if len(args) > 1 else "analyze"
        
        if self.orchestrator:
            result = self.orchestrator.execute_command(
                "cc-batch",
                {"input": input_file, "operation": operation}
            )
            return result
        else:
            return {
                "input": input_file,
                "operation": operation,
                "message": "批量处理需要ChemAgent核心模块"
            }
    
    def visualize_structure(self, args: List[str]) -> Dict[str, Any]:
        """可视化结构"""
        if not args:
            return {"error": "请提供分子结构"}
        
        molecule = args[0]
        style = args[1] if len(args) > 1 else "2D"
        
        return {
            "molecule": molecule,
            "style": style,
            "message": "可视化功能将在Gemini界面中显示"
        }
    
    def molecular_docking(self, args: List[str]) -> Dict[str, Any]:
        """分子对接"""
        if len(args) < 2:
            return {"error": "请提供配体和受体"}
        
        ligand = args[0]
        receptor = args[1]
        
        return {
            "ligand": ligand,
            "receptor": receptor,
            "message": "分子对接计算需要专业软件支持"
        }
    
    def assess_safety(self, args: List[str]) -> Dict[str, Any]:
        """安全性评估"""
        if not args:
            return {"error": "请提供分子结构"}
        
        molecule = args[0]
        
        if self.orchestrator:
            result = self.orchestrator._assess_safety(molecule)
            return result
        else:
            return {
                "molecule": molecule,
                "message": "安全性评估需要ChemAgent核心模块"
            }
    
    def check_patent(self, args: List[str]) -> Dict[str, Any]:
        """专利检查"""
        if not args:
            return {"error": "请提供分子结构或CAS号"}
        
        query = args[0]
        
        if self.orchestrator:
            result = self.orchestrator._check_patents(query)
            return result
        else:
            return {
                "query": query,
                "message": "专利检索需要ChemAgent核心模块"
            }


def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("用法: gemini chem:<command> [args]")
        print("可用命令:")
        print("  chem:analyze <SMILES>  - 分析分子")
        print("  chem:synthesize <target> - 规划合成")
        print("  chem:predict <reactants> - 预测反应")
        print("  chem:optimize <molecule> - 优化分子")
        print("  chem:search <query> - 搜索数据库")
        print("  chem:batch <file> - 批量处理")
        print("  chem:visualize <molecule> - 可视化")
        print("  chem:dock <ligand> <receptor> - 分子对接")
        print("  chem:safety <molecule> - 安全评估")
        print("  chem:patent <query> - 专利检查")
        return 0
    
    command = sys.argv[1]
    args = sys.argv[2:] if len(sys.argv) > 2 else []
    
    handler = GeminiChemHandler()
    
    try:
        result = handler.handle_command(command, args)
        
        # 格式化输出
        if "error" in result:
            print(f"错误: {result['error']}", file=sys.stderr)
            return 1
        else:
            print(json.dumps(result, indent=2, ensure_ascii=False))
            return 0
            
    except Exception as e:
        print(f"执行错误: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
'''
        
        main_file = self.chemagent_dir / "main.py"
        with open(main_file, "w", encoding="utf-8") as f:
            f.write(main_content)
        
        # 设置可执行权限
        os.chmod(main_file, 0o755)
        
        if not self.quiet:
            print("✅ 主入口脚本已安装")
    
    def install_prompts(self):
        """安装提示词模板"""
        if not self.quiet:
            print("\n📋 安装提示词模板...")
        
        prompts_dir = self.chemagent_dir / "prompts"
        
        prompts = {
            "analyze.txt": """分析以下分子结构：
{molecule}

请提供：
1. 基本性质（分子量、LogP、TPSA、HBA/HBD）
2. 药物相似性评估（Lipinski规则、QED分数）
3. ADMET预测
4. 毒性评估
5. 合成可及性
6. 可能的生物活性
7. 结构优化建议""",

            "synthesize.txt": """为以下目标分子设计合成路线：
目标: {target}

要求：
1. 逆合成分析（显示断键策略）
2. 前向合成路线（详细步骤）
3. 每步反应条件（温度、溶剂、催化剂）
4. 预期产率
5. 关键中间体的纯化方法
6. 替代路线（至少2条）
7. 成本和时间估算""",

            "predict.txt": """预测以下反应：
反应物: {reactants}
条件: {conditions}

分析内容：
1. 主要产物结构
2. 反应机理（逐步）
3. 副产物可能性
4. 产率估算
5. 反应选择性
6. 优化建议""",

            "optimize.txt": """优化以下分子以提高{property}：
起始分子: {molecule}

优化策略：
1. 识别可修饰位点
2. 提出修饰方案（至少5个）
3. 预测每个修饰的效果
4. 保持其他重要性质
5. 合成可行性评估
6. 推荐最优方案""",

            "batch.txt": """批量处理分子列表：
操作: {operation}
分子数: {count}

对每个分子执行：
1. {operation}分析
2. 结果汇总
3. 统计分析
4. 异常标记
5. 导出报告"""
        }
        
        for filename, content in prompts.items():
            prompt_file = prompts_dir / filename
            with open(prompt_file, "w", encoding="utf-8") as f:
                f.write(content)
        
        if not self.quiet:
            print(f"✅ 安装了{len(prompts)}个提示词模板")
    
    def configure_gemini_cli(self):
        """配置Gemini CLI"""
        if not self.quiet:
            print("\n⚙️  配置Gemini CLI...")
        
        # 确保配置目录存在
        self.gemini_home.mkdir(exist_ok=True)
        
        # 读取或创建配置
        config = {}
        if self.config_file.exists():
            try:
                with open(self.config_file, "r", encoding="utf-8") as f:
                    config = json.load(f)
            except:
                config = {}
        
        # 添加ChemAgent配置
        if "extensions" not in config:
            config["extensions"] = {}
        
        config["extensions"]["chemagent"] = {
            "enabled": True,
            "version": "1.0.0",
            "auto_load": True,
            "path": str(self.chemagent_dir),
            "entry_point": "main.py"
        }
        
        # 添加化学相关设置
        config["chemistry"] = {
            "default_format": "SMILES",
            "auto_validate": True,
            "safety_checks": True,
            "cache_results": True,
            "visualization": {
                "default_style": "2D",
                "color_scheme": "default",
                "show_hydrogens": False
            },
            "computation": {
                "max_atoms": 500,
                "timeout": 30,
                "parallel": True
            }
        }
        
        # 注册命令
        if "commands" not in config:
            config["commands"] = {}
        
        chem_commands = [
            "chem:analyze", "chem:synthesize", "chem:predict",
            "chem:optimize", "chem:search", "chem:batch",
            "chem:visualize", "chem:dock", "chem:safety", "chem:patent"
        ]
        
        for cmd in chem_commands:
            config["commands"][cmd] = {
                "extension": "chemagent",
                "handler": "main.py",
                "description": f"Chemistry command: {cmd.split(':')[1]}"
            }
        
        # 保存配置
        with open(self.config_file, "w", encoding="utf-8") as f:
            json.dump(config, f, indent=2, ensure_ascii=False)
        
        if not self.quiet:
            print("✅ Gemini CLI配置完成")
    
    def create_shell_integration(self):
        """创建Shell集成"""
        if not self.quiet:
            print("\n🐚 创建Shell集成...")
        
        # 创建别名文件
        aliases_content = """#!/bin/bash
# ChemAgent for Gemini CLI - Shell Integration
# Gemini CLI化学增强包 - Shell集成

# 基础别名
alias gchem='gemini chem:analyze'
alias gsynth='gemini chem:synthesize'
alias gpredict='gemini chem:predict'
alias gopt='gemini chem:optimize'
alias gsearch='gemini chem:search'
alias gbatch='gemini chem:batch'
alias gviz='gemini chem:visualize'
alias gdock='gemini chem:dock'
alias gsafety='gemini chem:safety'
alias gpatent='gemini chem:patent'

# 便捷函数

# 快速分子分析
chem() {
    if [ -z "$1" ]; then
        echo "用法: chem <SMILES或分子名称>"
        return 1
    fi
    gemini chem:analyze "$1" --full
}

# 合成路线设计
synthesize() {
    if [ -z "$1" ]; then
        echo "用法: synthesize <目标分子>"
        return 1
    fi
    gemini chem:synthesize "$1" --detailed
}

# 批量分析CSV文件
chembatch() {
    if [ -z "$1" ]; then
        echo "用法: chembatch <CSV文件>"
        return 1
    fi
    gemini chem:batch "$1" --operation analyze --output "${1%.csv}_results.csv"
}

# 分子优化
optimize() {
    if [ -z "$1" ]; then
        echo "用法: optimize <分子> [性质]"
        return 1
    fi
    local property="${2:-druglike}"
    gemini chem:optimize "$1" --target "$property"
}

# 安全性快速检查
safety() {
    if [ -z "$1" ]; then
        echo "用法: safety <分子>"
        return 1
    fi
    gemini chem:safety "$1" --comprehensive
}

# 显示帮助信息
chemhelp() {
    echo "ChemAgent for Gemini CLI - 可用命令："
    echo ""
    echo "基础命令："
    echo "  chem <molecule>      - 快速分子分析"
    echo "  synthesize <target>  - 合成路线设计"
    echo "  optimize <mol> [prop] - 分子优化"
    echo "  safety <molecule>    - 安全性评估"
    echo "  chembatch <file>     - 批量处理"
    echo ""
    echo "Gemini命令："
    echo "  gemini chem:analyze  - 详细分子分析"
    echo "  gemini chem:synthesize - 合成规划"
    echo "  gemini chem:predict  - 反应预测"
    echo "  gemini chem:search   - 数据库搜索"
    echo "  gemini chem:dock     - 分子对接"
    echo "  gemini chem:patent   - 专利检查"
    echo ""
    echo "别名："
    echo "  gchem, gsynth, gpredict, gopt, gsearch, gbatch, gviz, gdock, gsafety, gpatent"
}

# 启动消息
echo "✅ ChemAgent for Gemini CLI 已加载"
echo "输入 'chemhelp' 查看可用命令"
"""
        
        aliases_file = self.gemini_home / "chemagent_aliases.sh"
        with open(aliases_file, "w", encoding="utf-8") as f:
            f.write(aliases_content)
        
        os.chmod(aliases_file, 0o755)
        
        if not self.quiet:
            print(f"✅ Shell集成已创建: {aliases_file}")
            print(f"   请将以下行添加到您的shell配置文件（~/.bashrc或~/.zshrc）：")
            print(f"   source {aliases_file}")
    
    def install_examples(self):
        """安装示例文件"""
        if not self.quiet:
            print("\n📚 安装示例文件...")
        
        examples_dir = self.chemagent_dir / "examples"
        examples_dir.mkdir(exist_ok=True)
        
        # 示例分子列表
        molecules_csv = """SMILES,Name,Type
CC(=O)OC1=CC=CC=C1C(=O)O,Aspirin,Drug
CC1=CC=C(C=C1)C(C)CC(=O)O,Ibuprofen,Drug
CN1C=NC2=C1C(=O)N(C(=O)N2C)C,Caffeine,Natural
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O,Ibuprofen,Drug
CC1=CC(=O)CC(C1)(C)C,Camphor,Natural
"""
        
        with open(examples_dir / "molecules.csv", "w") as f:
            f.write(molecules_csv)
        
        # 示例脚本
        example_script = """#!/usr/bin/env python3
# ChemAgent Gemini CLI示例脚本

import subprocess
import json

def analyze_molecule(smiles):
    \"\"\"使用Gemini CLI分析分子\"\"\"
    result = subprocess.run(
        ["gemini", "chem:analyze", smiles],
        capture_output=True,
        text=True
    )
    
    if result.returncode == 0:
        return json.loads(result.stdout)
    else:
        return {"error": result.stderr}

# 分析阿司匹林
aspirin = "CC(=O)OC1=CC=CC=C1C(=O)O"
result = analyze_molecule(aspirin)
print("阿司匹林分析结果:")
print(json.dumps(result, indent=2, ensure_ascii=False))
"""
        
        with open(examples_dir / "example.py", "w") as f:
            f.write(example_script)
        
        os.chmod(examples_dir / "example.py", 0o755)
        
        if not self.quiet:
            print("✅ 示例文件已安装")
    
    def verify_installation(self) -> bool:
        """验证安装"""
        if not self.quiet:
            print("\n🔍 验证安装...")
        
        checks = {
            "扩展目录": self.chemagent_dir.exists(),
            "清单文件": (self.chemagent_dir / "manifest.json").exists(),
            "主脚本": (self.chemagent_dir / "main.py").exists(),
            "配置文件": self.config_file.exists(),
            "别名文件": (self.gemini_home / "chemagent_aliases.sh").exists()
        }
        
        all_ok = all(checks.values())
        
        if RICH_AVAILABLE and not self.quiet:
            table = Table(title="安装验证")
            table.add_column("组件", style="cyan")
            table.add_column("状态", style="green")
            
            for component, status in checks.items():
                status_text = "✅ 已安装" if status else "❌ 未找到"
                table.add_row(component, status_text)
            
            console.print(table)
        elif not self.quiet:
            for component, status in checks.items():
                status_text = "✅" if status else "❌"
                print(f"{status_text} {component}")
        
        return all_ok
    
    def test_command(self) -> bool:
        """测试命令是否工作"""
        if not self.quiet:
            print("\n🧪 测试化学命令...")
        
        try:
            # 测试基础命令
            test_cmd = [
                sys.executable,
                str(self.chemagent_dir / "main.py"),
                "chem:analyze",
                "CCO"  # 乙醇
            ]
            
            result = subprocess.run(
                test_cmd,
                capture_output=True,
                text=True,
                timeout=5
            )
            
            if result.returncode == 0:
                if not self.quiet:
                    print("✅ 命令测试成功")
                return True
            else:
                if not self.quiet:
                    print(f"❌ 命令测试失败: {result.stderr}")
                return False
                
        except Exception as e:
            if not self.quiet:
                print(f"❌ 测试出错: {e}")
            return False
    
    def print_summary(self):
        """打印安装摘要"""
        if self.quiet:
            return
        
        if RICH_AVAILABLE:
            console.print("\n[bold green]✅ ChemAgent for Gemini CLI 安装成功！[/bold green]")
            console.print("\n[yellow]快速开始：[/yellow]")
            console.print("1. 将以下行添加到您的shell配置文件：")
            console.print(f"   [cyan]source {self.gemini_home}/chemagent_aliases.sh[/cyan]")
            console.print("\n2. 重新加载shell配置或打开新终端")
            console.print("\n3. 尝试以下命令：")
            console.print("   [cyan]gemini chem:analyze CCO[/cyan]  # 分析乙醇")
            console.print("   [cyan]chem aspirin[/cyan]  # 快速分析阿司匹林")
            console.print("   [cyan]chemhelp[/cyan]  # 查看所有可用命令")
        else:
            print("\n✅ ChemAgent for Gemini CLI 安装成功！")
            print("\n快速开始：")
            print(f"1. 将以下行添加到您的shell配置文件：")
            print(f"   source {self.gemini_home}/chemagent_aliases.sh")
            print("\n2. 重新加载shell配置或打开新终端")
            print("\n3. 尝试以下命令：")
            print("   gemini chem:analyze CCO  # 分析乙醇")
            print("   chem aspirin  # 快速分析阿司匹林")
            print("   chemhelp  # 查看所有可用命令")
    
    def install(self) -> bool:
        """执行完整安装"""
        try:
            # 检查Gemini CLI
            if not self.check_gemini_cli():
                if not self.yes:
                    response = input("\n⚠️  未检测到Gemini CLI，是否继续安装？[y/N]: ")
                    if response.lower() != 'y':
                        print("安装已取消")
                        return False
                print("⚠️  警告：未检测到Gemini CLI，某些功能可能无法使用")
            else:
                version = self.detect_gemini_version()
                if version and not self.quiet:
                    print(f"✅ 检测到Gemini CLI: {version}")
            
            # 执行安装步骤
            self.create_extension_structure()
            self.install_manifest()
            self.install_main_script()
            self.install_prompts()
            self.configure_gemini_cli()
            self.create_shell_integration()
            self.install_examples()
            
            # 验证安装
            if self.verify_installation():
                self.test_command()
                self.print_summary()
                return True
            else:
                print("\n❌ 安装验证失败，请检查错误信息")
                return False
                
        except Exception as e:
            print(f"\n❌ 安装失败: {e}")
            return False
    
    def uninstall(self) -> bool:
        """卸载ChemAgent"""
        if not self.quiet:
            print("🗑️  卸载ChemAgent for Gemini CLI...")
        
        try:
            # 删除扩展目录
            if self.chemagent_dir.exists():
                shutil.rmtree(self.chemagent_dir)
                if not self.quiet:
                    print("✅ 删除扩展目录")
            
            # 清理配置
            if self.config_file.exists():
                with open(self.config_file, "r") as f:
                    config = json.load(f)
                
                # 移除ChemAgent相关配置
                if "extensions" in config and "chemagent" in config["extensions"]:
                    del config["extensions"]["chemagent"]
                
                if "chemistry" in config:
                    del config["chemistry"]
                
                if "commands" in config:
                    # 移除化学命令
                    chem_commands = [k for k in config["commands"] if k.startswith("chem:")]
                    for cmd in chem_commands:
                        del config["commands"][cmd]
                
                with open(self.config_file, "w") as f:
                    json.dump(config, f, indent=2)
                
                if not self.quiet:
                    print("✅ 清理配置文件")
            
            # 删除别名文件
            aliases_file = self.gemini_home / "chemagent_aliases.sh"
            if aliases_file.exists():
                aliases_file.unlink()
                if not self.quiet:
                    print("✅ 删除别名文件")
            
            if not self.quiet:
                print("\n✅ ChemAgent已完全卸载")
                print("   如果您已将source命令添加到shell配置，请手动移除")
            
            return True
            
        except Exception as e:
            print(f"❌ 卸载失败: {e}")
            return False


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="ChemAgent for Gemini CLI 安装器",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s                    # 交互式安装
  %(prog)s --yes             # 自动安装（跳过确认）
  %(prog)s --uninstall       # 卸载ChemAgent
  %(prog)s --quiet           # 静默模式
        """
    )
    
    parser.add_argument(
        "--yes", "-y",
        action="store_true",
        help="自动确认所有提示"
    )
    
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="静默模式，最小化输出"
    )
    
    parser.add_argument(
        "--uninstall",
        action="store_true",
        help="卸载ChemAgent"
    )
    
    parser.add_argument(
        "--verify",
        action="store_true",
        help="仅验证安装状态"
    )
    
    args = parser.parse_args()
    
    installer = GeminiCLIInstaller(
        quiet=args.quiet,
        yes=args.yes
    )
    
    try:
        if args.uninstall:
            success = installer.uninstall()
        elif args.verify:
            success = installer.verify_installation()
            if not args.quiet:
                if success:
                    print("✅ ChemAgent已正确安装")
                else:
                    print("❌ ChemAgent未安装或安装不完整")
        else:
            success = installer.install()
        
        sys.exit(0 if success else 1)
        
    except KeyboardInterrupt:
        print("\n\n⚠️  安装已被用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 发生错误: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
