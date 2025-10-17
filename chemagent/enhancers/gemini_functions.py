"""
Gemini CLI Function Handlers
处理Gemini CLI的化学命令
"""

import json
import sys
from typing import Dict, Any, List, Optional
from pathlib import Path

# 尝试导入ChemAgent核心模块
try:
    from chemagent.mcp_tools.orchestrator import MCPOrchestrator
    from chemagent.commands.loader import get_command_prompt, list_available_commands
    from chemagent.roles.loader import get_role, list_available_roles
    CHEMAGENT_CORE_AVAILABLE = True
except ImportError:
    CHEMAGENT_CORE_AVAILABLE = False


class GeminiChemistryHandler:
    """Gemini CLI化学命令处理器"""
    
    def __init__(self):
        """初始化处理器"""
        self.orchestrator = MCPOrchestrator() if CHEMAGENT_CORE_AVAILABLE else None
        self.commands = self._load_commands()
        self.current_role = "chemist"  # 默认角色
    
    def _load_commands(self) -> Dict[str, Any]:
        """加载可用命令映射"""
        commands = {
            "analyze": self.analyze_molecule,
            "synthesize": self.plan_synthesis, 
            "predict": self.predict_reaction,
            "optimize": self.optimize_molecule,
            "search": self.search_database,
            "design": self.design_molecule,
            "batch": self.batch_process,
            "workflow": self.run_workflow,
            "visualize": self.visualize_structure,
            "compare": self.compare_molecules,
            "simulate": self.run_simulation,
            "dock": self.molecular_docking,
            "safety": self.assess_safety,
            "patent": self.check_patent,
            "report": self.generate_report,
            "explain": self.explain_concept,
            "suggest": self.make_suggestions,
            "check": self.run_checks,
            "help": self.show_help
        }
        return commands
    
    def set_role(self, role: str):
        """设置当前角色"""
        available_roles = ["chemist", "drug-designer", "synthesist", "safety-expert", "data-analyst"]
        if role in available_roles:
            self.current_role = role
            return {"status": "success", "role": role}
        else:
            return {"status": "error", "message": f"Unknown role: {role}", "available": available_roles}
    
    def analyze_molecule(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """分析分子结构和性质"""
        molecule = args.get("molecule", args.get("structure", ""))
        
        if not molecule:
            return {"error": "No molecule provided"}
        
        if self.orchestrator:
            # 使用orchestrator执行命令
            result = self.orchestrator.execute_command("cc-analyze", {
                "molecule": molecule,
                "properties": args.get("properties", ["all"]),
                "visualize": args.get("visualize", False)
            })
            return result
        else:
            # 返回基础分析
            return {
                "molecule": molecule,
                "status": "basic_analysis",
                "properties": {
                    "input": molecule,
                    "format": "SMILES" if self._is_smiles(molecule) else "name"
                },
                "message": "Full analysis requires ChemAgent core"
            }
    
    def plan_synthesis(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """规划合成路线"""
        target = args.get("target", "")
        
        if not target:
            return {"error": "No target molecule provided"}
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-synthesize", {
                "target": target,
                "starting_materials": args.get("starting_materials", []),
                "max_steps": args.get("max_steps", 5),
                "constraints": args.get("constraints", {})
            })
            return result
        else:
            return {
                "target": target,
                "status": "requires_core",
                "message": "Synthesis planning requires ChemAgent core modules"
            }
    
    def predict_reaction(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """预测反应结果"""
        reactants = args.get("reactants", [])
        conditions = args.get("conditions", {})
        
        if not reactants:
            return {"error": "No reactants provided"}
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-predict", {
                "reactants": reactants,
                "conditions": conditions,
                "show_mechanism": args.get("show_mechanism", False)
            })
            return result
        else:
            return {
                "reactants": reactants,
                "conditions": conditions,
                "status": "requires_core",
                "message": "Reaction prediction requires ChemAgent core"
            }
    
    def optimize_molecule(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """优化分子结构"""
        molecule = args.get("molecule", "")
        target_property = args.get("target", "druglike")
        
        if not molecule:
            return {"error": "No molecule provided"}
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-optimize", {
                "molecule": molecule,
                "target": target_property,
                "constraints": args.get("constraints", {}),
                "n_suggestions": args.get("n_suggestions", 5)
            })
            return result
        else:
            return {
                "molecule": molecule,
                "target": target_property,
                "status": "requires_core",
                "message": "Molecule optimization requires ChemAgent core"
            }
    
    def search_database(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """搜索化学数据库"""
        query = args.get("query", "")
        database = args.get("database", "all")
        
        if not query:
            return {"error": "No search query provided"}
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-search", {
                "query": query,
                "database": database,
                "limit": args.get("limit", 10),
                "filters": args.get("filters", {})
            })
            return result
        else:
            return {
                "query": query,
                "database": database,
                "status": "requires_core",
                "message": "Database search requires ChemAgent core"
            }
    
    def design_molecule(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """设计新分子"""
        target_properties = args.get("properties", {})
        scaffold = args.get("scaffold", None)
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-design", {
                "target_properties": target_properties,
                "scaffold": scaffold,
                "n_designs": args.get("n_designs", 10)
            })
            return result
        else:
            return {
                "target_properties": target_properties,
                "status": "requires_core",
                "message": "Molecule design requires ChemAgent core"
            }
    
    def batch_process(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """批量处理分子"""
        input_file = args.get("input", "")
        operation = args.get("operation", "analyze")
        
        if not input_file:
            return {"error": "No input file provided"}
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-batch", {
                "input": input_file,
                "operation": operation,
                "output": args.get("output", "results.csv"),
                "parallel": args.get("parallel", True)
            })
            return result
        else:
            return {
                "input": input_file,
                "operation": operation,
                "status": "requires_core",
                "message": "Batch processing requires ChemAgent core"
            }
    
    def run_workflow(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """运行预定义工作流"""
        workflow_name = args.get("workflow", "")
        
        if not workflow_name:
            return {"error": "No workflow specified"}
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-workflow", {
                "workflow": workflow_name,
                "parameters": args.get("parameters", {})
            })
            return result
        else:
            # 列出可用工作流
            available_workflows = [
                "drug-discovery",
                "lead-optimization", 
                "synthesis-planning",
                "safety-assessment",
                "patent-analysis"
            ]
            return {
                "workflow": workflow_name,
                "available": available_workflows,
                "status": "requires_core",
                "message": "Workflow execution requires ChemAgent core"
            }
    
    def visualize_structure(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """可视化分子结构"""
        molecule = args.get("molecule", "")
        style = args.get("style", "2D")
        
        if not molecule:
            return {"error": "No molecule provided"}
        
        return {
            "molecule": molecule,
            "style": style,
            "status": "visualization",
            "message": "Structure will be rendered in Gemini interface",
            "options": {
                "show_hydrogens": args.get("show_hydrogens", False),
                "color_scheme": args.get("color_scheme", "default"),
                "labels": args.get("labels", True)
            }
        }
    
    def compare_molecules(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """比较多个分子"""
        molecules = args.get("molecules", [])
        
        if len(molecules) < 2:
            return {"error": "Need at least 2 molecules to compare"}
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-compare", {
                "molecules": molecules,
                "properties": args.get("properties", ["all"]),
                "similarity_metric": args.get("similarity", "tanimoto")
            })
            return result
        else:
            return {
                "molecules": molecules,
                "status": "requires_core",
                "message": "Molecule comparison requires ChemAgent core"
            }
    
    def run_simulation(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """运行分子模拟"""
        molecule = args.get("molecule", "")
        simulation_type = args.get("type", "dynamics")
        
        if not molecule:
            return {"error": "No molecule provided"}
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-simulate", {
                "molecule": molecule,
                "type": simulation_type,
                "duration": args.get("duration", 100),
                "temperature": args.get("temperature", 300)
            })
            return result
        else:
            return {
                "molecule": molecule,
                "simulation": simulation_type,
                "status": "requires_core",
                "message": "Molecular simulation requires ChemAgent core"
            }
    
    def molecular_docking(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """分子对接"""
        ligand = args.get("ligand", "")
        receptor = args.get("receptor", "")
        
        if not ligand or not receptor:
            return {"error": "Both ligand and receptor required"}
        
        return {
            "ligand": ligand,
            "receptor": receptor,
            "status": "computational",
            "message": "Docking calculation requires specialized software",
            "suggested_tools": ["AutoDock Vina", "Glide", "FlexX"]
        }
    
    def assess_safety(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """安全性评估"""
        molecule = args.get("molecule", "")
        
        if not molecule:
            return {"error": "No molecule provided"}
        
        if self.orchestrator:
            # 使用orchestrator的安全评估
            result = self.orchestrator._assess_safety(molecule)
            return result
        else:
            return {
                "molecule": molecule,
                "status": "requires_core",
                "message": "Safety assessment requires ChemAgent core"
            }
    
    def check_patent(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """专利检查"""
        query = args.get("query", "")
        
        if not query:
            return {"error": "No query provided"}
        
        if self.orchestrator:
            # 使用orchestrator的专利检查
            result = self.orchestrator._check_patents(query)
            return result
        else:
            return {
                "query": query,
                "status": "requires_core",
                "message": "Patent search requires ChemAgent core"
            }
    
    def generate_report(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """生成报告"""
        report_type = args.get("type", "summary")
        data = args.get("data", {})
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-report", {
                "type": report_type,
                "data": data,
                "format": args.get("format", "pdf")
            })
            return result
        else:
            return {
                "type": report_type,
                "status": "requires_core",
                "message": "Report generation requires ChemAgent core"
            }
    
    def explain_concept(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """解释化学概念"""
        concept = args.get("concept", "")
        level = args.get("level", "intermediate")
        
        if not concept:
            return {"error": "No concept provided"}
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-explain", {
                "concept": concept,
                "level": level
            })
            return result
        else:
            return {
                "concept": concept,
                "level": level,
                "status": "educational",
                "message": "Detailed explanation requires ChemAgent core"
            }
    
    def make_suggestions(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """提供建议"""
        context = args.get("context", "")
        goal = args.get("goal", "")
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-suggest", {
                "context": context,
                "goal": goal
            })
            return result
        else:
            return {
                "context": context,
                "goal": goal,
                "status": "requires_core",
                "message": "Suggestions require ChemAgent core"
            }
    
    def run_checks(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """运行各种检查"""
        check_type = args.get("type", "all")
        target = args.get("target", "")
        
        if not target:
            return {"error": "No target provided for checking"}
        
        if self.orchestrator:
            result = self.orchestrator.execute_command("cc-check", {
                "type": check_type,
                "target": target
            })
            return result
        else:
            available_checks = [
                "safety", "patent", "druglike",
                "synthetic_accessibility", "pains", "toxicity"
            ]
            return {
                "target": target,
                "type": check_type,
                "available_checks": available_checks,
                "status": "requires_core"
            }
    
    def show_help(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """显示帮助信息"""
        command = args.get("command", None)
        
        if command:
            # 显示特定命令的帮助
            if command in self.commands:
                return {
                    "command": command,
                    "description": self.commands[command].__doc__,
                    "status": "help"
                }
            else:
                return {"error": f"Unknown command: {command}"}
        else:
            # 显示所有命令
            commands_list = []
            for cmd_name, cmd_func in self.commands.items():
                commands_list.append({
                    "name": f"chem:{cmd_name}",
                    "description": cmd_func.__doc__.strip() if cmd_func.__doc__ else ""
                })
            
            return {
                "commands": commands_list,
                "roles": ["chemist", "drug-designer", "synthesist", "safety-expert", "data-analyst"],
                "version": "1.0.0",
                "status": "help"
            }
    
    def _is_smiles(self, text: str) -> bool:
        """简单判断是否为SMILES"""
        # 简单的SMILES特征检测
        smiles_chars = set("CNOSPFClBrI[]()=#@+-.0123456789")
        if len(text) > 0 and all(c in smiles_chars for c in text):
            return True
        return False


def handle_chemistry_command(command: str, args: List[str]) -> Dict[str, Any]:
    """处理化学命令的主入口"""
    handler = GeminiChemistryHandler()
    
    # 解析命令
    if command.startswith("chem:"):
        command = command[5:]
    
    # 将参数列表转换为字典
    args_dict = {}
    i = 0
    while i < len(args):
        if args[i].startswith("--"):
            key = args[i][2:]
            if i + 1 < len(args) and not args[i + 1].startswith("--"):
                args_dict[key] = args[i + 1]
                i += 2
            else:
                args_dict[key] = True
                i += 1
        else:
            # 位置参数
            if "molecule" not in args_dict and i == 0:
                args_dict["molecule"] = args[i]
            elif "target" not in args_dict and i == 0:
                args_dict["target"] = args[i]
            elif "query" not in args_dict and i == 0:
                args_dict["query"] = args[i]
            else:
                if "args" not in args_dict:
                    args_dict["args"] = []
                args_dict["args"].append(args[i])
            i += 1
    
    # 处理命令
    if command in handler.commands:
        return handler.commands[command](args_dict)
    else:
        return {
            "error": f"Unknown command: {command}",
            "available": list(handler.commands.keys())
        }


def main():
    """测试主函数"""
    if len(sys.argv) < 2:
        print("Usage: python gemini_functions.py <command> [args]")
        sys.exit(1)
    
    command = sys.argv[1]
    args = sys.argv[2:] if len(sys.argv) > 2 else []
    
    result = handle_chemistry_command(command, args)
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
