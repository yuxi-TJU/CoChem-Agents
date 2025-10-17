"""
Claude Code (Cursor) Enhancer
Provides chemistry enhancements for Claude Code in Cursor IDE
"""

from typing import Dict, Any, List, Optional
from pathlib import Path
import subprocess
import json
import shutil
import os
from .base import BaseEnhancer


class ClaudeCodeEnhancer(BaseEnhancer):
    """Enhancer for Claude Code in Cursor IDE"""
    
    def get_name(self) -> str:
        return "ChemAgent for Claude Code"
    
    def get_platform(self) -> str:
        return "claude-code"
    
    def get_default_config(self) -> Dict[str, Any]:
        """Get default configuration for Claude Code"""
        return {
            "platform": "claude-code",
            "version": "1.0.0",
            "cursor_workspace": str(Path.home() / ".cursor"),
            "rules_file": ".cursorrules",
            "commands_prefix": "cc-",
            "sub_agents_prefix": "@",
            "auto_load": True,
            "features": {
                "chemistry_commands": True,
                "sub_agents": True,
                "mcp_integration": True,
                "auto_complete": True,
                "inline_docs": True
            },
            "default_tools": [
                "rdkit",
                "pubchem",
                "chembl",
                "pdb"
            ]
        }
    

    def _setup_mcp_servers(self):
        """Install and configure MCP servers for Claude Desktop"""
        print("🔧 Setting up MCP servers...")
        
        # First, run the MCP installation script
        install_script = Path(__file__).parent.parent.parent / "install_all_mcp_servers.sh"
        if install_script.exists():
            try:
                print("  Installing MCP servers...")
                result = subprocess.run(
                    ["bash", str(install_script)], 
                    capture_output=True, 
                    text=True,
                    check=False
                )
                if result.returncode == 0:
                    print("  ✅ MCP servers installed")
                else:
                    print(f"  ⚠️  MCP installation had issues: {result.stderr}")
            except Exception as e:
                print(f"  ⚠️  Could not run MCP installer: {e}")
        
        # Configure Claude Desktop
        config_paths = [
            Path.home() / "Library" / "Application Support" / "Claude" / "claude_desktop_config.json",  # macOS
            Path.home() / ".config" / "Claude" / "claude_desktop_config.json",  # Linux
            Path(os.environ.get("APPDATA", "")) / "Claude" / "claude_desktop_config.json" if os.name == "nt" else None  # Windows
        ]
        
        config_file = None
        for path in config_paths:
            if path and path.parent.exists():
                config_file = path
                path.parent.mkdir(parents=True, exist_ok=True)
                break
        
        if not config_file:
            # Fallback to Cursor config
            config_file = Path.home() / ".cursor" / "mcp_config.json"
            config_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Load or create config
        if config_file.exists():
            with open(config_file, "r") as f:
                config = json.load(f)
        else:
            config = {}
        
        # Add MCP servers
        if "mcpServers" not in config:
            config["mcpServers"] = {}
        
        config["mcpServers"].update({
            "rdkit": {
                "command": "python",
                "args": ["-m", "mcp_rdkit"]
            },
            "pubchem": {
                "command": "python",
                "args": ["-m", "mcp_pubchem_server"]
            },
            "chembl": {
                "command": "python",
                "args": ["-m", "mcp_chembl_server"]
            },
            "openbabel": {
                "command": "python",
                "args": ["-m", "mcp_openbabel_server"]
            }
        })
        
        # Save config
        with open(config_file, "w") as f:
            json.dump(config, f, indent=2)
        
        print(f"  ✅ MCP configuration saved to {config_file}")

    def install(self) -> bool:
        """Install ChemAgent for Claude Code"""
        print("🚀 Installing ChemAgent for Claude Code...")
        
        # 1. Create .cursorrules file in project
        self._install_cursor_rules()
        
        # 2. Install chemistry commands
        self._install_commands()
        
        # 3. Setup sub-agents
        self._setup_sub_agents()
        
        # 4. Install and configure MCP servers
        self._setup_mcp_servers()
        
        # 5. Configure MCP if available (legacy)
        self._configure_mcp()
        
        # 6. Create workspace config
        self._create_workspace_config()
        
        print("✅ ChemAgent for Claude Code installed successfully!")
        print("\n📝 Quick Start:")
        print("1. Open Cursor in your chemistry project")
        print("2. Use cc- commands (e.g., cc-analyze, cc-synthesize)")
        print("3. Invoke sub-agents with @ (e.g., @organic-chemist)")
        print("4. Chemistry context is automatically loaded")
        
        return True
    
    def _install_cursor_rules(self):
        """Install .cursorrules file"""
        rules_content = """# ChemAgent Rules for Claude Code

## Chemistry Context
This project uses ChemAgent for chemistry-specific AI assistance.

## Available Commands (cc-series)
- cc-analyze: Molecular analysis
- cc-synthesize: Synthesis planning
- cc-predict: Reaction prediction
- cc-optimize: Property optimization
- cc-search: Database search
- cc-drug-design: Drug discovery workflow
- cc-materials: Materials design
- cc-retro: Retrosynthetic analysis
- cc-screen: Virtual screening
- cc-qsar: QSAR modeling

## Chemistry Sub-Agents
- @organic-chemist: Organic synthesis expert
- @drug-designer: Medicinal chemistry specialist
- @materials-scientist: Materials and polymers expert
- @comp-chemist: Computational chemistry specialist
- @analytical-chemist: Analytical methods expert

## Chemistry Tools
- RDKit: Molecular operations
- PubChem: Chemical database
- ChEMBL: Bioactivity data
- PDB: Protein structures
- ChemSpider: Chemical search

## Best Practices
1. Always validate SMILES/InChI before processing
2. Consider stereochemistry in reactions
3. Check drug-likeness for pharmaceuticals
4. Include safety warnings for hazardous compounds
5. Cite relevant literature and databases

## Response Format
- Use standard chemical notation (SMILES, InChI)
- Include 2D structure visualizations when helpful
- Provide step-by-step synthesis procedures
- Include reaction conditions and yields
- Add safety and disposal information

## Project-Specific Context
[Add your project-specific chemistry context here]
"""
        
        # Write to current directory
        with open(".cursorrules", "w") as f:
            f.write(rules_content)
        print("✅ Created .cursorrules file")
    
    def _install_commands(self):
        """Install chemistry commands"""
        from ..commands import command_loader
        
        # Commands are automatically loaded from Markdown files
        commands = command_loader.list_commands()
        
        print(f"✅ Loaded {len(commands)} chemistry commands from Markdown files")
    
    def _setup_sub_agents(self):
        """Setup chemistry sub-agents from Markdown files"""
        # Load roles from Markdown files
        roles_dir = Path(__file__).parent.parent.parent / "roles"
        
        sub_agents = []
        if roles_dir.exists():
            for role_file in roles_dir.glob("*.md"):
                sub_agents.append(role_file.stem)
        
        print(f"✅ Found {len(sub_agents)} chemistry roles: {', '.join(sub_agents)}")
    
    def _setup_mcp_servers(self):
        """Install and configure MCP servers for Claude Desktop"""
        print("🔧 Setting up MCP servers...")
        
        # First, run the MCP installation script
        install_script = Path(__file__).parent.parent.parent / "install_all_mcp_servers.sh"
        if install_script.exists():
            try:
                print("  Installing MCP servers...")
                result = subprocess.run(
                    ["bash", str(install_script)], 
                    capture_output=True, 
                    text=True,
                    check=False
                )
                if result.returncode == 0:
                    print("  ✅ MCP servers installed")
                else:
                    print(f"  ⚠️  MCP installation had issues")
            except Exception as e:
                print(f"  ⚠️  Could not run MCP installer: {e}")
        
        # Configure Claude Desktop
        config_paths = [
            Path.home() / "Library" / "Application Support" / "Claude" / "claude_desktop_config.json",  # macOS
            Path.home() / ".config" / "Claude" / "claude_desktop_config.json",  # Linux
            Path(os.environ.get("APPDATA", "")) / "Claude" / "claude_desktop_config.json" if os.name == "nt" else None  # Windows
        ]
        
        config_file = None
        for path in config_paths:
            if path and path.parent.exists():
                config_file = path
                path.parent.mkdir(parents=True, exist_ok=True)
                break
        
        if not config_file:
            # Fallback to Cursor config
            config_file = Path.home() / ".cursor" / "mcp_config.json"
            config_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Load or create config
        if config_file.exists():
            with open(config_file, "r") as f:
                config = json.load(f)
        else:
            config = {}
        
        # Add MCP servers
        if "mcpServers" not in config:
            config["mcpServers"] = {}
        
        config["mcpServers"].update({
            "rdkit": {
                "command": "python",
                "args": ["-m", "mcp_rdkit"]
            },
            "pubchem": {
                "command": "python",
                "args": ["-m", "mcp_pubchem_server"]
            },
            "chembl": {
                "command": "python",
                "args": ["-m", "mcp_chembl_server"]
            },
            "openbabel": {
                "command": "python",
                "args": ["-m", "mcp_openbabel_server"]
            }
        })
        
        # Save config
        with open(config_file, "w") as f:
            json.dump(config, f, indent=2)
        
        print(f"  ✅ MCP configuration saved to {config_file}")
    
    def _configure_mcp(self):
        """Configure MCP integration"""
        mcp_config = {
            "servers": [
                {
                    "name": "chemagent-mcp",
                    "command": "python",
                    "args": ["-m", "chemagent.mcp.server"],
                    "env": {},
                    "capabilities": {
                        "tools": True,
                        "context": True,
                        "resources": True
                    }
                }
            ]
        }
        
        # Write MCP config for Cursor
        config_dir = Path.home() / ".cursor" / "mcp"
        config_dir.mkdir(parents=True, exist_ok=True)
        
        config_file = config_dir / "chemagent.json"
        with open(config_file, "w") as f:
            json.dump(mcp_config, f, indent=2)
        
        print("✅ Configured MCP integration")
    
    def _create_workspace_config(self):
        """Create workspace configuration"""
        workspace_dir = Path.home() / ".chemagent" / "claude-code"
        workspace_dir.mkdir(parents=True, exist_ok=True)
        
        # Create command shortcuts
        shortcuts = {
            "cc-a": "cc-analyze",
            "cc-s": "cc-synthesize",
            "cc-p": "cc-predict",
            "cc-o": "cc-optimize",
            "cc-dd": "cc-drug-design",
        }
        
        shortcuts_file = workspace_dir / "shortcuts.json"
        with open(shortcuts_file, "w") as f:
            json.dump(shortcuts, f, indent=2)
        
        # Create templates
        templates_dir = workspace_dir / "templates"
        templates_dir.mkdir(exist_ok=True)
        
        # Save configuration
        config_file = workspace_dir / "config.json"
        with open(config_file, "w") as f:
            json.dump(self.config, f, indent=2)
        
        print("✅ Created workspace configuration")
    
    def uninstall(self) -> bool:
        """Uninstall ChemAgent from Claude Code"""
        print("🗑️  Uninstalling ChemAgent for Claude Code...")
        
        # Remove .cursorrules if it's ours
        if Path(".cursorrules").exists():
            with open(".cursorrules", "r") as f:
                content = f.read()
                if "ChemAgent Rules" in content:
                    Path(".cursorrules").unlink()
                    print("✅ Removed .cursorrules")
        
        # Remove MCP config
        mcp_config = Path.home() / ".cursor" / "mcp" / "chemagent.json"
        if mcp_config.exists():
            mcp_config.unlink()
            print("✅ Removed MCP configuration")
        
        # Remove workspace
        workspace = Path.home() / ".chemagent" / "claude-code"
        if workspace.exists():
            shutil.rmtree(workspace)
            print("✅ Removed workspace configuration")
        
        print("✅ ChemAgent for Claude Code uninstalled")
        return True
    
    def create_chemistry_context(self, project_type: str) -> str:
        """Create project-specific chemistry context"""
        contexts = {
            "drug-discovery": """
## Drug Discovery Project Context
- Focus on ADMET properties
- Consider synthetic accessibility
- Check for PAINS and toxicophores
- Optimize for target selectivity
- Include pharmacokinetic predictions
""",
            "materials": """
## Materials Science Project Context
- Consider polymer properties
- Focus on structure-property relationships
- Include processing conditions
- Consider scalability
- Add characterization methods
""",
            "synthesis": """
## Organic Synthesis Project Context
- Provide detailed reaction mechanisms
- Include reaction conditions
- Consider protecting group strategies
- Add purification methods
- Include yield estimations
""",
            "education": """
## Chemistry Education Project Context
- Explain concepts clearly
- Provide step-by-step solutions
- Include visual representations
- Add practice problems
- Reference learning resources
"""
        }
        
        return contexts.get(project_type, "")
    
    def enhance_prompt(self, original_prompt: str, context: Optional[Dict[str, Any]] = None) -> str:
        """Enhance a prompt with chemistry context"""
        enhanced = original_prompt
        
        # Add chemistry context if molecule is mentioned
        if any(term in original_prompt.lower() for term in ["molecule", "compound", "drug", "reaction"]):
            enhanced += "\n\n[Chemistry Context: Use appropriate chemical notation, validate structures, and consider safety]"
        
        # Add specific enhancements based on keywords
        if "synthesize" in original_prompt.lower():
            enhanced += "\n[Include: retrosynthetic analysis, reaction conditions, yields, and purification]"
        
        if "drug" in original_prompt.lower():
            enhanced += "\n[Consider: Lipinski's Rule of Five, ADMET properties, and toxicity]"
        
        if "material" in original_prompt.lower():
            enhanced += "\n[Include: structure-property relationships, processing conditions, and characterization]"
        
        return enhanced
