"""
Gemini CLI Enhancer
Provides chemistry enhancements for Google's Gemini CLI
"""

from typing import Dict, Any, List, Optional
from pathlib import Path
import json
import os
from .base import BaseEnhancer


class GeminiCLIEnhancer(BaseEnhancer):
    """Enhancer for Gemini CLI"""
    
    def get_name(self) -> str:
        return "ChemAgent for Gemini CLI"
    
    def get_platform(self) -> str:
        return "gemini-cli"
    
    def get_default_config(self) -> Dict[str, Any]:
        """Get default configuration for Gemini CLI"""
        return {
            "platform": "gemini-cli",
            "version": "1.0.0",
            "gemini_config_dir": str(Path.home() / ".gemini"),
            "extensions_dir": str(Path.home() / ".gemini" / "extensions"),
            "commands_prefix": "chem:",
            "auto_load": True,
            "features": {
                "chemistry_commands": True,
                "multimodal_molecules": True,  # Gemini's strength
                "batch_processing": True,
                "cloud_functions": True
            },
            "default_models": [
                "gemini-pro",
                "gemini-pro-vision"  # For molecule images
            ]
        }
    
    def install(self) -> bool:
        """Install ChemAgent for Gemini CLI"""
        print("🚀 Installing ChemAgent for Gemini CLI...")
        
        # 1. Create Gemini extension directory
        self._create_extension_structure()
        
        # 2. Install chemistry commands as Gemini functions
        self._install_gemini_functions()
        
        # 3. Setup chemistry prompts
        self._setup_chemistry_prompts()
        
        # 4. Configure Gemini CLI
        self._configure_gemini_cli()
        
        # 5. Create command aliases
        self._create_command_aliases()
        
        print("✅ ChemAgent for Gemini CLI installed successfully!")
        print("\n📝 Quick Start:")
        print("1. Use 'gemini chem:analyze <SMILES>' for molecular analysis")
        print("2. Use 'gemini chem:synthesize <target>' for synthesis planning")
        print("3. Upload molecule images for structure recognition")
        print("4. Batch process multiple molecules with 'gemini chem:batch'")
        
        return True
    
    def _create_extension_structure(self):
        """Create Gemini extension directory structure"""
        ext_dir = Path.home() / ".gemini" / "extensions" / "chemagent"
        ext_dir.mkdir(parents=True, exist_ok=True)
        
        # Create subdirectories
        (ext_dir / "functions").mkdir(exist_ok=True)
        (ext_dir / "prompts").mkdir(exist_ok=True)
        (ext_dir / "models").mkdir(exist_ok=True)
        (ext_dir / "data").mkdir(exist_ok=True)
        
        # Create extension manifest
        manifest = {
            "name": "ChemAgent",
            "version": "1.0.0",
            "description": "Chemistry enhancement for Gemini CLI",
            "author": "ChemAgent Team",
            "entry_point": "chemagent_main.py",
            "commands": [
                "chem:analyze",
                "chem:synthesize",
                "chem:predict",
                "chem:optimize",
                "chem:search",
                "chem:batch",
                "chem:visualize"
            ],
            "capabilities": [
                "molecular_analysis",
                "synthesis_planning",
                "reaction_prediction",
                "image_to_structure",
                "batch_processing"
            ]
        }
        
        manifest_file = ext_dir / "manifest.json"
        with open(manifest_file, "w") as f:
            json.dump(manifest, f, indent=2)
        
        print("✅ Created Gemini extension structure")
    
    def _install_gemini_functions(self):
        """Install chemistry functions for Gemini"""
        ext_dir = Path.home() / ".gemini" / "extensions" / "chemagent"
        
        # Create main entry point
        main_content = '''#!/usr/bin/env python3
"""
ChemAgent Main Entry Point for Gemini CLI
"""

import sys
import json
from pathlib import Path

# Add ChemAgent to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from chemagent.enhancers.gemini_functions import handle_chemistry_command

def main():
    """Main entry point for Gemini CLI"""
    if len(sys.argv) < 2:
        print("Usage: gemini chem:<command> [args]")
        return 1
    
    command = sys.argv[1]
    args = sys.argv[2:] if len(sys.argv) > 2 else []
    
    try:
        result = handle_chemistry_command(command, args)
        print(json.dumps(result, indent=2))
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())
'''
        
        main_file = ext_dir / "chemagent_main.py"
        with open(main_file, "w") as f:
            f.write(main_content)
        os.chmod(main_file, 0o755)
        
        # Create function definitions
        functions = {
            "analyze": {
                "name": "chem:analyze",
                "description": "Analyze molecular structure and properties",
                "parameters": {
                    "structure": "SMILES, InChI, or name",
                    "options": ["properties", "druglike", "toxicity", "visualize"]
                }
            },
            "synthesize": {
                "name": "chem:synthesize",
                "description": "Plan synthesis route for target molecule",
                "parameters": {
                    "target": "Target molecule SMILES or name",
                    "starting_materials": "Optional list of available starting materials",
                    "constraints": ["steps", "yield", "cost", "green"]
                }
            },
            "predict": {
                "name": "chem:predict",
                "description": "Predict reaction outcomes",
                "parameters": {
                    "reactants": "List of reactant SMILES",
                    "conditions": "Reaction conditions",
                    "mechanism": "Show mechanism details"
                }
            },
            "batch": {
                "name": "chem:batch",
                "description": "Process multiple molecules in batch",
                "parameters": {
                    "input_file": "CSV or SDF file with molecules",
                    "operation": "Operation to perform",
                    "output_file": "Output file path"
                }
            }
        }
        
        functions_dir = ext_dir / "functions"
        for name, func_def in functions.items():
            func_file = functions_dir / f"{name}.json"
            with open(func_file, "w") as f:
                json.dump(func_def, f, indent=2)
        
        print(f"✅ Installed {len(functions)} Gemini functions")
    
    def _setup_chemistry_prompts(self):
        """Setup chemistry-specific prompts for Gemini"""
        ext_dir = Path.home() / ".gemini" / "extensions" / "chemagent"
        prompts_dir = ext_dir / "prompts"
        
        prompts = {
            "molecule_analysis": """Analyze the following molecular structure:
{structure}

Provide:
1. Molecular properties (MW, LogP, TPSA, etc.)
2. Drug-likeness assessment
3. Potential applications
4. Safety considerations
5. Synthetic accessibility""",
            
            "synthesis_planning": """Design a synthesis route for:
Target: {target}
Starting materials: {starting_materials}
Constraints: {constraints}

Provide:
1. Retrosynthetic analysis
2. Step-by-step forward synthesis
3. Reaction conditions for each step
4. Expected yields
5. Alternative routes""",
            
            "image_to_structure": """Analyze the molecular structure in this image:
[Image provided]

Extract:
1. Chemical structure (SMILES/InChI)
2. Identify functional groups
3. Suggest the compound name
4. Predict key properties""",
            
            "batch_template": """Process the following molecules:
{molecules}

For each molecule, perform:
{operation}

Output format: {format}"""
        }
        
        for name, prompt in prompts.items():
            prompt_file = prompts_dir / f"{name}.txt"
            with open(prompt_file, "w") as f:
                f.write(prompt)
        
        print(f"✅ Created {len(prompts)} chemistry prompts")
    
    def _configure_gemini_cli(self):
        """Configure Gemini CLI for chemistry"""
        config_dir = Path.home() / ".gemini"
        config_dir.mkdir(exist_ok=True)
        
        # Load existing config or create new
        config_file = config_dir / "config.json"
        if config_file.exists():
            with open(config_file, "r") as f:
                config = json.load(f)
        else:
            config = {}
        
        # Add ChemAgent configuration
        config["extensions"] = config.get("extensions", {})
        config["extensions"]["chemagent"] = {
            "enabled": True,
            "auto_load": True,
            "default_model": "gemini-pro",
            "vision_model": "gemini-pro-vision",
            "batch_size": 10,
            "cache_results": True
        }
        
        # Add chemistry-specific settings
        config["chemistry"] = {
            "default_format": "SMILES",
            "auto_validate": True,
            "include_3d": False,
            "safety_checks": True,
            "citation_style": "ACS"
        }
        
        # Save updated config
        with open(config_file, "w") as f:
            json.dump(config, f, indent=2)
        
        print("✅ Configured Gemini CLI")
    
    def _create_command_aliases(self):
        """Create command aliases for easier access"""
        aliases_content = """# ChemAgent aliases for Gemini CLI
alias gchem='gemini chem:analyze'
alias gsynth='gemini chem:synthesize'
alias greact='gemini chem:predict'
alias gbatch='gemini chem:batch'
alias gmol='gemini chem:analyze --visualize'

# Function for quick molecule analysis
chem() {
    gemini chem:analyze "$1" --properties --druglike --visualize
}

# Function for synthesis planning
synthesize() {
    gemini chem:synthesize "$1" --steps 5 --green
}

# Function for batch processing
chembatch() {
    gemini chem:batch "$1" --operation analyze --output results.csv
}

echo "ChemAgent for Gemini CLI loaded ✅"
"""
        
        # Create aliases file
        aliases_file = Path.home() / ".gemini" / "chemagent_aliases.sh"
        with open(aliases_file, "w") as f:
            f.write(aliases_content)
        
        print("✅ Created command aliases")
        print(f"   Add 'source {aliases_file}' to your shell profile")
    
    def uninstall(self) -> bool:
        """Uninstall ChemAgent from Gemini CLI"""
        print("🗑️  Uninstalling ChemAgent for Gemini CLI...")
        
        # Remove extension directory
        ext_dir = Path.home() / ".gemini" / "extensions" / "chemagent"
        if ext_dir.exists():
            import shutil
            shutil.rmtree(ext_dir)
            print("✅ Removed extension directory")
        
        # Remove from config
        config_file = Path.home() / ".gemini" / "config.json"
        if config_file.exists():
            with open(config_file, "r") as f:
                config = json.load(f)
            
            if "extensions" in config and "chemagent" in config["extensions"]:
                del config["extensions"]["chemagent"]
            if "chemistry" in config:
                del config["chemistry"]
            
            with open(config_file, "w") as f:
                json.dump(config, f, indent=2)
            print("✅ Removed from configuration")
        
        # Remove aliases file
        aliases_file = Path.home() / ".gemini" / "chemagent_aliases.sh"
        if aliases_file.exists():
            aliases_file.unlink()
            print("✅ Removed aliases")
        
        print("✅ ChemAgent for Gemini CLI uninstalled")
        return True
    
    def handle_multimodal_input(self, text: str, image_path: Optional[str] = None) -> Dict[str, Any]:
        """Handle multimodal input (text + image) - Gemini's strength"""
        result = {
            "input_type": "multimodal" if image_path else "text",
            "text": text
        }
        
        if image_path:
            # Process image for molecular structure
            result["image"] = {
                "path": image_path,
                "analysis": "Structure recognition from image",
                "extracted_smiles": None  # Would use OCR/ML model
            }
        
        return result
