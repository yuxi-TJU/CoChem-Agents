#!/usr/bin/env python3
"""
ChemAgent Installation Script
Installs chemistry enhancements for Claude Code and Gemini CLI
"""

import os
import sys
import subprocess
import platform
import json
import shutil
from pathlib import Path
from typing import List, Dict, Optional
import argparse


class ChemAgentInstaller:
    """Installer for ChemAgent enhancements"""
    
    def __init__(self, platform: str = "auto", verbose: bool = True, dev: bool = False):
        self.verbose = verbose
        self.dev = dev
        self.target_platform = platform  # claude-code, gemini-cli, or auto
        self.system = platform.system().lower()
        self.python_version = sys.version_info
        self.errors = []
        self.warnings = []
        self.installed_platforms = []
        
    def log(self, message: str, level: str = "INFO"):
        """Log message with level"""
        if self.verbose or level in ["ERROR", "WARNING"]:
            prefix = {
                "INFO": "ℹ️ ",
                "SUCCESS": "✅",
                "WARNING": "⚠️ ",
                "ERROR": "❌"
            }.get(level, "")
            print(f"{prefix} {message}")
            
    def run_command(self, cmd: List[str], check: bool = True) -> Optional[str]:
        """Run shell command"""
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=check
            )
            return result.stdout
        except subprocess.CalledProcessError as e:
            self.log(f"Command failed: {' '.join(cmd)}", "ERROR")
            self.log(f"Error: {e.stderr}", "ERROR")
            if check:
                raise
            return None
            
    def check_python_version(self):
        """Check Python version compatibility"""
        self.log("Checking Python version...")
        if self.python_version < (3, 9):
            raise RuntimeError(f"Python 3.9+ required, found {sys.version}")
        self.log(f"Python {sys.version} ✓", "SUCCESS")
        
    def check_system_dependencies(self):
        """Check and install system dependencies"""
        self.log("Checking system dependencies...")
        
        dependencies = {
            "git": "Git version control",
            "curl": "HTTP client",
        }
        
        missing = []
        for dep, desc in dependencies.items():
            if shutil.which(dep) is None:
                missing.append(f"{dep} ({desc})")
                
        if missing:
            self.log(f"Missing system dependencies: {', '.join(missing)}", "WARNING")
            self.log("Please install missing dependencies manually", "WARNING")
            self.warnings.append("Some system dependencies are missing")
            
    def create_virtual_environment(self):
        """Create and activate virtual environment"""
        self.log("Setting up virtual environment...")
        
        venv_path = Path("venv")
        if not venv_path.exists():
            self.run_command([sys.executable, "-m", "venv", "venv"])
            self.log("Virtual environment created", "SUCCESS")
        else:
            self.log("Virtual environment already exists", "INFO")
            
        # Get pip path in venv
        if self.platform == "windows":
            self.pip_cmd = ["venv\\Scripts\\pip.exe"]
            self.python_cmd = ["venv\\Scripts\\python.exe"]
        else:
            self.pip_cmd = ["venv/bin/pip"]
            self.python_cmd = ["venv/bin/python"]
            
    def upgrade_pip(self):
        """Upgrade pip to latest version"""
        self.log("Upgrading pip...")
        self.run_command(self.pip_cmd + ["install", "--upgrade", "pip"])
        self.log("Pip upgraded", "SUCCESS")
        
    def install_base_dependencies(self):
        """Install base Python dependencies"""
        self.log("Installing base dependencies...")
        
        # Install from pyproject.toml
        if self.dev:
            self.run_command(self.pip_cmd + ["install", "-e", ".[dev,mcp,visualization]"])
        else:
            self.run_command(self.pip_cmd + ["install", "-e", "."])
            
        self.log("Base dependencies installed", "SUCCESS")
        
    def install_chemistry_tools(self):
        """Install specialized chemistry tools"""
        self.log("Installing chemistry tools...")
        
        # RDKit installation (special handling needed)
        try:
            # Try conda first if available
            if shutil.which("conda"):
                self.log("Installing RDKit via conda...")
                self.run_command(["conda", "install", "-c", "conda-forge", "rdkit", "-y"])
            else:
                # Fallback to pip
                self.log("Installing RDKit via pip...")
                self.run_command(self.pip_cmd + ["install", "rdkit"])
                
            self.log("RDKit installed", "SUCCESS")
        except Exception as e:
            self.log(f"RDKit installation failed: {e}", "WARNING")
            self.warnings.append("RDKit installation failed - some features may be limited")
            
    def setup_configuration(self):
        """Setup default configuration files"""
        self.log("Setting up configuration...")
        
        config_dir = Path.home() / ".chemagent"
        config_dir.mkdir(exist_ok=True)
        
        # Create default config
        default_config = {
            "model": {
                "provider": "anthropic",
                "model_name": "claude-3.5-sonnet",
                "temperature": 0.7,
                "max_tokens": 4096
            },
            "tools": [
                {"name": "rdkit_tool", "enabled": True},
                {"name": "pubchem", "enabled": True},
                {"name": "molecule_visualizer", "enabled": True}
            ],
            "mcp": {
                "enabled": True,
                "servers": []
            },
            "workspace_dir": str(config_dir),
            "cache_dir": str(config_dir / "cache"),
            "log_level": "INFO"
        }
        
        config_file = config_dir / "config.json"
        if not config_file.exists():
            with open(config_file, "w") as f:
                json.dump(default_config, f, indent=2)
            self.log(f"Configuration created at {config_file}", "SUCCESS")
        else:
            self.log("Configuration already exists", "INFO")
            
        # Create necessary directories
        (config_dir / "cache").mkdir(exist_ok=True)
        (config_dir / "logs").mkdir(exist_ok=True)
        (config_dir / "prompts").mkdir(exist_ok=True)
        (config_dir / "tools").mkdir(exist_ok=True)
        (config_dir / "roles").mkdir(exist_ok=True)
        
    def setup_environment_variables(self):
        """Setup environment variables"""
        self.log("Setting up environment variables...")
        
        env_file = Path(".env")
        if not env_file.exists():
            env_content = """# ChemAgent Environment Variables

# LLM API Keys (add your keys here)
ANTHROPIC_API_KEY=your_anthropic_api_key_here
OPENAI_API_KEY=your_openai_api_key_here
GEMINI_API_KEY=your_gemini_api_key_here

# ChemAgent Configuration
CHEMAGENT_LOG_LEVEL=INFO
CHEMAGENT_WORKSPACE_DIR=~/.chemagent

# MCP Server Configuration
CHEMAGENT_MCP_ENABLED=true
CHEMAGENT_MCP_PORT=8765
"""
            with open(env_file, "w") as f:
                f.write(env_content)
            self.log("Environment file created (.env)", "SUCCESS")
            self.log("⚠️  Please add your API keys to .env file", "WARNING")
        else:
            self.log("Environment file already exists", "INFO")
            
    def download_additional_resources(self):
        """Download additional resources and models"""
        self.log("Downloading additional resources...")
        
        resources_dir = Path("resources")
        resources_dir.mkdir(exist_ok=True)
        
        # Download example molecules database
        molecules_db = resources_dir / "molecules.json"
        if not molecules_db.exists():
            try:
                import requests
                # This would download from a real URL in production
                example_molecules = {
                    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
                    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                    "paracetamol": "CC(=O)NC1=CC=C(C=C1)O",
                    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
                }
                with open(molecules_db, "w") as f:
                    json.dump(example_molecules, f, indent=2)
                self.log("Example molecules database created", "SUCCESS")
            except Exception as e:
                self.log(f"Failed to create molecules database: {e}", "WARNING")
                
    def install_claude_code_enhancement(self):
        """Install Claude Code enhancement"""
        self.log("Installing Claude Code (Cursor) enhancement...")
        
        from chemagent.enhancers import ClaudeCodeEnhancer
        
        enhancer = ClaudeCodeEnhancer()
        if enhancer.install():
            self.log("Claude Code enhancement installed", "SUCCESS")
            self.installed_platforms.append("claude-code")
        else:
            self.log("Claude Code enhancement installation failed", "ERROR")
            self.errors.append("Claude Code installation failed")
    
    def install_gemini_cli_enhancement(self):
        """Install Gemini CLI enhancement"""
        self.log("Installing Gemini CLI enhancement...")
        
        from chemagent.enhancers import GeminiCLIEnhancer
        
        enhancer = GeminiCLIEnhancer()
        if enhancer.install():
            self.log("Gemini CLI enhancement installed", "SUCCESS")
            self.installed_platforms.append("gemini-cli")
        else:
            self.log("Gemini CLI enhancement installation failed", "ERROR")
            self.errors.append("Gemini CLI installation failed")
    
    def detect_platforms(self):
        """Detect available AI platforms"""
        self.log("Detecting available AI platforms...")
        
        platforms = []
        
        # Check for Cursor
        cursor_dir = Path.home() / ".cursor"
        if cursor_dir.exists() or shutil.which("cursor"):
            platforms.append("claude-code")
            self.log("Detected: Claude Code (Cursor)", "SUCCESS")
        
        # Check for Gemini CLI
        if shutil.which("gemini") or (Path.home() / ".gemini").exists():
            platforms.append("gemini-cli")
            self.log("Detected: Gemini CLI", "SUCCESS")
        
        if not platforms:
            self.log("No supported AI platforms detected", "WARNING")
            self.warnings.append("Install Cursor or Gemini CLI first")
        
        return platforms
        
    def run_tests(self):
        """Run basic tests to verify installation"""
        self.log("Running installation tests...")
        
        test_script = """
import sys
sys.path.insert(0, '.')

# Test imports
try:
    from chemagent import ChemAgent
    print("✓ ChemAgent import successful")
except ImportError as e:
    print(f"✗ ChemAgent import failed: {e}")
    sys.exit(1)

try:
    from chemagent.tools import RDKitTool
    print("✓ Tools import successful")
except ImportError as e:
    print(f"✗ Tools import failed: {e}")
    sys.exit(1)

try:
    from chemagent.roles import OrganicChemist
    print("✓ Roles import successful")
except ImportError as e:
    print(f"✗ Roles import failed: {e}")
    sys.exit(1)

print("All tests passed!")
"""
        
        # Write test script
        test_file = Path("test_install.py")
        with open(test_file, "w") as f:
            f.write(test_script)
            
        try:
            # Run test
            result = self.run_command(self.python_cmd + ["test_install.py"])
            if result:
                print(result)
            self.log("Installation tests passed", "SUCCESS")
        except Exception as e:
            self.log(f"Installation tests failed: {e}", "ERROR")
            self.errors.append("Installation tests failed")
        finally:
            # Clean up test file
            test_file.unlink(missing_ok=True)
            
    def print_summary(self):
        """Print installation summary"""
        print("\n" + "="*60)
        print("ChemAgent Enhancement Installation Summary")
        print("="*60)
        
        if self.errors:
            print("\n❌ Errors:")
            for error in self.errors:
                print(f"  - {error}")
                
        if self.warnings:
            print("\n⚠️  Warnings:")
            for warning in self.warnings:
                print(f"  - {warning}")
                
        if not self.errors:
            print("\n✅ Installation completed successfully!")
            print(f"\nInstalled for: {', '.join(self.installed_platforms)}")
            print("\nNext steps:")
            
            if "claude-code" in self.installed_platforms:
                print("\n📝 For Claude Code (Cursor):")
                print("   1. Open Cursor in your project")
                print("   2. Use cc- commands (e.g., cc-analyze)")
                print("   3. Invoke sub-agents with @ (e.g., @organic-chemist)")
            
            if "gemini-cli" in self.installed_platforms:
                print("\n📝 For Gemini CLI:")
                print("   1. Use 'gemini chem:analyze <SMILES>'")
                print("   2. Process batches with 'gemini chem:batch'")
                print("   3. Source aliases: source ~/.gemini/chemagent_aliases.sh")
            
            print("\n🔑 Don't forget to add your API keys to .env file")
            print("\nFor more information, see README.md")
        else:
            print("\n❌ Installation failed with errors")
            print("Please fix the errors and run the installer again")
            
    def install(self):
        """Run full installation"""
        try:
            self.log("Starting ChemAgent installation...", "INFO")
            
            self.check_python_version()
            self.check_system_dependencies()
            self.create_virtual_environment()
            self.upgrade_pip()
            self.install_base_dependencies()
            self.install_chemistry_tools()
            self.setup_configuration()
            self.setup_environment_variables()
            
            # Detect or use specified platforms
            if self.target_platform == "auto":
                platforms = self.detect_platforms()
            elif self.target_platform == "all":
                platforms = ["claude-code", "gemini-cli"]
            else:
                platforms = [self.target_platform]
            
            # Install enhancements for detected platforms
            for platform in platforms:
                if platform == "claude-code":
                    self.install_claude_code_enhancement()
                elif platform == "gemini-cli":
                    self.install_gemini_cli_enhancement()
            
            self.download_additional_resources()
            
            if not self.errors:
                self.run_tests()
                
        except Exception as e:
            self.log(f"Installation failed: {e}", "ERROR")
            self.errors.append(str(e))
            
        finally:
            self.print_summary()
            return len(self.errors) == 0


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description="ChemAgent Enhancement Installation")
    parser.add_argument("--platform", choices=["claude-code", "gemini-cli", "all", "auto"],
                       default="auto", help="Target platform (default: auto-detect)")
    parser.add_argument("--dev", action="store_true", help="Install development dependencies")
    parser.add_argument("--quiet", action="store_true", help="Minimal output")
    parser.add_argument("--skip-tests", action="store_true", help="Skip installation tests")
    
    args = parser.parse_args()
    
    installer = ChemAgentInstaller(
        platform=args.platform,
        verbose=not args.quiet,
        dev=args.dev
    )
    
    success = installer.install()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
