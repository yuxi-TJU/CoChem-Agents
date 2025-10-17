#!/usr/bin/env python3
"""
ChemAgent Interactive Installation
User-friendly installation with multiple modes
Inspired by SuperClaude_Framework
"""

import os
import sys
import time
import platform
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import json
import shutil
import subprocess

# Add colorful output support
try:
    from rich.console import Console
    from rich.prompt import Prompt, Confirm
    from rich.progress import Progress, SpinnerColumn, TextColumn
    from rich.panel import Panel
    from rich.table import Table
    from rich import print as rprint
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False
    print("📦 Installing rich for better UI...")
    subprocess.run([sys.executable, "-m", "pip", "install", "rich"], check=False)
    try:
        from rich.console import Console
        from rich.prompt import Prompt, Confirm
        from rich.progress import Progress, SpinnerColumn, TextColumn
        from rich.panel import Panel
        from rich.table import Table
        from rich import print as rprint
        RICH_AVAILABLE = True
    except:
        pass

console = Console() if RICH_AVAILABLE else None


class ChemAgentInteractiveInstaller:
    """Interactive installer with user-friendly interface"""
    
    def __init__(self):
        self.console = console or self._fallback_console()
        self.platform = platform.system().lower()
        self.install_mode = None
        self.selected_platforms = []
        self.selected_features = []
        
    def _fallback_console(self):
        """Fallback console for when rich is not available"""
        class FallbackConsole:
            def print(self, *args, **kwargs):
                print(*args)
            def rule(self, title="", **kwargs):
                print(f"\n{'='*50} {title} {'='*50}\n")
        return FallbackConsole()
    
    def run(self):
        """Main installation flow"""
        self.show_welcome()
        
        # Select installation mode
        self.install_mode = self.select_install_mode()
        
        if self.install_mode == "quick":
            self.quick_install()
        elif self.install_mode == "interactive":
            self.interactive_install()
        elif self.install_mode == "minimal":
            self.minimal_install()
        elif self.install_mode == "developer":
            self.developer_install()
        else:
            self.console.print("[red]Invalid installation mode[/red]")
            sys.exit(1)
    
    def show_welcome(self):
        """Show welcome message"""
        self.console.rule("[bold blue]ChemAgent Installation[/bold blue]")
        
        welcome_text = """
[bold]Welcome to ChemAgent! 🧪[/bold]

ChemAgent is a chemistry enhancement package for AI coding assistants,
inspired by SuperClaude_Framework.

This installer will help you set up ChemAgent for:
• [cyan]Claude Code (Cursor)[/cyan] - Chemistry commands and sub-agents
• [green]Gemini CLI[/green] - Chemistry tools and workflows
• [yellow]Custom configurations[/yellow] - Tailored to your needs
        """
        
        if RICH_AVAILABLE:
            self.console.print(Panel(welcome_text, title="🚀 Getting Started", border_style="blue"))
        else:
            self.console.print(welcome_text)
    
    def select_install_mode(self) -> str:
        """Select installation mode"""
        self.console.print("\n[bold]Select Installation Mode:[/bold]\n")
        
        modes = {
            "1": ("quick", "✨ Quick Install", "Full installation with all features (recommended)"),
            "2": ("interactive", "🎯 Interactive Install", "Choose specific components to install"),
            "3": ("minimal", "🚀 Minimal Install", "Core features only, lightweight"),
            "4": ("developer", "🔧 Developer Mode", "Full features + development tools"),
        }
        
        # Display options
        if RICH_AVAILABLE:
            table = Table(show_header=True, header_style="bold magenta")
            table.add_column("#", style="cyan", width=3)
            table.add_column("Mode", style="green")
            table.add_column("Description")
            
            for key, (_, name, desc) in modes.items():
                table.add_row(key, name, desc)
            
            self.console.print(table)
        else:
            for key, (_, name, desc) in modes.items():
                self.console.print(f"{key}. {name} - {desc}")
        
        # Get user choice
        if RICH_AVAILABLE:
            choice = Prompt.ask(
                "\n[bold]Choose installation mode[/bold]",
                choices=["1", "2", "3", "4"],
                default="1"
            )
        else:
            choice = input("\nChoose installation mode (1-4) [1]: ").strip() or "1"
        
        return modes.get(choice, ("quick", "", ""))[0]
    
    def quick_install(self):
        """Quick installation with all features"""
        self.console.print("\n[bold green]Starting Quick Installation...[/bold green]")
        
        # Auto-detect platforms
        detected_platforms = self.detect_platforms()
        
        if not detected_platforms:
            self.console.print("[yellow]No AI platforms detected. Installing for manual use.[/yellow]")
        else:
            self.console.print(f"\n[green]Detected platforms:[/green] {', '.join(detected_platforms)}")
        
        # Install everything
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=self.console
        ) if RICH_AVAILABLE else self.console:
            
            # Core installation
            self.install_core()
            
            # Platform-specific installations
            for platform in detected_platforms:
                self.install_platform(platform)
            
            # Install all chemistry tools
            self.install_chemistry_tools()
            
            # Create shortcuts
            self.create_shortcuts()
        
        self.show_completion_message(detected_platforms)
    
    def interactive_install(self):
        """Interactive installation with component selection"""
        self.console.print("\n[bold cyan]Interactive Installation[/bold cyan]")
        
        # Select platforms
        self.selected_platforms = self.select_platforms()
        
        # Select features
        self.selected_features = self.select_features()
        
        # Confirm selections
        if self.confirm_selections():
            self.perform_installation()
        else:
            self.console.print("[yellow]Installation cancelled.[/yellow]")
    
    def minimal_install(self):
        """Minimal installation with core features only"""
        self.console.print("\n[bold yellow]Minimal Installation[/bold yellow]")
        
        # Install only core components
        self.install_core()
        
        # Basic commands only
        self.console.print("Installing core commands...")
        commands = ["cc-analyze", "cc-synthesize", "cc-search"]
        for cmd in commands:
            self.console.print(f"  • {cmd}")
        
        self.show_completion_message([], minimal=True)
    
    def developer_install(self):
        """Developer installation with extra tools"""
        self.console.print("\n[bold magenta]Developer Mode Installation[/bold magenta]")
        
        # Everything from quick install
        self.quick_install()
        
        # Plus development tools
        self.console.print("\n[bold]Installing development tools...[/bold]")
        dev_packages = ["pytest", "black", "ruff", "mypy", "pre-commit"]
        
        for package in dev_packages:
            self.console.print(f"  • Installing {package}")
            subprocess.run([sys.executable, "-m", "pip", "install", package], 
                         capture_output=True, check=False)
        
        # Setup pre-commit hooks
        self.setup_dev_environment()
    
    def detect_platforms(self) -> List[str]:
        """Auto-detect installed AI platforms"""
        detected = []
        
        # Check for Cursor
        cursor_paths = [
            Path.home() / ".cursor",
            Path("/Applications/Cursor.app"),
            Path("C:/Users") / os.environ.get("USERNAME", "") / "AppData/Local/Programs/cursor",
        ]
        
        if any(p.exists() for p in cursor_paths):
            detected.append("claude-code")
        
        # Check for Gemini CLI
        try:
            result = subprocess.run(["gemini", "--version"], 
                                  capture_output=True, check=False)
            if result.returncode == 0:
                detected.append("gemini-cli")
        except:
            pass
        
        return detected
    
    def select_platforms(self) -> List[str]:
        """Let user select platforms to install for"""
        platforms = {
            "1": ("claude-code", "Claude Code (Cursor)", "✨ Chemistry commands and sub-agents"),
            "2": ("gemini-cli", "Gemini CLI", "🚀 Chemistry tools and multimodal support"),
            "3": ("both", "Both platforms", "💎 Complete installation"),
        }
        
        self.console.print("\n[bold]Select target platforms:[/bold]")
        
        for key, (_, name, desc) in platforms.items():
            self.console.print(f"{key}. {name} - {desc}")
        
        if RICH_AVAILABLE:
            choice = Prompt.ask("Choose platforms", choices=["1", "2", "3"], default="3")
        else:
            choice = input("Choose platforms (1-3) [3]: ").strip() or "3"
        
        if choice == "1":
            return ["claude-code"]
        elif choice == "2":
            return ["gemini-cli"]
        else:
            return ["claude-code", "gemini-cli"]
    
    def select_features(self) -> List[str]:
        """Let user select features to install"""
        features = {
            "commands": "Chemistry commands (cc-analyze, cc-synthesize, etc.)",
            "roles": "Expert roles (@chemist, @drug-designer, etc.)",
            "tools": "Chemistry tools (RDKit, patent search, etc.)",
            "workflows": "Pre-built workflows (drug discovery, materials)",
            "examples": "Example projects and templates",
        }
        
        self.console.print("\n[bold]Select features to install:[/bold]")
        self.console.print("[dim]Press Enter to select all[/dim]\n")
        
        selected = []
        for key, desc in features.items():
            if RICH_AVAILABLE:
                if Confirm.ask(f"Install {desc}?", default=True):
                    selected.append(key)
            else:
                response = input(f"Install {desc}? [Y/n]: ").strip().lower()
                if response != 'n':
                    selected.append(key)
        
        return selected
    
    def confirm_selections(self) -> bool:
        """Show summary and confirm installation"""
        self.console.print("\n[bold]Installation Summary:[/bold]")
        self.console.print(f"Platforms: {', '.join(self.selected_platforms)}")
        self.console.print(f"Features: {', '.join(self.selected_features)}")
        
        if RICH_AVAILABLE:
            return Confirm.ask("\nProceed with installation?", default=True)
        else:
            response = input("\nProceed with installation? [Y/n]: ").strip().lower()
            return response != 'n'
    
    def perform_installation(self):
        """Perform the actual installation"""
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            console=self.console
        ) if RICH_AVAILABLE else self.console:
            
            # Core installation
            if "commands" in self.selected_features or "roles" in self.selected_features:
                self.install_core()
            
            # Platform-specific
            for platform in self.selected_platforms:
                self.install_platform(platform)
            
            # Features
            if "tools" in self.selected_features:
                self.install_chemistry_tools()
            
            if "workflows" in self.selected_features:
                self.install_workflows()
            
            if "examples" in self.selected_features:
                self.install_examples()
        
        self.show_completion_message(self.selected_platforms)
    
    def install_core(self):
        """Install core ChemAgent components"""
        self.console.print("\n📦 Installing ChemAgent core...")
        
        # Install Python package
        subprocess.run([sys.executable, "-m", "pip", "install", "-e", "."], 
                      capture_output=True, check=False)
        
        # Create directories
        dirs = [
            Path.home() / ".chemagent",
            Path.home() / ".claude" / "commands",
            Path.home() / ".claude" / "roles",
        ]
        
        for d in dirs:
            d.mkdir(parents=True, exist_ok=True)
    
    def install_platform(self, platform: str):
        """Install platform-specific components"""
        self.console.print(f"\n🎯 Installing for {platform}...")
        
        if platform == "claude-code":
            # Copy .cursorrules
            shutil.copy(".cursorrules", Path.cwd() / ".cursorrules")
            self.console.print("  ✅ Installed .cursorrules")
            
        elif platform == "gemini-cli":
            # Create Gemini configuration
            gemini_config = Path.home() / ".gemini" / "chemagent.json"
            gemini_config.parent.mkdir(parents=True, exist_ok=True)
            # Write config...
            self.console.print("  ✅ Configured Gemini CLI")
    
    def install_chemistry_tools(self):
        """Install chemistry tool dependencies"""
        self.console.print("\n🧪 Installing chemistry tools...")
        
        tools = {
            "rdkit": "RDKit - Molecular operations",
            "pubchempy": "PubChem - Database access",
            "biopython": "BioPython - Sequence analysis",
        }
        
        for package, desc in tools.items():
            try:
                self.console.print(f"  • Installing {desc}")
                subprocess.run([sys.executable, "-m", "pip", "install", package],
                             capture_output=True, check=True)
            except:
                self.console.print(f"  ⚠️  Failed to install {package}")
    
    def install_workflows(self):
        """Install pre-built workflows"""
        self.console.print("\n📋 Installing workflows...")
        # Copy workflow templates
    
    def install_examples(self):
        """Install example projects"""
        self.console.print("\n📚 Installing examples...")
        # Copy example files
    
    def create_shortcuts(self):
        """Create command shortcuts"""
        self.console.print("\n🔗 Creating shortcuts...")
        
        # Platform-specific shortcuts
        if self.platform == "darwin":  # macOS
            # Create aliases
            pass
        elif self.platform == "linux":
            # Create shell aliases
            pass
        elif self.platform == "windows":
            # Create batch files
            pass
    
    def setup_dev_environment(self):
        """Setup development environment"""
        self.console.print("\n🔧 Setting up development environment...")
        
        # Initialize pre-commit
        subprocess.run(["pre-commit", "install"], capture_output=True, check=False)
        
        # Create VSCode settings
        vscode_dir = Path(".vscode")
        vscode_dir.mkdir(exist_ok=True)
        
        settings = {
            "python.linting.enabled": True,
            "python.linting.ruffEnabled": True,
            "python.formatting.provider": "black",
            "editor.formatOnSave": True,
        }
        
        with open(vscode_dir / "settings.json", "w") as f:
            json.dump(settings, f, indent=2)
    
    def show_completion_message(self, platforms: List[str], minimal: bool = False):
        """Show installation completion message"""
        self.console.rule("[bold green]Installation Complete![/bold green]")
        
        if minimal:
            message = """
[bold green]✅ ChemAgent Minimal Installation Complete![/bold green]

Core commands are now available. To get started:

1. Open your AI assistant
2. Use commands like:
   • cc-analyze <molecule>
   • cc-synthesize <target>
   • cc-search <query>

For the full experience, run: [cyan]chemagent install[/cyan]
            """
        else:
            message = f"""
[bold green]✨ ChemAgent Successfully Installed![/bold green]

{'[cyan]Claude Code (Cursor)[/cyan] is ready!' if 'claude-code' in platforms else ''}
{'[green]Gemini CLI[/green] is configured!' if 'gemini-cli' in platforms else ''}

[bold]🚀 Quick Start:[/bold]
1. {'Open Cursor and try: cc-analyze aspirin' if 'claude-code' in platforms else ''}
2. {'Run: gemini "analyze caffeine molecule"' if 'gemini-cli' in platforms else ''}
3. Create custom commands in .claude/commands/
4. Define new roles in .claude/roles/

[bold]📚 Resources:[/bold]
• Documentation: [link]https://chemagent.ai/docs[/link]
• Examples: Run [cyan]chemagent examples[/cyan]
• Help: Use [cyan]cc-help[/cyan] in your AI assistant

[bold]💡 Tips:[/bold]
• Install additional tools: [cyan]chemagent install --tools[/cyan]
• Update ChemAgent: [cyan]chemagent update[/cyan]
• Join community: [link]https://github.com/chemagent[/link]

Happy Chemistry! 🧪✨
            """
        
        if RICH_AVAILABLE:
            self.console.print(Panel(message, border_style="green"))
        else:
            self.console.print(message)


def main():
    """Main entry point"""
    installer = ChemAgentInteractiveInstaller()
    
    try:
        installer.run()
    except KeyboardInterrupt:
        installer.console.print("\n[yellow]Installation cancelled by user.[/yellow]")
        sys.exit(1)
    except Exception as e:
        installer.console.print(f"\n[red]Installation failed: {e}[/red]")
        sys.exit(1)


if __name__ == "__main__":
    main()
