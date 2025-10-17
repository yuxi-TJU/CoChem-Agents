#!/usr/bin/env python3
"""
ChemAgent Universal Installer
Supports multiple installation modes like SuperClaude_Framework
"""

import sys
import argparse
from pathlib import Path

# Import the interactive installer
sys.path.insert(0, str(Path(__file__).parent))
from install_interactive import ChemAgentInteractiveInstaller
from install import ChemAgentInstaller  # Original installer for backwards compatibility


def main():
    """Main entry point with argument parsing"""
    parser = argparse.ArgumentParser(
        description="ChemAgent Installation - Chemistry Enhancement for AI Assistants",
        epilog="Examples:\n"
               "  chemagent install              # Quick install (recommended)\n"
               "  chemagent install --interactive # Choose components\n"
               "  chemagent install --minimal     # Core features only\n"
               "  chemagent install --help        # Show all options",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Installation modes (mutually exclusive)
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument(
        "--interactive", "-i",
        action="store_true",
        help="Interactive installation with component selection"
    )
    mode_group.add_argument(
        "--minimal", "-m",
        action="store_true",
        help="Minimal installation (core features only)"
    )
    mode_group.add_argument(
        "--profile", "-p",
        choices=["developer", "researcher", "educator"],
        help="Install with a specific profile"
    )
    
    # Platform selection
    parser.add_argument(
        "--platform",
        choices=["claude-code", "gemini-cli", "all", "auto"],
        default="auto",
        help="Target platform (default: auto-detect)"
    )
    
    # Feature flags
    parser.add_argument(
        "--no-tools",
        action="store_true",
        help="Skip chemistry tools installation"
    )
    parser.add_argument(
        "--no-examples",
        action="store_true",
        help="Skip example files installation"
    )
    parser.add_argument(
        "--dev",
        action="store_true",
        help="Include development tools"
    )
    
    # Other options
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Quiet mode (minimal output)"
    )
    parser.add_argument(
        "--yes", "-y",
        action="store_true",
        help="Automatic yes to all prompts"
    )
    parser.add_argument(
        "--uninstall",
        action="store_true",
        help="Uninstall ChemAgent"
    )
    
    # Special commands
    parser.add_argument(
        "command",
        nargs="?",
        choices=["install", "update", "uninstall", "status", "examples", "mcp"],
        default="install",
        help="Command to execute (default: install)"
    )
    
    args = parser.parse_args()
    
    # Handle special commands
    if args.command == "status":
        show_status()
        return
    elif args.command == "examples":
        install_examples()
        return
    elif args.command == "mcp":
        install_mcp_servers()
        return
    elif args.command == "update":
        update_chemagent()
        return
    elif args.command == "uninstall" or args.uninstall:
        uninstall_chemagent()
        return
    
    # Determine installation mode
    if args.interactive:
        # Use interactive installer
        installer = ChemAgentInteractiveInstaller()
        installer.install_mode = "interactive"
    elif args.minimal:
        installer = ChemAgentInteractiveInstaller()
        installer.install_mode = "minimal"
    elif args.profile == "developer":
        installer = ChemAgentInteractiveInstaller()
        installer.install_mode = "developer"
    elif args.quiet or args.yes:
        # Use original non-interactive installer for automation
        installer = ChemAgentInstaller(
            platform=args.platform,
            verbose=not args.quiet,
            dev=args.dev
        )
        installer.install()
        return
    else:
        # Default to interactive installer with quick mode
        installer = ChemAgentInteractiveInstaller()
        installer.install_mode = "quick"
    
    # Run installation
    try:
        installer.run()
    except KeyboardInterrupt:
        print("\n⚠️  Installation cancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Installation failed: {e}")
        sys.exit(1)


def show_status():
    """Show ChemAgent installation status"""
    print("\n🔍 ChemAgent Installation Status\n")
    
    # Check if installed
    try:
        import chemagent
        print(f"✅ ChemAgent {chemagent.__version__} is installed")
    except ImportError:
        print("❌ ChemAgent is not installed")
        return
    
    # Check platforms
    from pathlib import Path
    
    if (Path.home() / ".cursor").exists() or Path(".cursorrules").exists():
        print("✅ Claude Code (Cursor) integration found")
    else:
        print("⚠️  Claude Code (Cursor) not configured")
    
    if (Path.home() / ".gemini").exists():
        gemini_ext = Path.home() / ".gemini" / "extensions" / "chemagent"
        if gemini_ext.exists():
            print("✅ Gemini CLI integration found")
        else:
            print("⚠️  Gemini CLI not configured")
    else:
        print("⚠️  Gemini CLI not installed")
    
    # Check commands
    try:
        from chemagent.commands import list_available_commands
        commands = list_available_commands()
        print(f"✅ {len(commands)} commands available")
    except:
        print("⚠️  Commands not properly loaded")
    
    # Check roles
    try:
        from chemagent.roles import list_available_roles
        roles = list_available_roles()
        print(f"✅ {len(roles)} roles available")
    except:
        print("⚠️  Roles not properly loaded")
    
    print("\nRun 'chemagent install' to fix any issues.")


def install_mcp_servers():
    """Install MCP servers for chemistry tools"""
    print("\n🔧 Installing MCP Servers...\n")
    
    import subprocess
    import os
    
    # Install all MCP servers
    if os.path.exists("install_all_mcp_servers.sh"):
        try:
            subprocess.run(["chmod", "+x", "install_all_mcp_servers.sh"], check=True)
            result = subprocess.run(["./install_all_mcp_servers.sh"], 
                                  capture_output=True, text=True, check=False)
            
            if result.returncode == 0:
                print("✅ All MCP servers installed successfully")
                
                # Try to start the servers
                print("\n🚀 Starting MCP servers...")
                mcp_manager = os.path.expanduser("~/.chemagent/mcp_servers/mcp_manager.py")
                if os.path.exists(mcp_manager):
                    subprocess.run([sys.executable, mcp_manager, "start-all"], 
                                 capture_output=True, check=False)
                    
                    # Show status
                    result = subprocess.run([sys.executable, mcp_manager, "status"],
                                          capture_output=True, text=True, check=False)
                    if result.stdout:
                        print(result.stdout)
            else:
                print("⚠️  Some MCP servers may have failed to install")
                if result.stderr:
                    print(f"Error: {result.stderr}")
                    
        except subprocess.CalledProcessError as e:
            print(f"⚠️  MCP installation error: {e}")
    else:
        # Fallback to individual installers
        print("Installing individual MCP servers...")
        
        # RDKit MCP
        if os.path.exists("install_rdkit_mcp.sh"):
            try:
                subprocess.run(["chmod", "+x", "install_rdkit_mcp.sh"], check=True)
                subprocess.run(["./install_rdkit_mcp.sh"], check=True)
                print("✅ RDKit MCP server installed")
            except:
                print("⚠️  RDKit MCP installation failed")
        
        # Install our custom MCPs from source
        for mcp_name in ["mcp_pubchem", "mcp_chembl", "mcp_openbabel"]:
            mcp_dir = f"mcp_servers/{mcp_name}"
            if os.path.exists(mcp_dir):
                try:
                    subprocess.run([sys.executable, "-m", "pip", "install", "-e", mcp_dir],
                                 capture_output=True, check=False)
                    print(f"✅ {mcp_name} installed from source")
                except:
                    print(f"⚠️  {mcp_name} installation failed")
    
    print("\n✨ MCP servers installation complete!")

def install_examples():
    """Install example files"""
    print("\n📚 Installing ChemAgent Examples...\n")
    
    examples_dir = Path("chemagent_examples")
    examples_dir.mkdir(exist_ok=True)
    
    # Create example command
    example_cmd = examples_dir / "drug-discovery.md"
    example_cmd.write_text("""---
description: Drug discovery workflow example
tools: [read_file, web_search]
---

# Drug Discovery Workflow

This example shows how to create a custom drug discovery workflow.

1. Load target protein structure
2. Search for known inhibitors
3. Generate analogs
4. Predict ADMET properties
5. Rank by drug-likeness
6. Plan synthesis routes
""")
    
    # Create example role
    example_role = examples_dir / "medicinal-chemist.md"
    example_role.write_text("""---
name: medicinal-chemist
description: Specialized medicinal chemistry expert
specialties: [drug-design, pharmacology, admet]
---

# Medicinal Chemistry Expert

You are a specialized medicinal chemist with expertise in:
- Structure-based drug design
- Lead optimization
- ADMET prediction
- Clinical candidate selection
""")
    
    print(f"✅ Examples installed in {examples_dir}")
    print("\nTry copying these to:")
    print("  • .claude/commands/ for commands")
    print("  • .claude/roles/ for roles")


def update_chemagent():
    """Update ChemAgent to latest version"""
    print("\n🔄 Updating ChemAgent...\n")
    
    import subprocess
    
    # Update from git
    try:
        subprocess.run(["git", "pull"], check=True)
        print("✅ Updated from repository")
    except:
        print("⚠️  Could not update from git")
    
    # Reinstall
    subprocess.run([sys.executable, "-m", "pip", "install", "-e", "."], check=False)
    print("✅ Reinstalled ChemAgent")
    
    print("\n✨ Update complete!")


def uninstall_chemagent():
    """Uninstall ChemAgent"""
    print("\n🗑️  Uninstalling ChemAgent...\n")
    
    from pathlib import Path
    import shutil
    
    # Remove installed files
    paths_to_remove = [
        Path.home() / ".chemagent",
        Path(".cursorrules"),
        Path.home() / ".gemini" / "extensions" / "chemagent",
        Path.home() / ".gemini" / "chemagent_aliases.sh",
    ]
    
    for path in paths_to_remove:
        if path.exists():
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()
            print(f"✅ Removed {path}")
    
    # Uninstall package
    import subprocess
    subprocess.run([sys.executable, "-m", "pip", "uninstall", "-y", "chemagent"], 
                   capture_output=True, check=False)
    
    print("\n✅ ChemAgent has been uninstalled")
    print("Note: Custom commands in .claude/ were preserved")


if __name__ == "__main__":
    main()
