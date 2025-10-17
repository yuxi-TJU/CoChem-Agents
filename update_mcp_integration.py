#!/usr/bin/env python3
"""
Update script to ensure MCP servers are properly integrated
"""

import json
import os
import subprocess
import sys
from pathlib import Path


def update_claude_enhancer():
    """Update Claude enhancer to properly install MCP servers"""
    
    enhancer_file = Path("chemagent/enhancers/claude_enhancer.py")
    
    # Read the current file
    with open(enhancer_file, "r") as f:
        content = f.read()
    
    # Check if MCP setup is already included in install method
    if "_setup_mcp_servers" not in content:
        # Add the new method
        new_method = '''
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
'''
        
        # Insert the new method before the install method
        insert_pos = content.find("    def install(self)")
        if insert_pos > 0:
            content = content[:insert_pos] + new_method + "\n" + content[insert_pos:]
        
        # Update the install method to call _setup_mcp_servers
        old_install = """        # Create workspace configuration
        self._create_workspace_config()"""
        
        new_install = """        # Install and configure MCP servers
        self._setup_mcp_servers()
        
        # Create workspace configuration
        self._create_workspace_config()"""
        
        content = content.replace(old_install, new_install)
        
        # Add imports if needed
        if "import subprocess" not in content:
            import_pos = content.find("import json")
            if import_pos > 0:
                content = content[:import_pos] + "import subprocess\n" + content[import_pos:]
        
        # Write back
        with open(enhancer_file, "w") as f:
            f.write(content)
        
        print("✅ Updated Claude enhancer with MCP integration")
    else:
        print("✅ Claude enhancer already has MCP integration")


def update_gemini_integration():
    """Update Gemini integration to include MCP references"""
    
    gemini_rules = Path("gemini_rules.md")
    
    if gemini_rules.exists():
        with open(gemini_rules, "r") as f:
            content = f.read()
        
        if "MCP Servers" not in content:
            # Add MCP section
            mcp_section = """
## Available MCP Servers

ChemAgent provides MCP servers for advanced chemistry operations:

1. **RDKit MCP** (port 8766)
   - Molecular analysis and descriptors
   - Fingerprint generation
   - 3D conformer generation

2. **PubChem MCP** (port 8767)
   - Compound search by name/SMILES/InChI
   - Property retrieval
   - Synonym lookup

3. **ChEMBL MCP** (port 8768)
   - Bioactivity data
   - Drug and target information
   - Clinical trial data

4. **OpenBabel MCP** (port 8769)
   - Format conversion
   - 3D structure generation
   - Geometry optimization

To start MCP servers:
```bash
mcp start-all
```

To use in commands:
```bash
# Example: Search PubChem
curl -X POST http://localhost:8767 -d '{"method":"search_compound","params":{"query":"aspirin"}}'
```
"""
            
            # Add before the "## Best Practices" section if it exists
            insert_pos = content.find("## Best Practices")
            if insert_pos > 0:
                content = content[:insert_pos] + mcp_section + "\n" + content[insert_pos:]
            else:
                content += "\n" + mcp_section
            
            with open(gemini_rules, "w") as f:
                f.write(content)
            
            print("✅ Updated Gemini rules with MCP information")
    else:
        print("⚠️  gemini_rules.md not found")


def update_main_installer():
    """Ensure main installer calls MCP installation"""
    
    installer = Path("chemagent_install.py")
    
    with open(installer, "r") as f:
        content = f.read()
    
    # Check if MCP is installed by default
    if "install_mcp_servers()" not in content:
        # Find the quick install section
        quick_install = content.find("def quick_install():")
        if quick_install > 0:
            # Find where to add MCP installation
            marker = "installer.install_platform(platform)"
            pos = content.find(marker, quick_install)
            if pos > 0:
                # Find the end of that line
                end_pos = content.find("\n", pos)
                # Add MCP installation after platform installation
                addition = """
        
        # Install MCP servers
        print("\\n📦 Installing MCP servers...")
        install_mcp_servers()"""
                
                content = content[:end_pos] + addition + content[end_pos:]
                
                with open(installer, "w") as f:
                    f.write(content)
                
                print("✅ Updated main installer to include MCP servers")
    else:
        print("✅ Main installer already includes MCP installation")


def main():
    print("🔧 Updating MCP Integration...")
    print("="*50)
    
    # Update Claude enhancer
    update_claude_enhancer()
    
    # Update Gemini integration
    update_gemini_integration()
    
    # Update main installer
    update_main_installer()
    
    print("\n" + "="*50)
    print("✅ MCP integration update complete!")
    print("\nNow MCP servers will be:")
    print("  1. Automatically installed during ChemAgent installation")
    print("  2. Configured in Claude Desktop config")
    print("  3. Referenced in Gemini rules")
    print("\nTo test the installation:")
    print("  python chemagent_install.py")


if __name__ == "__main__":
    main()
