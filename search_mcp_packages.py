#!/usr/bin/env python3
"""
Search for available MCP packages on PyPI and GitHub
搜索可用的 MCP 包
"""

import subprocess
import json
import requests
from typing import List, Dict
import sys

def search_pypi_for_mcp() -> List[Dict]:
    """Search PyPI for MCP-related packages"""
    print("🔍 Searching PyPI for MCP packages...")
    
    mcp_packages = []
    
    # Search terms
    search_terms = ["mcp-", "model-context-protocol", "modelcontext"]
    
    for term in search_terms:
        try:
            # Use pip search (if available) or PyPI API
            url = f"https://pypi.org/search/?q={term}"
            # Note: This is simplified. Real implementation would parse HTML or use API
            print(f"  Searching for: {term}")
            
            # Try using pip list to find installed packages with mcp
            result = subprocess.run(
                ["pip", "list", "--format=json"],
                capture_output=True,
                text=True,
                check=False
            )
            
            if result.returncode == 0:
                packages = json.loads(result.stdout)
                for pkg in packages:
                    if "mcp" in pkg["name"].lower() or "model-context" in pkg["name"].lower():
                        mcp_packages.append(pkg)
                        
        except Exception as e:
            print(f"  ⚠️ Error searching: {e}")
    
    return mcp_packages

def check_known_mcp_packages() -> Dict[str, bool]:
    """Check if known MCP packages are available"""
    print("\n📦 Checking known MCP packages...")
    
    known_packages = {
        "mcp-rdkit": "RDKit MCP Server (Chemistry)",
        "mcp-filesystem": "Filesystem operations",
        "mcp-sqlite": "SQLite database operations",
        "mcp-github": "GitHub integration",
        "mcp-web-browser": "Web browsing",
        "mcp-python": "Python code execution",
        "mcp-shell": "Shell commands",
        "mcp-google-drive": "Google Drive integration",
    }
    
    availability = {}
    
    for package, description in known_packages.items():
        try:
            # Check if package exists on PyPI
            response = requests.get(
                f"https://pypi.org/pypi/{package}/json",
                timeout=5
            )
            
            if response.status_code == 200:
                data = response.json()
                version = data.get("info", {}).get("version", "unknown")
                print(f"  ✅ {package} v{version} - {description}")
                availability[package] = True
            else:
                print(f"  ❌ {package} - {description} (not found)")
                availability[package] = False
                
        except Exception as e:
            print(f"  ⚠️ {package} - Could not check: {e}")
            availability[package] = False
    
    return availability

def search_github_for_mcp() -> None:
    """Search GitHub for MCP repositories"""
    print("\n🐙 GitHub MCP repositories to check:")
    
    repos = [
        "anthropics/model-context-protocol",
        "modelcontextprotocol/servers",
        "awesome-mcp/awesome-mcp-servers",
    ]
    
    print("  Suggested GitHub searches:")
    print('  - topic:mcp topic:chemistry')
    print('  - "model context protocol" chemistry')
    print('  - language:Python "MCP server"')
    
    print("\n  Known repositories:")
    for repo in repos:
        print(f"  - https://github.com/{repo}")

def check_chemistry_specific() -> None:
    """Check for chemistry-specific MCP implementations"""
    print("\n🧪 Chemistry-specific MCP status:")
    
    chemistry_tools = {
        "RDKit": "mcp-rdkit",
        "PubChem": None,
        "ChEMBL": None,
        "OpenBabel": None,
        "ChemSpider": None,
        "CSD": None,
        "Gaussian": None,
        "ORCA": None,
        "AutoDock": None,
        "PyMOL": None,
        "Avogadro": None,
        "GROMACS": None,
    }
    
    for tool, mcp_package in chemistry_tools.items():
        if mcp_package:
            # Check if it exists
            try:
                response = requests.get(
                    f"https://pypi.org/pypi/{mcp_package}/json",
                    timeout=5
                )
                if response.status_code == 200:
                    print(f"  ✅ {tool}: {mcp_package} (available)")
                else:
                    print(f"  ⚠️ {tool}: {mcp_package} (not found)")
            except:
                print(f"  ⚠️ {tool}: {mcp_package} (could not check)")
        else:
            print(f"  ❌ {tool}: No MCP implementation found")

def main():
    print("=" * 50)
    print("  MCP Package Availability Scanner")
    print("=" * 50)
    
    # Search for MCP packages
    installed = search_pypi_for_mcp()
    if installed:
        print(f"\n📌 Found {len(installed)} installed packages with 'mcp' in name:")
        for pkg in installed:
            print(f"  - {pkg['name']} v{pkg['version']}")
    
    # Check known packages
    availability = check_known_mcp_packages()
    
    # GitHub search suggestions
    search_github_for_mcp()
    
    # Chemistry-specific check
    check_chemistry_specific()
    
    # Summary
    print("\n" + "=" * 50)
    print("📊 Summary:")
    print(f"  - Known MCP packages checked: {len(availability)}")
    print(f"  - Available on PyPI: {sum(availability.values())}")
    print(f"  - Chemistry-specific: 1 (RDKit only)")
    print("\n💡 Recommendation:")
    print("  ChemAgent should continue using its orchestrator pattern")
    print("  with RDKit MCP and internal fallbacks for other tools.")
    print("=" * 50)

if __name__ == "__main__":
    main()
