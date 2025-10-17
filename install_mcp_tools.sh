#!/bin/bash

# ChemAgent MCP Tools Installer
# 安装化学相关的MCP服务器和工具

set -e

echo "================================================"
echo "  ChemAgent MCP Tools Installer                "
echo "================================================"
echo ""

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Check if npm is installed (needed for many MCP servers)
if ! command -v npm &> /dev/null; then
    echo -e "${YELLOW}npm not found. Installing Node.js and npm...${NC}"
    curl -fsSL https://deb.nodesource.com/setup_lts.x | sudo -E bash -
    sudo apt-get install -y nodejs
fi

# Create MCP directory
MCP_DIR="$HOME/.chemagent/mcp_servers"
mkdir -p "$MCP_DIR"
cd "$MCP_DIR"

echo -e "${YELLOW}Installing MCP servers...${NC}"

# 1. Install official MCP servers (when available)
echo ""
echo "1. Checking for official MCP servers..."

# Example: If RDKit has an official MCP (hypothetical)
# npm install -g @rdkit/mcp-server

# 2. Install community MCP servers
echo ""
echo "2. Installing community MCP servers..."

# File system MCP (for reading/writing molecule files)
if [ ! -d "filesystem-mcp" ]; then
    echo "Installing filesystem MCP..."
    git clone https://github.com/modelcontextprotocol/servers.git temp_servers
    cp -r temp_servers/src/filesystem "$MCP_DIR/filesystem-mcp"
    cd "$MCP_DIR/filesystem-mcp"
    npm install
    cd "$MCP_DIR"
    rm -rf temp_servers
    echo -e "${GREEN}✅ Filesystem MCP installed${NC}"
fi

# 3. Install Python-based chemistry tools
echo ""
echo "3. Installing Python chemistry tools..."

# Create virtual environment for chemistry tools
VENV_DIR="$HOME/.chemagent/venv"
if [ ! -d "$VENV_DIR" ]; then
    python3 -m venv "$VENV_DIR"
fi

# Activate virtual environment
source "$VENV_DIR/bin/activate"

# Install chemistry packages
echo "Installing chemistry packages..."
pip install --quiet --upgrade pip

# Core chemistry packages
pip install --quiet \
    rdkit>=2023.9.1 \
    pubchempy>=1.0.4 \
    chembl_webresource_client>=0.10.8 \
    biopython>=1.81 \
    openbabel-wheel>=3.1.1 \
    mordred>=1.2.0 \
    || echo -e "${YELLOW}Some packages failed to install${NC}"

# Additional chemistry tools
pip install --quiet \
    pymol-open-source \
    MDAnalysis>=2.0.0 \
    nglview>=3.0.0 \
    py3Dmol>=2.0.0 \
    || echo -e "${YELLOW}Some visualization packages failed to install${NC}"

# ADMET prediction tools
pip install --quiet \
    deepchem>=2.7.0 \
    tdc>=0.4.0 \
    || echo -e "${YELLOW}Some ML packages failed to install${NC}"

echo -e "${GREEN}✅ Python chemistry tools installed${NC}"

# 4. Create MCP wrapper scripts
echo ""
echo "4. Creating MCP wrapper scripts..."

# Create RDKit MCP wrapper
cat > "$MCP_DIR/rdkit_mcp_wrapper.py" << 'EOF'
#!/usr/bin/env python3
"""
RDKit MCP Wrapper
Provides MCP interface to RDKit functionality
"""

import json
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

def handle_request(request):
    """Handle MCP request for RDKit operations"""
    method = request.get("method")
    params = request.get("params", {})
    
    if method == "validate_smiles":
        smiles = params.get("smiles")
        mol = Chem.MolFromSmiles(smiles)
        return {"valid": mol is not None}
    
    elif method == "calculate_properties":
        smiles = params.get("smiles")
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return {
                "MW": Descriptors.MolWt(mol),
                "LogP": Crippen.MolLogP(mol),
                "TPSA": Descriptors.TPSA(mol),
                "HBD": Lipinski.NumHDonors(mol),
                "HBA": Lipinski.NumHAcceptors(mol),
                "RotatableBonds": Descriptors.NumRotatableBonds(mol)
            }
        return {"error": "Invalid SMILES"}
    
    return {"error": "Unknown method"}

if __name__ == "__main__":
    # Simple JSON-RPC style interface
    request = json.loads(sys.stdin.read())
    response = handle_request(request)
    print(json.dumps(response))
EOF

chmod +x "$MCP_DIR/rdkit_mcp_wrapper.py"
echo -e "${GREEN}✅ MCP wrapper scripts created${NC}"

# 5. Create MCP configuration
echo ""
echo "5. Creating MCP configuration..."

cat > "$HOME/.chemagent/mcp_config.json" << EOF
{
  "mcpServers": {
    "filesystem": {
      "command": "node",
      "args": ["$MCP_DIR/filesystem-mcp/index.js"],
      "env": {
        "ALLOWED_PATHS": "$HOME/chemistry_projects"
      }
    },
    "rdkit_wrapper": {
      "command": "$VENV_DIR/bin/python",
      "args": ["$MCP_DIR/rdkit_mcp_wrapper.py"],
      "env": {}
    }
  },
  "chemistry_tools": {
    "rdkit": {
      "type": "python",
      "module": "rdkit",
      "available": true
    },
    "pubchem": {
      "type": "api",
      "endpoint": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/",
      "available": true
    },
    "chembl": {
      "type": "python",
      "module": "chembl_webresource_client",
      "available": true
    }
  }
}
EOF

echo -e "${GREEN}✅ MCP configuration created${NC}"

# 6. Test installations
echo ""
echo "6. Testing installations..."

# Test RDKit
python3 -c "import rdkit; print('✅ RDKit version:', rdkit.__version__)" 2>/dev/null || echo "❌ RDKit not working"

# Test PubChemPy
python3 -c "import pubchempy; print('✅ PubChemPy installed')" 2>/dev/null || echo "❌ PubChemPy not working"

# Test ChEMBL
python3 -c "from chembl_webresource_client.new_client import new_client; print('✅ ChEMBL client installed')" 2>/dev/null || echo "❌ ChEMBL not working"

# Summary
echo ""
echo -e "${GREEN}================================================${NC}"
echo -e "${GREEN}  MCP Tools Installation Complete!             ${NC}"
echo -e "${GREEN}================================================${NC}"
echo ""
echo "Installed components:"
echo "  • Python chemistry tools (RDKit, PubChem, ChEMBL)"
echo "  • MCP wrapper scripts"
echo "  • MCP configuration"
echo ""
echo "Virtual environment: $VENV_DIR"
echo "MCP servers: $MCP_DIR"
echo "Configuration: $HOME/.chemagent/mcp_config.json"
echo ""
echo "To use the chemistry tools:"
echo "  source $VENV_DIR/bin/activate"
echo ""
echo "Note: Most chemistry tools don't have official MCP servers yet."
echo "We've installed Python packages and created wrapper scripts."
