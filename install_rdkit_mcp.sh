#!/bin/bash

# RDKit MCP Server Installation Script
# 安装RDKit的官方MCP服务器

set -e

echo "================================================"
echo "  RDKit MCP Server Installer                   "
echo "================================================"
echo ""

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# Check Python version
echo -e "${BLUE}Checking Python version...${NC}"
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version 2>&1 | grep -Po '(?<=Python )\d+\.\d+')
    echo -e "${GREEN}✅ Python $PYTHON_VERSION found${NC}"
else
    echo -e "${RED}❌ Python 3 is required${NC}"
    exit 1
fi

# Installation directory
MCP_DIR="$HOME/.chemagent/mcp_servers"
mkdir -p "$MCP_DIR"

echo ""
echo -e "${YELLOW}Installing RDKit MCP Server...${NC}"

# 1. Install mcp-rdkit from PyPI
echo ""
echo "1. Installing mcp-rdkit package..."
pip install --user mcp-rdkit || {
    echo -e "${YELLOW}Trying with pip3...${NC}"
    pip3 install --user mcp-rdkit
}
echo -e "${GREEN}✅ mcp-rdkit package installed${NC}"

# 2. Create virtual environment (optional but recommended)
echo ""
echo "2. Setting up dedicated environment..."
VENV_DIR="$MCP_DIR/rdkit_venv"
if [ ! -d "$VENV_DIR" ]; then
    python3 -m venv "$VENV_DIR"
    source "$VENV_DIR/bin/activate"
    pip install --upgrade pip
    pip install mcp-rdkit rdkit
    echo -e "${GREEN}✅ Virtual environment created${NC}"
else
    echo "Virtual environment already exists"
    source "$VENV_DIR/bin/activate"
fi

# 3. Create MCP server configuration
echo ""
echo "3. Creating MCP server configuration..."

cat > "$MCP_DIR/rdkit_mcp_config.json" << 'EOF'
{
  "name": "rdkit-mcp-server",
  "version": "1.0.0",
  "description": "RDKit MCP Server for chemistry operations",
  "capabilities": {
    "molecules": {
      "parse_smiles": true,
      "parse_mol": true,
      "generate_2d": true,
      "generate_3d": true
    },
    "descriptors": {
      "molecular_weight": true,
      "logp": true,
      "tpsa": true,
      "lipinski": true,
      "qed": true,
      "sa_score": true
    },
    "fingerprints": {
      "morgan": true,
      "rdkit": true,
      "maccs": true,
      "atom_pair": true
    },
    "similarity": {
      "tanimoto": true,
      "dice": true,
      "cosine": true
    },
    "reactions": {
      "reaction_smarts": true,
      "reaction_prediction": false
    },
    "visualization": {
      "2d_image": true,
      "3d_structure": true,
      "svg_output": true
    }
  },
  "server": {
    "host": "localhost",
    "port": 8766,
    "protocol": "mcp"
  }
}
EOF

echo -e "${GREEN}✅ Configuration created${NC}"

# 4. Create startup script
echo ""
echo "4. Creating startup script..."

cat > "$MCP_DIR/start_rdkit_mcp.sh" << 'EOF'
#!/bin/bash
# Start RDKit MCP Server

MCP_DIR="$HOME/.chemagent/mcp_servers"
VENV_DIR="$MCP_DIR/rdkit_venv"

# Activate virtual environment if it exists
if [ -d "$VENV_DIR" ]; then
    source "$VENV_DIR/bin/activate"
fi

# Start the MCP server
echo "Starting RDKit MCP Server..."
python -m mcp_rdkit --config "$MCP_DIR/rdkit_mcp_config.json"
EOF

chmod +x "$MCP_DIR/start_rdkit_mcp.sh"
echo -e "${GREEN}✅ Startup script created${NC}"

# 5. Create systemd service (optional)
echo ""
echo "5. Creating systemd service (optional)..."

cat > "$MCP_DIR/rdkit-mcp.service" << EOF
[Unit]
Description=RDKit MCP Server
After=network.target

[Service]
Type=simple
User=$USER
WorkingDirectory=$MCP_DIR
ExecStart=$MCP_DIR/start_rdkit_mcp.sh
Restart=on-failure
RestartSec=10

[Install]
WantedBy=multi-user.target
EOF

echo -e "${GREEN}✅ Systemd service file created${NC}"
echo "To install as a service, run:"
echo "  sudo cp $MCP_DIR/rdkit-mcp.service /etc/systemd/system/"
echo "  sudo systemctl enable rdkit-mcp"
echo "  sudo systemctl start rdkit-mcp"

# 6. Test the installation
echo ""
echo "6. Testing RDKit MCP installation..."

python3 -c "
try:
    import mcp_rdkit
    print('✅ mcp-rdkit module imported successfully')
except ImportError as e:
    print('❌ Failed to import mcp-rdkit:', e)
"

python3 -c "
try:
    from rdkit import Chem
    mol = Chem.MolFromSmiles('CCO')
    if mol:
        print('✅ RDKit working: Ethanol molecule created')
except ImportError as e:
    print('❌ RDKit not working:', e)
"

# 7. Create example usage script
echo ""
echo "7. Creating example usage script..."

cat > "$MCP_DIR/test_rdkit_mcp.py" << 'EOF'
#!/usr/bin/env python3
"""
Test script for RDKit MCP Server
"""

import json
import requests

# Example MCP request
def test_rdkit_mcp():
    # Server endpoint
    url = "http://localhost:8766/mcp"
    
    # Example: Calculate properties for aspirin
    request = {
        "method": "calculate_properties",
        "params": {
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "properties": ["MW", "LogP", "TPSA", "HBD", "HBA"]
        }
    }
    
    try:
        response = requests.post(url, json=request)
        result = response.json()
        print("Aspirin properties:")
        print(json.dumps(result, indent=2))
    except Exception as e:
        print(f"Server not running or error: {e}")
        print("Start the server with: $MCP_DIR/start_rdkit_mcp.sh")

if __name__ == "__main__":
    test_rdkit_mcp()
EOF

chmod +x "$MCP_DIR/test_rdkit_mcp.py"
echo -e "${GREEN}✅ Example script created${NC}"

# 8. Update ChemAgent configuration
echo ""
echo "8. Updating ChemAgent configuration..."

CHEMAGENT_CONFIG="$HOME/.chemagent/mcp_config.json"
if [ -f "$CHEMAGENT_CONFIG" ]; then
    # Backup existing config
    cp "$CHEMAGENT_CONFIG" "$CHEMAGENT_CONFIG.bak"
fi

cat > "$CHEMAGENT_CONFIG" << EOF
{
  "mcp_servers": {
    "rdkit": {
      "type": "official",
      "name": "RDKit MCP Server",
      "command": "$MCP_DIR/start_rdkit_mcp.sh",
      "endpoint": "http://localhost:8766/mcp",
      "capabilities": [
        "molecular_analysis",
        "descriptor_calculation",
        "fingerprint_generation",
        "similarity_search",
        "visualization"
      ],
      "status": "installed"
    }
  }
}
EOF

echo -e "${GREEN}✅ ChemAgent configuration updated${NC}"

# Summary
echo ""
echo -e "${GREEN}================================================${NC}"
echo -e "${GREEN}  RDKit MCP Server Installation Complete!      ${NC}"
echo -e "${GREEN}================================================${NC}"
echo ""
echo "Installation Summary:"
echo "  • mcp-rdkit package: ✅ Installed"
echo "  • Virtual environment: $VENV_DIR"
echo "  • Configuration: $MCP_DIR/rdkit_mcp_config.json"
echo "  • Startup script: $MCP_DIR/start_rdkit_mcp.sh"
echo ""
echo "To start the RDKit MCP server:"
echo -e "  ${BLUE}$MCP_DIR/start_rdkit_mcp.sh${NC}"
echo ""
echo "To test the server:"
echo -e "  ${BLUE}$MCP_DIR/test_rdkit_mcp.py${NC}"
echo ""
echo "To use with ChemAgent:"
echo "  The server is now registered in ChemAgent's MCP configuration"
echo "  It will be automatically used for molecular operations"
echo ""
echo -e "${YELLOW}Note: This is the official RDKit MCP integration${NC}"
echo "For more information: https://pypi.org/project/mcp-rdkit/"
