#!/bin/bash

# ChemAgent MCP Status Checker
# 检查各种化学工具的 MCP 支持状态

set -e

echo "================================================"
echo "  ChemAgent MCP Availability Status Check      "
echo "================================================"
echo ""

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# Function to check if a Python package is installed
check_python_package() {
    local package=$1
    if python3 -c "import $package" &> /dev/null; then
        return 0
    else
        return 1
    fi
}

# Function to check if a command exists
check_command() {
    local cmd=$1
    if command -v $cmd &> /dev/null; then
        return 0
    else
        return 1
    fi
}

echo -e "${BLUE}Checking MCP Support Status...${NC}"
echo ""

# 1. RDKit MCP
echo "1. RDKit MCP Server"
if check_python_package "mcp_rdkit"; then
    echo -e "   ${GREEN}✅ mcp-rdkit installed (OFFICIAL MCP)${NC}"
    echo "   Status: Full MCP support available"
    echo "   Run: python -m mcp_rdkit"
else
    echo -e "   ${YELLOW}⚠️  mcp-rdkit not installed${NC}"
    echo "   Install: pip install mcp-rdkit"
fi

if check_python_package "rdkit"; then
    echo -e "   ${GREEN}✅ RDKit library installed (fallback)${NC}"
else
    echo -e "   ${RED}❌ RDKit library not installed${NC}"
    echo "   Install: pip install rdkit-pypi"
fi
echo ""

# 2. PubChem
echo "2. PubChem MCP Server"
echo -e "   ${RED}❌ No official MCP available${NC}"
if check_python_package "pubchempy"; then
    echo -e "   ${GREEN}✅ PubChemPy installed (fallback)${NC}"
else
    echo -e "   ${YELLOW}⚠️  PubChemPy not installed${NC}"
    echo "   Install: pip install pubchempy"
fi
echo ""

# 3. ChEMBL
echo "3. ChEMBL MCP Server"
echo -e "   ${RED}❌ No official MCP available${NC}"
if check_python_package "chembl_webresource_client"; then
    echo -e "   ${GREEN}✅ ChEMBL client installed (fallback)${NC}"
else
    echo -e "   ${YELLOW}⚠️  ChEMBL client not installed${NC}"
    echo "   Install: pip install chembl-webresource-client"
fi
echo ""

# 4. OpenBabel
echo "4. OpenBabel MCP Server"
echo -e "   ${RED}❌ No official MCP available${NC}"
if check_python_package "openbabel"; then
    echo -e "   ${GREEN}✅ OpenBabel Python installed (fallback)${NC}"
elif check_command "obabel"; then
    echo -e "   ${GREEN}✅ OpenBabel CLI installed (fallback)${NC}"
else
    echo -e "   ${YELLOW}⚠️  OpenBabel not installed${NC}"
    echo "   Install: pip install openbabel-wheel"
fi
echo ""

# 5. USPTO Patents
echo "5. USPTO Patents MCP Server"
echo -e "   ${RED}❌ No official MCP available${NC}"
echo -e "   ${BLUE}ℹ️  Using USPTO API directly (fallback)${NC}"
echo ""

# 6. Chemical Suppliers
echo "6. Chemical Suppliers MCP"
echo -e "   ${RED}❌ No official MCP available${NC}"
echo -e "   ${BLUE}ℹ️  Using supplier APIs directly (fallback)${NC}"
echo ""

# 7. Literature Databases
echo "7. Literature MCP (PubMed, CrossRef)"
echo -e "   ${RED}❌ No official MCP available${NC}"
if check_python_package "biopython"; then
    echo -e "   ${GREEN}✅ BioPython installed (PubMed fallback)${NC}"
else
    echo -e "   ${YELLOW}⚠️  BioPython not installed${NC}"
    echo "   Install: pip install biopython"
fi
echo ""

# 8. Computational Chemistry
echo "8. Computational Chemistry MCP"
echo -e "   ${RED}❌ No official MCP for Gaussian/ORCA/GAMESS${NC}"
if check_command "orca"; then
    echo -e "   ${GREEN}✅ ORCA installed (CLI fallback)${NC}"
elif check_command "g16"; then
    echo -e "   ${GREEN}✅ Gaussian installed (CLI fallback)${NC}"
else
    echo -e "   ${YELLOW}⚠️  No QM software detected${NC}"
fi
echo ""

# 9. Molecular Docking
echo "9. Molecular Docking MCP"
echo -e "   ${RED}❌ No official MCP for AutoDock${NC}"
if check_command "vina"; then
    echo -e "   ${GREEN}✅ AutoDock Vina installed (CLI fallback)${NC}"
elif check_python_package "meeko"; then
    echo -e "   ${GREEN}✅ Meeko (AutoDock tools) installed${NC}"
else
    echo -e "   ${YELLOW}⚠️  No docking software detected${NC}"
fi
echo ""

# 10. Visualization
echo "10. Molecular Visualization MCP"
echo -e "   ${RED}❌ No official MCP for PyMOL${NC}"
if check_python_package "pymol"; then
    echo -e "   ${GREEN}✅ PyMOL Python installed (fallback)${NC}"
elif check_command "pymol"; then
    echo -e "   ${GREEN}✅ PyMOL installed (CLI fallback)${NC}"
else
    echo -e "   ${YELLOW}⚠️  PyMOL not detected${NC}"
fi
echo ""

# Summary
echo "================================================"
echo -e "${BLUE}Summary:${NC}"
echo ""
echo -e "${GREEN}Official MCP Support:${NC}"
echo "  • RDKit: ✅ (mcp-rdkit package available)"
echo ""
echo -e "${YELLOW}Using Fallback Implementations:${NC}"
echo "  • PubChem: PubChemPy"
echo "  • ChEMBL: Web Services API"
echo "  • OpenBabel: Python bindings"
echo "  • Patents: USPTO API"
echo "  • Suppliers: Direct APIs"
echo "  • Literature: PubMed/CrossRef APIs"
echo "  • QM Software: CLI interfaces"
echo "  • Docking: AutoDock CLI"
echo "  • Visualization: PyMOL API"
echo ""
echo -e "${BLUE}ChemAgent Status:${NC}"
echo "  ChemAgent's orchestrator intelligently routes requests:"
echo "  1. To official MCP servers when available (RDKit)"
echo "  2. To fallback implementations for other tools"
echo "  3. Users don't need to worry about the underlying implementation"
echo ""
echo "================================================"
echo ""
echo "To install missing components:"
echo -e "  ${BLUE}python chemagent_install.py --interactive${NC}"
echo ""
echo "To install RDKit MCP specifically:"
echo -e "  ${BLUE}python chemagent_install.py mcp${NC}"
echo ""
