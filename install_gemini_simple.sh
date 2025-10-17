#!/bin/bash

# ChemAgent for Gemini CLI - Simple Prompt-based Installation
# 基于提示词的简单安装脚本

set -e

echo "================================================"
echo "  ChemAgent for Gemini CLI - Prompt Installer  "
echo "================================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Gemini CLI paths
GEMINI_HOME="$HOME/.gemini"
GEMINI_PROMPTS="$GEMINI_HOME/prompts"
GEMINI_RULES="$GEMINI_HOME/rules"
GEMINI_CONFIG="$GEMINI_HOME/config"

# Create directories
echo -e "${YELLOW}Creating Gemini directories...${NC}"
mkdir -p "$GEMINI_PROMPTS"
mkdir -p "$GEMINI_RULES" 
mkdir -p "$GEMINI_CONFIG"

# Copy the main rules file
echo -e "${YELLOW}Installing chemistry rules...${NC}"
cp gemini_rules.md "$GEMINI_RULES/chemistry.md"

# Create chemistry prompts
echo -e "${YELLOW}Installing chemistry prompts...${NC}"

# Molecular Analysis Prompt
cat > "$GEMINI_PROMPTS/chem_analyze.md" << 'EOF'
When analyzing a molecule, provide:
1. Validated structure (SMILES/InChI)
2. Key molecular properties (MW, LogP, TPSA, etc.)
3. Drug-likeness assessment (Lipinski, QED)
4. Potential biological activities
5. Synthetic accessibility
6. Safety considerations
Format the output clearly with sections.
EOF

# Synthesis Planning Prompt
cat > "$GEMINI_PROMPTS/chem_synthesize.md" << 'EOF'
For synthesis planning:
1. Perform retrosynthetic analysis
2. Identify key disconnections
3. Suggest commercially available starting materials
4. Provide step-by-step forward synthesis
5. Include reaction conditions and expected yields
6. Consider green chemistry principles
7. Suggest alternative routes
EOF

# Drug Design Prompt
cat > "$GEMINI_PROMPTS/chem_drug_design.md" << 'EOF'
When designing drug molecules:
1. Consider the target and mechanism
2. Apply Lipinski's Rule of Five
3. Optimize ADMET properties
4. Check for structural alerts (PAINS)
5. Ensure synthetic feasibility
6. Consider intellectual property
7. Suggest analogs for SAR studies
EOF

# Safety Assessment Prompt
cat > "$GEMINI_PROMPTS/chem_safety.md" << 'EOF'
For safety assessment:
1. Check for known toxicophores
2. Predict acute and chronic toxicity
3. Assess environmental impact
4. Review regulatory status
5. Identify hazard classifications
6. Suggest safer alternatives
7. Provide handling precautions
EOF

# Create command aliases configuration
echo -e "${YELLOW}Creating command shortcuts...${NC}"
cat > "$GEMINI_CONFIG/chemistry_aliases.txt" << 'EOF'
# ChemAgent Command Shortcuts for Gemini CLI
# Use these shortcuts in your prompts:

chem:analyze <molecule> = Analyze molecular structure and properties
chem:synthesize <target> = Plan synthesis route
chem:predict <reaction> = Predict reaction outcomes
chem:optimize <molecule> <property> = Optimize for specific property
chem:search <query> = Search chemical databases
chem:design <requirements> = Design new molecules
chem:safety <molecule> = Assess safety and toxicity
chem:compare <mol1> <mol2> = Compare molecules
chem:explain <concept> = Explain chemistry concepts
chem:batch <file> = Process multiple molecules

# Role shortcuts:
@chemist = General chemistry expert
@drug-designer = Drug design specialist
@synthesist = Synthesis expert
@safety-expert = Safety and regulatory expert
@data-analyst = Chemical data analysis
EOF

# Create a Gemini CLI configuration snippet
echo -e "${YELLOW}Creating Gemini configuration...${NC}"
cat > "$GEMINI_CONFIG/chemistry_config.json" << 'EOF'
{
  "chemistry": {
    "enabled": true,
    "default_role": "chemist",
    "auto_validate_structures": true,
    "include_3d": false,
    "citation_style": "ACS",
    "safety_warnings": true,
    "commands": {
      "prefix": "chem:",
      "enabled": true
    },
    "prompts": {
      "auto_load": true,
      "directory": "~/.gemini/prompts"
    },
    "rules": {
      "auto_load": true,
      "directory": "~/.gemini/rules"
    }
  }
}
EOF

# Create integration script
echo -e "${YELLOW}Creating integration script...${NC}"
cat > "$GEMINI_HOME/load_chemistry.sh" << 'EOF'
#!/bin/bash
# Load ChemAgent for Gemini CLI

# Export chemistry rules path
export GEMINI_CHEMISTRY_RULES="$HOME/.gemini/rules/chemistry.md"
export GEMINI_CHEMISTRY_PROMPTS="$HOME/.gemini/prompts"

# Function to invoke Gemini with chemistry context
chem() {
    if [ -z "$1" ]; then
        echo "Usage: chem <command> [args]"
        echo "Commands: analyze, synthesize, predict, optimize, search, design, safety"
        return 1
    fi
    
    local cmd=$1
    shift
    
    # Build the prompt with chemistry context
    local prompt="Using chemistry expertise, $cmd: $*"
    
    # Call gemini with the chemistry rules loaded
    gemini --system-prompt "@$GEMINI_CHEMISTRY_RULES" "$prompt"
}

# Aliases for common chemistry tasks
alias chem-analyze='chem analyze'
alias chem-synth='chem synthesize'
alias chem-predict='chem predict'
alias chem-safety='chem safety'
alias chem-design='chem design'

echo "✅ ChemAgent for Gemini CLI loaded"
echo "Type 'chem help' for available commands"
EOF

chmod +x "$GEMINI_HOME/load_chemistry.sh"

# Create example usage file
echo -e "${YELLOW}Creating examples...${NC}"
cat > "$GEMINI_HOME/chemistry_examples.md" << 'EOF'
# ChemAgent for Gemini CLI - Examples

## Basic Usage

### Analyze a molecule
```bash
gemini "chem:analyze aspirin"
gemini "analyze the structure CC(=O)OC1=CC=CC=C1C(=O)O"
```

### Plan synthesis
```bash
gemini "chem:synthesize ibuprofen"
gemini "how do I synthesize paracetamol from phenol?"
```

### Predict properties
```bash
gemini "predict the logP of caffeine"
gemini "chem:predict reaction between benzene and bromine"
```

### Drug design
```bash
gemini "@drug-designer design a PDE5 inhibitor"
gemini "optimize this molecule for better solubility: [SMILES]"
```

### Safety assessment
```bash
gemini "@safety-expert assess the toxicity of benzene"
gemini "chem:safety check environmental impact of DDT"
```

## Advanced Usage

### Batch processing
```bash
gemini "analyze these molecules: aspirin, ibuprofen, paracetamol"
```

### Comparative analysis
```bash
gemini "chem:compare morphine vs codeine"
```

### Workflow
```bash
gemini "run a drug discovery workflow for COVID-19 protease inhibitors"
```

## Using with files
```bash
gemini "analyze the molecules in molecules.csv"
gemini "generate a safety report for compounds.sdf"
```

## Role-specific queries
```bash
gemini "@synthesist suggest a route to synthesize taxol"
gemini "@data-analyst analyze SAR for this series of compounds"
```
EOF

# Final setup message
echo ""
echo -e "${GREEN}✅ ChemAgent for Gemini CLI installed successfully!${NC}"
echo ""
echo "📝 Installation Summary:"
echo "  • Chemistry rules: $GEMINI_RULES/chemistry.md"
echo "  • Prompts: $GEMINI_PROMPTS/chem_*.md"
echo "  • Configuration: $GEMINI_CONFIG/chemistry_config.json"
echo "  • Examples: $GEMINI_HOME/chemistry_examples.md"
echo ""
echo "🚀 To activate ChemAgent:"
echo "  1. Add to your shell profile (~/.bashrc or ~/.zshrc):"
echo "     source $GEMINI_HOME/load_chemistry.sh"
echo ""
echo "  2. Reload your shell or run:"
echo "     source $GEMINI_HOME/load_chemistry.sh"
echo ""
echo "  3. Test with:"
echo "     gemini \"chem:analyze water\""
echo "     chem analyze aspirin"
echo ""
echo "📚 See $GEMINI_HOME/chemistry_examples.md for more examples"
echo ""

# Ask if user wants to add to shell profile
read -p "Would you like to add ChemAgent to your shell profile? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    SHELL_RC=""
    if [ -f "$HOME/.bashrc" ]; then
        SHELL_RC="$HOME/.bashrc"
    elif [ -f "$HOME/.zshrc" ]; then
        SHELL_RC="$HOME/.zshrc"
    fi
    
    if [ -n "$SHELL_RC" ]; then
        echo "" >> "$SHELL_RC"
        echo "# ChemAgent for Gemini CLI" >> "$SHELL_RC"
        echo "source $GEMINI_HOME/load_chemistry.sh" >> "$SHELL_RC"
        echo -e "${GREEN}✅ Added to $SHELL_RC${NC}"
        echo "   Reload your shell or run: source $SHELL_RC"
    fi
fi

echo ""
echo "🎉 Installation complete!"
