#!/bin/bash
# ChemAgent Command Aliases
# Inspired by SuperClaude_Framework's alias system
# Source this file in your shell: source chemagent_aliases.sh

# ============================================
# Core ChemAgent Commands - Short Aliases
# ============================================

# Analysis Commands
alias ca='chemagent cc-analyze'          # Analyze molecule
alias ccomp='chemagent cc-compare'       # Compare molecules
alias csearch='chemagent cc-search'      # Search databases

# Design Commands  
alias cdesign='chemagent cc-design'      # Design molecule
alias copt='chemagent cc-optimize'       # Optimize properties
alias csynth='chemagent cc-synthesize'   # Plan synthesis

# Prediction Commands
alias cpred='chemagent cc-predict'       # Predict properties/reactions
alias csim='chemagent cc-simulate'       # Run simulations

# Workflow Commands
alias cwork='chemagent cc-workflow'      # Run workflow
alias cbatch='chemagent cc-batch'        # Batch process
alias creport='chemagent cc-report'      # Generate report

# Helper Commands
alias cexplain='chemagent cc-explain'    # Explain concept
alias csuggest='chemagent cc-suggest'    # Get suggestions
alias ccheck='chemagent cc-check'        # Run checks
alias chelp='chemagent cc-help'          # Get help

# ============================================
# Convenience Functions
# ============================================

# Quick molecule analysis with visualization
mol() {
    if [ -z "$1" ]; then
        echo "Usage: mol <SMILES or molecule name>"
        echo "Example: mol aspirin"
        echo "Example: mol 'CC(=O)OC1=CC=CC=C1C(=O)O'"
        return 1
    fi
    chemagent cc-analyze "$1" --visualize
}

# Quick synthesis planning
synth() {
    if [ -z "$1" ]; then
        echo "Usage: synth <target molecule>"
        echo "Example: synth ibuprofen"
        return 1
    fi
    chemagent cc-synthesize "$1" --detailed
}

# Drug-likeness check
druglike() {
    if [ -z "$1" ]; then
        echo "Usage: druglike <molecule>"
        echo "Example: druglike caffeine"
        return 1
    fi
    chemagent cc-analyze "$1" --druglike --lipinski
}

# Safety assessment
safety() {
    if [ -z "$1" ]; then
        echo "Usage: safety <molecule>"
        echo "Example: safety benzene"
        return 1
    fi
    chemagent cc-check safety "$1" --comprehensive
}

# Patent check
patent() {
    if [ -z "$1" ]; then
        echo "Usage: patent <molecule or CAS>"
        echo "Example: patent 'aspirin'"
        echo "Example: patent '50-78-2'"
        return 1
    fi
    chemagent cc-check patent "$1"
}

# Batch analysis from file
batch_analyze() {
    if [ -z "$1" ]; then
        echo "Usage: batch_analyze <file>"
        echo "Example: batch_analyze molecules.csv"
        return 1
    fi
    chemagent cc-batch "$1" --operation analyze --output "${1%.csv}_results.csv"
}

# Compare two molecules
compare() {
    if [ -z "$1" ] || [ -z "$2" ]; then
        echo "Usage: compare <molecule1> <molecule2>"
        echo "Example: compare aspirin ibuprofen"
        return 1
    fi
    chemagent cc-compare "$1" "$2" --detailed
}

# Optimize for specific property
optimize() {
    if [ -z "$1" ] || [ -z "$2" ]; then
        echo "Usage: optimize <molecule> <property>"
        echo "Example: optimize 'CCO' solubility"
        echo "Properties: solubility, permeability, stability, potency"
        return 1
    fi
    chemagent cc-optimize "$1" --target "$2"
}

# ============================================
# Role-based Commands
# ============================================

# Drug Designer Role
drug_designer() {
    echo "Activating Drug Designer role..."
    chemagent --role drug-designer "$@"
}

# Synthesis Expert Role
synthesist() {
    echo "Activating Synthesis Expert role..."
    chemagent --role synthesist "$@"
}

# Safety Expert Role
safety_expert() {
    echo "Activating Safety Expert role..."
    chemagent --role safety-expert "$@"
}

# Data Analyst Role
data_analyst() {
    echo "Activating Data Analyst role..."
    chemagent --role data-analyst "$@"
}

# ============================================
# Workflow Shortcuts
# ============================================

# Drug Discovery Workflow
drug_discovery() {
    local target="${1:-kinase}"
    echo "Starting drug discovery workflow for: $target"
    chemagent cc-workflow drug-discovery --target "$target"
}

# Lead Optimization Workflow
lead_opt() {
    if [ -z "$1" ]; then
        echo "Usage: lead_opt <lead_molecule>"
        return 1
    fi
    chemagent cc-workflow lead-optimization --molecule "$1"
}

# ADMET Profiling
admet() {
    if [ -z "$1" ]; then
        echo "Usage: admet <molecule>"
        return 1
    fi
    chemagent cc-workflow admet-profile --molecule "$1"
}

# ============================================
# Gemini CLI Integration (if available)
# ============================================

if command -v gemini &> /dev/null; then
    # Gemini shortcuts
    alias gchem='gemini "chem:analyze"'
    alias gsynth='gemini "chem:synthesize"'
    alias gdrug='gemini "@drug-designer"'
    
    # Gemini chemistry function
    gem_chem() {
        if [ -z "$1" ]; then
            echo "Usage: gem_chem <query>"
            echo "Example: gem_chem 'analyze aspirin'"
            return 1
        fi
        gemini "chem: $1"
    }
fi

# ============================================
# Project Templates
# ============================================

# Create new drug discovery project
new_drug_project() {
    local name="${1:-drug_discovery_project}"
    echo "Creating new drug discovery project: $name"
    mkdir -p "$name"/{data,results,docs,structures}
    
    # Create README
    cat > "$name/README.md" << EOF
# $name

## Drug Discovery Project

### Structure
- data/       - Input data and datasets
- results/    - Analysis results
- docs/       - Documentation
- structures/ - Molecular structures

### Workflow
1. Target identification
2. Lead discovery
3. Lead optimization
4. ADMET profiling
5. Synthesis planning
EOF
    
    echo "Project created: $name/"
}

# Create materials design project
new_materials_project() {
    local name="${1:-materials_project}"
    echo "Creating new materials design project: $name"
    mkdir -p "$name"/{structures,properties,simulations,results}
    
    cat > "$name/README.md" << EOF
# $name

## Materials Design Project

### Structure
- structures/  - Material structures
- properties/  - Property calculations
- simulations/ - MD/DFT simulations
- results/     - Analysis results
EOF
    
    echo "Project created: $name/"
}

# ============================================
# Utility Functions
# ============================================

# Show all ChemAgent commands
chem_commands() {
    echo "ChemAgent Commands:"
    echo ""
    echo "Quick Aliases:"
    echo "  ca        - Analyze molecule"
    echo "  csynth    - Plan synthesis"
    echo "  cdesign   - Design molecule"
    echo "  copt      - Optimize properties"
    echo "  cpred     - Predict properties"
    echo "  csearch   - Search databases"
    echo "  ccheck    - Run checks"
    echo "  chelp     - Get help"
    echo ""
    echo "Functions:"
    echo "  mol       - Quick analysis with visualization"
    echo "  synth     - Detailed synthesis planning"
    echo "  druglike  - Check drug-likeness"
    echo "  safety    - Safety assessment"
    echo "  patent    - Patent check"
    echo "  compare   - Compare molecules"
    echo "  optimize  - Optimize for property"
    echo ""
    echo "Workflows:"
    echo "  drug_discovery - Start drug discovery"
    echo "  lead_opt      - Optimize lead compound"
    echo "  admet         - ADMET profiling"
    echo ""
    echo "Projects:"
    echo "  new_drug_project      - Create drug project"
    echo "  new_materials_project - Create materials project"
}

# Show ChemAgent status
chem_status() {
    chemagent status
}

# Update ChemAgent
chem_update() {
    chemagent update
}

# ============================================
# Welcome Message
# ============================================

echo "🧪 ChemAgent aliases loaded!"
echo "Type 'chem_commands' to see all available commands"
echo "Type 'chelp' for detailed help"
