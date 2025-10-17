#!/bin/bash

# ChemAgent One-Line Installer
# Inspired by SuperClaude_Framework
# Usage: curl -fsSL https://raw.githubusercontent.com/yourusername/chemagent/main/install.sh | bash

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
WHITE='\033[1;37m'
NC='\033[0m' # No Color

# ASCII Art Banner
echo -e "${CYAN}"
cat << "EOF"
   _____ _                    _                    _   
  / ____| |                  | |                  | |  
 | |    | |__   ___ _ __ ___ | |  __ _  ___ _ __ | |_ 
 | |    | '_ \ / _ \ '_ ` _ \| | / _` |/ _ \ '_ \| __|
 | |____| | | |  __/ | | | | | || (_| |  __/ | | | |_ 
  \_____|_| |_|\___|_| |_| |_|_| \__, |\___|_| |_|\__|
                                   __/ |               
  Chemistry Enhancement Package   |___/                
EOF
echo -e "${NC}"

echo -e "${WHITE}Welcome to ChemAgent Installer!${NC}"
echo -e "${YELLOW}Enhancing AI Assistants with Chemistry Expertise${NC}"
echo ""

# Detect OS
OS="Unknown"
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS="Linux"
elif [[ "$OSTYPE" == "darwin"* ]]; then
    OS="macOS"
elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    OS="Windows"
fi

echo -e "${BLUE}System Detection:${NC}"
echo "  • OS: $OS"
echo "  • User: $USER"
echo "  • Shell: $SHELL"
echo ""

# Check Python
echo -e "${BLUE}Checking Requirements:${NC}"
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version 2>&1 | grep -Po '(?<=Python )\d+\.\d+')
    echo -e "  ✅ Python $PYTHON_VERSION found"
    PYTHON_CMD="python3"
elif command -v python &> /dev/null; then
    PYTHON_VERSION=$(python --version 2>&1 | grep -Po '(?<=Python )\d+\.\d+')
    echo -e "  ✅ Python $PYTHON_VERSION found"
    PYTHON_CMD="python"
else
    echo -e "  ${RED}❌ Python not found. Please install Python 3.8+${NC}"
    exit 1
fi

# Check Git
if command -v git &> /dev/null; then
    echo -e "  ✅ Git found"
else
    echo -e "  ${YELLOW}⚠️  Git not found. Installing from archive...${NC}"
    USE_ARCHIVE=true
fi

# Detect AI Platforms
echo ""
echo -e "${BLUE}Detecting AI Platforms:${NC}"
PLATFORMS_FOUND=()

# Check Claude Code (Cursor)
if [ -d "$HOME/.cursor" ] || [ -d "$HOME/.config/Cursor" ] || [ -d "/Applications/Cursor.app" ]; then
    echo -e "  ✅ Claude Code (Cursor) detected"
    PLATFORMS_FOUND+=("claude-code")
else
    echo -e "  ⚠️  Claude Code (Cursor) not detected"
fi

# Check Gemini CLI
if [ -d "$HOME/.gemini" ] || command -v gemini &> /dev/null; then
    echo -e "  ✅ Gemini CLI detected"
    PLATFORMS_FOUND+=("gemini-cli")
else
    echo -e "  ⚠️  Gemini CLI not detected"
fi

# Check VS Code
if [ -d "$HOME/.vscode" ] || [ -d "$HOME/.config/Code" ] || [ -d "/Applications/Visual Studio Code.app" ]; then
    echo -e "  ℹ️  VS Code detected (can use with extensions)"
fi

# Installation directory
INSTALL_DIR="$HOME/.chemagent"
echo ""
echo -e "${BLUE}Installation Settings:${NC}"
echo "  • Install directory: $INSTALL_DIR"
echo "  • Platforms: ${PLATFORMS_FOUND[@]:-None detected}"
echo ""

# Ask for confirmation
read -p "$(echo -e ${YELLOW}Proceed with installation? [Y/n]:${NC} )" -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]] && [[ ! -z $REPLY ]]; then
    echo -e "${RED}Installation cancelled.${NC}"
    exit 1
fi

# Create installation directory
echo ""
echo -e "${BLUE}Installing ChemAgent...${NC}"
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

# Download ChemAgent
if [ "$USE_ARCHIVE" = true ]; then
    # Download as archive
    echo "  • Downloading ChemAgent archive..."
    curl -fsSL https://github.com/yourusername/chemagent/archive/main.tar.gz -o chemagent.tar.gz
    tar -xzf chemagent.tar.gz --strip-components=1
    rm chemagent.tar.gz
else
    # Clone repository
    echo "  • Cloning ChemAgent repository..."
    if [ -d ".git" ]; then
        git pull origin main
    else
        git clone https://github.com/yourusername/chemagent.git .
    fi
fi

# Install Python dependencies
echo "  • Installing Python dependencies..."
$PYTHON_CMD -m pip install --quiet --user -e . || {
    echo -e "${YELLOW}  ⚠️  Some dependencies failed to install${NC}"
}

# Platform-specific installation
echo ""
echo -e "${BLUE}Platform Configuration:${NC}"

# Claude Code installation
if [[ " ${PLATFORMS_FOUND[@]} " =~ " claude-code " ]]; then
    echo "  • Configuring Claude Code (Cursor)..."
    
    # Copy .cursorrules
    if [ -f ".cursorrules" ]; then
        cp .cursorrules "$HOME/"
        echo -e "    ${GREEN}✅ Cursor rules installed${NC}"
    fi
    
    # Create command symlinks
    mkdir -p "$HOME/.claude/commands"
    cp -r commands/*.md "$HOME/.claude/commands/" 2>/dev/null || true
    echo -e "    ${GREEN}✅ Commands installed${NC}"
    
    # Create role symlinks
    mkdir -p "$HOME/.claude/roles"
    cp -r roles/*.md "$HOME/.claude/roles/" 2>/dev/null || true
    echo -e "    ${GREEN}✅ Roles installed${NC}"
fi

# Gemini CLI installation
if [[ " ${PLATFORMS_FOUND[@]} " =~ " gemini-cli " ]]; then
    echo "  • Configuring Gemini CLI..."
    
    # Run Gemini installer
    if [ -f "install_gemini_simple.sh" ]; then
        chmod +x install_gemini_simple.sh
        ./install_gemini_simple.sh --quiet || {
            echo -e "    ${YELLOW}⚠️  Gemini configuration incomplete${NC}"
        }
    else
        # Manual configuration
        mkdir -p "$HOME/.gemini/rules"
        mkdir -p "$HOME/.gemini/prompts"
        [ -f "gemini_rules.md" ] && cp gemini_rules.md "$HOME/.gemini/rules/chemistry.md"
        echo -e "    ${GREEN}✅ Gemini rules installed${NC}"
    fi
fi

# Create command-line tool
echo ""
echo -e "${BLUE}Installing Command-line Tool:${NC}"
cat > "$HOME/.local/bin/chemagent" << 'SCRIPT_EOF'
#!/bin/bash
# ChemAgent CLI wrapper
CHEMAGENT_HOME="$HOME/.chemagent"
cd "$CHEMAGENT_HOME"
python3 chemagent_install.py "$@"
SCRIPT_EOF
chmod +x "$HOME/.local/bin/chemagent"
echo -e "  ${GREEN}✅ Command 'chemagent' installed${NC}"

# Add to PATH if needed
if [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
    echo ""
    echo -e "${YELLOW}Note: Add ~/.local/bin to your PATH:${NC}"
    echo '  export PATH="$HOME/.local/bin:$PATH"'
fi

# Shell integration
echo ""
echo -e "${BLUE}Shell Integration:${NC}"
SHELL_RC=""
if [ -f "$HOME/.bashrc" ]; then
    SHELL_RC="$HOME/.bashrc"
elif [ -f "$HOME/.zshrc" ]; then
    SHELL_RC="$HOME/.zshrc"
fi

if [ -n "$SHELL_RC" ]; then
    # Check if already added
    if ! grep -q "ChemAgent" "$SHELL_RC"; then
        echo "" >> "$SHELL_RC"
        echo "# ChemAgent" >> "$SHELL_RC"
        echo "export PATH=\"\$HOME/.local/bin:\$PATH\"" >> "$SHELL_RC"
        
        # Add Gemini integration if detected
        if [[ " ${PLATFORMS_FOUND[@]} " =~ " gemini-cli " ]]; then
            echo "[ -f \"\$HOME/.gemini/load_chemistry.sh\" ] && source \"\$HOME/.gemini/load_chemistry.sh\"" >> "$SHELL_RC"
        fi
        
        echo -e "  ${GREEN}✅ Added to $SHELL_RC${NC}"
    else
        echo -e "  ℹ️  ChemAgent already in $SHELL_RC"
    fi
fi

# Success message
echo ""
echo -e "${GREEN}╔══════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║         🎉 ChemAgent Installation Complete! 🎉          ║${NC}"
echo -e "${GREEN}╚══════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${WHITE}Quick Start:${NC}"
echo ""

if [[ " ${PLATFORMS_FOUND[@]} " =~ " claude-code " ]]; then
    echo -e "${CYAN}Claude Code (Cursor):${NC}"
    echo "  • Commands now available: cc-analyze, cc-synthesize, etc."
    echo "  • Roles available: @chemist, @drug-designer, etc."
    echo ""
fi

if [[ " ${PLATFORMS_FOUND[@]} " =~ " gemini-cli " ]]; then
    echo -e "${CYAN}Gemini CLI:${NC}"
    echo "  • Try: gemini \"chem:analyze aspirin\""
    echo "  • Or: gemini \"@drug-designer design a kinase inhibitor\""
    echo ""
fi

echo -e "${CYAN}Command-line:${NC}"
echo "  • chemagent --help        # Show all options"
echo "  • chemagent status        # Check installation"
echo "  • chemagent examples      # Install examples"
echo ""
echo -e "${YELLOW}Next steps:${NC}"
echo "  1. Reload your shell: source $SHELL_RC"
echo "  2. Test the installation: chemagent status"
echo "  3. Read the docs: https://github.com/yourusername/chemagent"
echo ""
echo -e "${MAGENTA}Thank you for installing ChemAgent!${NC}"
echo -e "${WHITE}Happy Chemistry Coding! 🧪⚗️🔬${NC}"