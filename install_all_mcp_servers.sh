#!/bin/bash

# Install and Configure All Chemistry MCP Servers
# 安装和配置所有化学MCP服务器

set -e

echo "================================================"
echo "  Installing Chemistry MCP Servers             "
echo "================================================"
echo ""

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m'

# MCP installation directory
MCP_DIR="$HOME/.chemagent/mcp_servers"
mkdir -p "$MCP_DIR"

# Log file
LOG_FILE="$MCP_DIR/installation.log"
echo "Installation started at $(date)" > "$LOG_FILE"

# Function to install a MCP server
install_mcp_server() {
    local name=$1
    local package=$2
    local port=$3
    local dependencies=$4
    
    echo ""
    echo -e "${BLUE}Installing $name MCP Server...${NC}"
    
    # Install dependencies
    if [ ! -z "$dependencies" ]; then
        echo "  Installing dependencies: $dependencies"
        pip install --user $dependencies >> "$LOG_FILE" 2>&1 || {
            echo -e "${YELLOW}  Warning: Some dependencies may have failed${NC}"
        }
    fi
    
    # Check if already installed
    if pip show $package &> /dev/null; then
        echo -e "${GREEN}  ✅ $package already installed${NC}"
    else
        # Try to install from PyPI
        echo "  Attempting to install $package from PyPI..."
        pip install --user $package >> "$LOG_FILE" 2>&1 || {
            echo -e "${YELLOW}  ⚠️  $package not on PyPI yet, installing from local${NC}"
            
            # Install from local if available
            if [ -d "mcp_servers/$package" ]; then
                echo "  Installing from local source..."
                cd "mcp_servers/$package"
                pip install --user -e . >> "$LOG_FILE" 2>&1
                cd - > /dev/null
                echo -e "${GREEN}  ✅ Installed from local source${NC}"
            else
                echo -e "${RED}  ❌ Could not install $package${NC}"
                return 1
            fi
        }
    fi
    
    # Create startup script
    cat > "$MCP_DIR/start_${name}_mcp.sh" << EOF
#!/bin/bash
# Start $name MCP Server
echo "Starting $name MCP Server on port $port..."
python -m ${package//-/_}_server --tcp $port
EOF
    chmod +x "$MCP_DIR/start_${name}_mcp.sh"
    
    # Create systemd service file
    cat > "$MCP_DIR/${name}-mcp.service" << EOF
[Unit]
Description=$name MCP Server
After=network.target

[Service]
Type=simple
User=$USER
WorkingDirectory=$MCP_DIR
ExecStart=/usr/bin/python3 -m ${package//-/_}_server
Restart=on-failure
RestartSec=10

[Install]
WantedBy=multi-user.target
EOF
    
    echo -e "${GREEN}  ✅ $name MCP Server configured${NC}"
    return 0
}

# Install base MCP server module
echo -e "${BLUE}Installing base MCP server module...${NC}"
if [ -f "mcp_servers/base_mcp_server.py" ]; then
    cp mcp_servers/base_mcp_server.py "$MCP_DIR/"
    echo -e "${GREEN}✅ Base MCP server module installed${NC}"
else
    echo -e "${YELLOW}⚠️  Base MCP server module not found${NC}"
fi

# 1. RDKit MCP (Official)
install_mcp_server "rdkit" "mcp-rdkit" "8766" "rdkit"

# 2. PubChem MCP (Our implementation)
install_mcp_server "pubchem" "mcp-pubchem" "8767" "pubchempy requests"

# 3. ChEMBL MCP (Our implementation)
install_mcp_server "chembl" "mcp-chembl" "8768" "chembl-webresource-client requests"

# 4. OpenBabel MCP (Our implementation)
install_mcp_server "openbabel" "mcp-openbabel" "8769" "openbabel-wheel"

# Create MCP manager script
echo ""
echo -e "${BLUE}Creating MCP Manager...${NC}"

cat > "$MCP_DIR/mcp_manager.py" << 'EOF'
#!/usr/bin/env python3
"""
MCP Server Manager for ChemAgent
Manages multiple MCP servers
"""

import subprocess
import psutil
import json
import sys
import time
from pathlib import Path

class MCPManager:
    def __init__(self):
        self.mcp_dir = Path.home() / ".chemagent" / "mcp_servers"
        self.config_file = self.mcp_dir / "mcp_config.json"
        self.servers = {
            "rdkit": {"port": 8766, "package": "mcp_rdkit"},
            "pubchem": {"port": 8767, "package": "mcp_pubchem_server"},
            "chembl": {"port": 8768, "package": "mcp_chembl_server"},
            "openbabel": {"port": 8769, "package": "mcp_openbabel_server"},
        }
    
    def start_server(self, name):
        """Start a specific MCP server"""
        if name not in self.servers:
            print(f"Unknown server: {name}")
            return False
        
        server = self.servers[name]
        
        # Check if already running
        if self.is_running(name):
            print(f"{name} MCP server is already running")
            return True
        
        # Start the server
        print(f"Starting {name} MCP server on port {server['port']}...")
        try:
            process = subprocess.Popen(
                [sys.executable, "-m", server['package'], "--tcp", str(server['port'])],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True
            )
            
            # Save PID
            pid_file = self.mcp_dir / f"{name}.pid"
            pid_file.write_text(str(process.pid))
            
            time.sleep(2)  # Wait for server to start
            
            if self.is_running(name):
                print(f"✅ {name} MCP server started successfully")
                return True
            else:
                print(f"❌ Failed to start {name} MCP server")
                return False
                
        except Exception as e:
            print(f"Error starting {name}: {e}")
            return False
    
    def stop_server(self, name):
        """Stop a specific MCP server"""
        pid_file = self.mcp_dir / f"{name}.pid"
        
        if not pid_file.exists():
            print(f"{name} MCP server is not running")
            return True
        
        try:
            pid = int(pid_file.read_text())
            process = psutil.Process(pid)
            process.terminate()
            process.wait(timeout=5)
            pid_file.unlink()
            print(f"✅ {name} MCP server stopped")
            return True
        except psutil.NoSuchProcess:
            pid_file.unlink()
            print(f"{name} MCP server was not running")
            return True
        except Exception as e:
            print(f"Error stopping {name}: {e}")
            return False
    
    def is_running(self, name):
        """Check if a server is running"""
        pid_file = self.mcp_dir / f"{name}.pid"
        
        if not pid_file.exists():
            return False
        
        try:
            pid = int(pid_file.read_text())
            process = psutil.Process(pid)
            return process.is_running()
        except:
            return False
    
    def status(self):
        """Show status of all servers"""
        print("\nMCP Server Status:")
        print("-" * 40)
        
        for name in self.servers:
            if self.is_running(name):
                port = self.servers[name]['port']
                print(f"  {name:12} ✅ Running (port {port})")
            else:
                print(f"  {name:12} ❌ Stopped")
        print("-" * 40)
    
    def start_all(self):
        """Start all MCP servers"""
        print("Starting all MCP servers...")
        for name in self.servers:
            self.start_server(name)
    
    def stop_all(self):
        """Stop all MCP servers"""
        print("Stopping all MCP servers...")
        for name in self.servers:
            self.stop_server(name)

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="MCP Server Manager")
    parser.add_argument("action", choices=["start", "stop", "restart", "status", "start-all", "stop-all"],
                       help="Action to perform")
    parser.add_argument("server", nargs="?", help="Server name (for start/stop/restart)")
    
    args = parser.parse_args()
    
    manager = MCPManager()
    
    if args.action == "status":
        manager.status()
    elif args.action == "start-all":
        manager.start_all()
    elif args.action == "stop-all":
        manager.stop_all()
    elif args.action in ["start", "stop", "restart"]:
        if not args.server:
            print("Server name required for this action")
            sys.exit(1)
        
        if args.action == "start":
            manager.start_server(args.server)
        elif args.action == "stop":
            manager.stop_server(args.server)
        elif args.action == "restart":
            manager.stop_server(args.server)
            time.sleep(1)
            manager.start_server(args.server)

if __name__ == "__main__":
    # Install psutil if needed
    try:
        import psutil
    except ImportError:
        print("Installing psutil...")
        subprocess.run([sys.executable, "-m", "pip", "install", "--user", "psutil"])
        import psutil
    
    main()
EOF

chmod +x "$MCP_DIR/mcp_manager.py"

# Create convenience script
cat > "$MCP_DIR/mcp" << 'EOF'
#!/bin/bash
python3 ~/.chemagent/mcp_servers/mcp_manager.py "$@"
EOF
chmod +x "$MCP_DIR/mcp"

# Add to PATH
echo ""
echo -e "${BLUE}Adding MCP manager to PATH...${NC}"
if ! grep -q "chemagent/mcp_servers" ~/.bashrc; then
    echo 'export PATH="$HOME/.chemagent/mcp_servers:$PATH"' >> ~/.bashrc
    echo -e "${GREEN}✅ Added to ~/.bashrc${NC}"
fi

# Create Claude Desktop configuration
echo ""
echo -e "${BLUE}Creating Claude Desktop configuration...${NC}"

CLAUDE_CONFIG="$HOME/.config/claude/mcp_config.json"
mkdir -p "$(dirname "$CLAUDE_CONFIG")"

cat > "$CLAUDE_CONFIG" << EOF
{
  "mcpServers": {
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
  }
}
EOF

echo -e "${GREEN}✅ Claude Desktop configuration created${NC}"

# Summary
echo ""
echo "================================================"
echo -e "${GREEN}  MCP Server Installation Complete!${NC}"
echo "================================================"
echo ""
echo "Installed MCP Servers:"
echo "  • RDKit MCP (port 8766)"
echo "  • PubChem MCP (port 8767)"
echo "  • ChEMBL MCP (port 8768)"
echo "  • OpenBabel MCP (port 8769)"
echo ""
echo "Usage:"
echo "  mcp status           # Check server status"
echo "  mcp start-all        # Start all servers"
echo "  mcp stop-all         # Stop all servers"
echo "  mcp start rdkit      # Start specific server"
echo "  mcp stop pubchem     # Stop specific server"
echo ""
echo "To apply PATH changes:"
echo -e "  ${BLUE}source ~/.bashrc${NC}"
echo ""
echo "Logs: $LOG_FILE"
echo ""
