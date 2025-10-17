#!/usr/bin/env python3
"""
MCP Orchestrator for ChemAgent
Orchestrates calls to official MCP servers and provides custom tools
"""

import asyncio
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from chemagent.mcp_tools.orchestrator import MCPOrchestrator


async def main():
    """Start the MCP orchestrator"""
    print("🚀 Starting ChemAgent MCP Orchestrator...")
    print("=" * 50)
    
    orchestrator = MCPOrchestrator()
    
    print("✅ Orchestrator initialized")
    print("\n📋 Orchestrating official MCP servers:")
    for name, config in orchestrator.official_mcps.items():
        print(f"  - {name}: {config['endpoint']}")
    
    print("\n🔧 Custom tools for missing functionality:")
    print("  - workflow_orchestration")
    print("  - best_practice_guidance")
    print("  - result_interpretation")
    print("  - literature_search (no official MCP)")
    print("  - safety_assessment (custom)")
    
    print("\n🎯 Available commands:")
    print("  - cc-analyze, cc-synthesize, cc-design")
    print("  - cc-optimize, cc-predict, cc-simulate")
    print("  - cc-workflow, cc-batch, cc-report")
    print("  - cc-explain, cc-suggest, cc-check, cc-help")
    
    print("\n🔌 Ready to orchestrate chemistry workflows")
    print("Press Ctrl+C to stop\n")
    
    # In production, this would start an actual MCP server
    # For now, just keep running
    try:
        while True:
            await asyncio.sleep(1)
    except KeyboardInterrupt:
        print("\n👋 Orchestrator stopped")
        

if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\nShutdown complete")
        sys.exit(0)
