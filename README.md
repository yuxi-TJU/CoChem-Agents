ChemAgent – Chemistry Enhancement Pack for Claude Code & Gemini CLI
Overview



ChemAgent is a chemistry-focused enhancement package designed for Claude Code (Cursor) and Gemini CLI. It extends these powerful AI programming assistants with pre-configured commands, sub-agents, specialized tools, and optimized prompt templates—empowering them to perform advanced tasks in chemical research, drug design, and materials science.




Inspired by the SuperClaude Framework, ChemAgent follows a modular “plug-and-play” model, offering seamless chemical intelligence augmentation for AI coding environments.

Core Features
Claude Code (Cursor) Enhancements
📝 cc-series commands: Specialized chemistry commands such as cc-analyze, cc-synthesize, cc-predict
🎭 Sub-agents: Domain experts like @organic-chemist, @drug-designer, etc.
🔧 Auto-context loading: .cursorrules auto-loads chemical rules and context
🧪 Integrated tools: Seamless access to RDKit, PubChem, ChEMBL, and more
Gemini CLI Enhancements
🌟 chem-series commands: e.g., gemini chem:analyze, gemini chem:synthesize
📋 Prompt-based integration: Lightweight rule and prompt system
🖼️ Multimodal capability: Use Gemini’s vision model to interpret molecular images
📦 Batch processing: Run bulk operations with gemini chem:batch
🎭 Chemical roles: Invoke experts like @drug-designer or @synthesist
General Features
🚀 One-click installation: Auto-detects and installs to Claude Code or Gemini CLI
🧪 Professional tools integration: RDKit, PubChem, ChemSpace, OpenBabel
📚 Optimized prompts: Ready-made templates for drug design, materials, synthesis
🔬 MCP support: Extend capabilities through the Model Context Protocol
Quick Start
Installation
🚀 Recommended: One-click Installation
curl -fsSL https://raw.githubusercontent.com/dazhaolang/ai-chemkit/main/install.sh | bash

Manual Installation
git clone https://github.com/dazhaolang/ai-chemkit.git
cd ai-chemkit
python chemagent_install.py

Install RDKit MCP Server (Official)
python chemagent_install.py mcp
# or
./install_rdkit_mcp.sh

Installation Modes



ChemAgent supports multiple installation profiles:

Quick Install – Auto-detects platform and installs all features
python chemagent_install.py
or chemagent install

Interactive Mode – Choose platforms, tools, and components
python chemagent_install.py --interactive

Minimal Mode – Lightweight setup for limited environments
python chemagent_install.py --minimal

Developer Mode – Includes dev utilities (testing, linting, type checks)
python chemagent_install.py --profile developer

Other Commands
python chemagent_install.py status        # Check installation status
python chemagent_install.py examples      # Install sample files
python chemagent_install.py update        # Update to latest version
python chemagent_install.py --yes --quiet # Silent install
python chemagent_install.py --help        # View all options

Command System



ChemAgent follows the Markdown-based command definition style of SuperClaude, offering flexibility and simplicity.




Command locations:

Project-level: .claude/commands/
User-level: ~/.claude/commands/
System-level: Built-in defaults



Example custom command:

mkdir -p .claude/commands
cat > .claude/commands/my-analysis.md << EOF
---
description: My molecule analysis workflow
tools: [read_file, web_search]
---
1. Validate molecular structure  
2. Compute basic properties  
3. Predict ADMET  
4. Generate a summary report
EOF

Usage Examples
Claude Code (Cursor)
cc-analyze aspirin
cc-synthesize "CC(=O)OC1=CC=CC=C1C(=O)O"
@chemist Design an eco-friendly aspirin synthesis route
cc-workflow drug-discovery target.pdb

Gemini CLI
gemini chem:analyze "CC(=O)OC1=CC=CC=C1C(=O)O" --properties
gemini chem:batch molecules.csv --operation analyze
gemini chem:analyze molecule.png --image-to-structure

# With aliases
chem aspirin
synthesize ibuprofen

Project Structure
chemagent/
├── commands/            # Markdown command definitions
│   ├── cc-analyze.md
│   ├── cc-synthesize.md
│   └── ...
├── chemagent/
│   ├── enhancers/
│   │   ├── claude_enhancer.py
│   │   └── gemini_enhancer.py
│   ├── commands/
│   │   └── loader.py
│   ├── roles/           # Expert sub-agents
│   │   ├── chemist.py
│   │   ├── drug_designer.py
│   │   ├── synthesist.py
│   │   ├── safety_expert.py
│   │   └── data_analyst.py
│   ├── tools/           # Chemistry utilities
│   │   ├── patent_search.py
│   │   ├── literature_enhanced.py
│   │   ├── safety_assessment.py
│   │   └── price_lookup.py
│   └── mcp_tools/
│       └── orchestrator.py
├── .cursorrules
└── install.py

Command Reference
Claude Code (cc- series)
Command	Description	Example
cc-analyze	Molecular analysis	cc-analyze aspirin --properties
cc-synthesize	Synthesis planning	cc-synthesize target.smi --retro
cc-predict	Reaction prediction	cc-predict "A + B" --mechanism
cc-optimize	Property optimization	cc-optimize mol.sdf --target logP=2.5
cc-drug-design	Drug design workflow	cc-drug-design lead.smi --admet
cc-screen	Virtual screening	cc-screen library.sdf --target protein.pdb
Gemini CLI (chem: series)
Command	Description	Example
chem	Molecular analysis	gemini chem:analyze "CCO"
chem	Synthesis planning	gemini chem:synthesize aspirin
chem	Batch operations	gemini chem:batch mols.csv
chem	Structure visualization	gemini chem:visualize "c1ccccc1"
Chemical Sub-Agents
Sub-Agent	Expertise	Invocation
@organic-chemist	Organic synthesis & mechanism	@organic-chemist
@drug-designer	Drug design & ADMET	@drug-designer
@materials-scientist	Materials & polymers	@materials-scientist
@comp-chemist	Quantum chemistry	@comp-chemist
@analytical-chemist	Spectroscopy & analysis	@analytical-chemist
Comparison with Related Projects
Feature	ChemAgent	ChemCrow	SuperClaude
Positioning	AI assistant enhancement pack	Standalone agent	Claude framework
Platform	Claude Code + Gemini CLI	LangChain	Claude Code
Chemistry tools	✅ Full integration	✅ Full integration	❌ Generic tools
Installation	One-click script	pip install	One-click script
Command system	cc-/chem: series	Python API	sc- series
Sub-agents	5 chemistry experts	None	15 general roles
Roadmap
Support more AI platforms (GitHub Copilot, Codeium)
Integrate additional databases (ChemSpider, ZINC)
Add quantum chemistry and docking simulation APIs
Support MD simulations and molecular docking
Develop a VS Code extension
Create a web-based UI
Contributing



We welcome contributions! Especially in:

New chemistry commands
Platform adapters
Tool integrations
Documentation & examples
Acknowledgments
SuperClaude Framework – Architectural inspiration
ChemCrow – Chemistry tool references
RDKit, PubChem, and other open-source chemistry libraries
License

MIT License
