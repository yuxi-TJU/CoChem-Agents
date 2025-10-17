CoChem Agents

A collaborative framework for building chemistry & materials agents with Gemini CLI + MCP.

CoChem Agents turns the “one-agent-per-domain” pattern on its head. Instead of crafting bespoke chemistry or materials agents, we use Gemini CLI as the general-purpose agent shell and plug in domain tools via the Model Context Protocol (MCP). Anyone can contribute an MCP server—RDKit, Materials Project, your internal pipeline—and it becomes instantly usable by the same agent. This creates an open, extensible ecosystem rather than a zoo of siloed agents. 


Why chemistry & materials agents matter (now)

AI is rapidly accelerating discovery across chemistry and materials—from structure/property prediction to polymer and crystal modeling—pushing research beyond static prediction toward agentic workflows that plan, act, and iterate. Surveys and community reports document both the momentum and the need for robust tooling to make these systems practical in the lab and in silico. 

What hasn’t worked

Most prior efforts ship a standalone agent per subfield (drug design, catalysis, crystals…), each with custom glue code, brittle integrations, and duplicated effort. Evaluations often emphasize reasoning but struggle with reproducibility and tool-generalization, so systems don’t travel well between tasks or labs. Meanwhile, tool access (APIs, DBs, codes) is fragmented and hard to standardize across agents. 


The core technical obstacles

Heterogeneous tools & schemas: cheminformatics libs, materials databases, simulation engines—all different call patterns and data models. 

Agent–tool wiring & maintenance: each agent re-implements connectors and auth, leading to drift and duplication. 

Security & governance: opening tools to agents raises questions around auth, data access, and isolation. 

Evaluation & provenance: agent benchmarks underweight reproducibility and end-to-end paper-to-protocol faithfulness. 

Our approach (what’s different)

One agent framework to rule them all
Use Gemini CLI as the generic agent runtime (chat + tools + prompts). No more domain-specific shells. 

Tools as MCP servers
Expose chemistry/materials capabilities as MCP tools (standardized names, schemas, metadata). Any MCP-compatible client (like Gemini CLI) can discover and call them—zero bespoke glue in the agent. 


Open ecosystem, not one-off agents

Cheminformatics: community MCP servers for RDKit provide descriptor calc, substructure search, rendering, and more. Plug and use. 


Materials data: connect to Materials Project via its public API (or an MCP wrapper) for structures, formation energies, and band gaps. 


Custom science: FastMCP + Gemini CLI make it straightforward to publish your lab’s pipeline as a reusable tool, not a bespoke agent. 

Result: we standardize the interface (MCP) and reuse the agent shell (Gemini CLI). New capabilities arrive as new MCP servers, instantly composable with existing ones. This reduces integration tax, improves reproducibility, and channels community energy into shared tooling rather than parallel agent rewrites. 


TL;DR (project intent)

Mission: build an open, multi-tool ecosystem for chemistry & materials agents by unifying on Gemini CLI + MCP.

Why it matters: agentic science needs interoperable tools, not more siloed agents. 

How you can help: contribute or refine an MCP server (RDKit, Materials Project, simulations, ELN/SDMS, robo-lab APIs). The agent comes for free. 
Model Context Protocol

“Integrate once, compose everywhere.”
