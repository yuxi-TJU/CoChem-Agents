# CoChem Agents

_A collaborative framework for building chemistry & materials agents with **Gemini CLI + MCP**._

> **CoChem Agents** turns the “one-agent-per-domain” pattern on its head. Instead of crafting bespoke chemistry or materials agents, we use **Gemini CLI** as the general-purpose agent shell and plug in domain tools via the **Model Context Protocol (MCP)**. Anyone can contribute an MCP server—RDKit, Materials Project, your internal pipeline—and it becomes instantly usable by the same agent. This creates an open, extensible ecosystem rather than a zoo of siloed agents.

---

## 🌍 Why chemistry & materials agents matter

AI is rapidly accelerating discovery across chemistry and materials—from structure/property prediction to polymer and crystal modeling—pushing research beyond static prediction toward **agentic** workflows that plan, act, and iterate. Surveys and community reports document both the momentum and the need for robust tooling to make these systems practical in the lab and in silico.

---

## 🧱 What hasn’t worked

Most prior efforts ship a **standalone agent per subfield** (drug design, catalysis, crystals…), each with custom glue code, brittle integrations, and duplicated effort.

- Agents are hard to reuse across domains.
- Tool interfaces are fragmented and not standardized.
- Reproducibility suffers due to ad-hoc wiring.

These systems often emphasize reasoning but struggle with **tool generalization, maintainability, and integration overhead**.

---

## 🧩 The core technical obstacles

- **Heterogeneous tools & schemas**  
  Cheminformatics libraries, materials databases, simulation engines—all use different APIs and data models.

- **Agent–tool wiring & maintenance**  
  Each agent re-implements connectors, leading to duplicated work and brittle interfaces.

- **Security & governance**  
  Agents often require open access to tools or keys, raising issues of control and reproducibility.

- **Evaluation & provenance**  
  Benchmarks often underweight reproducibility and end-to-end experimental fidelity.

---

## ✅ The CoChem Agents

### 1. One agent framework to rule them all  
Use **Gemini CLI** as the **unified agent runtime**: chat, ReAct loop, tools, streaming — no need to build your own shell.

### 2. Tools as MCP servers  
Expose chemistry or materials tools (e.g., RDKit, Materials Project) via the **Model Context Protocol (MCP)**.  
Gemini CLI can auto-discover and invoke them, with **no extra glue code.**

### 3. Open ecosystem, not one-off agents  
- **Cheminformatics:** use community MCP servers for **RDKit** — compute descriptors, search substructures, render molecules, etc.  
- **Materials data:** plug in **Materials Project** for band gaps, formation energies, structures.  
- **Custom tools:** package your lab's own script into a **FastMCP** server — instantly usable by the agent.

---

## 🌱 Result

We standardize the _interface_ (MCP) and reuse the _agent shell_ (Gemini CLI). New capabilities arrive as new MCP servers, instantly composable with existing ones.

This:
- reduces integration burden  
- improves reproducibility  
- encourages community tool sharing instead of siloed agent development

> _“Integrate once, compose everywhere.”_

---

## 🎯 TL;DR (project intent)

- **Mission:** build an open, multi-tool ecosystem for chemistry & materials agents using **Gemini CLI + MCP**
- **Why it matters:** agentic science needs interoperable tools, not a zoo of siloed agents
- **How you can help:** contribute or refine an **MCP server** — for RDKit, Materials Project, simulations, ELNs, robo-lab APIs, etc.  
  The agent comes for free.
