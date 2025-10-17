# mcp-openbabel

OpenBabel MCP Server - Model Context Protocol for molecular format conversion and 3D generation

## Overview

This MCP server provides access to OpenBabel functionality, enabling:
- Molecular format conversion (SMILES, MOL, SDF, PDB, XYZ, etc.)
- 3D coordinate generation and optimization
- Property calculations
- Drug-likeness assessment
- Canonical SMILES generation

## Installation

### Prerequisites

OpenBabel can be installed via:

**Option 1: Python bindings (recommended)**
```bash
pip install openbabel-wheel
```

**Option 2: System package**
```bash
# Ubuntu/Debian
sudo apt-get install openbabel python3-openbabel

# macOS
brew install open-babel

# Windows
# Download from http://openbabel.org/wiki/Category:Installation
```

### MCP Server Installation

**From PyPI (when available):**
```bash
pip install mcp-openbabel
```

**From source:**
```bash
git clone https://github.com/dazhaolang/ai-chemkit.git
cd ai-chemkit/mcp_servers/mcp_openbabel
pip install -e .
```

## Configuration

### For Claude Desktop

Add to your configuration file:

**macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
**Windows**: `%APPDATA%\Claude\claude_desktop_config.json`

```json
{
  "mcpServers": {
    "openbabel": {
      "command": "python",
      "args": ["-m", "mcp_openbabel_server"]
    }
  }
}
```

## Available Tools

### 1. convert_molecule
Convert molecules between different formats

**Parameters:**
- `molecule` (string, required): Input molecule data
- `input_format` (string, required): Input format (smiles, mol, sdf, pdb, etc.)
- `output_format` (string, required): Output format
- `options` (object): Additional options
  - `add_hydrogens` (boolean): Add explicit hydrogens
  - `remove_hydrogens` (boolean): Remove hydrogens
  - `gen3d` (boolean): Generate 3D coordinates

**Example:**
```json
{
  "method": "tools/call",
  "params": {
    "name": "convert_molecule",
    "arguments": {
      "molecule": "CCO",
      "input_format": "smiles",
      "output_format": "mol",
      "options": {
        "gen3d": true,
        "add_hydrogens": true
      }
    }
  }
}
```

### 2. generate_3d
Generate 3D coordinates from SMILES

**Parameters:**
- `smiles` (string, required): Input SMILES string
- `force_field` (string): Force field for optimization - "mmff94", "uff", "ghemical" (default: "mmff94")
- `steps` (integer): Optimization steps (default: 500)

**Example:**
```json
{
  "method": "tools/call",
  "params": {
    "name": "generate_3d",
    "arguments": {
      "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
      "force_field": "mmff94",
      "steps": 1000
    }
  }
}
```

**Response:**
```json
{
  "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
  "mol_block": "...[3D MOL data]...",
  "xyz": "...[XYZ coordinates]...",
  "energy": -45.678
}
```

### 3. optimize_geometry
Optimize molecular geometry

**Parameters:**
- `molecule` (string, required): Input molecule
- `format` (string): Input format (default: "mol")
- `force_field` (string): Force field (default: "mmff94")
- `steps` (integer): Optimization steps (default: 2500)

**Example:**
```json
{
  "method": "tools/call",
  "params": {
    "name": "optimize_geometry",
    "arguments": {
      "molecule": "...[MOL data]...",
      "format": "mol",
      "force_field": "uff",
      "steps": 5000
    }
  }
}
```

### 4. canonicalize_smiles
Generate canonical SMILES

**Parameters:**
- `smiles` (string, required): Input SMILES

**Example:**
```json
{
  "method": "tools/call",
  "params": {
    "name": "canonicalize_smiles",
    "arguments": {
      "smiles": "C(O)C"
    }
  }
}
```

**Response:**
```json
{
  "input_smiles": "C(O)C",
  "canonical_smiles": "CCO",
  "is_valid": true
}
```

### 5. calculate_properties
Calculate molecular properties

**Parameters:**
- `molecule` (string, required): Input molecule
- `format` (string): Input format (default: "smiles")

**Example:**
```json
{
  "method": "tools/call",
  "params": {
    "name": "calculate_properties",
    "arguments": {
      "molecule": "CC(=O)OC1=CC=CC=C1C(=O)O",
      "format": "smiles"
    }
  }
}
```

**Response:**
```json
{
  "molecular_weight": 180.16,
  "formula": "C9H8O4",
  "num_atoms": 21,
  "num_bonds": 21,
  "num_rotatable_bonds": 3,
  "charge": 0,
  "inchi": "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
  "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
}
```

### 6. check_drug_like
Check if molecule passes Lipinski's Rule of Five

**Parameters:**
- `smiles` (string, required): Input SMILES

**Example:**
```json
{
  "method": "tools/call",
  "params": {
    "name": "check_drug_like",
    "arguments": {
      "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
    }
  }
}
```

**Response:**
```json
{
  "molecular_weight": 180.16,
  "logp": 1.31,
  "h_bond_donors": 1,
  "h_bond_acceptors": 4,
  "lipinski_violations": 0,
  "is_drug_like": true
}
```

## Usage Examples

### Python Client

```python
import json
import asyncio
from mcp_openbabel_server import OpenBabelMCPServer

async def convert_smiles_to_3d():
    server = OpenBabelMCPServer()
    
    # Convert SMILES to 3D MOL
    request = json.dumps({
        "jsonrpc": "2.0",
        "method": "tools/call",
        "params": {
            "name": "generate_3d",
            "arguments": {
                "smiles": "CCO"
            }
        },
        "id": 1
    })
    
    response = await server.process_message(request)
    result = json.loads(response)
    print(f"3D coordinates generated: {result['result']['xyz']}")

asyncio.run(convert_smiles_to_3d())
```

### Command Line

```bash
# Start TCP server
python -m mcp_openbabel_server --tcp 8769

# Test with curl
curl -X POST http://localhost:8769 \
  -H "Content-Type: application/json" \
  -d '{"jsonrpc":"2.0","method":"tools/call","params":{"name":"canonicalize_smiles","arguments":{"smiles":"C(C)O"}},"id":1}'
```

### Batch Conversion Script

```python
import asyncio
import json
from pathlib import Path

async def batch_convert(input_file, from_format, to_format):
    """Convert all molecules in a file"""
    server = OpenBabelMCPServer()
    
    with open(input_file) as f:
        molecules = f.readlines()
    
    results = []
    for i, mol in enumerate(molecules):
        request = json.dumps({
            "jsonrpc": "2.0",
            "method": "tools/call",
            "params": {
                "name": "convert_molecule",
                "arguments": {
                    "molecule": mol.strip(),
                    "input_format": from_format,
                    "output_format": to_format
                }
            },
            "id": i
        })
        
        response = await server.process_message(request)
        results.append(json.loads(response))
    
    return results
```

## Supported Formats

### Input/Output Formats
- **SMILES**: Simplified molecular input line entry system
- **MOL/SDF**: MDL molecule/structure data file
- **PDB**: Protein Data Bank format
- **XYZ**: Cartesian coordinates
- **InChI**: International Chemical Identifier
- **CML**: Chemical Markup Language
- **MOL2**: Tripos molecule format
- **PDBQT**: AutoDock format

### Force Fields for Optimization
- **MMFF94**: Merck Molecular Force Field
- **UFF**: Universal Force Field
- **Ghemical**: Ghemical force field

## Performance Tips

1. **Batch Operations**: Process multiple molecules in one session
2. **Format Selection**: Use binary formats (MOL, SDF) for 3D data
3. **Optimization Steps**: Balance accuracy vs speed (500-2500 steps)
4. **Caching**: Cache canonical SMILES for repeated conversions

## Troubleshooting

### OpenBabel Not Found

If you get "OpenBabel not available" errors:

1. Check Python bindings:
```python
import openbabel
print(openbabel.__version__)
```

2. Check CLI:
```bash
obabel -V
```

3. Reinstall:
```bash
pip uninstall openbabel-wheel
pip install --no-cache-dir openbabel-wheel
```

### Slow 3D Generation

For faster 3D generation:
- Reduce optimization steps
- Use UFF instead of MMFF94 for initial coordinates
- Pre-filter invalid SMILES

## License

MIT License

## Support

- GitHub: https://github.com/dazhaolang/ai-chemkit
- OpenBabel Docs: http://openbabel.org/docs/current/

## Citation

```bibtex
@software{mcp_openbabel,
  title = {OpenBabel MCP Server},
  author = {ChemAgent Team},
  year = {2024},
  url = {https://github.com/dazhaolang/ai-chemkit}
}
```
