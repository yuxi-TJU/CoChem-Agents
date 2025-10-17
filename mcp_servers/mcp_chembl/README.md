# mcp-chembl

ChEMBL MCP Server - Model Context Protocol for ChEMBL bioactivity database

## Overview

This MCP server provides access to the ChEMBL database, a manually curated database of bioactive molecules with drug-like properties. It contains 2D structures, calculated properties and abstracted bioactivities for millions of compounds.

## Installation

### From PyPI (when available)
```bash
pip install mcp-chembl
```

### From source
```bash
git clone https://github.com/dazhaolang/ai-chemkit.git
cd ai-chemkit/mcp_servers/mcp_chembl
pip install -e .
```

## Configuration

### For Claude Desktop

Add to your Claude Desktop configuration file:

**macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
**Windows**: `%APPDATA%\Claude\claude_desktop_config.json`

```json
{
  "mcpServers": {
    "chembl": {
      "command": "python",
      "args": ["-m", "mcp_chembl_server"]
    }
  }
}
```

### For other MCP clients

```json
{
  "servers": {
    "chembl": {
      "type": "stdio",
      "command": "python -m mcp_chembl_server"
    }
  }
}
```

## Available Tools

### 1. search_molecule
Search for molecules in ChEMBL database

**Parameters:**
- `query` (string, required): Search query
- `search_type` (string): Type of search - "name", "smiles", "chembl_id" (default: "name")
- `max_results` (integer): Maximum number of results (default: 10)

**Example:**
```json
{
  "method": "tools/call",
  "params": {
    "name": "search_molecule",
    "arguments": {
      "query": "imatinib",
      "search_type": "name",
      "max_results": 5
    }
  }
}
```

**Response:**
```json
{
  "molecules": [
    {
      "chembl_id": "CHEMBL941",
      "name": "IMATINIB",
      "smiles": "CN1CCN(Cc2ccc(cc2)C(=O)Nc2ccc(C)c(Nc3nccc(n3)c3cccnc3)c2)CC1",
      "molecular_weight": 493.603,
      "max_phase": 4
    }
  ]
}
```

### 2. search_target
Search for protein targets

**Parameters:**
- `query` (string, required): Search query
- `organism` (string): Organism filter (default: "Homo sapiens")
- `max_results` (integer): Maximum results (default: 10)

**Example:**
```json
{
  "method": "tools/call",
  "params": {
    "name": "search_target",
    "arguments": {
      "query": "EGFR",
      "organism": "Homo sapiens"
    }
  }
}
```

### 3. get_bioactivity
Get bioactivity data for molecules or targets

**Parameters:**
- `molecule_chembl_id` (string, optional): Molecule ChEMBL ID
- `target_chembl_id` (string, optional): Target ChEMBL ID
- `min_pchembl` (number): Minimum pChEMBL value (default: 5)
- `max_results` (integer): Maximum results (default: 50)

**Example:**
```json
{
  "method": "tools/call",
  "params": {
    "name": "get_bioactivity",
    "arguments": {
      "molecule_chembl_id": "CHEMBL941",
      "min_pchembl": 6
    }
  }
}
```

### 4. get_drug_info
Get detailed drug information

**Parameters:**
- `chembl_id` (string, required): ChEMBL ID of the drug

**Example:**
```json
{
  "method": "tools/call",
  "params": {
    "name": "get_drug_info",
    "arguments": {
      "chembl_id": "CHEMBL941"
    }
  }
}
```

### 5. find_drugs_for_target
Find approved drugs for a specific target

**Parameters:**
- `target_chembl_id` (string, required): Target ChEMBL ID
- `min_phase` (integer): Minimum clinical phase (default: 3)

**Example:**
```json
{
  "method": "tools/call",
  "params": {
    "name": "find_drugs_for_target",
    "arguments": {
      "target_chembl_id": "CHEMBL203",
      "min_phase": 4
    }
  }
}
```

## Usage Examples

### Python Client Example

```python
import json
import asyncio
from mcp_chembl_server import ChEMBLMCPServer

async def example():
    server = ChEMBLMCPServer()
    
    # Search for a drug
    request = json.dumps({
        "jsonrpc": "2.0",
        "method": "tools/call",
        "params": {
            "name": "search_molecule",
            "arguments": {
                "query": "aspirin",
                "search_type": "name"
            }
        },
        "id": 1
    })
    
    response = await server.process_message(request)
    print(json.loads(response))

asyncio.run(example())
```

### Command Line Usage

```bash
# Start as TCP server for testing
python -m mcp_chembl_server --tcp 8768

# In another terminal, test with netcat
echo '{"jsonrpc":"2.0","method":"initialize","params":{},"id":1}' | nc localhost 8768
```

### Integration with LLMs

When integrated with Claude or other LLMs, you can use natural language:

```
"Search ChEMBL for EGFR inhibitors"
"What drugs target CHEMBL203?"
"Find bioactivity data for imatinib"
"Show me approved drugs for treating leukemia targets"
```

## Data Sources

- **ChEMBL Database**: Version 33 (or latest)
- **API**: ChEMBL Web Services
- **Update Frequency**: ChEMBL releases new versions approximately every 3-4 months

## Limitations

1. **API Rate Limits**: ChEMBL web services may have rate limits
2. **Data Size**: Large queries may be slow or timeout
3. **Network Required**: Requires internet connection to access ChEMBL
4. **Data Currency**: Data is as current as the latest ChEMBL release

## Error Handling

The server returns standard JSON-RPC 2.0 errors:

```json
{
  "jsonrpc": "2.0",
  "error": {
    "code": -32603,
    "message": "No target data found for CHEMBL999999"
  },
  "id": 1
}
```

## Development

### Running Tests

```bash
cd mcp_servers/mcp_chembl
python -m pytest tests/
```

### Debug Mode

```bash
# Run with debug logging
LOGLEVEL=DEBUG python -m mcp_chembl_server
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## License

MIT License - see LICENSE file

## Support

- GitHub Issues: https://github.com/dazhaolang/ai-chemkit/issues
- Documentation: https://github.com/dazhaolang/ai-chemkit/wiki
- ChEMBL Documentation: https://chembl.gitbook.io/chembl-interface-documentation/

## Citation

If you use this MCP server in your research, please cite:

```bibtex
@software{mcp_chembl,
  title = {ChEMBL MCP Server},
  author = {ChemAgent Team},
  year = {2024},
  url = {https://github.com/dazhaolang/ai-chemkit}
}
```

And the ChEMBL database:

```bibtex
@article{chembl_2023,
  title = {ChEMBL Database},
  author = {EMBL-EBI},
  journal = {Nucleic Acids Research},
  year = {2023}
}
```
