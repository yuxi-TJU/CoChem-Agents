# mcp-pubchem

PubChem MCP Server - Model Context Protocol for PubChem database access

## Installation

```bash
pip install mcp-pubchem
```

## Usage

### As MCP Server (stdio mode)

```bash
python -m mcp_pubchem
```

### As TCP Server (for testing)

```bash
python -m mcp_pubchem --tcp 8766
```

### In Claude Desktop

Add to your Claude Desktop configuration:

```json
{
  "mcpServers": {
    "pubchem": {
      "command": "python",
      "args": ["-m", "mcp_pubchem"]
    }
  }
}
```

## Available Methods

### search_compound
Search for compounds in PubChem database

Parameters:
- `query` (string, required): Search query
- `search_type` (string): Type of search - "name", "smiles", "inchi", "formula"
- `max_results` (integer): Maximum number of results

### get_properties
Get properties for a compound by CID

Parameters:
- `cid` (integer, required): PubChem Compound ID
- `properties` (array): List of properties to retrieve

### get_synonyms
Get synonyms for a compound

Parameters:
- `cid` (integer, required): PubChem Compound ID

### find_similar
Find similar compounds

Parameters:
- `cid` (integer, required): PubChem Compound ID
- `threshold` (number): Similarity threshold (0-1)

## Example

```python
# Search for aspirin
{
  "method": "search_compound",
  "params": {
    "query": "aspirin",
    "search_type": "name",
    "max_results": 3
  }
}

# Get properties
{
  "method": "get_properties",
  "params": {
    "cid": 2244,
    "properties": ["molecular_weight", "canonical_smiles", "xlogp"]
  }
}
```

## License

MIT
