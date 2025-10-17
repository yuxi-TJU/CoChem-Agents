# MCP Integration Guide

Complete guide for integrating Chemistry MCP Servers with AI assistants and applications.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Claude Desktop Integration](#claude-desktop-integration)
3. [Programmatic Integration](#programmatic-integration)
4. [Best Practices](#best-practices)
5. [Troubleshooting](#troubleshooting)

## Quick Start

### 1. Install MCP Servers

```bash
# Clone repository
git clone https://github.com/dazhaolang/ai-chemkit.git
cd ai-chemkit

# Install all MCP servers
./install_all_mcp_servers.sh

# Or install individually
pip install mcp-rdkit
pip install mcp-pubchem
pip install mcp-chembl
pip install mcp-openbabel
```

### 2. Start Servers

```bash
# Start all servers
mcp start-all

# Check status
mcp status

# Start individual server
mcp start pubchem
```

### 3. Test Connection

```bash
# Run client examples
python mcp_servers/examples/mcp_client_example.py
```

## Claude Desktop Integration

### Configuration File Location

- **macOS**: `~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows**: `%APPDATA%\Claude\claude_desktop_config.json`
- **Linux**: `~/.config/Claude/claude_desktop_config.json`

### Complete Configuration

```json
{
  "mcpServers": {
    "rdkit": {
      "command": "python",
      "args": ["-m", "mcp_rdkit"],
      "env": {
        "PYTHONPATH": "/path/to/ai-chemkit"
      }
    },
    "pubchem": {
      "command": "python",
      "args": ["-m", "mcp_pubchem_server"],
      "env": {
        "PYTHONPATH": "/path/to/ai-chemkit/mcp_servers"
      }
    },
    "chembl": {
      "command": "python",
      "args": ["-m", "mcp_chembl_server"],
      "env": {
        "PYTHONPATH": "/path/to/ai-chemkit/mcp_servers"
      }
    },
    "openbabel": {
      "command": "python",
      "args": ["-m", "mcp_openbabel_server"],
      "env": {
        "PYTHONPATH": "/path/to/ai-chemkit/mcp_servers"
      }
    }
  }
}
```

### Using in Claude

Once configured, you can use natural language:

```
"Search PubChem for aspirin and show me its properties"
"Find EGFR inhibitors in ChEMBL with pIC50 > 7"
"Convert this SMILES to 3D structure: CC(=O)OC1=CC=CC=C1C(=O)O"
"Check if this molecule is drug-like: CCN(CC)C(=O)C1=CC=CC=C1"
```

## Programmatic Integration

### Python Integration

```python
import asyncio
import json
from pathlib import Path
import sys

# Add MCP servers to path
sys.path.append(str(Path(__file__).parent / "mcp_servers"))

from mcp_pubchem.mcp_pubchem_server import PubChemMCPServer
from mcp_chembl.mcp_chembl_server import ChEMBLMCPServer
from mcp_openbabel.mcp_openbabel_server import OpenBabelMCPServer

async def use_mcp_servers():
    # Initialize servers
    pubchem = PubChemMCPServer()
    chembl = ChEMBLMCPServer()
    openbabel = OpenBabelMCPServer()
    
    # Search PubChem
    request = json.dumps({
        "jsonrpc": "2.0",
        "method": "tools/call",
        "params": {
            "name": "search_compound",
            "arguments": {
                "query": "caffeine",
                "search_type": "name"
            }
        },
        "id": 1
    })
    
    response = await pubchem.process_message(request)
    result = json.loads(response)
    print(f"PubChem result: {result}")
    
    # Use other servers similarly...

asyncio.run(use_mcp_servers())
```

### Node.js Integration

```javascript
const { spawn } = require('child_process');
const readline = require('readline');

class MCPClient {
  constructor(command, args) {
    this.process = spawn(command, args);
    this.rl = readline.createInterface({
      input: this.process.stdout,
      output: this.process.stdin
    });
    this.requestId = 0;
  }
  
  async sendRequest(method, params) {
    const request = {
      jsonrpc: "2.0",
      method: method,
      params: params,
      id: ++this.requestId
    };
    
    return new Promise((resolve, reject) => {
      this.rl.once('line', (line) => {
        try {
          const response = JSON.parse(line);
          if (response.error) {
            reject(new Error(response.error.message));
          } else {
            resolve(response.result);
          }
        } catch (e) {
          reject(e);
        }
      });
      
      this.process.stdin.write(JSON.stringify(request) + '\n');
    });
  }
  
  async initialize() {
    return await this.sendRequest('initialize', {
      clientInfo: {
        name: 'Node.js Client',
        version: '1.0.0'
      }
    });
  }
}

// Usage
async function main() {
  const client = new MCPClient('python', ['-m', 'mcp_pubchem_server']);
  
  await client.initialize();
  
  const result = await client.sendRequest('tools/call', {
    name: 'search_compound',
    arguments: {
      query: 'aspirin',
      search_type: 'name'
    }
  });
  
  console.log('Result:', result);
}

main().catch(console.error);
```

### REST API Wrapper

Create a REST API wrapper for HTTP clients:

```python
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import asyncio
import json

app = FastAPI()

# Initialize MCP servers
servers = {
    "pubchem": PubChemMCPServer(),
    "chembl": ChEMBLMCPServer(),
    "openbabel": OpenBabelMCPServer()
}

class MCPRequest(BaseModel):
    server: str
    method: str
    params: dict

@app.post("/mcp/call")
async def call_mcp(request: MCPRequest):
    if request.server not in servers:
        raise HTTPException(status_code=404, detail="Server not found")
    
    server = servers[request.server]
    
    mcp_request = json.dumps({
        "jsonrpc": "2.0",
        "method": request.method,
        "params": request.params,
        "id": 1
    })
    
    response = await server.process_message(mcp_request)
    result = json.loads(response)
    
    if "error" in result:
        raise HTTPException(status_code=400, detail=result["error"])
    
    return result["result"]

# Run with: uvicorn api_wrapper:app --reload
```

## Best Practices

### 1. Error Handling

Always handle MCP errors gracefully:

```python
async def safe_mcp_call(server, method, params):
    try:
        request = json.dumps({
            "jsonrpc": "2.0",
            "method": method,
            "params": params,
            "id": 1
        })
        
        response = await server.process_message(request)
        result = json.loads(response)
        
        if "error" in result:
            logger.error(f"MCP error: {result['error']}")
            return None
        
        return result.get("result")
        
    except Exception as e:
        logger.error(f"Failed to call MCP: {e}")
        return None
```

### 2. Connection Pooling

For high-throughput applications:

```python
import asyncio
from concurrent.futures import ThreadPoolExecutor

class MCPPool:
    def __init__(self, server_class, pool_size=5):
        self.servers = [server_class() for _ in range(pool_size)]
        self.semaphore = asyncio.Semaphore(pool_size)
    
    async def call(self, method, params):
        async with self.semaphore:
            server = self.servers.pop(0)
            try:
                result = await server.call_tool(method, params)
                return result
            finally:
                self.servers.append(server)

# Usage
pool = MCPPool(PubChemMCPServer, pool_size=10)
results = await asyncio.gather(*[
    pool.call("search_compound", {"query": name})
    for name in compound_names
])
```

### 3. Caching Results

Implement caching for frequently accessed data:

```python
from functools import lru_cache
import hashlib

class CachedMCPClient:
    def __init__(self, server):
        self.server = server
        self.cache = {}
    
    def _cache_key(self, method, params):
        """Generate cache key from method and params"""
        key_str = f"{method}:{json.dumps(params, sort_keys=True)}"
        return hashlib.md5(key_str.encode()).hexdigest()
    
    async def call(self, method, params, cache_ttl=3600):
        cache_key = self._cache_key(method, params)
        
        # Check cache
        if cache_key in self.cache:
            cached_time, cached_result = self.cache[cache_key]
            if time.time() - cached_time < cache_ttl:
                return cached_result
        
        # Call server
        result = await self.server.call_tool(method, params)
        
        # Update cache
        self.cache[cache_key] = (time.time(), result)
        
        return result
```

### 4. Batch Processing

Process multiple requests efficiently:

```python
async def batch_process(server, requests):
    """Process multiple requests as a batch"""
    
    # Create batch request
    batch = [
        {
            "jsonrpc": "2.0",
            "method": req["method"],
            "params": req["params"],
            "id": i
        }
        for i, req in enumerate(requests)
    ]
    
    # Send as batch
    response = await server.process_message(json.dumps(batch))
    results = json.loads(response)
    
    # Map results back
    return {r["id"]: r.get("result") for r in results}
```

## Troubleshooting

### Common Issues

#### 1. Server Not Starting

```bash
# Check if port is in use
lsof -i :8767

# Kill existing process
kill -9 <PID>

# Restart server
mcp restart pubchem
```

#### 2. Import Errors

```bash
# Check Python path
python -c "import sys; print(sys.path)"

# Add to PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:/path/to/ai-chemkit/mcp_servers"
```

#### 3. Connection Timeouts

```python
# Increase timeout in client
async with aiohttp.ClientSession(
    timeout=aiohttp.ClientTimeout(total=60)
) as session:
    # ... make requests
```

#### 4. Memory Issues

```python
# Limit concurrent requests
semaphore = asyncio.Semaphore(10)

async def limited_call(server, method, params):
    async with semaphore:
        return await server.call_tool(method, params)
```

### Debug Mode

Enable debug logging:

```python
import logging

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

# Or set environment variable
os.environ['LOGLEVEL'] = 'DEBUG'
```

### Performance Monitoring

```python
import time
from contextlib import asynccontextmanager

@asynccontextmanager
async def timed_operation(name):
    start = time.time()
    try:
        yield
    finally:
        duration = time.time() - start
        print(f"{name} took {duration:.2f} seconds")

# Usage
async with timed_operation("PubChem search"):
    result = await client.call_tool("search_compound", {...})
```

## Advanced Topics

### Custom MCP Server

Create your own MCP server:

```python
from base_mcp_server_v2 import BaseMCPServerV2, mcp_tool, ToolAnnotation, RiskLevel

class CustomChemServer(BaseMCPServerV2):
    def __init__(self):
        super().__init__("custom-chem", "1.0.0")
    
    def _setup_capabilities(self):
        self.capabilities = {
            "custom_analysis": True
        }
    
    def _register_tools(self):
        self.register_tool(
            "analyze_custom",
            self.analyze_custom,
            ToolAnnotation(
                risk_level=RiskLevel.LOW,
                description="Custom analysis"
            )
        )
    
    async def analyze_custom(self, data: str) -> dict:
        """Perform custom analysis"""
        # Your implementation
        return {"result": "analyzed"}
```

### Security Considerations

1. **Input Validation**: Always validate SMILES/InChI inputs
2. **Rate Limiting**: Implement rate limiting for public APIs
3. **Authentication**: Add API keys for production use
4. **Sandboxing**: Run MCP servers in containers for isolation

## Resources

- [MCP Specification](https://modelcontextprotocol.io)
- [ChemAgent GitHub](https://github.com/dazhaolang/ai-chemkit)
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [PubChem API](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)
- [ChEMBL Web Services](https://chembl.gitbook.io/chembl-interface-documentation/)

## Support

- GitHub Issues: https://github.com/dazhaolang/ai-chemkit/issues
- Discord: [Join our community](https://discord.gg/chemagent)
- Email: support@chemagent.org
