#!/usr/bin/env python3
"""
Test script for Chemistry MCP Servers
测试化学MCP服务器
"""

import asyncio
import json
import sys
from pathlib import Path

# Add mcp_servers to path
sys.path.insert(0, str(Path(__file__).parent / "mcp_servers"))

async def test_mcp_server(server_class, test_requests):
    """Test a MCP server with sample requests"""
    server = server_class()
    
    print(f"\n{'='*50}")
    print(f"Testing {server.name} v{server.version}")
    print('='*50)
    
    # Test initialize
    init_request = json.dumps({
        "method": "initialize",
        "params": {},
        "id": "init-1"
    })
    
    response = await server.process_message(init_request)
    result = json.loads(response)
    
    print("\n1. Server Info:")
    if "result" in result:
        print(f"   Name: {result['result']['name']}")
        print(f"   Version: {result['result']['version']}")
        print(f"   Protocol: {result['result']['protocol']}")
    
    # Test list_tools
    list_request = json.dumps({
        "method": "list_tools",
        "params": {},
        "id": "list-1"
    })
    
    response = await server.process_message(list_request)
    result = json.loads(response)
    
    print("\n2. Available Tools:")
    if "result" in result:
        for tool in result['result'][:5]:  # Show first 5 tools
            print(f"   • {tool['name']}: {tool['description'][:50]}...")
    
    # Test specific methods
    print("\n3. Testing Methods:")
    for i, request_data in enumerate(test_requests, 1):
        request = json.dumps({
            "method": request_data["method"],
            "params": request_data["params"],
            "id": f"test-{i}"
        })
        
        print(f"\n   Test {i}: {request_data['method']}")
        print(f"   Params: {request_data['params']}")
        
        response = await server.process_message(request)
        result = json.loads(response)
        
        if "result" in result:
            print(f"   ✅ Success:")
            # Show partial result
            result_str = json.dumps(result['result'], indent=2)
            lines = result_str.split('\n')
            for line in lines[:10]:  # Show first 10 lines
                print(f"      {line}")
            if len(lines) > 10:
                print("      ...")
        elif "error" in result:
            print(f"   ❌ Error: {result['error']['message']}")


async def main():
    print("="*60)
    print("   Chemistry MCP Servers Test Suite")
    print("="*60)
    
    # Test PubChem MCP
    try:
        from mcp_pubchem.mcp_pubchem_server import PubChemMCPServer
        
        pubchem_tests = [
            {
                "method": "search_compound",
                "params": {
                    "query": "aspirin",
                    "search_type": "name",
                    "max_results": 2
                }
            },
            {
                "method": "get_properties",
                "params": {
                    "cid": 2244,
                    "properties": ["molecular_weight", "canonical_smiles", "xlogp"]
                }
            }
        ]
        
        await test_mcp_server(PubChemMCPServer, pubchem_tests)
    except ImportError as e:
        print(f"\n⚠️  PubChem MCP not available: {e}")
    except Exception as e:
        print(f"\n❌ PubChem MCP test failed: {e}")
    
    # Test ChEMBL MCP
    try:
        from mcp_chembl.mcp_chembl_server import ChEMBLMCPServer
        
        chembl_tests = [
            {
                "method": "search_molecule",
                "params": {
                    "query": "imatinib",
                    "search_type": "name",
                    "max_results": 2
                }
            },
            {
                "method": "search_target",
                "params": {
                    "query": "kinase",
                    "organism": "Homo sapiens",
                    "max_results": 2
                }
            }
        ]
        
        await test_mcp_server(ChEMBLMCPServer, chembl_tests)
    except ImportError as e:
        print(f"\n⚠️  ChEMBL MCP not available: {e}")
    except Exception as e:
        print(f"\n❌ ChEMBL MCP test failed: {e}")
    
    # Test OpenBabel MCP
    try:
        from mcp_openbabel.mcp_openbabel_server import OpenBabelMCPServer
        
        openbabel_tests = [
            {
                "method": "convert_molecule",
                "params": {
                    "molecule": "CCO",
                    "input_format": "smiles",
                    "output_format": "inchi"
                }
            },
            {
                "method": "canonicalize_smiles",
                "params": {
                    "smiles": "C(C)O"
                }
            },
            {
                "method": "calculate_properties",
                "params": {
                    "molecule": "CC(=O)OC1=CC=CC=C1C(=O)O",
                    "format": "smiles"
                }
            }
        ]
        
        await test_mcp_server(OpenBabelMCPServer, openbabel_tests)
    except ImportError as e:
        print(f"\n⚠️  OpenBabel MCP not available: {e}")
    except Exception as e:
        print(f"\n❌ OpenBabel MCP test failed: {e}")
    
    print("\n" + "="*60)
    print("   Test Suite Complete")
    print("="*60)


if __name__ == "__main__":
    # Check dependencies
    print("\nChecking dependencies...")
    
    dependencies = {
        "pubchempy": "PubChem MCP",
        "chembl_webresource_client": "ChEMBL MCP",
        "openbabel": "OpenBabel MCP (Python)",
    }
    
    missing = []
    for module, name in dependencies.items():
        try:
            __import__(module.replace("-", "_"))
            print(f"✅ {name}: {module} installed")
        except ImportError:
            print(f"⚠️  {name}: {module} not installed")
            missing.append(module)
    
    if missing:
        print(f"\nTo install missing dependencies:")
        print(f"pip install {' '.join(missing)}")
        print("\nContinuing with available servers...\n")
    
    # Run tests
    asyncio.run(main())
