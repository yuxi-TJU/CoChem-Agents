#!/usr/bin/env python3
"""
MCP Client Example
Demonstrates how to interact with Chemistry MCP Servers
"""

import asyncio
import json
import sys
from typing import Any, Dict, List, Optional
from dataclasses import dataclass
import aiohttp


@dataclass
class MCPClient:
    """Simple MCP client for testing and demonstration"""
    
    def __init__(self, server_url: str = "http://localhost:8766"):
        self.server_url = server_url
        self.session = None
        self.request_id = 0
    
    async def __aenter__(self):
        self.session = aiohttp.ClientSession()
        await self.initialize()
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        if self.session:
            await self.session.close()
    
    def _next_id(self) -> int:
        """Generate next request ID"""
        self.request_id += 1
        return self.request_id
    
    async def send_request(self, method: str, params: Dict[str, Any] = None) -> Dict[str, Any]:
        """Send a JSON-RPC request to the server"""
        request = {
            "jsonrpc": "2.0",
            "method": method,
            "params": params or {},
            "id": self._next_id()
        }
        
        async with self.session.post(self.server_url, json=request) as response:
            result = await response.json()
            
            if "error" in result:
                raise Exception(f"Server error: {result['error']}")
            
            return result.get("result")
    
    async def initialize(self) -> Dict[str, Any]:
        """Initialize connection with the server"""
        return await self.send_request("initialize", {
            "clientInfo": {
                "name": "MCP Python Client",
                "version": "1.0.0"
            },
            "protocolVersion": "1.0.0"
        })
    
    async def list_tools(self) -> List[Dict[str, Any]]:
        """Get list of available tools"""
        return await self.send_request("tools/list")
    
    async def call_tool(self, tool_name: str, arguments: Dict[str, Any]) -> Any:
        """Call a specific tool"""
        return await self.send_request("tools/call", {
            "name": tool_name,
            "arguments": arguments
        })


# Example workflows for each MCP server

async def pubchem_example():
    """Example using PubChem MCP Server"""
    print("\n" + "="*50)
    print("PubChem MCP Server Example")
    print("="*50)
    
    async with MCPClient("http://localhost:8767") as client:
        # Search for aspirin
        print("\n1. Searching for aspirin...")
        result = await client.call_tool("search_compound", {
            "query": "aspirin",
            "search_type": "name",
            "max_results": 3
        })
        
        print(f"Found {result['count']} compounds")
        for compound in result['compounds']:
            print(f"  - CID: {compound['cid']}")
            print(f"    Formula: {compound['formula']}")
            print(f"    Weight: {compound['weight']}")
            print(f"    SMILES: {compound['smiles']}")
        
        # Get properties for aspirin (CID: 2244)
        print("\n2. Getting properties for aspirin (CID: 2244)...")
        props = await client.call_tool("get_properties", {
            "cid": 2244,
            "properties": ["molecular_weight", "xlogp", "tpsa", "h_bond_donor_count"]
        })
        
        print("Properties:")
        for key, value in props.items():
            if key != "cid":
                print(f"  - {key}: {value}")


async def chembl_example():
    """Example using ChEMBL MCP Server"""
    print("\n" + "="*50)
    print("ChEMBL MCP Server Example")
    print("="*50)
    
    async with MCPClient("http://localhost:8768") as client:
        # Search for EGFR inhibitors
        print("\n1. Searching for EGFR target...")
        result = await client.call_tool("search_target", {
            "query": "EGFR",
            "organism": "Homo sapiens",
            "max_results": 2
        })
        
        print(f"Found {result['count']} targets")
        for target in result['targets']:
            print(f"  - {target['chembl_id']}: {target['name']}")
            print(f"    Type: {target['type']}")
            print(f"    Organism: {target['organism']}")
        
        # Find drugs for EGFR
        if result['targets']:
            target_id = result['targets'][0]['chembl_id']
            print(f"\n2. Finding drugs for {target_id}...")
            
            drugs = await client.call_tool("find_drugs_for_target", {
                "target_chembl_id": target_id,
                "min_phase": 3
            })
            
            print(f"Found {drugs['count']} drugs")
            for drug in drugs['drugs'][:5]:
                print(f"  - {drug['chembl_id']}: {drug.get('name', 'N/A')}")
                print(f"    Max Phase: {drug['max_phase']}")


async def openbabel_example():
    """Example using OpenBabel MCP Server"""
    print("\n" + "="*50)
    print("OpenBabel MCP Server Example")
    print("="*50)
    
    async with MCPClient("http://localhost:8769") as client:
        # Convert SMILES to 3D
        smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
        
        print(f"\n1. Converting SMILES to 3D: {smiles}")
        result = await client.call_tool("generate_3d", {
            "smiles": smiles,
            "force_field": "mmff94",
            "steps": 1000
        })
        
        if "error" not in result:
            print("3D generation successful!")
            print(f"  Energy: {result.get('energy', 'N/A')}")
            # Show first few lines of XYZ
            xyz_lines = result['xyz'].split('\n')[:5]
            print("  XYZ coordinates (first 5 lines):")
            for line in xyz_lines:
                print(f"    {line}")
        
        # Check drug-likeness
        print("\n2. Checking drug-likeness...")
        drug_check = await client.call_tool("check_drug_like", {
            "smiles": smiles
        })
        
        if "error" not in drug_check:
            print(f"  Molecular Weight: {drug_check['molecular_weight']}")
            print(f"  LogP: {drug_check['logp']}")
            print(f"  H-Bond Donors: {drug_check['h_bond_donors']}")
            print(f"  H-Bond Acceptors: {drug_check['h_bond_acceptors']}")
            print(f"  Lipinski Violations: {drug_check['lipinski_violations']}")
            print(f"  Is Drug-like: {'✅' if drug_check['is_drug_like'] else '❌'}")
        
        # Canonicalize SMILES
        print("\n3. Canonicalizing SMILES...")
        test_smiles = ["C(C)O", "OCC", "C(O)C"]
        for smi in test_smiles:
            result = await client.call_tool("canonicalize_smiles", {
                "smiles": smi
            })
            if "error" not in result:
                print(f"  {smi} → {result['canonical_smiles']}")


async def batch_example():
    """Example of batch processing with MCP servers"""
    print("\n" + "="*50)
    print("Batch Processing Example")
    print("="*50)
    
    # List of drug molecules to analyze
    drugs = [
        ("Aspirin", "CC(=O)OC1=CC=CC=C1C(=O)O"),
        ("Ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"),
        ("Paracetamol", "CC(=O)NC1=CC=C(O)C=C1"),
    ]
    
    async with MCPClient("http://localhost:8769") as client:
        print("\nAnalyzing drug molecules:")
        
        for name, smiles in drugs:
            print(f"\n{name}:")
            
            # Calculate properties
            props = await client.call_tool("calculate_properties", {
                "molecule": smiles,
                "format": "smiles"
            })
            
            if "error" not in props:
                print(f"  Formula: {props['formula']}")
                print(f"  MW: {props['molecular_weight']:.2f}")
                print(f"  Rotatable bonds: {props['num_rotatable_bonds']}")
            
            # Check drug-likeness
            drug_check = await client.call_tool("check_drug_like", {
                "smiles": smiles
            })
            
            if "error" not in drug_check:
                print(f"  Drug-like: {'✅' if drug_check['is_drug_like'] else '❌'}")
                if drug_check['lipinski_violations'] > 0:
                    print(f"  Lipinski violations: {drug_check['lipinski_violations']}")


async def integrated_workflow():
    """
    Integrated workflow using multiple MCP servers
    Example: Find a drug, get its structure, and analyze it
    """
    print("\n" + "="*50)
    print("Integrated Workflow: Drug Discovery Pipeline")
    print("="*50)
    
    drug_name = "imatinib"
    
    # Step 1: Search in PubChem
    print(f"\n1. Searching PubChem for {drug_name}...")
    async with MCPClient("http://localhost:8767") as pubchem:
        search_result = await pubchem.call_tool("search_compound", {
            "query": drug_name,
            "search_type": "name",
            "max_results": 1
        })
        
        if search_result['count'] > 0:
            compound = search_result['compounds'][0]
            smiles = compound['smiles']
            print(f"  Found: {compound['iupac']}")
            print(f"  SMILES: {smiles}")
        else:
            print("  Not found in PubChem")
            return
    
    # Step 2: Search in ChEMBL for bioactivity
    print(f"\n2. Searching ChEMBL for bioactivity data...")
    async with MCPClient("http://localhost:8768") as chembl:
        mol_result = await chembl.call_tool("search_molecule", {
            "query": drug_name,
            "search_type": "name",
            "max_results": 1
        })
        
        if mol_result['count'] > 0:
            chembl_id = mol_result['molecules'][0]['chembl_id']
            print(f"  ChEMBL ID: {chembl_id}")
            
            # Get bioactivity data
            bioactivity = await chembl.call_tool("get_bioactivity", {
                "molecule_chembl_id": chembl_id,
                "min_pchembl": 6,
                "max_results": 5
            })
            
            print(f"  Found {bioactivity['count']} bioactivities")
            for activity in bioactivity['bioactivities'][:3]:
                print(f"    - Target: {activity['target_chembl_id']}")
                print(f"      Type: {activity['standard_type']}")
                print(f"      Value: {activity['standard_value']} {activity['standard_units']}")
    
    # Step 3: Generate 3D structure and analyze
    print(f"\n3. Generating 3D structure and analyzing...")
    async with MCPClient("http://localhost:8769") as openbabel:
        # Generate 3D
        structure_3d = await openbabel.call_tool("generate_3d", {
            "smiles": smiles,
            "force_field": "mmff94"
        })
        
        if "error" not in structure_3d:
            print(f"  3D structure generated")
            print(f"  Energy: {structure_3d.get('energy', 'N/A')}")
        
        # Check drug properties
        drug_check = await openbabel.call_tool("check_drug_like", {
            "smiles": smiles
        })
        
        if "error" not in drug_check:
            print(f"\n4. Drug-likeness assessment:")
            print(f"  Molecular Weight: {drug_check['molecular_weight']:.2f}")
            print(f"  LogP: {drug_check['logp']:.2f}")
            print(f"  Lipinski's Rule of Five: {'✅ Pass' if drug_check['is_drug_like'] else '❌ Fail'}")


async def main():
    """Run all examples"""
    
    print("\n" + "="*60)
    print("   Chemistry MCP Servers - Client Examples")
    print("="*60)
    
    print("\nMake sure MCP servers are running:")
    print("  mcp start-all")
    print("\nOr start individually:")
    print("  python -m mcp_pubchem_server --tcp 8767")
    print("  python -m mcp_chembl_server --tcp 8768")
    print("  python -m mcp_openbabel_server --tcp 8769")
    
    try:
        # Run examples
        await pubchem_example()
        await chembl_example()
        await openbabel_example()
        await batch_example()
        await integrated_workflow()
        
    except aiohttp.ClientConnectorError:
        print("\n❌ Could not connect to MCP server")
        print("Make sure the servers are running")
    except Exception as e:
        print(f"\n❌ Error: {e}")


if __name__ == "__main__":
    # Install aiohttp if needed
    try:
        import aiohttp
    except ImportError:
        print("Installing aiohttp...")
        import subprocess
        subprocess.run([sys.executable, "-m", "pip", "install", "aiohttp"])
        import aiohttp
    
    asyncio.run(main())
