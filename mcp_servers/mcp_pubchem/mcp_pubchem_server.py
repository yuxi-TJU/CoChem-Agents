#!/usr/bin/env python3
"""
PubChem MCP Server
Provides MCP interface to PubChem database
"""

import asyncio
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from base_mcp_server import BaseMCPServer, MCPRequest, MCPResponse, mcp_tool
import pubchempy as pcp
from typing import Any, Dict, List, Optional
import logging

logger = logging.getLogger(__name__)


class PubChemMCPServer(BaseMCPServer):
    """
    MCP Server for PubChem database access
    """
    
    def __init__(self):
        super().__init__("pubchem-mcp", "0.1.0")
        
    def _setup_capabilities(self):
        """Setup PubChem-specific capabilities"""
        self.capabilities = {
            "search": {
                "by_name": True,
                "by_smiles": True,
                "by_inchi": True,
                "by_formula": True,
                "by_cid": True
            },
            "properties": {
                "basic": ["molecular_weight", "molecular_formula", "smiles", "inchi"],
                "computed": ["xlogp", "tpsa", "complexity", "h_bond_donor_count"],
                "identifiers": ["iupac_name", "synonyms", "cas"]
            },
            "bioassays": True,
            "patents": True,
            "literature": True,
            "similar_compounds": True
        }
    
    async def handle_request(self, request: MCPRequest) -> MCPResponse:
        """Handle PubChem-specific requests"""
        
        method_map = {
            "search_compound": self.tool_search_compound,
            "get_properties": self.tool_get_properties,
            "get_synonyms": self.tool_get_synonyms,
            "get_bioassays": self.tool_get_bioassays,
            "find_similar": self.tool_find_similar,
            "get_patents": self.tool_get_patents,
        }
        
        if request.method in method_map:
            try:
                result = await method_map[request.method](**request.params)
                return MCPResponse(result=result, id=request.id)
            except Exception as e:
                logger.error(f"Error handling {request.method}: {e}")
                return MCPResponse(
                    error={
                        "code": -32603,
                        "message": str(e)
                    },
                    id=request.id
                )
        else:
            return MCPResponse(
                error={
                    "code": -32601,
                    "message": f"Method not found: {request.method}"
                },
                id=request.id
            )
    
    @mcp_tool({
        "query": {"type": "string", "required": True},
        "search_type": {"type": "string", "enum": ["name", "smiles", "inchi", "formula"], "default": "name"},
        "max_results": {"type": "integer", "default": 5}
    })
    async def tool_search_compound(self, query: str, search_type: str = "name", max_results: int = 5) -> Dict[str, Any]:
        """Search for compounds in PubChem"""
        try:
            compounds = []
            
            if search_type == "name":
                results = pcp.get_compounds(query, 'name')
            elif search_type == "smiles":
                results = pcp.get_compounds(query, 'smiles')
            elif search_type == "inchi":
                results = pcp.get_compounds(query, 'inchi')
            elif search_type == "formula":
                results = pcp.get_compounds(query, 'formula')
            else:
                return {"error": f"Invalid search type: {search_type}"}
            
            for compound in results[:max_results]:
                compounds.append({
                    "cid": compound.cid,
                    "smiles": compound.canonical_smiles,
                    "formula": compound.molecular_formula,
                    "weight": compound.molecular_weight,
                    "iupac": compound.iupac_name
                })
            
            return {
                "query": query,
                "search_type": search_type,
                "count": len(compounds),
                "compounds": compounds
            }
            
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "cid": {"type": "integer", "required": True},
        "properties": {"type": "array", "items": {"type": "string"}}
    })
    async def tool_get_properties(self, cid: int, properties: List[str] = None) -> Dict[str, Any]:
        """Get properties for a compound by CID"""
        try:
            compound = pcp.Compound.from_cid(cid)
            
            if not compound:
                return {"error": f"Compound with CID {cid} not found"}
            
            # Default properties
            if not properties:
                properties = ["molecular_weight", "molecular_formula", "canonical_smiles", 
                            "xlogp", "tpsa", "complexity", "h_bond_donor_count", 
                            "h_bond_acceptor_count", "rotatable_bond_count"]
            
            result = {"cid": cid}
            
            for prop in properties:
                try:
                    value = getattr(compound, prop, None)
                    if value is not None:
                        result[prop] = value
                except:
                    result[prop] = None
            
            return result
            
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "cid": {"type": "integer", "required": True}
    })
    async def tool_get_synonyms(self, cid: int) -> Dict[str, Any]:
        """Get synonyms for a compound"""
        try:
            compound = pcp.Compound.from_cid(cid)
            
            if not compound:
                return {"error": f"Compound with CID {cid} not found"}
            
            synonyms = compound.synonyms or []
            
            return {
                "cid": cid,
                "iupac_name": compound.iupac_name,
                "synonyms": synonyms[:20],  # Limit to 20 synonyms
                "total_synonyms": len(synonyms)
            }
            
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "cid": {"type": "integer", "required": True}
    })
    async def tool_get_bioassays(self, cid: int) -> Dict[str, Any]:
        """Get bioassay data for a compound"""
        try:
            # Note: This is a simplified version
            # Real implementation would query PubChem bioassay database
            
            return {
                "cid": cid,
                "message": "Bioassay data retrieval requires additional API setup",
                "suggestion": "Use PubChem REST API for detailed bioassay data"
            }
            
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "cid": {"type": "integer", "required": True},
        "threshold": {"type": "number", "default": 0.9}
    })
    async def tool_find_similar(self, cid: int, threshold: float = 0.9) -> Dict[str, Any]:
        """Find similar compounds"""
        try:
            # Note: This requires additional setup for similarity search
            # Using PubChem's similarity search API
            
            return {
                "cid": cid,
                "threshold": threshold,
                "message": "Similarity search requires PubChem REST API",
                "api_endpoint": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/{cid}/cids/JSON"
            }
            
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "cid": {"type": "integer", "required": True}
    })
    async def tool_get_patents(self, cid: int) -> Dict[str, Any]:
        """Get patent information for a compound"""
        try:
            # Note: Patent data requires additional API setup
            
            return {
                "cid": cid,
                "message": "Patent data requires PubChem patent API",
                "api_endpoint": f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
            }
            
        except Exception as e:
            return {"error": str(e)}


async def main():
    """Run the PubChem MCP server"""
    server = PubChemMCPServer()
    
    # Check command line arguments
    if len(sys.argv) > 1 and sys.argv[1] == "--tcp":
        # Run as TCP server for testing
        port = int(sys.argv[2]) if len(sys.argv) > 2 else 8766
        await server.run_tcp(port=port)
    else:
        # Run as stdio server (standard MCP mode)
        await server.run_stdio()


if __name__ == "__main__":
    asyncio.run(main())
