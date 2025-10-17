#!/usr/bin/env python3
"""
ChEMBL MCP Server
Provides MCP interface to ChEMBL bioactivity database
"""

import asyncio
import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from base_mcp_server import BaseMCPServer, MCPRequest, MCPResponse, mcp_tool
from chembl_webresource_client.new_client import new_client
from typing import Any, Dict, List, Optional
import logging

logger = logging.getLogger(__name__)


class ChEMBLMCPServer(BaseMCPServer):
    """
    MCP Server for ChEMBL database access
    """
    
    def __init__(self):
        super().__init__("chembl-mcp", "0.1.0")
        self.molecule = new_client.molecule
        self.target = new_client.target
        self.activity = new_client.activity
        self.assay = new_client.assay
        self.drug = new_client.drug
        
    def _setup_capabilities(self):
        """Setup ChEMBL-specific capabilities"""
        self.capabilities = {
            "search": {
                "molecules": True,
                "targets": True,
                "drugs": True,
                "assays": True
            },
            "bioactivity": {
                "by_molecule": True,
                "by_target": True,
                "ic50_values": True,
                "ki_values": True
            },
            "drug_data": {
                "approved_drugs": True,
                "clinical_candidates": True,
                "drug_mechanisms": True
            },
            "target_data": {
                "protein_targets": True,
                "target_families": True,
                "organism_filter": True
            }
        }
    
    async def handle_request(self, request: MCPRequest) -> MCPResponse:
        """Handle ChEMBL-specific requests"""
        
        method_map = {
            "search_molecule": self.tool_search_molecule,
            "search_target": self.tool_search_target,
            "get_bioactivity": self.tool_get_bioactivity,
            "get_drug_info": self.tool_get_drug_info,
            "get_target_info": self.tool_get_target_info,
            "find_drugs_for_target": self.tool_find_drugs_for_target,
            "get_assay_data": self.tool_get_assay_data,
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
        "search_type": {"type": "string", "enum": ["smiles", "name", "chembl_id"], "default": "name"},
        "max_results": {"type": "integer", "default": 10}
    })
    async def tool_search_molecule(self, query: str, search_type: str = "name", max_results: int = 10) -> Dict[str, Any]:
        """Search for molecules in ChEMBL"""
        try:
            results = []
            
            if search_type == "smiles":
                molecules = self.molecule.filter(molecule_structures__canonical_smiles=query)
            elif search_type == "chembl_id":
                molecules = self.molecule.filter(molecule_chembl_id=query)
            else:  # name search
                molecules = self.molecule.filter(pref_name__icontains=query)
            
            for mol in molecules[:max_results]:
                results.append({
                    "chembl_id": mol['molecule_chembl_id'],
                    "name": mol.get('pref_name', 'N/A'),
                    "smiles": mol['molecule_structures']['canonical_smiles'] if 'molecule_structures' in mol else None,
                    "molecular_weight": mol['molecule_properties']['mw_freebase'] if 'molecule_properties' in mol else None,
                    "max_phase": mol.get('max_phase', 0)
                })
            
            return {
                "query": query,
                "search_type": search_type,
                "count": len(results),
                "molecules": results
            }
            
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "query": {"type": "string", "required": True},
        "organism": {"type": "string", "default": "Homo sapiens"},
        "max_results": {"type": "integer", "default": 10}
    })
    async def tool_search_target(self, query: str, organism: str = "Homo sapiens", max_results: int = 10) -> Dict[str, Any]:
        """Search for targets in ChEMBL"""
        try:
            results = []
            
            targets = self.target.filter(
                pref_name__icontains=query,
                organism__icontains=organism
            )
            
            for target in targets[:max_results]:
                results.append({
                    "chembl_id": target['target_chembl_id'],
                    "name": target.get('pref_name', 'N/A'),
                    "type": target.get('target_type', 'Unknown'),
                    "organism": target.get('organism', 'N/A'),
                    "gene_names": target.get('target_components', [])
                })
            
            return {
                "query": query,
                "organism": organism,
                "count": len(results),
                "targets": results
            }
            
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "molecule_chembl_id": {"type": "string", "required": False},
        "target_chembl_id": {"type": "string", "required": False},
        "min_pchembl": {"type": "number", "default": 5},
        "max_results": {"type": "integer", "default": 50}
    })
    async def tool_get_bioactivity(self, molecule_chembl_id: str = None, target_chembl_id: str = None, 
                                   min_pchembl: float = 5, max_results: int = 50) -> Dict[str, Any]:
        """Get bioactivity data from ChEMBL"""
        try:
            filters = {}
            
            if molecule_chembl_id:
                filters['molecule_chembl_id'] = molecule_chembl_id
            if target_chembl_id:
                filters['target_chembl_id'] = target_chembl_id
            if min_pchembl:
                filters['pchembl_value__gte'] = min_pchembl
            
            if not filters:
                return {"error": "Either molecule_chembl_id or target_chembl_id must be provided"}
            
            activities = self.activity.filter(**filters)
            
            results = []
            for act in activities[:max_results]:
                results.append({
                    "molecule_chembl_id": act.get('molecule_chembl_id'),
                    "target_chembl_id": act.get('target_chembl_id'),
                    "standard_type": act.get('standard_type'),
                    "standard_value": act.get('standard_value'),
                    "standard_units": act.get('standard_units'),
                    "pchembl_value": act.get('pchembl_value'),
                    "assay_description": act.get('assay_description')
                })
            
            return {
                "filters": filters,
                "count": len(results),
                "bioactivities": results
            }
            
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "chembl_id": {"type": "string", "required": True}
    })
    async def tool_get_drug_info(self, chembl_id: str) -> Dict[str, Any]:
        """Get drug information from ChEMBL"""
        try:
            drug_data = self.drug.filter(molecule_chembl_id=chembl_id)
            
            if not drug_data:
                return {"error": f"No drug data found for {chembl_id}"}
            
            drug = drug_data[0]
            
            return {
                "chembl_id": chembl_id,
                "name": drug.get('pref_name'),
                "synonyms": drug.get('synonyms', []),
                "max_phase": drug.get('max_phase'),
                "indication": drug.get('indication_class'),
                "mechanism": drug.get('drug_mechanisms', []),
                "atc_classifications": drug.get('atc_classifications', [])
            }
            
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "chembl_id": {"type": "string", "required": True}
    })
    async def tool_get_target_info(self, chembl_id: str) -> Dict[str, Any]:
        """Get target information from ChEMBL"""
        try:
            target_data = self.target.filter(target_chembl_id=chembl_id)
            
            if not target_data:
                return {"error": f"No target data found for {chembl_id}"}
            
            target = target_data[0]
            
            return {
                "chembl_id": chembl_id,
                "name": target.get('pref_name'),
                "type": target.get('target_type'),
                "organism": target.get('organism'),
                "gene_names": target.get('target_components', []),
                "protein_accessions": [comp.get('accession') for comp in target.get('target_components', [])]
            }
            
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "target_chembl_id": {"type": "string", "required": True},
        "min_phase": {"type": "integer", "default": 3}
    })
    async def tool_find_drugs_for_target(self, target_chembl_id: str, min_phase: int = 3) -> Dict[str, Any]:
        """Find drugs targeting a specific target"""
        try:
            # Get bioactivity data for the target
            activities = self.activity.filter(
                target_chembl_id=target_chembl_id,
                pchembl_value__gte=6  # Good potency
            )
            
            molecule_ids = set()
            for act in activities[:100]:  # Limit to prevent too many queries
                molecule_ids.add(act['molecule_chembl_id'])
            
            drugs = []
            for mol_id in list(molecule_ids)[:20]:  # Limit results
                mol = self.molecule.filter(molecule_chembl_id=mol_id)
                if mol and mol[0].get('max_phase', 0) >= min_phase:
                    drugs.append({
                        "chembl_id": mol_id,
                        "name": mol[0].get('pref_name'),
                        "max_phase": mol[0].get('max_phase'),
                        "first_approval": mol[0].get('first_approval')
                    })
            
            return {
                "target_chembl_id": target_chembl_id,
                "min_phase": min_phase,
                "count": len(drugs),
                "drugs": drugs
            }
            
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "assay_chembl_id": {"type": "string", "required": True}
    })
    async def tool_get_assay_data(self, assay_chembl_id: str) -> Dict[str, Any]:
        """Get assay information from ChEMBL"""
        try:
            assay_data = self.assay.filter(assay_chembl_id=assay_chembl_id)
            
            if not assay_data:
                return {"error": f"No assay data found for {assay_chembl_id}"}
            
            assay = assay_data[0]
            
            return {
                "chembl_id": assay_chembl_id,
                "description": assay.get('description'),
                "assay_type": assay.get('assay_type'),
                "target": assay.get('target_chembl_id'),
                "organism": assay.get('assay_organism'),
                "confidence_score": assay.get('confidence_score')
            }
            
        except Exception as e:
            return {"error": str(e)}


async def main():
    """Run the ChEMBL MCP server"""
    server = ChEMBLMCPServer()
    
    # Check command line arguments
    if len(sys.argv) > 1 and sys.argv[1] == "--tcp":
        # Run as TCP server for testing
        port = int(sys.argv[2]) if len(sys.argv) > 2 else 8767
        await server.run_tcp(port=port)
    else:
        # Run as stdio server (standard MCP mode)
        await server.run_stdio()


if __name__ == "__main__":
    asyncio.run(main())
