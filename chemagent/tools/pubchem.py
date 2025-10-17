"""PubChem database integration"""

from typing import Dict, List, Any, Optional
import logging
import asyncio
import pubchempy as pcp
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

from ..core.registry import BaseTool, ToolMetadata, ToolCategory


logger = logging.getLogger(__name__)


class PubChemTool(BaseTool):
    """Tool for querying PubChem database"""
    
    def __init__(self, config: Optional[Dict[str, Any]] = None):
        super().__init__(config)
        self.session = self._create_session()
        
    def _create_session(self):
        """Create HTTP session with retry logic"""
        session = requests.Session()
        retry = Retry(
            total=3,
            read=3,
            connect=3,
            backoff_factor=0.3,
            status_forcelist=(500, 502, 504)
        )
        adapter = HTTPAdapter(max_retries=retry)
        session.mount('http://', adapter)
        session.mount('https://', adapter)
        return session
        
    def get_metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="pubchem",
            description="Query PubChem database for chemical information",
            category=ToolCategory.DATABASE,
            version="1.0.0",
            author="ChemAgent",
            dependencies=["pubchempy", "requests"],
            parameters={
                "required": ["query", "query_type"],
                "properties": {
                    "query": {"type": "string", "description": "Search query"},
                    "query_type": {
                        "type": "string",
                        "enum": ["name", "smiles", "inchi", "cid", "formula"],
                        "description": "Type of query"
                    },
                    "properties": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Properties to retrieve",
                        "default": ["MolecularFormula", "MolecularWeight", "CanonicalSMILES"]
                    },
                    "limit": {
                        "type": "integer",
                        "description": "Maximum number of results",
                        "default": 5
                    }
                }
            },
            examples=[
                {"query": "aspirin", "query_type": "name"},
                {"query": "CC(=O)OC1=CC=CC=C1C(=O)O", "query_type": "smiles"},
                {"query": "C9H8O4", "query_type": "formula"}
            ]
        )
        
    async def execute(self,
                      query: str,
                      query_type: str,
                      properties: Optional[List[str]] = None,
                      limit: int = 5,
                      **kwargs) -> Dict[str, Any]:
        """Query PubChem database"""
        try:
            # Default properties if not specified
            if not properties:
                properties = [
                    "MolecularFormula",
                    "MolecularWeight", 
                    "CanonicalSMILES",
                    "IUPACName",
                    "XLogP",
                    "TPSA",
                    "Complexity"
                ]
                
            # Run query in executor to avoid blocking
            loop = asyncio.get_event_loop()
            compounds = await loop.run_in_executor(
                None,
                self._search_compounds,
                query,
                query_type,
                limit
            )
            
            if not compounds:
                return {
                    "query": query,
                    "query_type": query_type,
                    "results": [],
                    "message": "No compounds found"
                }
                
            # Extract properties
            results = []
            for compound in compounds[:limit]:
                comp_data = await loop.run_in_executor(
                    None,
                    self._extract_properties,
                    compound,
                    properties
                )
                results.append(comp_data)
                
            return {
                "query": query,
                "query_type": query_type,
                "num_results": len(results),
                "results": results
            }
            
        except Exception as e:
            logger.error(f"PubChem query failed: {e}")
            return {"error": str(e)}
            
    def _search_compounds(self, query: str, query_type: str, limit: int):
        """Search for compounds in PubChem"""
        try:
            if query_type == "name":
                compounds = pcp.get_compounds(query, 'name', listkey_count=limit)
            elif query_type == "smiles":
                compounds = pcp.get_compounds(query, 'smiles', listkey_count=limit)
            elif query_type == "inchi":
                compounds = pcp.get_compounds(query, 'inchi', listkey_count=limit)
            elif query_type == "cid":
                compounds = [pcp.Compound.from_cid(int(query))]
            elif query_type == "formula":
                compounds = pcp.get_compounds(query, 'formula', listkey_count=limit)
            else:
                compounds = []
                
            return compounds
            
        except Exception as e:
            logger.error(f"PubChem search failed: {e}")
            return []
            
    def _extract_properties(self, compound, properties: List[str]) -> Dict[str, Any]:
        """Extract properties from PubChem compound"""
        data = {
            "cid": compound.cid,
            "name": compound.iupac_name or compound.synonyms[0] if compound.synonyms else "Unknown"
        }
        
        # Map property names to PubChemPy attributes
        property_map = {
            "MolecularFormula": "molecular_formula",
            "MolecularWeight": "molecular_weight",
            "CanonicalSMILES": "canonical_smiles",
            "IsomericSMILES": "isomeric_smiles",
            "IUPACName": "iupac_name",
            "XLogP": "xlogp",
            "TPSA": "tpsa",
            "Complexity": "complexity",
            "HBondDonorCount": "h_bond_donor_count",
            "HBondAcceptorCount": "h_bond_acceptor_count",
            "RotatableBondCount": "rotatable_bond_count",
            "HeavyAtomCount": "heavy_atom_count",
            "AtomStereoCount": "atom_stereo_count",
            "DefinedAtomStereoCount": "defined_atom_stereo_count",
            "UndefinedAtomStereoCount": "undefined_atom_stereo_count",
            "BondStereoCount": "bond_stereo_count",
            "CovalentUnitCount": "covalent_unit_count"
        }
        
        for prop in properties:
            if prop in property_map:
                attr_name = property_map[prop]
                try:
                    value = getattr(compound, attr_name, None)
                    if value is not None:
                        data[prop] = value
                except:
                    pass
                    
        # Add synonyms if available
        if compound.synonyms:
            data["synonyms"] = compound.synonyms[:5]  # Limit to 5 synonyms
            
        return data
        
    async def get_bioactivity(self, cid: int, **kwargs) -> Dict[str, Any]:
        """Get bioactivity data for a compound"""
        try:
            loop = asyncio.get_event_loop()
            
            # Get assay data
            assays = await loop.run_in_executor(
                None,
                pcp.get_assays,
                cid,
                'cid'
            )
            
            bioactivity_data = {
                "cid": cid,
                "num_assays": len(assays) if assays else 0,
                "assays": []
            }
            
            if assays:
                for assay in assays[:10]:  # Limit to 10 assays
                    assay_info = {
                        "aid": assay.aid,
                        "name": assay.name,
                        "description": assay.description[:200] if assay.description else None,
                        "activity_outcome": assay.activity_outcome
                    }
                    bioactivity_data["assays"].append(assay_info)
                    
            return bioactivity_data
            
        except Exception as e:
            logger.error(f"Failed to get bioactivity data: {e}")
            return {"error": str(e)}
            
    async def get_similar_compounds(self,
                                   smiles: str,
                                   threshold: float = 0.9,
                                   limit: int = 10,
                                   **kwargs) -> Dict[str, Any]:
        """Find similar compounds using structure similarity"""
        try:
            loop = asyncio.get_event_loop()
            
            # Search for similar compounds
            similar = await loop.run_in_executor(
                None,
                lambda: pcp.get_compounds(
                    smiles,
                    'smiles',
                    searchtype='similarity',
                    threshold=threshold,
                    listkey_count=limit
                )
            )
            
            results = []
            for compound in similar:
                comp_data = {
                    "cid": compound.cid,
                    "name": compound.iupac_name or (compound.synonyms[0] if compound.synonyms else "Unknown"),
                    "smiles": compound.canonical_smiles,
                    "molecular_weight": compound.molecular_weight,
                    "molecular_formula": compound.molecular_formula
                }
                results.append(comp_data)
                
            return {
                "query_smiles": smiles,
                "threshold": threshold,
                "num_results": len(results),
                "similar_compounds": results
            }
            
        except Exception as e:
            logger.error(f"Similarity search failed: {e}")
            return {"error": str(e)}
