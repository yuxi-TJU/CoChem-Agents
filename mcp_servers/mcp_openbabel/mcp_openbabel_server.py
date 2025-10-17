#!/usr/bin/env python3
"""
OpenBabel MCP Server
Provides MCP interface for molecular format conversion and 3D generation
"""

import asyncio
import sys
import os
import tempfile
import subprocess
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from base_mcp_server import BaseMCPServer, MCPRequest, MCPResponse, mcp_tool
from typing import Any, Dict, List, Optional
import logging

logger = logging.getLogger(__name__)

# Try to import openbabel
try:
    from openbabel import openbabel as ob
    from openbabel import pybel
    OPENBABEL_AVAILABLE = True
except ImportError:
    OPENBABEL_AVAILABLE = False
    logger.warning("OpenBabel Python bindings not available, will use CLI if available")


class OpenBabelMCPServer(BaseMCPServer):
    """
    MCP Server for OpenBabel molecular conversion
    """
    
    def __init__(self):
        super().__init__("openbabel-mcp", "0.1.0")
        self.check_openbabel()
        
    def check_openbabel(self):
        """Check OpenBabel availability"""
        self.has_python = OPENBABEL_AVAILABLE
        self.has_cli = self._check_cli()
        
        if not self.has_python and not self.has_cli:
            logger.error("Neither OpenBabel Python bindings nor CLI found!")
            logger.error("Install with: pip install openbabel-wheel or apt-get install openbabel")
    
    def _check_cli(self) -> bool:
        """Check if obabel CLI is available"""
        try:
            result = subprocess.run(["obabel", "--help"], 
                                  capture_output=True, 
                                  text=True, 
                                  timeout=5)
            return result.returncode == 0
        except:
            return False
    
    def _setup_capabilities(self):
        """Setup OpenBabel-specific capabilities"""
        self.capabilities = {
            "conversion": {
                "formats": ["smiles", "mol", "mol2", "sdf", "pdb", "xyz", "inchi", "can"],
                "2d_to_3d": True,
                "3d_optimization": True,
                "add_hydrogens": True,
                "remove_hydrogens": True
            },
            "generation": {
                "conformers": True,
                "3d_coordinates": True,
                "canonical_smiles": True,
                "inchi_key": True
            },
            "properties": {
                "molecular_weight": True,
                "formula": True,
                "charge": True,
                "spin_multiplicity": True,
                "descriptors": True
            },
            "filters": {
                "lipinski": True,
                "drug_like": True,
                "lead_like": True
            }
        }
    
    async def handle_request(self, request: MCPRequest) -> MCPResponse:
        """Handle OpenBabel-specific requests"""
        
        method_map = {
            "convert_molecule": self.tool_convert_molecule,
            "generate_3d": self.tool_generate_3d,
            "optimize_geometry": self.tool_optimize_geometry,
            "generate_conformers": self.tool_generate_conformers,
            "calculate_properties": self.tool_calculate_properties,
            "canonicalize_smiles": self.tool_canonicalize_smiles,
            "add_hydrogens": self.tool_add_hydrogens,
            "check_drug_like": self.tool_check_drug_like,
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
        "molecule": {"type": "string", "required": True},
        "input_format": {"type": "string", "required": True},
        "output_format": {"type": "string", "required": True},
        "options": {"type": "object", "default": {}}
    })
    async def tool_convert_molecule(self, molecule: str, input_format: str, 
                                   output_format: str, options: Dict = None) -> Dict[str, Any]:
        """Convert molecule between different formats"""
        try:
            if self.has_python:
                # Use Python bindings
                mol = pybel.readstring(input_format, molecule)
                
                # Apply options
                if options:
                    if options.get("add_hydrogens"):
                        mol.addh()
                    if options.get("remove_hydrogens"):
                        mol.removeh()
                    if options.get("gen3d"):
                        mol.make3D()
                
                output = mol.write(output_format)
                
                return {
                    "input_format": input_format,
                    "output_format": output_format,
                    "converted_molecule": output,
                    "method": "python"
                }
                
            elif self.has_cli:
                # Use CLI as fallback
                with tempfile.NamedTemporaryFile(mode='w', suffix=f'.{input_format}', delete=False) as f:
                    f.write(molecule)
                    input_file = f.name
                
                output_file = input_file.replace(f'.{input_format}', f'.{output_format}')
                
                cmd = ["obabel", input_file, "-O", output_file]
                if options:
                    if options.get("add_hydrogens"):
                        cmd.append("-h")
                    if options.get("gen3d"):
                        cmd.append("--gen3d")
                
                subprocess.run(cmd, check=True, capture_output=True)
                
                with open(output_file, 'r') as f:
                    output = f.read()
                
                # Cleanup
                os.unlink(input_file)
                os.unlink(output_file)
                
                return {
                    "input_format": input_format,
                    "output_format": output_format,
                    "converted_molecule": output,
                    "method": "cli"
                }
            else:
                return {"error": "OpenBabel not available"}
                
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "smiles": {"type": "string", "required": True},
        "force_field": {"type": "string", "enum": ["mmff94", "uff", "ghemical"], "default": "mmff94"},
        "steps": {"type": "integer", "default": 500}
    })
    async def tool_generate_3d(self, smiles: str, force_field: str = "mmff94", steps: int = 500) -> Dict[str, Any]:
        """Generate 3D coordinates from SMILES"""
        try:
            if self.has_python:
                mol = pybel.readstring("smiles", smiles)
                mol.make3D(forcefield=force_field, steps=steps)
                
                return {
                    "smiles": smiles,
                    "mol_block": mol.write("mol"),
                    "xyz": mol.write("xyz"),
                    "energy": mol.energy if hasattr(mol, 'energy') else None
                }
            else:
                return {"error": "3D generation requires OpenBabel Python bindings"}
                
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "molecule": {"type": "string", "required": True},
        "format": {"type": "string", "default": "mol"},
        "force_field": {"type": "string", "default": "mmff94"},
        "steps": {"type": "integer", "default": 2500}
    })
    async def tool_optimize_geometry(self, molecule: str, format: str = "mol", 
                                    force_field: str = "mmff94", steps: int = 2500) -> Dict[str, Any]:
        """Optimize molecular geometry"""
        try:
            if self.has_python:
                mol = pybel.readstring(format, molecule)
                mol.localopt(forcefield=force_field, steps=steps)
                
                return {
                    "optimized_molecule": mol.write("mol"),
                    "energy": mol.energy if hasattr(mol, 'energy') else None,
                    "convergence": "completed"
                }
            else:
                return {"error": "Geometry optimization requires OpenBabel Python bindings"}
                
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "smiles": {"type": "string", "required": True},
        "num_conformers": {"type": "integer", "default": 10},
        "force_field": {"type": "string", "default": "mmff94"}
    })
    async def tool_generate_conformers(self, smiles: str, num_conformers: int = 10, 
                                      force_field: str = "mmff94") -> Dict[str, Any]:
        """Generate multiple conformers"""
        try:
            if self.has_python:
                mol = pybel.readstring("smiles", smiles)
                mol.make3D()
                
                # Note: Simplified conformer generation
                conformers = []
                for i in range(min(num_conformers, 10)):
                    mol.localopt(forcefield=force_field, steps=500)
                    conformers.append({
                        "id": i + 1,
                        "mol_block": mol.write("mol"),
                        "energy": mol.energy if hasattr(mol, 'energy') else None
                    })
                
                return {
                    "smiles": smiles,
                    "num_conformers": len(conformers),
                    "conformers": conformers
                }
            else:
                return {"error": "Conformer generation requires OpenBabel Python bindings"}
                
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "molecule": {"type": "string", "required": True},
        "format": {"type": "string", "default": "smiles"}
    })
    async def tool_calculate_properties(self, molecule: str, format: str = "smiles") -> Dict[str, Any]:
        """Calculate molecular properties"""
        try:
            if self.has_python:
                mol = pybel.readstring(format, molecule)
                
                return {
                    "molecular_weight": mol.molwt,
                    "formula": mol.formula,
                    "num_atoms": len(mol.atoms),
                    "num_bonds": len(mol.OBMol.GetBonds()) if hasattr(mol.OBMol, 'GetBonds') else None,
                    "num_rotatable_bonds": mol.OBMol.NumRotors(),
                    "charge": mol.charge,
                    "spin": mol.spin,
                    "inchi": mol.write("inchi").strip(),
                    "inchikey": mol.write("inchikey").strip()
                }
            else:
                return {"error": "Property calculation requires OpenBabel Python bindings"}
                
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "smiles": {"type": "string", "required": True}
    })
    async def tool_canonicalize_smiles(self, smiles: str) -> Dict[str, Any]:
        """Generate canonical SMILES"""
        try:
            if self.has_python:
                mol = pybel.readstring("smiles", smiles)
                canonical = mol.write("can").strip()
                
                return {
                    "input_smiles": smiles,
                    "canonical_smiles": canonical,
                    "is_valid": True
                }
            elif self.has_cli:
                result = subprocess.run(
                    ["obabel", "-ismi", "-ocan"],
                    input=smiles,
                    capture_output=True,
                    text=True
                )
                
                return {
                    "input_smiles": smiles,
                    "canonical_smiles": result.stdout.strip(),
                    "is_valid": result.returncode == 0
                }
            else:
                return {"error": "OpenBabel not available"}
                
        except Exception as e:
            return {"error": str(e), "is_valid": False}
    
    @mcp_tool({
        "molecule": {"type": "string", "required": True},
        "format": {"type": "string", "default": "smiles"},
        "ph": {"type": "number", "default": 7.4}
    })
    async def tool_add_hydrogens(self, molecule: str, format: str = "smiles", ph: float = 7.4) -> Dict[str, Any]:
        """Add hydrogens to molecule"""
        try:
            if self.has_python:
                mol = pybel.readstring(format, molecule)
                mol.OBMol.AddHydrogens(False, True, ph)
                
                return {
                    "molecule_with_h": mol.write(format),
                    "num_hydrogens": sum(1 for atom in mol.atoms if atom.atomicnum == 1),
                    "ph": ph
                }
            else:
                return {"error": "Adding hydrogens requires OpenBabel Python bindings"}
                
        except Exception as e:
            return {"error": str(e)}
    
    @mcp_tool({
        "smiles": {"type": "string", "required": True}
    })
    async def tool_check_drug_like(self, smiles: str) -> Dict[str, Any]:
        """Check if molecule is drug-like (Lipinski's Rule of Five)"""
        try:
            if self.has_python:
                mol = pybel.readstring("smiles", smiles)
                
                mw = mol.molwt
                logp = mol.calcdesc(['logP'])['logP']
                hbd = mol.calcdesc(['HBD'])['HBD']
                hba = mol.calcdesc(['HBA1'])['HBA1']
                
                lipinski_violations = 0
                if mw > 500:
                    lipinski_violations += 1
                if logp > 5:
                    lipinski_violations += 1
                if hbd > 5:
                    lipinski_violations += 1
                if hba > 10:
                    lipinski_violations += 1
                
                return {
                    "smiles": smiles,
                    "molecular_weight": mw,
                    "logp": logp,
                    "h_bond_donors": hbd,
                    "h_bond_acceptors": hba,
                    "lipinski_violations": lipinski_violations,
                    "is_drug_like": lipinski_violations <= 1
                }
            else:
                return {"error": "Drug-likeness check requires OpenBabel Python bindings"}
                
        except Exception as e:
            return {"error": str(e)}


async def main():
    """Run the OpenBabel MCP server"""
    server = OpenBabelMCPServer()
    
    # Check command line arguments
    if len(sys.argv) > 1 and sys.argv[1] == "--tcp":
        # Run as TCP server for testing
        port = int(sys.argv[2]) if len(sys.argv) > 2 else 8768
        await server.run_tcp(port=port)
    else:
        # Run as stdio server (standard MCP mode)
        await server.run_stdio()


if __name__ == "__main__":
    asyncio.run(main())
