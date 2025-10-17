"""RDKit-based chemistry tools"""

from typing import Dict, List, Any, Optional, Union
import logging
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import IPythonConsole
import numpy as np

from ..core.registry import BaseTool, ToolMetadata, ToolCategory


logger = logging.getLogger(__name__)


class RDKitTool(BaseTool):
    """General RDKit tool for molecular operations"""
    
    def get_metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="rdkit_tool",
            description="General RDKit tool for molecular operations",
            category=ToolCategory.MOLECULAR,
            version="1.0.0",
            author="ChemAgent",
            dependencies=["rdkit"],
            parameters={
                "required": ["smiles", "operation"],
                "properties": {
                    "smiles": {"type": "string", "description": "SMILES string of molecule"},
                    "operation": {
                        "type": "string",
                        "enum": ["validate", "canonicalize", "properties", "fingerprint"],
                        "description": "Operation to perform"
                    }
                }
            },
            examples=[
                {"smiles": "CCO", "operation": "properties"},
                {"smiles": "c1ccccc1", "operation": "canonicalize"}
            ]
        )
        
    async def execute(self, smiles: str, operation: str, **kwargs) -> Dict[str, Any]:
        """Execute RDKit operation"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"error": f"Invalid SMILES: {smiles}"}
                
            if operation == "validate":
                return {"valid": True, "smiles": smiles}
                
            elif operation == "canonicalize":
                canonical = Chem.MolToSmiles(mol)
                return {"canonical_smiles": canonical}
                
            elif operation == "properties":
                return self._calculate_properties(mol)
                
            elif operation == "fingerprint":
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                return {"fingerprint": fp.ToBitString()}
                
            else:
                return {"error": f"Unknown operation: {operation}"}
                
        except Exception as e:
            logger.error(f"RDKit operation failed: {e}")
            return {"error": str(e)}
            
    def _calculate_properties(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Calculate molecular properties"""
        return {
            "molecular_weight": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "tpsa": Descriptors.TPSA(mol),
            "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "aromatic_rings": Descriptors.NumAromaticRings(mol),
            "heavy_atoms": Descriptors.HeavyAtomCount(mol),
            "formal_charge": Chem.rdmolops.GetFormalCharge(mol),
            "num_rings": rdMolDescriptors.CalcNumRings(mol),
            "lipinski_violations": self._check_lipinski(mol)
        }
        
    def _check_lipinski(self, mol: Chem.Mol) -> int:
        """Check Lipinski's Rule of Five violations"""
        violations = 0
        if Descriptors.MolWt(mol) > 500:
            violations += 1
        if Descriptors.MolLogP(mol) > 5:
            violations += 1
        if Descriptors.NumHDonors(mol) > 5:
            violations += 1
        if Descriptors.NumHAcceptors(mol) > 10:
            violations += 1
        return violations


class MoleculeParser(BaseTool):
    """Tool for parsing different molecular formats"""
    
    def get_metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="molecule_parser",
            description="Parse and convert between different molecular formats",
            category=ToolCategory.MOLECULAR,
            version="1.0.0",
            author="ChemAgent",
            dependencies=["rdkit"],
            parameters={
                "required": ["input", "input_format"],
                "properties": {
                    "input": {"type": "string", "description": "Molecular input string"},
                    "input_format": {
                        "type": "string",
                        "enum": ["smiles", "inchi", "mol", "pdb"],
                        "description": "Input format"
                    },
                    "output_format": {
                        "type": "string",
                        "enum": ["smiles", "inchi", "mol", "pdb"],
                        "description": "Output format (default: smiles)"
                    }
                }
            },
            examples=[
                {"input": "CCO", "input_format": "smiles", "output_format": "inchi"},
                {"input": "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3", "input_format": "inchi", "output_format": "smiles"}
            ]
        )
        
    async def execute(self, 
                      input: str, 
                      input_format: str,
                      output_format: str = "smiles",
                      **kwargs) -> Dict[str, Any]:
        """Parse and convert molecular formats"""
        try:
            # Parse input
            mol = None
            if input_format == "smiles":
                mol = Chem.MolFromSmiles(input)
            elif input_format == "inchi":
                mol = Chem.MolFromInchi(input)
            elif input_format == "mol":
                mol = Chem.MolFromMolBlock(input)
            elif input_format == "pdb":
                mol = Chem.MolFromPDBBlock(input)
            else:
                return {"error": f"Unsupported input format: {input_format}"}
                
            if mol is None:
                return {"error": f"Failed to parse {input_format} input"}
                
            # Convert to output format
            output = None
            if output_format == "smiles":
                output = Chem.MolToSmiles(mol)
            elif output_format == "inchi":
                output = Chem.MolToInchi(mol)
            elif output_format == "mol":
                output = Chem.MolToMolBlock(mol)
            elif output_format == "pdb":
                output = Chem.MolToPDBBlock(mol)
            else:
                return {"error": f"Unsupported output format: {output_format}"}
                
            return {
                "output": output,
                "format": output_format,
                "num_atoms": mol.GetNumAtoms(),
                "num_bonds": mol.GetNumBonds()
            }
            
        except Exception as e:
            logger.error(f"Molecule parsing failed: {e}")
            return {"error": str(e)}


class ReactionPredictor(BaseTool):
    """Tool for predicting chemical reactions"""
    
    def get_metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="reaction_predictor",
            description="Predict products of chemical reactions",
            category=ToolCategory.REACTION,
            version="1.0.0",
            author="ChemAgent",
            dependencies=["rdkit"],
            parameters={
                "required": ["reactants", "reaction_type"],
                "properties": {
                    "reactants": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of reactant SMILES"
                    },
                    "reaction_type": {
                        "type": "string",
                        "enum": ["addition", "elimination", "substitution", "oxidation", "reduction"],
                        "description": "Type of reaction"
                    },
                    "conditions": {
                        "type": "object",
                        "description": "Reaction conditions (temperature, solvent, catalyst)"
                    }
                }
            },
            examples=[
                {
                    "reactants": ["CC=O", "CCO"],
                    "reaction_type": "addition",
                    "conditions": {"catalyst": "H+"}
                }
            ]
        )
        
    async def execute(self,
                      reactants: List[str],
                      reaction_type: str,
                      conditions: Optional[Dict[str, Any]] = None,
                      **kwargs) -> Dict[str, Any]:
        """Predict reaction products"""
        try:
            # Parse reactants
            mols = []
            for smiles in reactants:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return {"error": f"Invalid reactant SMILES: {smiles}"}
                mols.append(mol)
                
            # This is a simplified example - in production, you'd use
            # more sophisticated reaction prediction methods
            products = self._predict_reaction(mols, reaction_type, conditions)
            
            return {
                "reactants": reactants,
                "reaction_type": reaction_type,
                "products": products,
                "conditions": conditions or {},
                "mechanism": self._get_mechanism(reaction_type)
            }
            
        except Exception as e:
            logger.error(f"Reaction prediction failed: {e}")
            return {"error": str(e)}
            
    def _predict_reaction(self,
                         mols: List[Chem.Mol],
                         reaction_type: str,
                         conditions: Optional[Dict[str, Any]]) -> List[str]:
        """Predict reaction products (simplified)"""
        # This is a placeholder - implement actual reaction prediction
        # You could integrate with reaction prediction models or databases
        
        if reaction_type == "addition" and len(mols) == 2:
            # Simple addition example
            combined = Chem.CombineMols(mols[0], mols[1])
            return [Chem.MolToSmiles(combined)]
            
        # Return original molecules as fallback
        return [Chem.MolToSmiles(mol) for mol in mols]
        
    def _get_mechanism(self, reaction_type: str) -> str:
        """Get reaction mechanism description"""
        mechanisms = {
            "addition": "Nucleophilic or electrophilic addition to unsaturated bond",
            "elimination": "Removal of atoms/groups to form double/triple bond",
            "substitution": "Replacement of one functional group with another",
            "oxidation": "Loss of electrons or increase in oxidation state",
            "reduction": "Gain of electrons or decrease in oxidation state"
        }
        return mechanisms.get(reaction_type, "Unknown mechanism")
