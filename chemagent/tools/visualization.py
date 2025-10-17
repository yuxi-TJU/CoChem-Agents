"""Molecular visualization tools"""

from typing import Dict, List, Any, Optional
import logging
import base64
from io import BytesIO
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib.pyplot as plt
import numpy as np

from ..core.registry import BaseTool, ToolMetadata, ToolCategory


logger = logging.getLogger(__name__)


class MoleculeVisualizer(BaseTool):
    """Tool for visualizing molecular structures"""
    
    def get_metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="molecule_visualizer",
            description="Visualize molecular structures in various formats",
            category=ToolCategory.VISUALIZATION,
            version="1.0.0",
            author="ChemAgent",
            dependencies=["rdkit", "matplotlib", "pillow"],
            parameters={
                "required": ["smiles"],
                "properties": {
                    "smiles": {
                        "type": "string",
                        "description": "SMILES string or list of SMILES"
                    },
                    "format": {
                        "type": "string",
                        "enum": ["png", "svg", "ascii", "3d"],
                        "description": "Output format",
                        "default": "png"
                    },
                    "size": {
                        "type": "array",
                        "items": {"type": "integer"},
                        "description": "Image size [width, height]",
                        "default": [300, 300]
                    },
                    "highlight_atoms": {
                        "type": "array",
                        "items": {"type": "integer"},
                        "description": "Atom indices to highlight"
                    },
                    "highlight_bonds": {
                        "type": "array",
                        "items": {"type": "integer"},
                        "description": "Bond indices to highlight"
                    }
                }
            },
            examples=[
                {"smiles": "CC(=O)OC1=CC=CC=C1C(=O)O", "format": "png"},
                {"smiles": "CCO", "format": "ascii"},
                {"smiles": "c1ccccc1", "format": "svg", "size": [400, 400]}
            ]
        )
        
    async def execute(self,
                      smiles: str,
                      format: str = "png",
                      size: List[int] = None,
                      highlight_atoms: List[int] = None,
                      highlight_bonds: List[int] = None,
                      **kwargs) -> Dict[str, Any]:
        """Visualize molecular structure"""
        try:
            if size is None:
                size = [300, 300]
                
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"error": f"Invalid SMILES: {smiles}"}
                
            # Generate 2D coordinates
            AllChem.Compute2DCoords(mol)
            
            if format == "png":
                return self._generate_png(mol, size, highlight_atoms, highlight_bonds)
            elif format == "svg":
                return self._generate_svg(mol, size, highlight_atoms, highlight_bonds)
            elif format == "ascii":
                return self._generate_ascii(mol)
            elif format == "3d":
                return self._generate_3d_coords(mol)
            else:
                return {"error": f"Unsupported format: {format}"}
                
        except Exception as e:
            logger.error(f"Visualization failed: {e}")
            return {"error": str(e)}
            
    def _generate_png(self,
                     mol: Chem.Mol,
                     size: List[int],
                     highlight_atoms: Optional[List[int]] = None,
                     highlight_bonds: Optional[List[int]] = None) -> Dict[str, Any]:
        """Generate PNG image"""
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        
        # Set highlighting
        if highlight_atoms or highlight_bonds:
            drawer.DrawMolecule(
                mol,
                highlightAtoms=highlight_atoms or [],
                highlightBonds=highlight_bonds or []
            )
        else:
            drawer.DrawMolecule(mol)
            
        drawer.FinishDrawing()
        
        # Convert to base64
        img_data = drawer.GetDrawingText()
        img_base64 = base64.b64encode(img_data).decode('utf-8')
        
        return {
            "format": "png",
            "size": size,
            "image_base64": img_base64,
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds()
        }
        
    def _generate_svg(self,
                     mol: Chem.Mol,
                     size: List[int],
                     highlight_atoms: Optional[List[int]] = None,
                     highlight_bonds: Optional[List[int]] = None) -> Dict[str, Any]:
        """Generate SVG image"""
        drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
        
        # Set highlighting
        if highlight_atoms or highlight_bonds:
            drawer.DrawMolecule(
                mol,
                highlightAtoms=highlight_atoms or [],
                highlightBonds=highlight_bonds or []
            )
        else:
            drawer.DrawMolecule(mol)
            
        drawer.FinishDrawing()
        
        svg_string = drawer.GetDrawingText()
        
        return {
            "format": "svg",
            "size": size,
            "svg": svg_string,
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds()
        }
        
    def _generate_ascii(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Generate ASCII representation"""
        # Simple ASCII representation
        ascii_art = []
        ascii_art.append(f"Molecule: {Chem.MolToSmiles(mol)}")
        ascii_art.append(f"Atoms: {mol.GetNumAtoms()}")
        ascii_art.append(f"Bonds: {mol.GetNumBonds()}")
        ascii_art.append("\nAtom List:")
        
        for atom in mol.GetAtoms():
            ascii_art.append(f"  {atom.GetIdx()}: {atom.GetSymbol()} (charge: {atom.GetFormalCharge()})")
            
        ascii_art.append("\nBond List:")
        for bond in mol.GetBonds():
            ascii_art.append(
                f"  {bond.GetBeginAtomIdx()}-{bond.GetEndAtomIdx()}: {bond.GetBondType()}"
            )
            
        return {
            "format": "ascii",
            "ascii": "\n".join(ascii_art),
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds()
        }
        
    def _generate_3d_coords(self, mol: Chem.Mol) -> Dict[str, Any]:
        """Generate 3D coordinates"""
        # Add hydrogens and generate 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Extract coordinates
        conf = mol.GetConformer()
        coords = []
        atoms = []
        
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append([pos.x, pos.y, pos.z])
            atoms.append(mol.GetAtomWithIdx(i).GetSymbol())
            
        # Extract bonds
        bonds = []
        for bond in mol.GetBonds():
            bonds.append([
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                str(bond.GetBondType())
            ])
            
        return {
            "format": "3d",
            "atoms": atoms,
            "coordinates": coords,
            "bonds": bonds,
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds(),
            "energy": self._calculate_energy(mol)
        }
        
    def _calculate_energy(self, mol: Chem.Mol) -> float:
        """Calculate molecular energy using MMFF"""
        try:
            ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
            if ff:
                return ff.CalcEnergy()
        except:
            pass
        return 0.0


class ReactionVisualizer(BaseTool):
    """Tool for visualizing chemical reactions"""
    
    def get_metadata(self) -> ToolMetadata:
        return ToolMetadata(
            name="reaction_visualizer",
            description="Visualize chemical reactions",
            category=ToolCategory.VISUALIZATION,
            version="1.0.0",
            author="ChemAgent",
            dependencies=["rdkit", "matplotlib"],
            parameters={
                "required": ["reaction_smarts"],
                "properties": {
                    "reaction_smarts": {
                        "type": "string",
                        "description": "Reaction SMARTS string"
                    },
                    "reactants": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "Reactant SMILES (optional)"
                    },
                    "size": {
                        "type": "array",
                        "items": {"type": "integer"},
                        "description": "Image size [width, height]",
                        "default": [600, 300]
                    }
                }
            },
            examples=[
                {"reaction_smarts": "[C:1]=[O:2]>>[C:1]-[O:2]"},
                {"reaction_smarts": "[C:1](=[O:2])-[OH].[N:3]>>[C:1](=[O:2])-[N:3]", "reactants": ["CC(=O)O", "CN"]}
            ]
        )
        
    async def execute(self,
                      reaction_smarts: str,
                      reactants: Optional[List[str]] = None,
                      size: List[int] = None,
                      **kwargs) -> Dict[str, Any]:
        """Visualize chemical reaction"""
        try:
            if size is None:
                size = [600, 300]
                
            # Parse reaction
            rxn = AllChem.ReactionFromSmarts(reaction_smarts)
            if rxn is None:
                return {"error": f"Invalid reaction SMARTS: {reaction_smarts}"}
                
            # Apply to reactants if provided
            products = []
            if reactants:
                reactant_mols = []
                for smiles in reactants:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        reactant_mols.append(mol)
                        
                if len(reactant_mols) == rxn.GetNumReactantTemplates():
                    products_tuples = rxn.RunReactants(tuple(reactant_mols))
                    if products_tuples:
                        products = [Chem.MolToSmiles(p) for p in products_tuples[0]]
                        
            # Generate visualization
            drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
            drawer.DrawReaction(rxn)
            drawer.FinishDrawing()
            
            img_data = drawer.GetDrawingText()
            img_base64 = base64.b64encode(img_data).decode('utf-8')
            
            return {
                "reaction_smarts": reaction_smarts,
                "num_reactants": rxn.GetNumReactantTemplates(),
                "num_products": rxn.GetNumProductTemplates(),
                "image_base64": img_base64,
                "size": size,
                "predicted_products": products if products else []
            }
            
        except Exception as e:
            logger.error(f"Reaction visualization failed: {e}")
            return {"error": str(e)}
