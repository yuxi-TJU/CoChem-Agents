"""
Chemical tools - Reference implementations and fallbacks
These are used when official MCP servers are not available
In production, prefer official MCP servers over these implementations
"""

# Import only if needed as fallback
# These serve as reference implementations
from .rdkit_tools import RDKitTool, MoleculeParser, ReactionPredictor
from .pubchem import PubChemTool
from .visualization import MoleculeVisualizer
from .literature import LiteratureSearchTool

# New enhanced tools (ChemCrow-like features)
from .patent_search import PatentSearchTool, CASLookupTool
from .literature_enhanced import EnhancedLiteratureSearch
from .safety_assessment import SafetyAssessmentTool
from .price_lookup import PriceLookupTool, ChemicalEconomicsAnalyzer

__all__ = [
    # Core chemistry tools
    "RDKitTool",
    "MoleculeParser", 
    "ReactionPredictor",
    "PubChemTool",
    "MoleculeVisualizer",
    
    # Enhanced tools (ChemCrow-inspired)
    "PatentSearchTool",      # Patent checking
    "CASLookupTool",         # CAS number conversion
    "EnhancedLiteratureSearch",  # Advanced literature search
    "SafetyAssessmentTool",  # Comprehensive safety evaluation
    "PriceLookupTool",       # Price and availability
    "ChemicalEconomicsAnalyzer",  # Make vs buy analysis
    
    # Legacy
    "LiteratureSearchTool",  # Basic version
]
