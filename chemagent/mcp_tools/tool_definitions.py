"""
MCP Tool Definitions for Chemistry
定义通过 MCP 协议暴露给 Claude Code 和 Gemini CLI 的工具
"""

def get_chemistry_tool_definitions():
    """
    返回所有化学工具的 MCP 定义
    这些定义会被 Claude Code 和 Gemini CLI 识别和调用
    """
    return [
        {
            "name": "analyze_molecule",
            "description": "Analyze molecular structure and calculate properties",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "structure": {
                        "type": "string",
                        "description": "SMILES, InChI, or molecule name"
                    },
                    "calculations": {
                        "type": "array",
                        "items": {
                            "type": "string",
                            "enum": ["properties", "descriptors", "druglike", "toxicity", "3d"]
                        },
                        "description": "Types of calculations to perform"
                    }
                },
                "required": ["structure"]
            }
        },
        {
            "name": "predict_synthesis",
            "description": "Predict synthesis routes for a target molecule",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "target": {
                        "type": "string",
                        "description": "Target molecule SMILES or name"
                    },
                    "max_steps": {
                        "type": "integer",
                        "description": "Maximum number of synthetic steps",
                        "default": 5
                    },
                    "constraints": {
                        "type": "object",
                        "properties": {
                            "available_reagents": {"type": "array", "items": {"type": "string"}},
                            "avoid_reagents": {"type": "array", "items": {"type": "string"}},
                            "green_chemistry": {"type": "boolean", "default": false}
                        }
                    }
                },
                "required": ["target"]
            }
        },
        {
            "name": "predict_reaction",
            "description": "Predict products of a chemical reaction",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "reactants": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of reactant SMILES"
                    },
                    "conditions": {
                        "type": "object",
                        "properties": {
                            "temperature": {"type": "string"},
                            "solvent": {"type": "string"},
                            "catalyst": {"type": "string"},
                            "time": {"type": "string"}
                        }
                    },
                    "mechanism": {
                        "type": "boolean",
                        "description": "Include reaction mechanism",
                        "default": false
                    }
                },
                "required": ["reactants"]
            }
        },
        {
            "name": "search_chemical_database",
            "description": "Search chemical databases (PubChem, ChEMBL, etc.)",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search query"
                    },
                    "database": {
                        "type": "string",
                        "enum": ["pubchem", "chembl", "chemspider", "all"],
                        "default": "pubchem"
                    },
                    "search_type": {
                        "type": "string",
                        "enum": ["name", "smiles", "inchi", "similarity", "substructure"],
                        "default": "name"
                    },
                    "limit": {
                        "type": "integer",
                        "default": 10
                    }
                },
                "required": ["query"]
            }
        },
        {
            "name": "optimize_molecule",
            "description": "Optimize molecular properties for specific targets",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "molecule": {
                        "type": "string",
                        "description": "Starting molecule SMILES"
                    },
                    "optimization_targets": {
                        "type": "object",
                        "properties": {
                            "logP": {"type": "number"},
                            "molecular_weight": {"type": "number"},
                            "tpsa": {"type": "number"},
                            "hbd": {"type": "integer"},
                            "hba": {"type": "integer"},
                            "qed": {"type": "number"}
                        }
                    },
                    "constraints": {
                        "type": "object",
                        "properties": {
                            "maintain_scaffold": {"type": "boolean", "default": true},
                            "max_heavy_atoms": {"type": "integer"},
                            "allowed_elements": {"type": "array", "items": {"type": "string"}}
                        }
                    }
                },
                "required": ["molecule", "optimization_targets"]
            }
        },
        {
            "name": "visualize_molecule",
            "description": "Generate molecular visualizations",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "molecule": {
                        "type": "string",
                        "description": "Molecule SMILES or InChI"
                    },
                    "visualization_type": {
                        "type": "string",
                        "enum": ["2d", "3d", "surface", "orbital"],
                        "default": "2d"
                    },
                    "highlight": {
                        "type": "object",
                        "properties": {
                            "atoms": {"type": "array", "items": {"type": "integer"}},
                            "bonds": {"type": "array", "items": {"type": "integer"}},
                            "substructure": {"type": "string"}
                        }
                    },
                    "format": {
                        "type": "string",
                        "enum": ["png", "svg", "mol", "pdb"],
                        "default": "png"
                    }
                },
                "required": ["molecule"]
            }
        },
        {
            "name": "batch_process_molecules",
            "description": "Process multiple molecules in batch",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "molecules": {
                        "type": "array",
                        "items": {"type": "string"},
                        "description": "List of molecules (SMILES)"
                    },
                    "operation": {
                        "type": "string",
                        "enum": ["analyze", "filter", "cluster", "scaffold_analysis"],
                        "description": "Operation to perform"
                    },
                    "filters": {
                        "type": "object",
                        "properties": {
                            "lipinski": {"type": "boolean"},
                            "mw_range": {"type": "array", "items": {"type": "number"}},
                            "logp_range": {"type": "array", "items": {"type": "number"}},
                            "similarity_threshold": {"type": "number"}
                        }
                    }
                },
                "required": ["molecules", "operation"]
            }
        },
        {
            "name": "calculate_descriptors",
            "description": "Calculate molecular descriptors for QSAR/QSPR",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "molecule": {
                        "type": "string",
                        "description": "Molecule SMILES"
                    },
                    "descriptor_sets": {
                        "type": "array",
                        "items": {
                            "type": "string",
                            "enum": ["rdkit", "mordred", "constitutional", "topological", "electronic", "3d"]
                        },
                        "default": ["rdkit"]
                    }
                },
                "required": ["molecule"]
            }
        },
        {
            "name": "dock_molecule",
            "description": "Perform molecular docking",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "ligand": {
                        "type": "string",
                        "description": "Ligand SMILES"
                    },
                    "protein": {
                        "type": "string",
                        "description": "Protein PDB ID or structure"
                    },
                    "docking_params": {
                        "type": "object",
                        "properties": {
                            "exhaustiveness": {"type": "integer", "default": 8},
                            "num_poses": {"type": "integer", "default": 10},
                            "box_center": {"type": "array", "items": {"type": "number"}},
                            "box_size": {"type": "array", "items": {"type": "number"}}
                        }
                    }
                },
                "required": ["ligand", "protein"]
            }
        },
        {
            "name": "retrosynthetic_analysis",
            "description": "Perform retrosynthetic analysis",
            "inputSchema": {
                "type": "object",
                "properties": {
                    "target": {
                        "type": "string",
                        "description": "Target molecule SMILES"
                    },
                    "depth": {
                        "type": "integer",
                        "description": "Maximum depth of analysis",
                        "default": 3
                    },
                    "strategies": {
                        "type": "array",
                        "items": {
                            "type": "string",
                            "enum": ["functional_group", "ring_analysis", "strategic_bond", "all"]
                        },
                        "default": ["all"]
                    }
                },
                "required": ["target"]
            }
        }
    ]
