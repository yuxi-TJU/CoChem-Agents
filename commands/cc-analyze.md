---
description: Analyze molecular structure and properties - 分析分子结构和性质
tools: [read_file, codebase_search, grep]
---

# Molecular Analysis - 分子分析

Please perform a comprehensive molecular analysis using the following steps:

## 1. Structure Validation
- Validate the SMILES/MOL file structure
- Check for common structural errors
- Identify functional groups and key features

## 2. Property Calculation
Calculate the following molecular properties:
- Molecular weight and formula
- LogP, LogD, and LogS (solubility)
- Topological polar surface area (TPSA)
- Number of rotatable bonds
- Hydrogen bond donors and acceptors
- Aromatic rings and heteroatoms

## 3. Drug-likeness Assessment
Evaluate drug-like properties:
- Lipinski's Rule of Five
- Veber's Rules
- QED (Quantitative Estimate of Drug-likeness)
- PAINS (Pan-Assay Interference Compounds) alerts
- Synthetic accessibility score

## 4. ADMET Prediction
Predict ADMET properties:
- Absorption: Caco-2 permeability, HIA
- Distribution: BBB penetration, Plasma protein binding
- Metabolism: CYP inhibition/substrate
- Excretion: Clearance predictions
- Toxicity: hERG, hepatotoxicity, mutagenicity

## 5. Visualization
Generate molecular visualizations:
- 2D structure depiction
- 3D conformation (if applicable)
- Electrostatic potential map
- Property heatmaps

## 6. Safety Assessment
Check for safety concerns:
- Structural alerts for toxicity
- Reactive functional groups
- Environmental hazards
- Controlled substance similarity

## 7. Report Generation
Provide a comprehensive report including:
- Summary of key findings
- Detailed property table
- Visual representations
- Recommendations for optimization
- Comparison with known drugs (if applicable)

**Note**: Use official RDKit MCP for calculations when available, fallback to custom implementations if needed.
