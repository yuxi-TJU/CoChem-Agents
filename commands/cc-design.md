---
description: Design molecules based on objectives - 基于目标设计分子
tools: [codebase_search]
parameters:
  objective:
    description: Design objective (e.g., 'improve solubility', 'reduce toxicity')
    required: true
  scaffold:
    description: Starting scaffold or reference molecule
    required: false
---

# Molecular Design - 分子设计

Design new molecules based on the specified objective: ${objective}

## Design Strategy

### 1. Objective Analysis
- Understand the specific goals
- Identify key molecular properties to optimize
- Define success criteria and constraints

### 2. Starting Point
${scaffold ? `Starting from scaffold: ${scaffold}` : 'De novo design approach'}

### 3. Design Approaches

#### Structure-Based Design
- Analyze target binding site (if applicable)
- Identify key interactions
- Design molecules to enhance binding
- Consider selectivity requirements

#### Property-Based Design
- Current property profile analysis
- Property optimization strategies:
  - Solubility: Add polar groups, reduce LogP
  - Permeability: Balance polarity and lipophilicity
  - Metabolic stability: Block metabolic hot spots
  - Toxicity: Remove structural alerts

#### Fragment-Based Design
- Identify key fragments
- Fragment growing/linking/merging
- Maintain ligand efficiency

### 4. Chemical Modifications

Suggest specific modifications:
- Bioisosteric replacements
- Scaffold hopping options
- Side chain variations
- Stereochemical considerations

### 5. ADMET Optimization

For each design:
- Predict ADMET properties
- Compare with target profile
- Identify potential issues
- Suggest solutions

### 6. Synthetic Feasibility

Evaluate each design for:
- Synthetic accessibility score
- Number of synthetic steps
- Commercial starting materials
- Key transformations required

## Output

Provide 3-5 designed molecules with:
1. **Structure** (SMILES and 2D depiction)
2. **Rationale** for design choices
3. **Predicted properties** vs objectives
4. **Synthetic route** (brief outline)
5. **Prioritization** based on overall profile

Rank designs by likelihood of success.
