---
description: Predict chemical properties and outcomes - 预测化学性质和结果
tools: [read_file]
parameters:
  prediction_type:
    description: Type of prediction (properties, reactions, activity, metabolism, toxicity)
    required: true
  input:
    description: Input molecule or reaction
    required: true
---

# Chemical Prediction - 化学预测

Perform ${prediction_type} prediction for: ${input}

## Prediction Types

### Property Prediction
Calculate and predict molecular properties:
- **Physicochemical**: MW, LogP, LogD, LogS, PSA, HBA/HBD
- **Quantum mechanical**: HOMO/LUMO, dipole, polarizability
- **Thermodynamic**: melting point, boiling point, vapor pressure
- **Spectroscopic**: NMR shifts, UV-Vis λmax, IR frequencies

### Reaction Prediction
Predict reaction outcomes:
- **Products**: Major and minor products
- **Mechanism**: Step-by-step pathway
- **Selectivity**: Regio-, stereo-, chemoselectivity
- **Conditions**: Optimal temperature, solvent, catalyst
- **Yield**: Expected yield range
- **Side reactions**: Potential byproducts

### Activity Prediction
Predict biological activity:
- **Target binding**: Affinity predictions (Ki, IC50)
- **Selectivity profile**: Off-target predictions
- **Mechanism of action**: Binding mode
- **Pharmacophore matching**: Key features
- **QSAR predictions**: Activity from structure

### Metabolism Prediction
Predict metabolic transformations:
- **Phase I metabolism**: CYP-mediated oxidations
- **Phase II metabolism**: Conjugation reactions
- **Metabolic hotspots**: Sites of metabolism
- **Major metabolites**: Structure and properties
- **Species differences**: Human vs rat/mouse
- **Enzyme involvement**: CYP isoforms

### Toxicity Prediction
Predict safety issues:
- **Acute toxicity**: LD50 estimates
- **Organ toxicity**: Hepato-, cardio-, neurotoxicity
- **Genetic toxicity**: Mutagenicity, clastogenicity
- **Carcinogenicity**: Structural alerts
- **Environmental**: Bioaccumulation, persistence

## Prediction Methods

### 1. Knowledge-Based
- Rule-based predictions
- Expert systems
- Structural alerts
- Reaction databases

### 2. Statistical Models
- QSAR/QSPR models
- Machine learning (RF, SVM, XGBoost)
- Deep learning (Graph neural networks)
- Ensemble predictions

### 3. Physics-Based
- Quantum mechanical calculations
- Molecular dynamics
- Docking simulations
- Free energy perturbation

## Output Format

### Predictions
Present results with:
- **Predicted values** with units
- **Confidence scores** or intervals
- **Model applicability** domain check
- **Similar compounds** for comparison

### Visualization
Include relevant plots:
- Property distributions
- Reaction energy profiles
- Metabolic trees
- Toxicity alerts highlighted

### Interpretation
Provide context:
- Compare to typical ranges
- Highlight unusual predictions
- Suggest experimental validation
- Note limitations and assumptions

### Recommendations
Based on predictions:
- Flag potential issues
- Suggest modifications
- Prioritize experiments
- Guide decision making
