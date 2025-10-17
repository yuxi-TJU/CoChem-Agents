---
description: Optimize molecular properties - 优化分子性质
tools: [read_file, codebase_search]
parameters:
  molecule:
    description: Input molecule (SMILES, name, or file)
    required: true
  objectives:
    description: Properties to optimize (e.g., 'solubility,potency,selectivity')
    required: true
  constraints:
    description: Constraints to maintain (e.g., 'MW<500,LogP<5')
    required: false
---

# Molecular Optimization - 分子优化

Optimize the molecule ${molecule} for: ${objectives}
${constraints ? `While maintaining constraints: ${constraints}` : ''}

## Multi-Objective Optimization

### 1. Current Profile Assessment
Analyze the starting molecule:
- Calculate all relevant properties
- Identify strengths and weaknesses
- Set optimization targets

### 2. Optimization Strategies

#### For Solubility Enhancement
- Introduce polar functional groups
- Disrupt crystalline packing
- Add solubilizing groups (PEG, amines)
- Consider prodrug approaches

#### For Potency Improvement
- Enhance key interactions
- Optimize binding geometry
- Increase complementarity
- Consider covalent strategies

#### For Selectivity Enhancement
- Exploit unique binding site features
- Add bulky groups for steric clash
- Optimize electrostatic interactions
- Design out off-target binding

#### For ADMET Optimization
- **Absorption**: Optimize LogP, PSA, MW
- **Distribution**: Manage plasma protein binding
- **Metabolism**: Block/redirect metabolic sites
- **Excretion**: Optimize renal clearance
- **Toxicity**: Remove problematic substructures

### 3. Chemical Modifications

Generate analogs using:
- **Bioisosteric replacement**
- **Homologation/chain extension**
- **Ring opening/closing**
- **Heteroatom walking**
- **Stereochemical modification**

### 4. Machine Learning Predictions

For each optimization:
- Predict property changes
- Calculate optimization scores
- Assess synthetic accessibility
- Evaluate drug-likeness

### 5. Pareto Optimization

When optimizing multiple properties:
- Identify Pareto optimal solutions
- Visualize property trade-offs
- Rank by weighted objectives
- Highlight best compromises

## Deliverables

### Optimized Molecules (5-10 analogs)
For each analog provide:
1. **Structure and SMILES**
2. **Property improvements** (% change)
3. **Radar chart** comparing to original
4. **Synthetic accessibility**
5. **Key modifications** highlighted

### Summary Table
| Analog | ${objectives.split(',').join(' | ')} | Overall Score |
|--------|${objectives.split(',').map(() => '---').join('|')}|-------------|
| Original | baseline values | - |
| Analog 1 | improved values | score |
| ... | ... | ... |

### Recommendations
- Top 3 candidates with rationale
- Synthesis priority order
- Further optimization suggestions
- Potential risks/considerations
