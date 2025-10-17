---
description: Execute predefined chemistry workflows - 执行预定义的化学工作流
tools: [read_file, web_search, codebase_search]
parameters:
  workflow_type:
    description: Type of workflow (drug-discovery, materials, synthesis, qsar)
    required: true
  input:
    description: Input data (molecules, targets, objectives)
    required: true
---

# Chemistry Workflow Execution - 化学工作流执行

Execute ${workflow_type} workflow with input: ${input}

## Available Workflows

### Drug Discovery Workflow
Complete drug discovery pipeline:
1. **Target Analysis**: Understand binding site, known ligands
2. **Virtual Screening**: Search chemical databases
3. **Hit Identification**: Filter by drug-likeness, ADMET
4. **Lead Optimization**: Improve potency, selectivity, properties
5. **ADMET Profiling**: Detailed predictions
6. **Synthesis Planning**: Routes for top candidates
7. **Report Generation**: Executive summary with recommendations

### Materials Design Workflow
Materials discovery pipeline:
1. **Property Requirements**: Define target properties
2. **Structure Generation**: Design candidates
3. **Property Prediction**: Calculate material properties
4. **Stability Analysis**: Thermodynamic and kinetic stability
5. **Synthesis Feasibility**: Synthetic routes
6. **Performance Ranking**: Multi-objective optimization
7. **Scale-up Considerations**: Manufacturing analysis

### Synthesis Optimization Workflow
Optimize synthetic routes:
1. **Retrosynthetic Analysis**: Multiple disconnection strategies
2. **Route Comparison**: Evaluate different pathways
3. **Condition Optimization**: Best reagents and conditions
4. **Cost Analysis**: Material and process costs
5. **Green Chemistry**: Environmental impact assessment
6. **Risk Assessment**: Identify potential issues
7. **Experimental Protocol**: Detailed procedures

### QSAR Modeling Workflow
Build predictive models:
1. **Data Collection**: Gather training data
2. **Molecular Descriptors**: Calculate relevant features
3. **Feature Selection**: Identify key descriptors
4. **Model Building**: Multiple algorithms
5. **Validation**: Cross-validation, test set
6. **Applicability Domain**: Define model scope
7. **Predictions**: Apply to new compounds

### High-Throughput Screening Workflow
Process large compound libraries:
1. **Library Preparation**: Clean and standardize structures
2. **Property Filtering**: Remove undesirable compounds
3. **Similarity Clustering**: Group related structures
4. **Batch Predictions**: Parallel processing
5. **Hit Selection**: Multi-criteria ranking
6. **Diversity Analysis**: Ensure chemical diversity
7. **Visualization**: Chemical space plots

## Workflow Execution

### Step-by-Step Process
For each workflow step:
- ✅ Show completion status
- 📊 Display intermediate results
- ⚠️ Flag any issues or warnings
- 💡 Provide insights and observations

### Parallel Processing
When possible:
- Run independent steps in parallel
- Optimize for efficiency
- Report progress in real-time

### Error Handling
- Gracefully handle failures
- Provide alternative approaches
- Continue with remaining steps
- Summarize issues at end

## Output Format

### Executive Summary
- Workflow completed: ✅/❌
- Key findings (3-5 bullet points)
- Recommended next steps
- Critical decisions needed

### Detailed Results
- Step-by-step outcomes
- All generated data
- Visualizations and plots
- Supporting evidence

### Exportable Deliverables
- CSV/Excel files with data
- PDF report with figures
- SDF/MOL files for structures
- JSON for programmatic access

### Follow-up Actions
- Suggested experiments
- Additional analyses
- Validation requirements
- Timeline recommendations
