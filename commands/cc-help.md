---
description: Get help on ChemAgent commands and features - 获取ChemAgent命令和功能的帮助
tools: []
parameters:
  topic:
    description: Specific command or topic to get help for
    required: false
---

# ChemAgent Help - ChemAgent 帮助

${topic ? `Getting help for: ${topic}` : 'ChemAgent Command Reference'}

## Available Commands

### 🔬 Analysis Commands
- **cc-analyze**: Analyze molecular structure and properties
- **cc-compare**: Compare multiple molecules
- **cc-explain**: Explain chemistry concepts and results

### 🎯 Design Commands  
- **cc-design**: Design molecules based on objectives
- **cc-optimize**: Optimize molecular properties
- **cc-synthesize**: Plan synthesis routes

### 🔮 Prediction Commands
- **cc-predict**: Predict properties and outcomes
- **cc-simulate**: Run molecular simulations

### 🔍 Search Commands
- **cc-search**: Search databases and literature
- **cc-check**: Check safety, patents, and compliance

### ⚙️ Workflow Commands
- **cc-workflow**: Execute chemistry workflows
- **cc-batch**: Process multiple molecules
- **cc-report**: Generate comprehensive reports

### 💡 Utility Commands
- **cc-suggest**: Get suggestions for next steps
- **cc-help**: Show this help information

## Command Usage

All commands follow the pattern:
```
cc-[command] [arguments] [options]
```

### Examples:
- `cc-analyze aspirin.mol`
- `cc-search literature "kinase inhibitors" --limit 20`
- `cc-check safety "CCO" --detailed`
- `cc-workflow drug-discovery target.pdb`

## Common Options

Most commands support:
- `--detailed`: Show comprehensive results
- `--export`: Save results to file
- `--format`: Output format (json, csv, pdf)
- `--limit`: Limit number of results

## Sub-Agents

You can invoke specialized chemistry experts:
- **@chemist**: General chemistry expertise
- **@drug-designer**: Medicinal chemistry and ADMET
- **@synthesist**: Synthesis and process chemistry
- **@safety-expert**: Safety and regulatory compliance
- **@data-analyst**: Chemical data analysis and ML

## Getting Started

1. **Analyze a molecule**: Start with `cc-analyze` to understand properties
2. **Search for information**: Use `cc-search` to find related data
3. **Check feasibility**: Run `cc-check` for safety and IP
4. **Design improvements**: Use `cc-optimize` or `cc-design`
5. **Plan synthesis**: Execute `cc-synthesize` for routes

## Best Practices

### For Drug Discovery:
1. Start with target analysis
2. Search for known ligands
3. Design new molecules
4. Check ADMET and safety
5. Plan synthesis routes

### For Chemical Safety:
1. Always run `cc-check safety` first
2. Review GHS classification
3. Check regulatory status
4. Assess environmental impact
5. Document all findings

### For Literature Research:
1. Use specific search terms
2. Filter by recent years
3. Check both papers and patents
4. Follow citation networks
5. Summarize key findings

## Additional Resources

- **Documentation**: Full guide at `/docs`
- **Examples**: See `/examples` directory  
- **Custom Commands**: Create your own in `.claude/commands/`
- **Configuration**: Adjust settings in `config.yaml`

## Need More Help?

- For specific command help: `cc-help [command-name]`
- Report issues: GitHub repository
- Request features: Open an issue
- Contribute: Pull requests welcome!

**Remember**: ChemAgent enhances your AI assistant's chemistry capabilities. Use it to accelerate research, ensure safety, and make informed decisions.
