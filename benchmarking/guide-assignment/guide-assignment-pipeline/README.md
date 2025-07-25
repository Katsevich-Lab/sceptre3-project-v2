# Benchmarking Pipeline for Guide Assignment Methods

A Nextflow pipeline to benchmark guide assignment methods on CRISPR screening datasets, starting with crispat.

## Quick Start

1. **Configure datasets**: Edit `configs/datasets.csv` to specify which datasets to benchmark
2. **Set dataset path**: Update `dataset_base_dir` in `nextflow.config` 
3. **Choose methods**: Set `params.methods` in `nextflow.config` or command line
4. **Run pipeline**: `nextflow run main.nf`

## Current Status

🚧 **This is scaffolding for initial development** 🚧

### What's Ready
- [x] Basic pipeline structure
- [x] Module template for crispat
- [x] Dataset configuration system
- [x] Output organization

### What Needs Implementation (TODOs)
- [ ] Actual crispat implementation in modules/crispat.nf
- [ ] Conda environment specification
- [ ] Real dataset configurations
- [ ] Additional methods (cleanser, etc.)
- [ ] Result comparison framework

## Structure

```
benchmarking/
├── main.nf              # Main pipeline workflow
├── nextflow.config      # Configuration and parameters
├── workflows/           # Subworkflows
│   └── guide_assignment.nf  # Method routing logic
├── modules/             # Method implementations
│   └── crispat.nf      # crispat guide assignment  
├── configs/             # Configuration
│   └── config.csv      # Methods, datasets, and resources
└── results/            # Output directory (created on run)
```

## Adding New Methods

1. Create new module in `modules/new_method.nf`
2. Add import and branch to `workflows/guide_assignment.nf`
3. Add rows to `configs/config.csv` for new method
4. Configure conda environment in `nextflow.config`

## Development Notes

- Start with toy dataset for testing
- Currently focused on crispat only for simplicity
- Easy to add more methods later by following the same pattern
- Each method outputs standardized CSV format for comparison