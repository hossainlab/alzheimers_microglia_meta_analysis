# Alzheimer's Disease Microglia Meta-Analysis

## Overview

This project conducts a comprehensive meta-analysis of microglial activation states across multiple Alzheimer's disease single-cell RNA-sequencing studies. The analysis aims to identify conserved activation states that are reproducible across different datasets, species, and experimental conditions.

## Research Question

**Are there conserved microglial activation states across Alzheimer's single-cell studies?**

## Key Features

- **Multi-dataset Integration**: Harmonizes data from 5+ major Alzheimer's scRNA-seq studies
- **Cross-species Analysis**: Includes both human and mouse datasets
- **Advanced Integration**: Uses scVI and Scanorama for batch correction
- **Statistical Validation**: Meta-analysis with effect size calculation and conservation testing
- **Comprehensive Annotation**: Automated cell type identification using CellTypist
- **Publication-ready Output**: High-quality figures and detailed reports

## Project Structure

```
alzheimers_microglia_meta_analysis/
├── data/
│   ├── raw/                    # Original downloaded datasets
│   ├── processed/              # Processed and integrated data
│   └── metadata/               # Dataset metadata and catalogs
├── scripts/
│   ├── preprocessing/          # Data download and initial processing
│   ├── integration/            # Dataset integration scripts
│   ├── analysis/              # Activation state and meta-analysis
│   └── visualization/         # Additional plotting scripts
├── results/
│   ├── figures/               # All generated plots
│   ├── tables/                # Statistical results and gene lists
│   ├── reports/               # Analysis summary reports
│   └── meta_analysis/         # Meta-analysis specific results
├── notebooks/                 # Jupyter notebooks for exploration
├── docs/                      # Documentation and protocols
└── run_pipeline.py           # Main pipeline execution script
```

## Datasets Included

### Primary GEO Datasets:
1. **GSE98969** - Keren-Shaul et al. (2017) - Mouse AD model microglia states (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98969)
2. **GSE103334** - Mathys et al. (2019) - Human AD hippocampus (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103334)
3. **GSE135437** - Sankowski et al. (2019) - Human AD cortex (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135437)
4. **GSE157827** - Leng et al. (2021) - Human AD prefrontal cortex (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827)
5. **GSE129788** - Aging mouse brain microglia (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129788)

### Specialized Databases:
- **Single Cell Portal** - Additional AD datasets
- **Allen Brain Map** - Reference brain cell types

## Quick Start

### Requirements

```bash
# Core dependencies
scanpy>=1.9.0
pandas>=1.5.0
numpy>=1.21.0
scipy>=1.9.0
matplotlib>=3.6.0
seaborn>=0.12.0

# Integration and analysis
scvi-tools>=1.0.0
scanpy[scanorama]>=1.9.0
decoupler>=1.4.0
celltypist>=1.6.0

# Statistics
statsmodels>=0.14.0
scikit-learn>=1.1.0

# Data access
GEOparse>=2.0.3
pandas>=1.5.0
```

### Installation

```bash
# Clone or extract the project
cd alzheimers_microglia_meta_analysis

# Install dependencies
pip install -r requirements.txt

# Run complete pipeline
python run_pipeline.py
```

### Alternative: Step-by-step Execution

```bash
# 1. Download data
python scripts/preprocessing/01_download_data.py

# 2. Initial processing
python scripts/preprocessing/02_initial_processing.py

# 3. Cell annotation  
python scripts/preprocessing/03_cell_annotation.py

# 4. Dataset integration
python scripts/integration/01_integrate_datasets.py

# 5. Activation state analysis
python scripts/analysis/02_activation_states.py

# 6. Meta-analysis validation
python scripts/analysis/03_meta_analysis.py
```

## Analysis Pipeline

### 1. Data Acquisition & Processing
- Downloads datasets from GEO using GEOparse
- Quality control filtering (cells: >200 genes, genes: >3 cells)
- Normalization to 10,000 counts per cell + log1p transformation
- Mitochondrial and ribosomal gene annotation

### 2. Cell Type Annotation
- Uses CellTypist with Immune_All_High.pkl model
- Filters for high-confidence microglial cells
- Validates annotations using known microglial markers (CX3CR1, P2RY12, TMEM119)

### 3. Dataset Integration  
- **Method 1**: scVI - Deep generative model for integration
- **Method 2**: Scanorama - Mutual nearest neighbors approach
- Evaluates integration quality using silhouette scores
- Preserves biological signal while removing batch effects

### 4. Activation State Discovery
- Multi-resolution Leiden clustering (0.1-1.0)
- Optimal resolution selection via silhouette score
- Differential expression analysis between states
- Pathway enrichment using MSigDB via decoupler

### 5. Conservation Analysis
- Identifies gene signatures conserved across datasets
- Calculates conservation scores for each activation state
- Cross-species validation between human and mouse data

### 6. Meta-Analysis
- Fixed-effects meta-analysis of activation state proportions  
- Effect size calculation with Wilson score confidence intervals
- Heterogeneity assessment using I² statistic
- Statistical significance testing for conservation

## Key Outputs

### Statistical Results
- `results/meta_analysis/meta_analysis_results.csv` - Pooled effect sizes
- `results/tables/gene_signatures/conserved_signatures.csv` - Conserved gene lists
- `results/tables/gene_signatures/activation_state_markers.csv` - DE genes

### Visualizations
- `results/figures/activation_states/comprehensive_analysis.pdf` - Main analysis
- `results/meta_analysis/forest_plots.pdf` - Meta-analysis forest plots
- `results/meta_analysis/conservation_significance.pdf` - Conservation testing

### Reports
- `results/reports/activation_states_summary.md` - Analysis summary
- `results/meta_analysis/meta_analysis_report.md` - Meta-analysis results
- `docs/analysis_protocol.md` - Detailed methodology

## Expected Results

Based on the analysis design, you should expect to find:

1. **3-5 conserved microglial activation states** across datasets
2. **Gene signatures** for each state (50-200 genes per state)
3. **Cross-species conservation** for major activation programs
4. **Statistical validation** of conservation (p<0.05 for conserved states)
5. **Pathway enrichment** showing functional relevance (neuroinflammation, phagocytosis, etc.)

## Methodology References

### Integration Methods
- **scVI**: Lopez et al., Nature Methods 2018 [1]
- **Scanorama**: Hie et al., Nature Biotechnology 2019 [2]

### Cell Annotation
- **CellTypist**: Domínguez Conde et al., Science 2022 [3]

### Meta-analysis
- **Fixed-effects model**: Borenstein et al., Introduction to Meta-Analysis 2009 [4]
- **Effect size calculation**: Wilson, JASA 1927 [5]

## Troubleshooting

### Common Issues

1. **Download failures**: Check internet connection and GEO access
2. **Memory errors**: Reduce dataset size or increase system memory
3. **Integration failures**: Try alternative integration method
4. **Annotation errors**: Verify CellTypist model installation

### Performance Optimization

- Use GPU for scVI integration if available
- Enable parallel processing for differential expression
- Consider subsampling for initial testing

## Contributing

This is a research analysis pipeline. For questions or improvements:

1. Check existing documentation in `docs/`
2. Review analysis protocol for methodology details  
3. Examine log outputs for debugging information

## License

This project is for research purposes. Dataset usage subject to original publication terms.

## Citation

If you use this analysis pipeline, please cite the relevant methods papers and original datasets as listed in the results reports.

---

## References

[1] López, R., Regier, J., Cole, M.B., Jordan, M.I. & Yosef, N. Deep generative modeling for single-cell transcriptomics. Nat. Methods 15, 1053–1058 (2018).

[2] Hie, B., Bryson, B. & Berger, B. Efficient integration of heterogeneous single-cell transcriptomes using Scanorama. Nat. Biotechnol. 37, 685–691 (2019).

[3] Domínguez Conde, C. et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science 376, eabl5197 (2022).

[4] Borenstein, M., Hedges, L.V., Higgins, J.P.T. & Rothstein, H.R. Introduction to Meta-Analysis (John Wiley & Sons, 2009).

[5] Wilson, E.B. Probable inference, the law of succession, and statistical inference. J. Am. Stat. Assoc. 22, 209–212 (1927).
