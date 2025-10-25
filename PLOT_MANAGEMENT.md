# Plot Management System

## Overview

The plot management system provides centralized, organized storage of all plots generated during the Alzheimer's Microglia Meta-Analysis. This system ensures:

- **No plots appear during analysis** - All plots are automatically saved without showing on screen
- **Organized directory structure** - Plots are categorized by analysis type
- **Multiple formats** - Saves plots in PDF and PNG formats by default
- **Manifest tracking** - Keeps a log of all generated plots with timestamps and descriptions
- **Consistent styling** - Unified plot parameters across all analyses

## Directory Structure

```
plots/
├── plot_manifest.txt           # Log of all generated plots
├── quality_control/           # QC metrics and validation plots
├── preprocessing/             # Data preprocessing visualization
├── integration/              # Dataset integration plots
├── cell_annotation/          # Cell type annotation validation
├── activation_states/        # Microglia activation analysis
├── meta_analysis/           # Meta-analysis results
├── comparative_analysis/    # Cross-dataset comparisons
└── supplementary/          # Additional plots and figures
```

## Usage

### Basic Usage

```python
from utils.plot_manager import PlotManager

# Initialize plot manager
plot_manager = PlotManager()

# Create and save a plot
fig, ax = plt.subplots()
ax.plot(data)
plot_manager.save_plot(fig, 'my_plot', 'qc', 'Description of plot')
```

### Context Manager (Recommended)

```python
from utils.plot_manager import PlotContext

# Automatic plot saving with context manager
with PlotContext('dataset_overview', 'qc', 'Dataset size comparison') as (fig, ax):
    ax.bar(datasets, counts)
    ax.set_title('Dataset Sizes')
    # Plot automatically saved when context exits
```

### Integration in Analysis Scripts

All analysis scripts have been updated to use the plot management system:

1. **Preprocessing** (`scripts/preprocessing/02_preprocess_data.py`)
   - QC plots saved to `plots/quality_control/`
   
2. **Cell Annotation** (`scripts/preprocessing/03_cell_annotation.py`)
   - Annotation validation plots saved to `plots/cell_annotation/`
   
3. **Analysis Scripts** (to be updated)
   - Activation state plots → `plots/activation_states/`
   - Meta-analysis plots → `plots/meta_analysis/`

## Configuration

### Plot Manager Options

```python
plot_manager = PlotManager(
    base_dir='plots',    # Base directory for plots
    dpi=300,            # Plot resolution
    figsize=(8, 6)      # Default figure size
)
```

### Supported File Formats

- PDF (vector graphics, publication quality)
- PNG (raster graphics, web-friendly)
- SVG (optional, vector graphics)

## Features

### Automatic Plot Settings

- **No interactive display**: `plt.ioff()` ensures plots don't show during analysis
- **High DPI**: 300 DPI for publication-quality figures
- **White background**: Consistent styling across all plots
- **Tight layout**: Automatic bbox_inches='tight' for clean borders

### Plot Manifest

Every saved plot is logged in `plots/plot_manifest.txt`:

```
2025-09-04 16:05:49 | qc | GSE98969_qc_metrics | Quality control metrics for GSE98969
    -> plots/quality_control/GSE98969_qc_metrics.pdf
    -> plots/quality_control/GSE98969_qc_metrics.png
```

### Utility Functions

```python
# List plots in a category
qc_plots = plot_manager.list_plots('qc')

# Get path for a specific plot
plot_path = plot_manager.get_plot_path('qc', 'dataset_overview', 'pdf')

# Clean up empty directories
plot_manager.cleanup_empty_dirs()
```

## Benefits

1. **Organization**: Clear categorization of all plots
2. **Reproducibility**: Consistent plot generation across runs
3. **Documentation**: Automatic logging of all generated figures
4. **Flexibility**: Easy to change output formats and directories
5. **Clean Analysis**: No plots interrupting analysis workflow
6. **Publication Ready**: High-quality output suitable for papers

## Migration from Old System

The old system scattered plots across multiple directories:
- `data/processed/*.pdf` (QC plots)
- `results/figures/*.png` (annotation plots)  
- `results/meta_analysis/*.pdf` (meta-analysis plots)

The new system consolidates everything under `plots/` with proper categorization.

## Troubleshooting

### Common Issues

1. **Import Error**: Make sure `scripts` directory is in Python path
2. **Permission Error**: Check write permissions in the plots directory
3. **Missing Plots**: Check the manifest file for save locations

### Matplotlib Configuration

The plot manager automatically configures matplotlib for non-interactive use. If you need to show plots during development, temporarily use:

```python
plt.ion()  # Turn on interactive mode
plt.show()  # Show specific plots
```