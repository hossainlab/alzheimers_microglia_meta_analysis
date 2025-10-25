# Installation Guide

## Prerequisites

- Python >= 3.8
- pip >= 20.0
- Git

## Installation

### Option 1: Standard Installation (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/alzheimers_microglia_meta_analysis.git
cd alzheimers_microglia_meta_analysis

# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install the package in editable mode
pip install -e .
```

### Option 2: Development Installation

For contributors or developers:

```bash
# Clone the repository
git clone https://github.com/yourusername/alzheimers_microglia_meta_analysis.git
cd alzheimers_microglia_meta_analysis

# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install with development dependencies
pip install -e ".[dev]"
```

## Verify Installation

```bash
# Test imports
python -c "from ad_microglia import PlotManager; print('Installation successful!')"

# Run tests
pytest tests/
```

## Configuration

1. Copy the example environment file:
```bash
cp .env.example .env
```

2. Edit `.env` to customize parameters (optional):
```bash
# Edit with your preferred editor
nano .env
```

## GPU Support (Optional)

For faster scVI integration, install GPU-enabled PyTorch:

```bash
# For CUDA 11.8
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# For CUDA 12.1
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
```

## Troubleshooting

### ImportError: No module named 'ad_microglia'

Make sure you installed the package in editable mode:
```bash
pip install -e .
```

### scvi-tools installation issues

Try installing with conda:
```bash
conda install scvi-tools -c conda-forge
```

### CellTypist model download fails

Manually download the model:
```bash
python -c "import celltypist; celltypist.models.download_models()"
```

## Next Steps

After installation, see:
- [Quick Start Guide](quickstart.md)
- [Usage Examples](usage.md)
- [API Documentation](api/)
