# Migration Guide: Old Structure → New Structure

This guide helps you understand the changes made during the restructuring.

## Summary of Changes

The project has been reorganized from a flat script-based structure to a proper Python package with best practices.

## Directory Mapping

### Old → New

```
OLD STRUCTURE                        NEW STRUCTURE
════════════════════════════════════════════════════════════════

run_pipeline.py                  →   scripts/run_full_pipeline.py
test_pipeline.py                 →   tests/test_pipeline.py
quick_test.py                    →   tests/test_quick.py

scripts/preprocessing/
├── 01_download_data.py          →   src/ad_microglia/data/download.py
├── 02_preprocess_data.py        →   src/ad_microglia/data/preprocessing.py
└── 03_cell_annotation.py        →   src/ad_microglia/analysis/cell_annotation.py

scripts/integration/
└── 01_integrate_datasets.py     →   src/ad_microglia/integration/integration.py

scripts/analysis/
├── 02_activation_states.py      →   src/ad_microglia/analysis/activation_states.py
└── 03_meta_analysis.py          →   src/ad_microglia/analysis/meta_analysis.py

scripts/utils/
└── plot_manager.py               →   src/ad_microglia/utils/plot_manager.py

(NEW) Executable scripts         →   scripts/
                                     ├── 01_download_data.py
                                     ├── 02_preprocess_data.py
                                     ├── 03_annotate_cells.py
                                     ├── 04_integrate_datasets.py
                                     ├── 05_analyze_activation_states.py
                                     ├── 06_run_meta_analysis.py
                                     └── run_full_pipeline.py
```

## Import Changes

### Before (Bad)
```python
# Old way with sys.path hacks
import sys
sys.path.append(str(Path(__file__).parent.parent))
from utils.plot_manager import PlotManager
```

### After (Good)
```python
# New way with proper imports
from ad_microglia.utils import PlotManager
from ad_microglia.config import get_settings
```

## Running the Pipeline

### Before
```bash
python run_pipeline.py
```

### After
```bash
# Option 1: Run individual scripts
python scripts/01_download_data.py
python scripts/02_preprocess_data.py
# ... etc

# Option 2: Run full pipeline
python scripts/run_full_pipeline.py

# Option 3: After installing package
pip install -e .
# Then you can import anywhere
```

## Configuration

### Before
Hardcoded values in scripts:
```python
raw_data_dir = 'data/raw'
processed_data_dir = 'data/processed'
min_genes = 200
```

### After
Centralized configuration:
```python
from ad_microglia.config import get_settings

settings = get_settings()
raw_data_dir = settings.RAW_DATA_DIR
min_genes = settings.QC_MIN_GENES
```

Or use environment variables in `.env`:
```bash
QC_MIN_GENES=200
QC_MAX_MT_PERCENT=20
```

Or edit `config/analysis_params.yaml`:
```yaml
quality_control:
  min_genes: 200
  max_mt_percent: 20
```

## New Features

### 1. Package Installation
```bash
# Install as editable package
pip install -e .

# Now you can import from anywhere
python -c "from ad_microglia import PlotManager"
```

### 2. Centralized Settings
```python
from ad_microglia.config import get_settings

settings = get_settings()
print(settings.RAW_DATA_DIR)
print(settings.QC_MIN_GENES)
```

### 3. Environment Variables
```bash
# Copy template
cp .env.example .env

# Edit .env
QC_MIN_GENES=300
INTEGRATION_METHOD=scanorama
```

### 4. Proper Testing
```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_pipeline.py

# With coverage
pytest --cov=src/ad_microglia
```

### 5. Code Quality Tools
```bash
# Format code
black src/ scripts/ tests/

# Check style
flake8 src/

# Type checking
mypy src/
```

## File Locations

### Data Files
- **Raw data**: `data/raw/` (same)
- **Processed data**: `data/processed/` (same)
- **Metadata**: `data/metadata/` (same)

### Results
- **Figures**: `results/figures/` (same)
- **Tables**: `results/tables/` (same)
- **Reports**: `results/reports/` (same)
- **Meta-analysis**: `results/meta_analysis/` (same)

### Plots (New centralized system)
- **All plots**: `plots/` with subdirectories by category

### Configuration
- **YAML configs**: `config/`
- **Environment**: `.env` (create from `.env.example`)

## Breaking Changes

### Import Paths
All imports must now use the `ad_microglia` package name:

❌ **Old (broken)**:
```python
from utils.plot_manager import PlotManager
```

✅ **New (correct)**:
```python
from ad_microglia.utils import PlotManager
```

### Script Locations
Direct script execution moved to `scripts/`:

❌ **Old (broken)**:
```bash
python run_pipeline.py
```

✅ **New (correct)**:
```bash
python scripts/run_full_pipeline.py
```

## Compatibility Notes

### Old Scripts
Old scripts in `scripts/preprocessing/`, `scripts/analysis/`, etc. are preserved for reference but should not be used directly. Use the new executable scripts in `scripts/` instead.

### Data and Results
All existing data and results are compatible. No changes needed.

## Recommended Migration Steps

1. **Install the package**:
   ```bash
   pip install -e .
   ```

2. **Copy environment template**:
   ```bash
   cp .env.example .env
   ```

3. **Update any custom scripts**:
   - Change imports to use `ad_microglia`
   - Use `get_settings()` for configuration

4. **Run tests**:
   ```bash
   pytest
   ```

5. **Run pipeline**:
   ```bash
   python scripts/run_full_pipeline.py
   ```

## Getting Help

If you encounter issues:

1. Check [installation.md](installation.md) for setup help
2. Review [usage.md](usage.md) for examples
3. See [troubleshooting.md](troubleshooting.md) for common issues
4. Check the API documentation in `docs/api/`

## Benefits of New Structure

✅ **Proper Python package** - Installable with `pip install -e .`
✅ **Clean imports** - No more `sys.path` hacks
✅ **Centralized config** - Settings, YAML, environment variables
✅ **Better testing** - Proper test infrastructure
✅ **Code quality** - Black, flake8, mypy support
✅ **Professional** - Follows Python best practices
✅ **Scalable** - Easy to add new features
✅ **Maintainable** - Clear separation of concerns
