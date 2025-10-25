#!/usr/bin/env python3
"""
Preprocess datasets with quality control and normalization.
"""

import sys
from pathlib import Path
import json

# Ensure package is importable
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from ad_microglia.data.preprocessing import DataPreprocessor
from ad_microglia.config.settings import get_settings


def main():
    """Preprocess all datasets."""
    settings = get_settings()
    preprocessor = DataPreprocessor()

    # Dataset list
    datasets = [
        (name, settings.get_raw_data_path(name))
        for name in settings.DATASETS.keys()
    ]

    print(f"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                 DATA PREPROCESSING - AD MICROGLIA META-ANALYSIS              ║
╚══════════════════════════════════════════════════════════════════════════════╝

Processing {len(datasets)} datasets...
""")

    success_count = 0

    for dataset_name, dataset_path in datasets:
        if preprocessor.process_dataset(dataset_path, dataset_name):
            success_count += 1

    print(f"""
\n{'='*80}
PREPROCESSING COMPLETE
{'='*80}

Successfully processed: {success_count}/{len(datasets)} datasets

Results saved to: {settings.PROCESSED_DATA_DIR}
Plots saved to: plots/quality_control/

Next step: Run cell annotation
  python scripts/03_annotate_cells.py
""")

    # Generate summary report
    summary = {
        'total_datasets': len(datasets),
        'successful': success_count,
        'failed': len(datasets) - success_count,
        'datasets': [name for name, _ in datasets]
    }

    settings.METADATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(settings.METADATA_DIR / 'preprocessing_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)


if __name__ == "__main__":
    main()
