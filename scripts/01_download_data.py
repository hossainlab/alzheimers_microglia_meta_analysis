#!/usr/bin/env python3
"""
Download datasets from GEO for Alzheimer's Microglia Meta-Analysis.
"""

import sys
from pathlib import Path

# Ensure package is importable
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from ad_microglia.data.download import download_geo_dataset
from ad_microglia.config.settings import get_settings


def main():
    """Download all required datasets."""
    settings = get_settings()

    # Dataset list
    datasets = list(settings.DATASETS.keys())

    print(f"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                      DATA DOWNLOAD - AD MICROGLIA META-ANALYSIS              ║
╚══════════════════════════════════════════════════════════════════════════════╝

Downloading {len(datasets)} datasets from GEO...
""")

    success_count = 0

    for geo_id in datasets:
        print(f"\n{'='*60}")
        print(f"Downloading {geo_id}: {settings.DATASETS[geo_id]}")
        print(f"{'='*60}")

        try:
            dataset_dir = download_geo_dataset(geo_id, settings.RAW_DATA_DIR)
            if dataset_dir:
                success_count += 1
                print(f"✓ {geo_id} downloaded successfully")
        except Exception as e:
            print(f"✗ Error downloading {geo_id}: {e}")

    print(f"""
\n{'='*80}
DOWNLOAD SUMMARY
{'='*80}

Successfully downloaded: {success_count}/{len(datasets)} datasets

Next step: Run preprocessing
  python scripts/02_preprocess_data.py
""")


if __name__ == "__main__":
    main()
