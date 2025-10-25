#!/usr/bin/env python3
"""
Main pipeline runner for Alzheimer's Microglia Meta-Analysis.
Executes the complete analysis workflow.
"""

import subprocess
import sys
from pathlib import Path
import argparse

# Ensure package is importable
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from ad_microglia.config.settings import get_settings


def run_script(script_path, description):
    """Run a script and handle errors."""
    print(f"\n{'='*60}")
    print(f"RUNNING: {description}")
    print(f"{'='*60}")

    try:
        result = subprocess.run([sys.executable, script_path],
                              capture_output=True, text=True)

        if result.returncode == 0:
            print(f"✓ {description} completed successfully")
            if result.stdout:
                print("Output:")
                print(result.stdout[-1000:])  # Last 1000 chars
        else:
            print(f"✗ {description} failed with return code {result.returncode}")
            print("Error output:")
            print(result.stderr)
            return False

    except Exception as e:
        print(f"✗ {description} failed with exception: {e}")
        return False

    return True


def main():
    """Execute the full analysis pipeline."""
    parser = argparse.ArgumentParser(
        description='Alzheimer Microglia Meta-Analysis Pipeline'
    )
    parser.add_argument('--skip-download', action='store_true',
                       help='Skip data download step')
    parser.add_argument('--integration-method', default='scvi',
                       choices=['scvi', 'scanorama'],
                       help='Integration method to use')

    args = parser.parse_args()

    # Define script paths
    scripts_dir = Path(__file__).parent
    scripts = [
        (scripts_dir / '01_download_data.py', 'Data Download'),
        (scripts_dir / '02_preprocess_data.py', 'Data Preprocessing'),
        (scripts_dir / '03_annotate_cells.py', 'Cell Type Annotation'),
        (scripts_dir / '04_integrate_datasets.py', 'Dataset Integration'),
        (scripts_dir / '05_analyze_activation_states.py', 'Activation State Analysis'),
        (scripts_dir / '06_run_meta_analysis.py', 'Meta-Analysis Validation')
    ]

    # Skip download if requested
    if args.skip_download:
        scripts = scripts[1:]

    print(f"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                ALZHEIMER'S MICROGLIA META-ANALYSIS PIPELINE                  ║
╚══════════════════════════════════════════════════════════════════════════════╝

Pipeline Overview:
1. Data Download - Retrieve datasets from GEO
2. Initial Processing - QC, filtering, and normalization
3. Cell Annotation - Identify microglial cells
4. Dataset Integration - Harmonize datasets using {args.integration_method}
5. Activation State Analysis - Identify conserved activation states
6. Meta-Analysis - Statistical validation across studies

Starting pipeline execution...
""")

    success_count = 0
    total_steps = len(scripts)

    for script_path, description in scripts:
        if not script_path.exists():
            print(f"✗ Script not found: {script_path}")
            continue

        success = run_script(script_path, description)

        if success:
            success_count += 1
        else:
            print(f"\n⚠️  Pipeline stopped due to failure in: {description}")
            print("You can fix the issue and rerun from this step.")
            break

    settings = get_settings()

    print(f"""
\n{'='*80}
PIPELINE EXECUTION SUMMARY
{'='*80}

Completed steps: {success_count}/{total_steps}

""" + ('✓ All steps completed successfully!' if success_count == total_steps else
     f'⚠️  Stopped at step {success_count + 1}') + f"""

Results are saved in:
- {settings.PROCESSED_DATA_DIR}  - Processed datasets
- {settings.FIGURES_DIR}          - Analysis plots
- {settings.TABLES_DIR}           - Statistical results
- {settings.REPORTS_DIR}          - Summary reports
- {settings.META_ANALYSIS_DIR}    - Meta-analysis results

For questions or issues, check the documentation in docs/
""")


if __name__ == "__main__":
    main()
