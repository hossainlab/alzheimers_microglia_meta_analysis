#!/usr/bin/env python3
"""
Main Pipeline for Alzheimer's Microglia Meta-Analysis
Executes complete analysis workflow from data download to final results
"""

import subprocess
import sys
from pathlib import Path
import argparse

def run_script(script_path, description):
    """Run a script and handle errors"""
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
    parser = argparse.ArgumentParser(description='Alzheimer Microglia Meta-Analysis Pipeline')
    parser.add_argument('--skip-download', action='store_true',
                       help='Skip data download step')
    parser.add_argument('--integration-method', default='scvi',
                       choices=['scvi', 'scanorama'], 
                       help='Integration method to use')
    
    args = parser.parse_args()
    
    # Define script paths
    base_dir = Path('.')
    scripts = [
        ('scripts/preprocessing/01_download_data.py', 'Data Download'),
        ('scripts/preprocessing/02_preprocess_data.py', 'Data Preprocessing'),
        ('scripts/preprocessing/03_cell_annotation.py', 'Cell Type Annotation'),
        ('scripts/integration/01_integrate_datasets.py', 'Dataset Integration'),
        ('scripts/analysis/01_activation_states.py', 'Activation State Analysis'),
        ('scripts/analysis/02_meta_analysis.py', 'Meta-Analysis Validation')
    ]
    
    # Skip download if requested
    if args.skip_download:
        scripts = scripts[1:]
        
    print(f"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                ALZHEIMER'S MICROGLIA META-ANALYSIS PIPELINE                  ║
╚══════════════════════════════════════════════════════════════════════════════╝

Pipeline Overview:
1. Data Download - Retrieve datasets from GEO and other sources
2. Initial Processing - QC, filtering, and normalization
3. Cell Annotation - Identify microglial cells using CellTypist
4. Dataset Integration - Harmonize datasets using {args.integration_method}
5. Activation State Analysis - Identify conserved activation states
6. Meta-Analysis - Statistical validation across studies

Starting pipeline execution...
""")
    
    success_count = 0
    total_steps = len(scripts)
    
    for script_path, description in scripts:
        full_path = base_dir / script_path
        
        if not full_path.exists():
            print(f"✗ Script not found: {full_path}")
            continue
            
        success = run_script(full_path, description)
        
        if success:
            success_count += 1
        else:
            print(f"\n⚠️  Pipeline stopped due to failure in: {description}")
            print("You can fix the issue and rerun from this step.")
            break
            
    print(f"""
\n{'='*80}
PIPELINE EXECUTION SUMMARY
{'='*80}

Completed steps: {success_count}/{total_steps}

""" + ('✓ All steps completed successfully!' if success_count == total_steps else 
     f'⚠️  Stopped at step {success_count + 1}: {scripts[success_count][1]}') + """

Results are saved in:
- data/processed/      - Processed datasets  
- results/figures/     - Analysis plots
- results/tables/      - Statistical results
- results/reports/     - Summary reports
- results/meta_analysis/ - Meta-analysis results

Next steps:
1. Review results in results/reports/
2. Examine key figures in results/figures/
3. Check meta-analysis validation in results/meta_analysis/

For questions or issues, check the documentation in docs/
""")

if __name__ == "__main__":
    main()
