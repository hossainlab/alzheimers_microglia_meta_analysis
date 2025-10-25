#!/usr/bin/env python3
"""
Quick test to verify pipeline works end-to-end on one small dataset
"""

import sys
sys.path.append('.')

import importlib.util
spec = importlib.util.spec_from_file_location("preprocess_data", "scripts/preprocessing/02_preprocess_data.py")
preprocess_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preprocess_module)
DataPreprocessor = preprocess_module.DataPreprocessor
from pathlib import Path
import traceback

def main():
    """Test preprocessing on one dataset"""
    print("Testing preprocessing pipeline...")
    
    preprocessor = DataPreprocessor()
    
    # Test on GSE98969 (first dataset)
    dataset_name = 'GSE98969'
    dataset_path = 'data/raw/GSE98969'
    
    try:
        success = preprocessor.process_dataset(dataset_path, dataset_name)
        
        if success:
            print(f"✓ {dataset_name} processed successfully!")
            
            # Check if output file exists
            output_path = Path(f'data/processed/{dataset_name}_processed.h5ad')
            if output_path.exists():
                print(f"✓ Output file created: {output_path}")
                import scanpy as sc
                adata = sc.read_h5ad(output_path)
                print(f"✓ Processed data: {adata.n_obs} cells, {adata.n_vars} genes")
            else:
                print("✗ Output file not found")
                
        else:
            print(f"✗ {dataset_name} processing failed")
            
    except Exception as e:
        print(f"✗ Error: {e}")
        traceback.print_exc()

if __name__ == "__main__":
    main()