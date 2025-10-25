#!/usr/bin/env python3
"""
Quick test to verify pipeline components work
"""

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

def test_basic_imports():
    """Test that all required packages can be imported"""
    try:
        print("Testing basic imports...")
        
        # Core packages
        import scanpy as sc
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns
        print("✓ Core packages imported")
        
        # Check if data directory exists
        if Path('data/raw').exists():
            print("✓ Data directory found")
        else:
            print("✗ Data directory not found")
            
        # Test basic scanpy functionality
        adata = sc.AnnData(np.random.randn(100, 50))
        sc.pp.normalize_total(adata)
        print("✓ Basic scanpy operations work")
        
        return True
        
    except Exception as e:
        print(f"✗ Import test failed: {e}")
        return False

def test_data_loading():
    """Test data loading from one dataset"""
    try:
        print("\nTesting data loading...")
        
        # Check if GSE98969 data exists
        dataset_path = Path('data/raw/GSE98969')
        if not dataset_path.exists():
            print("✗ GSE98969 dataset not found")
            return False
            
        # Look for data files
        txt_files = list(dataset_path.glob('**/*.txt.gz'))
        if txt_files:
            print(f"✓ Found {len(txt_files)} data files in GSE98969")
            
            # Try loading one file
            test_file = txt_files[0]
            print(f"Testing load of: {test_file.name}")
            
            df = pd.read_csv(test_file, index_col=0, compression='gzip', 
                           sep=None, engine='python', nrows=10)
            print(f"✓ Successfully loaded sample data: {df.shape}")
            
        return True
        
    except Exception as e:
        print(f"✗ Data loading test failed: {e}")
        return False

def main():
    """Run all tests"""
    print("="*50)
    print("PIPELINE COMPONENT TESTS")
    print("="*50)
    
    tests = [
        test_basic_imports,
        test_data_loading
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
    
    print(f"\n{'='*50}")
    print(f"TEST RESULTS: {passed}/{total} tests passed")
    print(f"{'='*50}")
    
    if passed == total:
        print("✓ All tests passed! Pipeline components are working.")
        return True
    else:
        print("✗ Some tests failed. Check the output above.")
        return False

if __name__ == "__main__":
    main()