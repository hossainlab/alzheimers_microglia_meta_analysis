#!/usr/bin/env python3
"""
Data Download Script for Alzheimer's Microglia Meta-Analysis
Downloads datasets from GEO and processes them for analysis
"""

import os
import pandas as pd
import scanpy as sc
import numpy as np
from pathlib import Path
import urllib.request
import gzip
import tarfile

# Configure scanpy
sc.settings.verbosity = 3  # verbosity level
sc.settings.set_figure_params(dpi=80, facecolor='white')

def download_geo_dataset(geo_id, output_dir):
    """Download dataset from GEO using GEOparse"""
    print(f"Downloading {geo_id}...")
    
    try:
        import GEOparse
    except ImportError:
        print("GEOparse not available. Install with: pip install GEOparse")
        return None
    
    dataset_dir = os.path.join(output_dir, geo_id)
    os.makedirs(dataset_dir, exist_ok=True)
    
    try:
        # Download GEO dataset
        gse = GEOparse.get_GEO(geo=geo_id, destdir=dataset_dir)
        print(f"Successfully downloaded {geo_id}")
        print(f"Title: {gse.metadata.get('title', ['Unknown'])[0]}")
        print(f"Samples: {len(gse.gsms)}")
        
        # Look for supplementary files
        if hasattr(gse, 'download_supplementary_files'):
            print("Downloading supplementary files...")
            gse.download_supplementary_files(directory=dataset_dir)
        
        return dataset_dir
        
    except Exception as e:
        print(f"Error downloading {geo_id}: {e}")
        print("Creating placeholder directory for manual download")
        
        # Create instructions file
        with open(os.path.join(dataset_dir, 'DOWNLOAD_INSTRUCTIONS.txt'), 'w') as f:
            f.write(f"""Manual Download Instructions for {geo_id}

1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={geo_id}
2. Download the following files:
   - Raw count matrix (usually .txt, .csv, or .h5 format)
   - Cell metadata (if available)
   - Gene annotation (if available)
3. Place files in this directory: {dataset_dir}
4. Ensure files follow naming convention: {geo_id}_expression.csv, {geo_id}_metadata.csv

Common file patterns to look for:
- *_matrix.mtx (Matrix Market format)
- *_barcodes.tsv (cell identifiers) 
- *_features.tsv or *_genes.tsv (gene information)
- *.h5 (HDF5 format)
- *.csv or *.txt (comma/tab separated)
""")
        
        return dataset_dir

def process_count_matrix(dataset_dir, geo_id):
    """Process raw count matrix into AnnData object"""
    print(f"Processing {geo_id}...")
    
    # Look for common file formats in the dataset directory
    dataset_path = Path(dataset_dir)
    
    # Handle specific dataset formats
    if geo_id == 'GSE103334':
        return process_geo103334(dataset_path, geo_id)
    elif geo_id == 'GSE135437':
        return process_geo135437(dataset_path, geo_id)
    elif geo_id == 'GSE157827':
        return process_geo157827(dataset_path, geo_id)
    
    # Try different file formats for other datasets
    file_patterns = [
        '*.h5',
        '*.h5ad', 
        '*matrix.mtx*',
        '*expression*.csv',
        '*expression*.txt',
        '*counts*.csv',
        '*counts*.txt',
        '*.csv',
        '*.txt'
    ]
    
    adata = None
    
    for pattern in file_patterns:
        files = list(dataset_path.glob(pattern))
        if files:
            file_path = files[0]  # Take first match
            print(f"Found data file: {file_path}")
            
            try:
                if file_path.suffix == '.h5':
                    adata = sc.read_10x_h5(file_path)
                elif file_path.suffix == '.h5ad':
                    adata = sc.read_h5ad(file_path)
                elif 'mtx' in file_path.name:
                    # Look for barcodes and features files
                    barcodes_files = list(dataset_path.glob('*barcode*'))
                    features_files = list(dataset_path.glob('*feature*')) or list(dataset_path.glob('*gene*'))
                    
                    if barcodes_files and features_files:
                        adata = sc.read_10x_mtx(
                            dataset_path,
                            var_names='gene_symbols',
                            cache=True
                        )
                    else:
                        print("MTX format requires barcodes and features files")
                        continue
                elif file_path.suffix in ['.csv', '.txt']:
                    # Try to determine if genes are rows or columns
                    df = pd.read_csv(file_path, index_col=0, nrows=5)
                    
                    # Heuristic: if first column looks like gene names, transpose
                    if df.index[0].startswith(('ENS', 'ENSG')) or len(df.index[0]) > 10:
                        adata = sc.read_csv(file_path, first_column_names=True).T
                    else:
                        adata = sc.read_csv(file_path, first_column_names=True)
                else:
                    print(f"Unsupported file format: {file_path}")
                    continue
                    
                if adata is not None:
                    print(f"Successfully loaded: {adata.n_obs} cells, {adata.n_vars} genes")
                    break
                    
            except Exception as e:
                print(f"Error loading {file_path}: {e}")
                continue
    
    if adata is None:
        print(f"No compatible data files found in {dataset_dir}")
        return None
    
    # Basic preprocessing
    adata.var_names_make_unique()
    adata.obs['dataset'] = geo_id
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    
    # Add mitochondrial gene info
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    print(f"Processed {geo_id}: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata

def process_geo103334(dataset_path, geo_id):
    """Process GSE103334 dataset - single cell RNA-seq from macrophages"""
    print("Processing GSE103334 (macrophage single cell data)...")
    
    # This dataset has individual sample directories
    sample_dirs = [d for d in dataset_path.iterdir() if d.is_dir() and d.name.startswith('Supp_')]
    
    if not sample_dirs:
        print("No sample directories found")
        return None
    
    print(f"Found {len(sample_dirs)} sample directories")
    
    # For this dataset, we need to check the GEO page for proper processing
    # Creating placeholder since this requires specialized handling
    print("GSE103334 requires manual processing - creating placeholder")
    return None

def process_geo135437(dataset_path, geo_id):
    """Process GSE135437 dataset - MS brain single cell data"""
    print("Processing GSE135437 (MS brain single cell data)...")
    
    # Find all .coutt.tsv.gz files (count matrices)
    count_files = list(dataset_path.glob('**/*.coutt.tsv.gz'))
    
    if not count_files:
        print("No count matrix files found")
        return None
    
    print(f"Found {len(count_files)} count matrix files")
    
    # Combine all samples into single AnnData object
    all_data = []
    sample_names = []
    
    for count_file in count_files[:5]:  # Process first 5 samples as test
        try:
            print(f"Processing {count_file.name}...")
            
            # Read count matrix
            df = pd.read_csv(count_file, sep='\t', index_col=0, compression='gzip')
            
            # Create AnnData object
            sample_adata = sc.AnnData(df.T)  # Transpose so cells are rows
            sample_adata.obs['sample'] = count_file.stem.replace('.coutt.tsv', '')
            sample_adata.var_names_make_unique()
            
            all_data.append(sample_adata)
            sample_names.append(count_file.stem)
            
        except Exception as e:
            print(f"Error processing {count_file}: {e}")
            continue
    
    if not all_data:
        print("No samples successfully processed")
        return None
    
    # Concatenate all samples
    adata = sc.concat(all_data, join='outer', index_unique='-')
    adata.obs['dataset'] = geo_id
    
    print(f"Combined dataset: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata

def process_geo157827(dataset_path, geo_id):
    """Process GSE157827 dataset - 10X format single cell data"""
    print("Processing GSE157827 (10X format single cell data)...")
    
    # Find sample directories with 10X format files
    sample_dirs = []
    for subdir in dataset_path.rglob('*'):
        if subdir.is_dir():
            mtx_files = list(subdir.glob('*matrix.mtx.gz'))
            if mtx_files:
                sample_dirs.append(subdir)
    
    if not sample_dirs:
        print("No 10X format directories found")
        return None
    
    print(f"Found {len(sample_dirs)} 10X format directories")
    
    # Process each sample and combine
    all_data = []
    
    for sample_dir in sample_dirs[:3]:  # Process first 3 samples as test
        try:
            print(f"Processing {sample_dir.name}...")
            
            # Read 10X format data
            sample_adata = sc.read_10x_mtx(
                sample_dir,
                var_names='gene_symbols',
                cache=True
            )
            
            sample_adata.var_names_make_unique()
            sample_adata.obs['sample'] = sample_dir.name
            all_data.append(sample_adata)
            
        except Exception as e:
            print(f"Error processing {sample_dir}: {e}")
            continue
    
    if not all_data:
        print("No samples successfully processed")
        return None
    
    # Concatenate all samples
    adata = sc.concat(all_data, join='outer', index_unique='-')
    adata.obs['dataset'] = geo_id
    
    print(f"Combined dataset: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata

def main():
    # Dataset list from catalog
    datasets = [
        'GSE98969',
        'GSE103334', 
        'GSE135437',
        'GSE157827',
        'GSE129788'
    ]
    
    raw_data_dir = 'data/raw'
    processed_data_dir = 'data/processed'
    
    os.makedirs(raw_data_dir, exist_ok=True)
    os.makedirs(processed_data_dir, exist_ok=True)
    
    # Download and process each dataset
    processed_datasets = {}
    
    for geo_id in datasets:
        print(f"\n{'='*50}")
        print(f"PROCESSING {geo_id}")
        print(f"{'='*50}")
        
        try:
            # Download dataset
            dataset_dir = download_geo_dataset(geo_id, raw_data_dir)
            
            if dataset_dir is None:
                print(f"Failed to download {geo_id}")
                continue
            
            # Process count matrix
            adata = process_count_matrix(dataset_dir, geo_id)
            
            if adata is not None:
                # Save processed data
                output_file = os.path.join(processed_data_dir, f"{geo_id}_raw.h5ad")
                adata.write(output_file)
                processed_datasets[geo_id] = {
                    'file_path': output_file,
                    'n_cells': adata.n_obs,
                    'n_genes': adata.n_vars
                }
                print(f"Saved to {output_file}")
            else:
                print(f"Failed to process {geo_id}")
            
        except Exception as e:
            print(f"Failed to process {geo_id}: {e}")
            continue
    
    # Create summary report
    if processed_datasets:
        summary_df = pd.DataFrame.from_dict(processed_datasets, orient='index')
        summary_df.index.name = 'dataset'
        summary_file = os.path.join(processed_data_dir, 'download_summary.csv')
        summary_df.to_csv(summary_file)
        
        print(f"\n{'='*60}")
        print("DOWNLOAD SUMMARY")
        print(f"{'='*60}")
        print(f"Successfully processed: {len(processed_datasets)} datasets")
        print(f"Total cells: {summary_df['n_cells'].sum():,}")
        print(f"Summary saved to: {summary_file}")
        print(f"\nDatasets processed:")
        for dataset, info in processed_datasets.items():
            print(f"  {dataset}: {info['n_cells']:,} cells, {info['n_genes']:,} genes")
    else:
        print("No datasets were successfully processed!")
        print("Check download instructions in data/raw/ directories")
    
    print("\nNext step: Run preprocessing pipeline")

if __name__ == "__main__":
    main()
