"""
Data preprocessing module for Alzheimer's Microglia Meta-Analysis.
Handles quality control, normalization, and dataset loading.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from ad_microglia.utils.plot_manager import PlotManager, PlotContext
from ad_microglia.config.settings import get_settings

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white', figsize=(6, 6))


class DataPreprocessor:
    """Preprocess single-cell RNA-seq datasets with quality control and normalization."""

    def __init__(self, raw_data_dir=None, processed_data_dir=None):
        """
        Initialize the preprocessor.

        Parameters
        ----------
        raw_data_dir : str or Path, optional
            Directory containing raw data. Uses settings default if None.
        processed_data_dir : str or Path, optional
            Directory for processed data. Uses settings default if None.
        """
        settings = get_settings()

        self.raw_data_dir = Path(raw_data_dir) if raw_data_dir else settings.RAW_DATA_DIR
        self.processed_data_dir = Path(processed_data_dir) if processed_data_dir else settings.PROCESSED_DATA_DIR
        self.processed_data_dir.mkdir(exist_ok=True, parents=True)

        # QC parameters from settings
        self.qc_min_genes = settings.QC_MIN_GENES
        self.qc_min_cells = settings.QC_MIN_CELLS
        self.qc_max_mt_percent = settings.QC_MAX_MT_PERCENT
        self.qc_max_genes = settings.QC_MAX_GENES
        self.qc_min_counts = settings.QC_MIN_COUNTS
        self.norm_target_sum = settings.NORM_TARGET_SUM

        # Initialize plot manager
        self.plot_manager = PlotManager()

    def quality_control(self, adata, dataset_name):
        """Apply quality control filters."""
        print(f"\nQuality control for {dataset_name}")
        print(f"Initial: {adata.n_obs} cells, {adata.n_vars} genes")

        # Calculate QC metrics
        adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'Mt-', 'mt-'))
        adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL', 'Rps', 'Rpl'))
        adata.var['hb'] = adata.var_names.str.contains(r'HB[AB]|Hb[ab]', regex=True)

        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=['mt', 'ribo', 'hb'],
            percent_top=None,
            log1p=False,
            inplace=True
        )

        # Store pre-filter counts
        adata.obs['n_genes_pre_filter'] = adata.obs['n_genes_by_counts']
        adata.obs['total_counts_pre_filter'] = adata.obs['total_counts']

        # Quality filtering
        sc.pp.filter_cells(adata, min_genes=self.qc_min_genes)
        sc.pp.filter_genes(adata, min_cells=self.qc_min_cells)

        # Remove cells with too many mitochondrial genes (likely dying cells)
        adata = adata[adata.obs.pct_counts_mt < self.qc_max_mt_percent, :]

        # Remove cells with unusually high gene counts (potential doublets)
        adata = adata[adata.obs.n_genes_by_counts < self.qc_max_genes, :]

        # Remove cells with very low counts
        adata = adata[adata.obs.total_counts > self.qc_min_counts, :]

        print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes")

        return adata

    def normalize_and_scale(self, adata):
        """Normalize and scale data."""
        # Store raw counts
        adata.raw = adata

        # Normalize to target sum per cell
        sc.pp.normalize_total(adata, target_sum=self.norm_target_sum)

        # Log transform
        sc.pp.log1p(adata)

        # Find highly variable genes
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

        # Keep highly variable genes
        adata.raw = adata  # Save full data
        adata = adata[:, adata.var.highly_variable]

        # Scale to unit variance
        sc.pp.scale(adata, max_value=10)

        return adata

    def generate_qc_plots(self, adata, dataset_name, output_dir):
        """Generate quality control plots."""
        with PlotContext(
            f'{dataset_name}_qc_metrics',
            'qc',
            f'Quality control metrics for {dataset_name}',
            ['pdf', 'png'],
            figsize=(15, 10)
        ) as (fig, ax):

            # Create subplots manually since PlotContext gives us only one axis
            plt.close(fig)  # Close the single axis figure
            fig, axes = plt.subplots(2, 3, figsize=(15, 10))
            fig.suptitle(f'Quality Control Metrics - {dataset_name}', fontsize=16)

            # Total counts
            sc.pl.violin(adata, ['total_counts'],
                        jitter=0.4, multi_panel=True, ax=axes[0,0], show=False)
            axes[0,0].set_title('Total Counts')

            # Gene counts
            sc.pl.violin(adata, ['n_genes_by_counts'],
                        jitter=0.4, multi_panel=True, ax=axes[0,1], show=False)
            axes[0,1].set_title('Gene Counts')

            # Mitochondrial percentage
            sc.pl.violin(adata, ['pct_counts_mt'],
                        jitter=0.4, multi_panel=True, ax=axes[0,2], show=False)
            axes[0,2].set_title('Mitochondrial %')

            # Scatter plots
            sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=axes[1,0], show=False)
            sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=axes[1,1], show=False)
            sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_mt', ax=axes[1,2], show=False)

            plt.tight_layout()

            # Save using plot manager
            self.plot_manager.save_plot(
                fig,
                f'{dataset_name}_qc_metrics',
                'qc',
                f'Quality control metrics for {dataset_name}',
                ['pdf', 'png']
            )

    def load_raw_data(self, dataset_path, dataset_name):
        """Load raw data from various formats."""
        print(f"Loading data from {dataset_path}")

        # Look for common file formats in the dataset directory
        dataset_path = Path(dataset_path)

        # Try different file formats (including nested directories and compressed files)
        file_patterns = [
            '*.h5',
            '*.h5ad',
            '*matrix.mtx*',
            '*expression*.csv',
            '*expression*.txt',
            '*counts*.csv',
            '*counts*.txt',
            '*.csv',
            '*.txt',
            '**/*.h5',
            '**/*.h5ad'
        ]

        adata = None

        # Try standard formats first
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
                        # Handle uncompressed files
                        df = pd.read_csv(file_path, index_col=0, nrows=5, sep=None, engine='python')

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

        # If no standard formats found, try to combine multiple files (like GSE98969 and GSE129788)
        if adata is None:
            print("No standard format found, checking for multiple sample files...")

            # Find all compressed txt/csv files including in subdirectories
            txt_files = list(dataset_path.glob('**/*.txt.gz')) + list(dataset_path.glob('**/*.csv.gz'))
            txt_files += list(dataset_path.glob('**/*.txt')) + list(dataset_path.glob('**/*.csv'))

            if len(txt_files) >= 1:
                print(f"Found {len(txt_files)} sample files. Attempting to combine them...")
                adata = self.combine_sample_files(txt_files, dataset_name)

        if adata is None:
            print(f"No compatible data files found in {dataset_path}")
            return None

        # Basic preprocessing
        adata.var_names_make_unique()
        adata.obs['dataset'] = dataset_name

        return adata

    def combine_sample_files(self, file_paths, dataset_name):
        """Combine multiple sample files into one AnnData object."""
        print(f"Combining {len(file_paths)} sample files for {dataset_name}")

        combined_data = []
        sample_names = []

        # For GSE129788, process all files but limit for other datasets during testing
        max_files = len(file_paths) if 'GSE129788' in dataset_name else 5

        for i, file_path in enumerate(file_paths[:max_files]):
            try:
                compression = 'gzip' if str(file_path).endswith('.gz') else None

                # Load the data - handle tab-separated format properly
                try:
                    # Read with explicit tab separator and handle the first column as index
                    df = pd.read_csv(file_path, index_col=0, compression=compression, sep='\t', low_memory=False)
                except Exception as e1:
                    try:
                        # Fallback to auto-detect separator
                        df = pd.read_csv(file_path, index_col=0, compression=compression, sep=None, engine='python', low_memory=False)
                    except Exception as e2:
                        print(f"  Error parsing {file_path}: {e1}, {e2}")
                        continue

                # Extract sample name from file path (handle GSE129788 format)
                if 'GSM' in str(file_path):
                    # Extract GSM ID from directory name
                    parts = file_path.parent.name.split('_')
                    if len(parts) > 1:
                        gsm_id = parts[1]
                        # Extract condition from GSM ID or directory name
                        if 'YX' in gsm_id or 'Young' in file_path.parent.name:
                            condition = 'Young'
                        elif 'OX' in gsm_id or 'Old' in file_path.parent.name:
                            condition = 'Old'
                        else:
                            condition = 'Unknown'
                        sample_name = f"{gsm_id}_{condition}"
                    else:
                        sample_name = file_path.stem.replace('.txt', '').replace('.csv', '').replace('.gz', '')
                else:
                    sample_name = file_path.stem.replace('.txt', '').replace('.csv', '').replace('.gz', '')

                print(f"  Loading sample {sample_name}: {df.shape}")

                # Convert to numeric data, handling any string/mixed columns
                for col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)

                # For 10X format (genes as rows, which is standard), keep as is
                # Add sample prefix to cell barcodes
                df.columns = [f"{sample_name}_{col}" for col in df.columns]
                combined_data.append(df)
                sample_names.extend([sample_name] * df.shape[1])

            except Exception as e:
                print(f"  Error loading {file_path}: {e}")
                continue

        if not combined_data:
            return None

        # Concatenate all samples
        print("Concatenating samples...")
        try:
            # Use outer join to include all genes from all samples
            combined_df = pd.concat(combined_data, axis=1, join='outer', sort=True)
            # Fill NaN values with 0 (genes not detected in some samples)
            combined_df = combined_df.fillna(0)
        except Exception as e:
            print(f"Error concatenating samples: {e}")
            return None

        print(f"Combined data shape: {combined_df.shape}")

        # Ensure all data is numeric
        combined_df = combined_df.astype('float32')

        # Convert to AnnData (genes as vars, cells as obs)
        # combined_df has genes as rows and cells as columns, so transpose for AnnData
        adata = sc.AnnData(X=combined_df.T.values,
                          obs=pd.DataFrame(index=combined_df.columns),
                          var=pd.DataFrame(index=combined_df.index))

        # Add sample information
        adata.obs['sample'] = sample_names

        # Extract condition information for GSE129788
        if 'GSE129788' in dataset_name:
            adata.obs['condition'] = [name.split('_')[-1] if '_' in name else 'Unknown' for name in sample_names]
            adata.obs['age_group'] = adata.obs['condition'].map({'Young': 'Young', 'Old': 'Old'}).fillna('Unknown')

        return adata

    def process_dataset(self, dataset_path, dataset_name):
        """
        Process a single dataset through the full pipeline.

        Parameters
        ----------
        dataset_path : str or Path
            Path to raw dataset directory
        dataset_name : str
            Name of the dataset

        Returns
        -------
        bool
            True if successful, False otherwise
        """
        print(f"\n{'='*50}")
        print(f"Processing {dataset_name}")
        print(f"{'='*50}")

        try:
            # Load raw data
            adata = self.load_raw_data(dataset_path, dataset_name)

            if adata is None:
                print(f"Failed to load {dataset_name}")
                return False

            print(f"Dataset {dataset_name} loaded successfully")

            # Apply preprocessing pipeline
            adata = self.quality_control(adata, dataset_name)
            adata = self.normalize_and_scale(adata)

            # Generate QC plots
            self.generate_qc_plots(adata, dataset_name, self.processed_data_dir)

            # Save processed data
            output_path = self.processed_data_dir / f'{dataset_name}_processed.h5ad'
            adata.write(output_path)
            print(f"Processed data saved to {output_path}")

            return True

        except Exception as e:
            print(f"Error processing {dataset_name}: {e}")
            import traceback
            traceback.print_exc()
            return False
