#!/usr/bin/env python3
"""
Cell Type Annotation Script for Alzheimer's Microglia Meta-Analysis
Uses CellTypist to identify microglial cells across datasets
"""

import scanpy as sc
import pandas as pd
import numpy as np
import celltypist
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import warnings
warnings.filterwarnings('ignore')

# Add utils to path
sys.path.append(str(Path(__file__).parent.parent))
from utils.plot_manager import PlotManager

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')

class CellAnnotator:
    def __init__(self, processed_data_dir='data/processed'):
        self.processed_data_dir = Path(processed_data_dir)
        self.results_dir = Path('results/figures')
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize plot manager
        self.plot_manager = PlotManager()
        
    def load_celltypist_model(self):
        """Load CellTypist model for immune cells"""
        print("Loading CellTypist model...")
        try:
            # Download and load the Immune_All_High model
            celltypist.models.download_models(force_update=True)
            model = celltypist.models.Model.load(model='Immune_All_High.pkl')
            return model
        except Exception as e:
            print(f"Error loading CellTypist model: {e}")
            print("Using fallback marker-based annotation")
            return None
    
    def annotate_with_celltypist(self, adata, model):
        """Annotate cells using CellTypist"""
        print("Running CellTypist annotation...")
        
        # Ensure proper gene names
        adata.var_names_unique()
        
        # Run prediction
        predictions = celltypist.annotate(adata, model=model, majority_voting=True)
        
        # Add predictions to adata
        adata.obs['celltypist_predicted_labels'] = predictions.predicted_labels.predicted_labels
        adata.obs['celltypist_conf_score'] = predictions.predicted_labels.conf_score
        
        return adata
    
    def annotate_with_markers(self, adata):
        """Fallback marker-based annotation for microglia"""
        print("Using marker-based annotation...")
        
        # Define microglia markers (both human and mouse versions)
        microglia_markers = [
            'CX3CR1', 'Cx3cr1',    # CX3C chemokine receptor 1
            'P2RY12', 'P2ry12',    # Purinergic receptor P2Y12
            'TMEM119', 'Tmem119',  # Transmembrane protein 119
            'CSF1R', 'Csf1r',     # Colony stimulating factor 1 receptor
            'AIF1', 'Aif1',       # Allograft inflammatory factor 1
            'ITGAM', 'Itgam'      # Integrin alpha M
        ]
        
        # Check which markers are available
        available_markers = [m for m in microglia_markers if m in adata.var_names]
        
        if not available_markers:
            print("Warning: No microglia markers found in data")
            adata.obs['is_microglia'] = False
            return adata
        
        print(f"Using markers: {available_markers}")
        
        # Calculate marker expression scores
        sc.tl.score_genes(adata, available_markers, score_name='microglia_score')
        
        # Set threshold for microglia identification
        threshold = np.percentile(adata.obs['microglia_score'], 75)
        adata.obs['is_microglia'] = adata.obs['microglia_score'] > threshold
        
        # Create cell type labels
        adata.obs['cell_type'] = 'Other'
        adata.obs.loc[adata.obs['is_microglia'], 'cell_type'] = 'Microglia'
        
        return adata
    
    def filter_microglia(self, adata):
        """Filter for high-confidence microglial cells"""
        print("Filtering for microglial cells...")
        
        if 'celltypist_predicted_labels' in adata.obs:
            # Use CellTypist predictions
            microglia_mask = (
                adata.obs['celltypist_predicted_labels'].str.contains('Microglia', case=False) |
                adata.obs['celltypist_predicted_labels'].str.contains('Macrophage', case=False)
            ) & (adata.obs['celltypist_conf_score'] > 0.5)
        else:
            # Use marker-based annotation
            microglia_mask = adata.obs['is_microglia']
        
        print(f"Found {microglia_mask.sum()} microglial cells out of {len(adata)} total cells")
        
        if microglia_mask.sum() < 50:
            print("Warning: Very few microglia found. Check annotation parameters.")
        
        return adata[microglia_mask].copy()
    
    def validate_annotation(self, adata, dataset_name):
        """Validate microglia annotation using known markers"""
        print(f"Validating annotation for {dataset_name}")
        
        # Microglia validation markers (both human and mouse versions)
        validation_markers = {
            'core': ['CX3CR1', 'Cx3cr1', 'P2RY12', 'P2ry12', 'TMEM119', 'Tmem119'],
            'extended': ['CSF1R', 'Csf1r', 'AIF1', 'Aif1', 'ITGAM', 'Itgam', 'CD68', 'Cd68', 'TREM2', 'Trem2']
        }
        
        validation_results = {}
        
        for marker_set, markers in validation_markers.items():
            available_markers = [m for m in markers if m in adata.var_names]
            if available_markers:
                # Calculate mean expression
                mean_expr = adata[:, available_markers].X.mean(axis=1)
                validation_results[f'{marker_set}_expression'] = np.array(mean_expr).flatten()
        
        # Add validation results to obs
        for key, values in validation_results.items():
            adata.obs[key] = values
        
        # Create validation plots
        self.plot_validation(adata, dataset_name, validation_results)
        
        return adata
    
    def plot_validation(self, adata, dataset_name, validation_results):
        """Create validation plots"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f'Microglia Annotation Validation - {dataset_name}')
        
        # Plot 1: Cell type distribution
        if 'celltypist_predicted_labels' in adata.obs:
            cell_counts = adata.obs['celltypist_predicted_labels'].value_counts().head(10)
            axes[0, 0].bar(range(len(cell_counts)), cell_counts.values)
            axes[0, 0].set_xticks(range(len(cell_counts)))
            axes[0, 0].set_xticklabels(cell_counts.index, rotation=45, ha='right')
            axes[0, 0].set_title('Predicted Cell Types')
            axes[0, 0].set_ylabel('Cell Count')
        
        # Plot 2: Confidence score distribution
        if 'celltypist_conf_score' in adata.obs:
            axes[0, 1].hist(adata.obs['celltypist_conf_score'], bins=30, alpha=0.7)
            axes[0, 1].axvline(0.5, color='red', linestyle='--', label='Threshold')
            axes[0, 1].set_xlabel('Confidence Score')
            axes[0, 1].set_ylabel('Cell Count')
            axes[0, 1].set_title('CellTypist Confidence')
            axes[0, 1].legend()
        
        # Plot 3: Marker expression
        if 'core_expression' in validation_results:
            axes[1, 0].hist(validation_results['core_expression'], bins=30, alpha=0.7)
            axes[1, 0].set_xlabel('Core Marker Expression')
            axes[1, 0].set_ylabel('Cell Count')
            axes[1, 0].set_title('Core Microglia Markers')
        
        # Plot 4: QC metrics
        if 'total_counts' in adata.obs:
            axes[1, 1].scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'], 
                             alpha=0.5, s=1)
            axes[1, 1].set_xlabel('Total UMI Counts')
            axes[1, 1].set_ylabel('Number of Genes')
            axes[1, 1].set_title('Quality Metrics')
        
        plt.tight_layout()
        
        # Save using plot manager
        self.plot_manager.save_plot(
            fig, 
            f'{dataset_name}_annotation_validation',
            'annotation',
            f'Cell type annotation validation for {dataset_name}',
            ['png', 'pdf']
        )
    
    def process_dataset(self, dataset_path):
        """Process single dataset for cell annotation"""
        dataset_name = dataset_path.stem
        print(f"\nProcessing {dataset_name} for cell annotation...")
        
        # Load processed data
        adata = sc.read_h5ad(dataset_path)
        print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
        
        # For now, skip CellTypist and use marker-based annotation only
        print("Using marker-based annotation (skipping CellTypist for faster processing)")
        adata = self.annotate_with_markers(adata)
        
        # Filter for microglia
        microglia_adata = self.filter_microglia(adata)
        
        # Validate annotation
        microglia_adata = self.validate_annotation(microglia_adata, dataset_name)
        
        # Save annotated microglia data
        output_path = self.processed_data_dir / f"{dataset_name}_microglia.h5ad"
        microglia_adata.write(output_path)
        print(f"Saved {microglia_adata.n_obs} microglia to {output_path}")
        
        return microglia_adata
    
    def run_annotation(self):
        """Run cell annotation on all processed datasets"""
        print("Starting cell type annotation...")
        
        # Find all processed datasets
        processed_files = list(self.processed_data_dir.glob("*_processed.h5ad"))
        
        if not processed_files:
            print("No processed datasets found. Run preprocessing first.")
            return
        
        annotated_datasets = []
        
        for dataset_path in processed_files:
            try:
                microglia_data = self.process_dataset(dataset_path)
                annotated_datasets.append({
                    'dataset': dataset_path.stem.replace('_processed', ''),
                    'n_microglia': microglia_data.n_obs,
                    'file_path': str(dataset_path).replace('_processed', '_microglia')
                })
            except Exception as e:
                print(f"Error processing {dataset_path}: {e}")
                continue
        
        # Create summary report
        summary_df = pd.DataFrame(annotated_datasets)
        summary_path = Path('results/tables') / 'microglia_annotation_summary.csv'
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(summary_path, index=False)
        
        print(f"\nAnnotation complete! Summary saved to {summary_path}")
        print(f"Total datasets processed: {len(annotated_datasets)}")
        print(f"Total microglia identified: {summary_df['n_microglia'].sum()}")

def main():
    annotator = CellAnnotator()
    annotator.run_annotation()

if __name__ == "__main__":
    main()