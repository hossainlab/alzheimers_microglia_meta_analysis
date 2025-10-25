#!/usr/bin/env python3
"""
Cell Type Annotation Script for Alzheimer's Microglia Meta-Analysis
Uses CellTypist for automated annotation and manual validation
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import celltypist
from celltypist import models

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')

class CellTypeAnnotator:
    def __init__(self, processed_data_dir='data/processed'):
        self.processed_data_dir = Path(processed_data_dir)
        
        # Download appropriate models
        print("Checking CellTypist models...")
        self.setup_models()
        
    def setup_models(self):
        """Download and setup CellTypist models"""
        # Check available models
        models.download_models(force_update=False)
        available_models = models.models_description()
        print("Available models:")
        print(available_models[['Model', 'Tissue', 'Description']].head(10))
        
        # Select appropriate models for brain tissue
        self.brain_models = [
            'Human_Brain_Atlas.pkl',  # Human brain atlas
            'Adult_Mouse_Brain.pkl',  # Mouse brain atlas  
            'Immune_All_Low.pkl',     # Immune cells (includes microglia)
            'Immune_All_High.pkl'     # High resolution immune cells
        ]
        
    def annotate_with_celltypist(self, adata, model_name, dataset_name):
        """Annotate cells using CellTypist"""
        print(f"\nAnnotating {dataset_name} with {model_name}")
        
        try:
            # Load model
            model = models.Model.load(model=model_name)
            
            # Predict cell types
            predictions = celltypist.annotate(
                adata, 
                model=model_name,
                majority_voting=True
            )
            
            # Add predictions to adata
            adata.obs[f'celltypist_{model_name}'] = predictions.predicted_labels.predicted_labels
            adata.obs[f'celltypist_{model_name}_conf'] = predictions.predicted_labels.conf_score
            
            # Add probability matrix
            prob_matrix = predictions.probability_matrix
            for cell_type in prob_matrix.columns:
                adata.obs[f'{model_name}_prob_{cell_type}'] = prob_matrix[cell_type].values
                
            print(f"Annotation completed for {dataset_name}")
            return adata, predictions
            
        except Exception as e:
            print(f"Error in annotation: {e}")
            return adata, None
            
    def identify_microglia(self, adata, dataset_name):
        """Identify and extract microglial cells"""
        print(f"\nIdentifying microglia in {dataset_name}")
        
        # Define microglial markers
        microglia_markers = ['CX3CR1', 'P2RY12', 'TMEM119', 'AIF1', 'CSF1R', 'TREM2']
        
        # Calculate microglial signature score
        if any(marker in adata.var_names for marker in microglia_markers):
            available_markers = [m for m in microglia_markers if m in adata.var_names]
            print(f"Using markers: {available_markers}")
            
            # Score cells for microglial signature
            sc.tl.score_genes(adata, available_markers, ctrl_size=50, 
                             gene_pool=None, n_bins=25, score_name='microglia_score')
        
        # Identify microglia based on:
        # 1. CellTypist predictions
        # 2. Gene expression signatures
        # 3. Manual curation
        
        microglia_mask = pd.Series(False, index=adata.obs_names)
        
        # From CellTypist predictions
        for col in adata.obs.columns:
            if 'celltypist' in col and 'conf' not in col and 'prob' not in col:
                microglia_predictions = adata.obs[col].str.contains(
                    'Microglia|microglia|Macro|macro', case=False, na=False
                )
                microglia_mask |= microglia_predictions
                
        # From gene signature (if available)
        if 'microglia_score' in adata.obs.columns:
            high_microglia_score = adata.obs['microglia_score'] > adata.obs['microglia_score'].quantile(0.8)
            microglia_mask |= high_microglia_score
            
        # Create microglia subset
        adata.obs['is_microglia'] = microglia_mask
        microglia_adata = adata[microglia_mask].copy()
        
        print(f"Identified {microglia_adata.n_obs} microglial cells ({microglia_adata.n_obs/adata.n_obs*100:.1f}%)")
        
        return microglia_adata
        
    def validate_annotations(self, adata, predictions, dataset_name, output_dir):
        """Generate validation plots for annotations"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f'Cell Type Annotations - {dataset_name}', fontsize=16)
        
        # Cell type composition
        cell_counts = adata.obs['celltypist_Immune_All_High.pkl'].value_counts()
        cell_counts.head(10).plot(kind='bar', ax=axes[0,0])
        axes[0,0].set_title('Cell Type Distribution')
        axes[0,0].tick_params(axis='x', rotation=45)
        
        # Confidence scores
        adata.obs['celltypist_Immune_All_High.pkl_conf'].hist(bins=50, ax=axes[0,1])
        axes[0,1].set_title('Annotation Confidence')
        axes[0,1].set_xlabel('Confidence Score')
        
        # UMAP of cell types (if coordinates available)
        if 'X_umap' in adata.obsm:
            sc.pl.umap(adata, color='celltypist_Immune_All_High.pkl', 
                      ax=axes[1,0], show=False, frameon=False)
            axes[1,0].set_title('Cell Types (UMAP)')
            
        # Microglial signature
        if 'microglia_score' in adata.obs.columns:
            sc.pl.umap(adata, color='microglia_score', 
                      ax=axes[1,1], show=False, frameon=False)
            axes[1,1].set_title('Microglial Signature')
            
        plt.tight_layout()
        plt.savefig(output_dir / f'{dataset_name}_annotations.pdf', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
    def process_dataset(self, dataset_file):
        """Process annotations for a single dataset"""
        dataset_name = dataset_file.stem.replace('_processed', '')
        print(f"\nProcessing annotations for {dataset_name}")
        
        try:
            # Load processed data
            adata = sc.read_h5ad(dataset_file)
            
            # Ensure data is properly formatted for CellTypist
            # Need log-normalized data with 10k scaling
            if adata.raw is not None:
                adata = adata.raw.to_adata()
                sc.pp.normalize_total(adata, target_sum=1e4)
                sc.pp.log1p(adata)
            
            # Run multiple annotation models
            for model_name in self.brain_models:
                try:
                    adata, predictions = self.annotate_with_celltypist(
                        adata, model_name, dataset_name
                    )
                except Exception as e:
                    print(f"Failed to annotate with {model_name}: {e}")
                    continue
                    
            # Identify microglial populations
            microglia_adata = self.identify_microglia(adata, dataset_name)
            
            # Generate validation plots
            if predictions is not None:
                self.validate_annotations(
                    adata, predictions, dataset_name, 
                    Path('results/figures')
                )
            
            # Save annotated data
            output_path = self.processed_data_dir / f'{dataset_name}_annotated.h5ad'
            adata.write(output_path)
            
            # Save microglia subset
            microglia_path = self.processed_data_dir / f'{dataset_name}_microglia.h5ad'
            microglia_adata.write(microglia_path)
            
            print(f"Annotations saved to {output_path}")
            print(f"Microglia subset saved to {microglia_path}")
            
            return True, microglia_adata.n_obs
            
        except Exception as e:
            print(f"Error processing {dataset_name}: {e}")
            return False, 0

def main():
    """Main annotation pipeline"""
    annotator = CellTypeAnnotator()
    
    # Find processed datasets
    processed_files = list(Path('data/processed').glob('*_processed.h5ad'))
    
    if not processed_files:
        print("No processed datasets found. Run preprocessing first.")
        return
        
    print(f"Found {len(processed_files)} datasets to annotate")
    
    # Process annotations
    results = {}
    total_microglia = 0
    
    for dataset_file in processed_files:
        success, n_microglia = annotator.process_dataset(dataset_file)
        dataset_name = dataset_file.stem.replace('_processed', '')
        results[dataset_name] = {
            'success': success,
            'n_microglia': n_microglia
        }
        total_microglia += n_microglia
        
    # Generate summary
    print(f"\n{'='*50}")
    print("ANNOTATION SUMMARY")
    print(f"{'='*50}")
    
    successful_datasets = sum(1 for r in results.values() if r['success'])
    print(f"Successfully annotated: {successful_datasets}/{len(processed_files)} datasets")
    print(f"Total microglia identified: {total_microglia}")
    
    for dataset, result in results.items():
        status = "✓" if result['success'] else "✗"
        print(f"{status} {dataset}: {result['n_microglia']} microglia")
        
    # Save summary
    import json
    with open('data/metadata/annotation_summary.json', 'w') as f:
        json.dump(results, f, indent=2)
        
    print("\nNext step: Run integration pipeline")

if __name__ == "__main__":
    main()
