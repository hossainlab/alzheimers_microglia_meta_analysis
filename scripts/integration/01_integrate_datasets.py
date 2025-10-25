#!/usr/bin/env python3
"""
Integration Pipeline for Alzheimer's Microglia Meta-Analysis
Integrates microglial populations across datasets using scANVI and Harmony
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import scvi
import scanorama
import warnings
warnings.filterwarnings('ignore')

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')

class MicrogliaIntegrator:
    def __init__(self, processed_data_dir='data/processed', results_dir='results'):
        self.processed_data_dir = Path(processed_data_dir)
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(exist_ok=True)
        
    def load_microglial_datasets(self):
        """Load all microglial datasets"""
        print("Loading microglial datasets...")
        
        microglia_files = list(self.processed_data_dir.glob('*_microglia.h5ad'))
        
        if not microglia_files:
            print("No microglial datasets found. Run annotation first.")
            return None
            
        datasets = []
        dataset_names = []
        
        for file_path in microglia_files:
            dataset_name = file_path.stem.replace('_microglia', '')
            print(f"Loading {dataset_name}...")
            
            try:
                adata = sc.read_h5ad(file_path)
                adata.obs['dataset'] = dataset_name
                
                # Add species information based on dataset
                if dataset_name in ['GSE98969', 'GSE103334', 'GSE129788']:
                    adata.obs['species'] = 'mouse'
                else:
                    adata.obs['species'] = 'human'
                    
                # Add study metadata
                study_info = self.get_study_metadata(dataset_name)
                for key, value in study_info.items():
                    adata.obs[key] = value
                    
                datasets.append(adata)
                dataset_names.append(dataset_name)
                
                print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
                
            except Exception as e:
                print(f"Error loading {dataset_name}: {e}")
                continue
                
        print(f"\nLoaded {len(datasets)} datasets with microglial cells")
        return datasets, dataset_names
        
    def get_study_metadata(self, dataset_name):
        """Get metadata for each study"""
        metadata = {
            'GSE98969': {
                'study': 'Keren-Shaul_2017',
                'tissue': 'whole_brain',
                'condition': 'AD_model',
                'age_group': 'mixed'
            },
            'GSE103334': {
                'study': 'Mathys_2019', 
                'tissue': 'hippocampus',
                'condition': 'AD_model',
                'age_group': 'multiple_timepoints'
            },
            'GSE135437': {
                'study': 'Sankowski_2019',
                'tissue': 'cortex',
                'condition': 'AD_human',
                'age_group': 'elderly'
            },
            'GSE157827': {
                'study': 'Leng_2021',
                'tissue': 'prefrontal_cortex', 
                'condition': 'AD_human',
                'age_group': 'elderly'
            },
            'GSE129788': {
                'study': 'Aging_Mouse',
                'tissue': 'whole_brain',
                'condition': 'aging',
                'age_group': 'multiple'
            }
        }
        
        return metadata.get(dataset_name, {})
        
    def preprocess_for_integration(self, datasets):
        """Preprocess datasets for integration"""
        print("\nPreprocessing for integration...")
        
        processed_datasets = []
        
        for i, adata in enumerate(datasets):
            print(f"Processing dataset {i+1}/{len(datasets)}")
            
            # Ensure unique gene names
            adata.var_names_make_unique()
            
            # Basic filtering
            sc.pp.filter_genes(adata, min_cells=3)
            sc.pp.filter_cells(adata, min_genes=200)
            
            # Store raw counts
            adata.raw = adata
            
            # Normalize and log transform
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            
            processed_datasets.append(adata)
            
        return processed_datasets
        
    def concatenate_datasets(self, datasets, dataset_names):
        """Concatenate all datasets"""
        print("\nConcatenating datasets...")
        
        # Find common genes across all datasets
        common_genes = None
        for adata in datasets:
            if common_genes is None:
                common_genes = set(adata.var_names)
            else:
                common_genes = common_genes.intersection(set(adata.var_names))
                
        common_genes = sorted(list(common_genes))
        print(f"Found {len(common_genes)} common genes across datasets")
        
        if len(common_genes) < 100:
            print("Warning: Very few common genes found. This may indicate data format issues.")
        
        # Subset to common genes and clean data
        for adata in datasets:
            adata._inplace_subset_var(common_genes)
            
            # Clean up any remaining NaN or inf values
            if hasattr(adata.X, 'toarray'):
                X_data = adata.X.toarray()
            else:
                X_data = adata.X
            
            # Replace NaN and inf with 0
            import numpy as np
            X_data = np.where(np.isnan(X_data), 0, X_data)
            X_data = np.where(np.isinf(X_data), 0, X_data)
            
            # Update the data
            adata.X = X_data
            
        # Concatenate
        adata_concat = datasets[0].concatenate(
            datasets[1:], 
            batch_categories=dataset_names,
            index_unique='-'
        )
        
        # Final cleanup of concatenated data
        if hasattr(adata_concat.X, 'toarray'):
            X_data = adata_concat.X.toarray()
        else:
            X_data = adata_concat.X
            
        # Clean any remaining NaN/inf values
        import numpy as np
        X_data = np.where(np.isnan(X_data), 0, X_data)
        X_data = np.where(np.isinf(X_data), 0, X_data)
        adata_concat.X = X_data
        
        print(f"Combined dataset: {adata_concat.n_obs} cells, {adata_concat.n_vars} genes")
        
        # Final validation
        nan_count = np.isnan(adata_concat.X).sum()
        inf_count = np.isinf(adata_concat.X).sum()
        if nan_count > 0 or inf_count > 0:
            print(f"Warning: Still have {nan_count} NaN and {inf_count} inf values")
        else:
            print("Data cleaning successful - no NaN or inf values remaining")
            
        return adata_concat
        
    def integrate_with_scanorama(self, adata_concat):
        """Integrate datasets using Scanorama"""
        print("\nRunning Scanorama integration...")
        
        try:
            # Prepare data for Scanorama
            datasets = []
            genes_list = []
            dataset_names = adata_concat.obs['batch'].unique()
            
            for batch in dataset_names:
                batch_data = adata_concat[adata_concat.obs['batch'] == batch]
                if hasattr(batch_data.X, 'toarray'):
                    X_data = batch_data.X.toarray()
                else:
                    X_data = batch_data.X
                
                # Ensure no NaN values
                import numpy as np
                X_data = np.where(np.isnan(X_data), 0, X_data)
                datasets.append(X_data)
                genes_list.append(list(batch_data.var_names))
                
            # Run integration
            integrated, genes = scanorama.integrate(datasets, genes_list)
            
            # Create integrated AnnData object
            adata_integrated = adata_concat.copy()
            adata_integrated.X = np.concatenate(integrated)
            
            # Compute embeddings
            sc.pp.highly_variable_genes(adata_integrated, min_mean=0.0125, max_mean=3, min_disp=0.5)
            sc.pp.scale(adata_integrated, max_value=10)
            sc.tl.pca(adata_integrated, svd_solver='arpack', n_comps=50)
            sc.pp.neighbors(adata_integrated, n_neighbors=15, n_pcs=50)
            sc.tl.umap(adata_integrated)
            
            return adata_integrated
            
        except Exception as e:
            print(f"Scanorama integration failed: {e}")
            return None
            
    def integrate_with_scvi(self, adata_concat):
        """Integrate datasets using scVI"""
        print("\nRunning scVI integration...")
        
        try:
            # Setup scVI model
            scvi.model.SCVI.setup_anndata(
                adata_concat,
                layer=None,
                batch_key='batch'
            )
            
            # Train model
            model = scvi.model.SCVI(adata_concat, n_latent=30, n_hidden=128)
            model.train(max_epochs=200, early_stopping=True, plan_kwargs={'lr':1e-3})
            
            # Get latent representation
            adata_integrated = adata_concat.copy()
            adata_integrated.obsm['X_scVI'] = model.get_latent_representation()
            
            # Compute neighborhood graph and UMAP
            sc.pp.neighbors(adata_integrated, use_rep='X_scVI', n_neighbors=15)
            sc.tl.umap(adata_integrated)
            
            return adata_integrated, model
            
        except Exception as e:
            print(f"scVI integration failed: {e}")
            return None, None
            
    def evaluate_integration(self, adata_integrated, method_name):
        """Evaluate integration quality"""
        print(f"\nEvaluating {method_name} integration...")
        
        # Calculate mixing metrics
        try:
            # Silhouette score for batch mixing
            from sklearn.metrics import silhouette_score
            
            if 'X_scVI' in adata_integrated.obsm:
                embedding = adata_integrated.obsm['X_scVI']
            else:
                embedding = adata_integrated.obsm['X_pca']
                
            # Batch silhouette (lower is better for mixing)
            batch_sil = silhouette_score(embedding, adata_integrated.obs['batch'])
            
            # Cell type silhouette (higher is better for preservation)  
            if 'celltypist_Immune_All_High.pkl' in adata_integrated.obs.columns:
                celltype_sil = silhouette_score(
                    embedding, 
                    adata_integrated.obs['celltypist_Immune_All_High.pkl']
                )
            else:
                celltype_sil = np.nan
                
            metrics = {
                'method': method_name,
                'batch_silhouette': batch_sil,
                'celltype_silhouette': celltype_sil,
                'n_cells': adata_integrated.n_obs,
                'n_genes': adata_integrated.n_vars,
                'n_batches': len(adata_integrated.obs['batch'].unique())
            }
            
            print(f"Batch silhouette: {batch_sil:.3f}")
            print(f"Cell type silhouette: {celltype_sil:.3f}")
            
            return metrics
            
        except Exception as e:
            print(f"Error evaluating integration: {e}")
            return None
            
    def generate_integration_plots(self, adata_integrated, method_name):
        """Generate integration visualization plots"""
        print(f"\nGenerating {method_name} visualization plots...")
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle(f'Integration Results - {method_name}', fontsize=16)
        
        # UMAP by batch
        sc.pl.umap(adata_integrated, color='batch', ax=axes[0,0], 
                  show=False, frameon=False, legend_fontsize=8)
        axes[0,0].set_title('Datasets')
        
        # UMAP by species
        sc.pl.umap(adata_integrated, color='species', ax=axes[0,1],
                  show=False, frameon=False)  
        axes[0,1].set_title('Species')
        
        # UMAP by study condition
        sc.pl.umap(adata_integrated, color='condition', ax=axes[0,2],
                  show=False, frameon=False, legend_fontsize=8)
        axes[0,2].set_title('Condition')
        
        # Cell type annotations
        if 'celltypist_Immune_All_High.pkl' in adata_integrated.obs.columns:
            sc.pl.umap(adata_integrated, color='celltypist_Immune_All_High.pkl', 
                      ax=axes[1,0], show=False, frameon=False, legend_fontsize=6)
            axes[1,0].set_title('Cell Types')
            
        # Dataset composition
        dataset_counts = adata_integrated.obs['batch'].value_counts()
        dataset_counts.plot(kind='bar', ax=axes[1,1])
        axes[1,1].set_title('Dataset Composition')
        axes[1,1].tick_params(axis='x', rotation=45)
        
        # Species composition by dataset
        comp_df = pd.crosstab(adata_integrated.obs['batch'], 
                             adata_integrated.obs['species'])
        comp_df.plot(kind='bar', stacked=True, ax=axes[1,2])
        axes[1,2].set_title('Species by Dataset')  
        axes[1,2].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(self.results_dir / f'figures/{method_name}_integration.pdf',
                   dpi=300, bbox_inches='tight')
        plt.close()

def main():
    """Main integration pipeline"""
    integrator = MicrogliaIntegrator()
    
    # Load microglial datasets
    result = integrator.load_microglial_datasets()
    
    if result is None:
        print("No microglial datasets found. Please run cell annotation first.")
        return
    
    datasets, dataset_names = result
        
    # Preprocess for integration
    datasets = integrator.preprocess_for_integration(datasets)
    
    # Concatenate datasets
    adata_concat = integrator.concatenate_datasets(datasets, dataset_names)
    
    # Save concatenated data
    adata_concat.write(integrator.processed_data_dir / 'microglia_concatenated.h5ad')
    
    integration_results = {}
    
    # Try multiple integration methods
    methods = ['scanorama', 'scvi']
    
    for method in methods:
        print(f"\n{'='*50}")
        print(f"INTEGRATION METHOD: {method.upper()}")
        print(f"{'='*50}")
        
        try:
            if method == 'scanorama':
                adata_integrated = integrator.integrate_with_scanorama(adata_concat)
                model = None
            elif method == 'scvi':
                adata_integrated, model = integrator.integrate_with_scvi(adata_concat)
                
            if adata_integrated is not None:
                # Evaluate integration
                metrics = integrator.evaluate_integration(adata_integrated, method)
                integration_results[method] = metrics
                
                # Generate plots
                integrator.generate_integration_plots(adata_integrated, method)
                
                # Save integrated data
                output_path = integrator.processed_data_dir / f'microglia_integrated_{method}.h5ad'
                adata_integrated.write(output_path)
                print(f"Integrated data saved to {output_path}")
                
                # Save model if applicable
                if model is not None:
                    model.save(integrator.processed_data_dir / f'scvi_model_{method}')
                    
        except Exception as e:
            print(f"Integration with {method} failed: {e}")
            continue
            
    # Save integration summary
    import json
    with open(integrator.results_dir / 'integration_summary.json', 'w') as f:
        json.dump(integration_results, f, indent=2, default=str)
        
    print(f"\n{'='*50}")
    print("INTEGRATION COMPLETE")
    print(f"{'='*50}")
    print("Next step: Run activation state analysis")

if __name__ == "__main__":
    main()
