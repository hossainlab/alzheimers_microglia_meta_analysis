#!/usr/bin/env python3
"""
Activation States Analysis Script for Alzheimer's Microglia Meta-Analysis
Identifies and characterizes conserved microglial activation states
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Try importing additional analysis tools
try:
    import decoupler as dc
    DECOUPLER_AVAILABLE = True
except ImportError:
    print("Warning: decoupler not available. Install with: pip install decoupler")
    DECOUPLER_AVAILABLE = False

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')

class ActivationStatesAnalyzer:
    def __init__(self, processed_data_dir='data/processed', results_dir='results'):
        self.processed_data_dir = Path(processed_data_dir)
        self.results_dir = Path(results_dir)
        self.tables_dir = self.results_dir / 'tables' / 'gene_signatures'
        self.figures_dir = self.results_dir / 'figures' / 'activation_states'
        
        # Create directories
        self.tables_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        
    def load_integrated_data(self):
        """Load integrated microglia data"""
        print("Loading integrated microglia data...")
        
        # Try to find integrated data - prefer scvi then scanorama
        for method in ['scvi', 'scanorama']:
            integrated_file = self.processed_data_dir / f'microglia_integrated_{method}.h5ad'
            if integrated_file.exists():
                print(f"Loading {method} integrated data...")
                adata = sc.read_h5ad(integrated_file)
                print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
                return adata, method
                
        # Fallback to concatenated data
        concat_file = self.processed_data_dir / 'microglia_concatenated.h5ad'
        if concat_file.exists():
            print("Loading concatenated data (no integration found)...")
            adata = sc.read_h5ad(concat_file)
            return adata, 'concatenated'
            
        print("No integrated data found. Run integration first.")
        return None, None
    
    def identify_activation_states(self, adata):
        """Identify microglial activation states using clustering"""
        print("\nIdentifying activation states...")
        
        # Use appropriate representation for clustering
        if 'X_scVI' in adata.obsm:
            use_rep = 'X_scVI'
            print("Using scVI representation for clustering")
        elif 'X_pca' in adata.obsm:
            use_rep = 'X_pca'
            print("Using PCA representation for clustering")
        else:
            # Clean data before computing PCA
            import numpy as np
            
            # Clean expression data
            if hasattr(adata.X, 'toarray'):
                X_data = adata.X.toarray()
            else:
                X_data = adata.X.copy()
            
            # Replace NaN and inf with 0
            X_data = np.where(np.isnan(X_data), 0, X_data)
            X_data = np.where(np.isinf(X_data), 0, X_data)
            adata.X = X_data
            
            # Clean obs data
            for col in adata.obs.columns:
                if adata.obs[col].dtype in ['float64', 'float32']:
                    adata.obs[col] = adata.obs[col].fillna(0)
            
            print("Cleaned data before PCA computation")
            
            # Compute PCA if not available
            try:
                sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
                sc.pp.scale(adata, max_value=10)
                sc.tl.pca(adata, svd_solver='arpack')
                use_rep = 'X_pca'
                print("Computed PCA for clustering")
            except Exception as e:
                print(f"PCA computation failed: {e}")
                print("Using scaled expression data for clustering")
                sc.pp.scale(adata, max_value=10)
                use_rep = None
        
        # Test multiple resolutions
        resolutions = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0]
        clustering_metrics = []
        
        for resolution in resolutions:
            # Compute neighbors if not done
            if 'neighbors' not in adata.uns:
                sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15)
            
            # Leiden clustering
            sc.tl.leiden(adata, resolution=resolution, key_added=f'leiden_{resolution}')
            
            # Calculate silhouette score
            from sklearn.metrics import silhouette_score
            
            if use_rep == 'X_scVI':
                embedding = adata.obsm['X_scVI']
            else:
                embedding = adata.obsm['X_pca']
                
            sil_score = silhouette_score(embedding, adata.obs[f'leiden_{resolution}'])
            n_clusters = len(adata.obs[f'leiden_{resolution}'].unique())
            
            clustering_metrics.append({
                'resolution': resolution,
                'n_clusters': n_clusters,
                'silhouette_score': sil_score
            })
            
            print(f"Resolution {resolution}: {n_clusters} clusters, silhouette = {sil_score:.3f}")
        
        # Select optimal resolution
        metrics_df = pd.DataFrame(clustering_metrics)
        optimal_idx = metrics_df['silhouette_score'].idxmax()
        optimal_resolution = metrics_df.loc[optimal_idx, 'resolution']
        
        print(f"\nOptimal resolution: {optimal_resolution}")
        adata.obs['activation_state'] = adata.obs[f'leiden_{optimal_resolution}']
        
        # Save clustering metrics
        metrics_df.to_csv(self.tables_dir / 'clustering_metrics.csv', index=False)
        
        return adata, optimal_resolution
    
    def characterize_activation_states(self, adata):
        """Characterize activation states through differential expression"""
        print("\nCharacterizing activation states...")
        
        # Differential expression between states
        sc.tl.rank_genes_groups(
            adata, 
            'activation_state', 
            method='wilcoxon',
            use_raw=True,
            n_genes=100
        )
        
        # Get results as DataFrame
        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        
        de_results = []
        for group in groups:
            group_df = pd.DataFrame({
                'gene': result['names'][group],
                'pvalue': result['pvals'][group],
                'pvalue_adj': result['pvals_adj'][group], 
                'logfoldchange': result['logfoldchanges'][group],
                'activation_state': group
            })
            de_results.append(group_df)
        
        de_df = pd.concat(de_results, ignore_index=True)
        
        # Filter for significant genes
        sig_genes = de_df[
            (de_df['pvalue_adj'] < 0.05) & 
            (de_df['logfoldchange'] > 0.5)
        ]
        
        # Save results
        de_df.to_csv(self.tables_dir / 'activation_state_markers.csv', index=False)
        sig_genes.to_csv(self.tables_dir / 'significant_markers.csv', index=False)
        
        print(f"Found {len(sig_genes)} significant marker genes")
        
        return de_df, sig_genes
    
    def pathway_enrichment_analysis(self, adata, sig_genes):
        """Perform pathway enrichment analysis for activation states"""
        print("\nPerforming pathway enrichment analysis...")
        
        if not DECOUPLER_AVAILABLE:
            print("Decoupler not available, skipping pathway analysis")
            return None
        
        try:
            # Get MSigDB gene sets
            msigdb = dc.get_resource('MSigDB')
            
            # Filter for relevant collections
            relevant_collections = ['hallmark', 'kegg', 'reactome', 'go_bp']
            msigdb = msigdb[msigdb['collection'].isin(relevant_collections)]
            
            # Prepare gene expression matrix
            X = adata.to_df()
            
            # Run enrichment analysis for each activation state
            enrichment_results = []
            
            for state in adata.obs['activation_state'].unique():
                print(f"Analyzing pathways for activation state {state}...")
                
                # Get marker genes for this state
                state_markers = sig_genes[sig_genes['activation_state'] == state]['gene'].tolist()
                
                if len(state_markers) < 5:
                    print(f"Skipping state {state} - too few markers")
                    continue
                
                # Subset to state markers
                X_state = X[state_markers]
                
                # Run Over Representation Analysis (ORA)
                ora_results = dc.get_ora_df(
                    df=X_state.mean().to_frame('score'),
                    net=msigdb,
                    source='source',
                    target='target'
                )
                
                ora_results['activation_state'] = state
                enrichment_results.append(ora_results)
            
            if enrichment_results:
                enrichment_df = pd.concat(enrichment_results, ignore_index=True)
                
                # Filter significant pathways
                sig_pathways = enrichment_df[enrichment_df['FDR'] < 0.05]
                
                # Save results
                enrichment_df.to_csv(self.tables_dir / 'pathway_enrichment.csv', index=False)
                sig_pathways.to_csv(self.tables_dir / 'significant_pathways.csv', index=False)
                
                print(f"Found {len(sig_pathways)} significant pathways")
                return enrichment_df
            
        except Exception as e:
            print(f"Pathway analysis failed: {e}")
            return None
    
    def calculate_conservation_scores(self, adata, de_df):
        """Calculate conservation scores across datasets"""
        print("\nCalculating conservation scores...")
        
        conservation_results = []
        
        for state in adata.obs['activation_state'].unique():
            state_markers = de_df[
                (de_df['activation_state'] == state) &
                (de_df['pvalue_adj'] < 0.05) &
                (de_df['logfoldchange'] > 0.5)
            ]['gene'].tolist()
            
            if len(state_markers) < 10:
                continue
                
            # Calculate expression consistency across datasets
            dataset_consistency = []
            
            for dataset in adata.obs['batch'].unique():
                dataset_mask = adata.obs['batch'] == dataset
                state_mask = adata.obs['activation_state'] == state
                combined_mask = dataset_mask & state_mask
                
                if combined_mask.sum() < 5:  # Need at least 5 cells
                    continue
                
                # Mean expression in this dataset/state combination
                dataset_state_expr = adata[combined_mask, state_markers].X.mean(axis=0)
                dataset_other_expr = adata[dataset_mask & ~state_mask, state_markers].X.mean(axis=0)
                
                # Calculate fold change
                fold_change = np.array(dataset_state_expr).flatten() / (np.array(dataset_other_expr).flatten() + 1e-8)
                mean_fold_change = np.mean(np.log2(fold_change + 1e-8))
                
                dataset_consistency.append({
                    'activation_state': state,
                    'dataset': dataset,
                    'mean_log2fc': mean_fold_change,
                    'n_cells': combined_mask.sum(),
                    'n_markers': len(state_markers)
                })
            
            conservation_results.extend(dataset_consistency)
        
        conservation_df = pd.DataFrame(conservation_results)
        
        # Calculate overall conservation score per state
        if not conservation_df.empty:
            conservation_summary = conservation_df.groupby('activation_state').agg({
                'mean_log2fc': ['mean', 'std'],
                'dataset': 'count',
                'n_cells': 'sum'
            }).round(3)
            
            conservation_summary.columns = ['mean_log2fc', 'std_log2fc', 'n_datasets', 'total_cells']
            conservation_summary['conservation_score'] = conservation_summary['mean_log2fc'] / (conservation_summary['std_log2fc'] + 0.1)
        else:
            # Create empty summary if no conservation data
            conservation_summary = pd.DataFrame(columns=['mean_log2fc', 'std_log2fc', 'n_datasets', 'total_cells', 'conservation_score'])
            print("No conservation data available - insufficient marker genes")
        
        # Save results
        conservation_df.to_csv(self.tables_dir / 'conservation_analysis.csv', index=False)
        conservation_summary.to_csv(self.tables_dir / 'conservation_summary.csv')
        
        print(f"Conservation analysis complete for {len(conservation_summary)} states")
        
        return conservation_df, conservation_summary
    
    def create_comprehensive_plots(self, adata, de_df, conservation_df):
        """Create comprehensive visualization plots"""
        print("\nCreating comprehensive plots...")
        
        try:
            # Create a simple but effective plot
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            fig.suptitle('Microglia Activation States Analysis', fontsize=16, fontweight='bold')
            
            # Plot 1: Cell counts per state
            if 'activation_state' in adata.obs.columns:
                state_counts = adata.obs['activation_state'].value_counts().sort_index()
                state_counts.plot(kind='bar', ax=axes[0,0], color='skyblue')
                axes[0,0].set_title('Cells per Activation State')
                axes[0,0].set_xlabel('Activation State')
                axes[0,0].set_ylabel('Number of Cells')
                axes[0,0].tick_params(axis='x', rotation=45)
            
            # Plot 2: Dataset composition
            if 'batch' in adata.obs.columns:
                dataset_counts = adata.obs['batch'].value_counts()
                dataset_counts.plot(kind='pie', ax=axes[0,1], autopct='%1.1f%%')
                axes[0,1].set_title('Dataset Composition')
                axes[0,1].set_ylabel('')
            
            # Plot 3: State composition by dataset
            if 'activation_state' in adata.obs.columns and 'batch' in adata.obs.columns:
                comp_df = pd.crosstab(adata.obs['batch'], adata.obs['activation_state'], normalize='index')
                sns.heatmap(comp_df, ax=axes[1,0], cmap='viridis', annot=True, fmt='.2f')
                axes[1,0].set_title('State Composition by Dataset')
                axes[1,0].set_xlabel('Activation State')
                axes[1,0].set_ylabel('Dataset')
            
            # Plot 4: Marker gene counts
            if not de_df.empty:
                marker_counts = de_df[de_df['pvalue_adj'] < 0.05].groupby('activation_state').size()
                if len(marker_counts) > 0:
                    marker_counts.plot(kind='bar', ax=axes[1,1], color='green')
                    axes[1,1].set_title('Significant Markers per State')
                    axes[1,1].set_xlabel('Activation State')
                    axes[1,1].set_ylabel('Number of Markers')
                    axes[1,1].tick_params(axis='x', rotation=45)
                else:
                    axes[1,1].text(0.5, 0.5, 'No significant markers found', ha='center', va='center')
            else:
                axes[1,1].text(0.5, 0.5, 'No differential expression results', ha='center', va='center')
            
            plt.tight_layout()
            plt.savefig(self.figures_dir / 'comprehensive_analysis.pdf', 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
            print(f"Comprehensive plot saved to {self.figures_dir / 'comprehensive_analysis.pdf'}")
            
        except Exception as e:
            print(f"Error creating comprehensive plots: {e}")
            print("Creating basic fallback plot...")
            
            # Create a simple fallback plot
            fig, ax = plt.subplots(1, 1, figsize=(8, 6))
            
            try:
                if 'activation_state' in adata.obs.columns:
                    state_counts = adata.obs['activation_state'].value_counts()
                    state_counts.plot(kind='bar', ax=ax)
                    ax.set_title('Microglia Activation States')
                    ax.set_xlabel('State')
                    ax.set_ylabel('Number of Cells')
                else:
                    ax.text(0.5, 0.5, 'No activation states found', ha='center', va='center')
                    ax.set_title('Analysis Results')
                    
                plt.tight_layout()
                plt.savefig(self.figures_dir / 'comprehensive_analysis.pdf', 
                           dpi=300, bbox_inches='tight')
                plt.close()
                print("Fallback plot saved")
                
            except Exception as e2:
                print(f"Fallback plot also failed: {e2}")
                plt.close('all')  # Close any remaining figures
    
    def generate_summary_report(self, adata, de_df, conservation_df):
        """Generate summary report"""
        print("\nGenerating summary report...")
        
        # Calculate summary statistics
        n_states = len(adata.obs['activation_state'].unique())
        n_cells = adata.n_obs
        n_datasets = len(adata.obs['batch'].unique())
        n_markers = len(de_df[de_df['pvalue_adj'] < 0.05])
        
        report = f"""# Microglial Activation States Analysis Report

## Overview
- **Total Cells**: {n_cells:,}
- **Activation States Identified**: {n_states}
- **Datasets Integrated**: {n_datasets}
- **Significant Marker Genes**: {n_markers}

## Activation States Summary
"""
        
        # Add state-specific information
        for state in sorted(adata.obs['activation_state'].unique()):
            state_cells = (adata.obs['activation_state'] == state).sum()
            state_markers = de_df[
                (de_df['activation_state'] == state) & 
                (de_df['pvalue_adj'] < 0.05)
            ]
            
            report += f"""
### State {state}
- **Cells**: {state_cells} ({state_cells/n_cells*100:.1f}%)
- **Marker Genes**: {len(state_markers)}
- **Top Markers**: {', '.join(state_markers.head(5)['gene'].tolist())}
"""
        
        # Conservation analysis
        if not conservation_df.empty:
            conservation_summary = conservation_df.groupby('activation_state')['mean_log2fc'].mean().sort_values(ascending=False)
            report += f"""
## Conservation Analysis
Most conserved states (by expression consistency):
"""
            for state, score in conservation_summary.head(3).items():
                report += f"- **State {state}**: {score:.3f}\n"
        
        report += f"""
## Files Generated
- **Marker genes**: `{self.tables_dir / 'activation_state_markers.csv'}`
- **Conservation analysis**: `{self.tables_dir / 'conservation_summary.csv'}`
- **Comprehensive plot**: `{self.figures_dir / 'comprehensive_analysis.pdf'}`

## Next Steps
1. Validate activation states with known microglia biology
2. Perform meta-analysis across studies
3. Compare with published microglia states
"""
        
        # Save report
        report_path = self.results_dir / 'reports' / 'activation_states_summary.md'
        report_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(report_path, 'w') as f:
            f.write(report)
        
        print(f"Report saved to {report_path}")
    
    def run_analysis(self):
        """Run complete activation states analysis"""
        print("Starting activation states analysis...")
        
        # Load integrated data
        adata, method = self.load_integrated_data()
        if adata is None:
            return
        
        print(f"Using {method} integrated data")
        
        # Clean data at the start
        print("Cleaning data...")
        import numpy as np
        
        # Clean expression data
        if hasattr(adata.X, 'toarray'):
            X_data = adata.X.toarray()
        else:
            X_data = adata.X.copy()
        
        # Replace NaN and inf with 0
        nan_count = np.isnan(X_data).sum()
        inf_count = np.isinf(X_data).sum()
        print(f"Found {nan_count} NaN and {inf_count} inf values - replacing with 0")
        
        X_data = np.where(np.isnan(X_data), 0, X_data)
        X_data = np.where(np.isinf(X_data), 0, X_data)
        adata.X = X_data
        
        # Clean obs data
        for col in adata.obs.columns:
            if adata.obs[col].dtype in ['float64', 'float32']:
                nan_count = adata.obs[col].isna().sum()
                if nan_count > 0:
                    print(f"Cleaning {nan_count} NaN values in obs column '{col}'")
                    adata.obs[col] = adata.obs[col].fillna(0)
        
        print("Data cleaning completed")
        
        # Identify activation states
        adata, optimal_resolution = self.identify_activation_states(adata)
        
        # Characterize states
        de_df, sig_genes = self.characterize_activation_states(adata)
        
        # Pathway enrichment
        enrichment_df = self.pathway_enrichment_analysis(adata, sig_genes)
        
        # Conservation analysis
        conservation_df, conservation_summary = self.calculate_conservation_scores(adata, de_df)
        
        # Create comprehensive plots
        self.create_comprehensive_plots(adata, de_df, conservation_df)
        
        # Generate report
        self.generate_summary_report(adata, de_df, conservation_df)
        
        # Save analyzed data
        output_path = self.processed_data_dir / 'microglia_activation_states.h5ad'
        adata.write(output_path)
        
        print(f"\nActivation states analysis complete!")
        print(f"Results saved in: {self.results_dir}")
        print(f"Analyzed data: {output_path}")

def main():
    analyzer = ActivationStatesAnalyzer()
    analyzer.run_analysis()

if __name__ == "__main__":
    main()