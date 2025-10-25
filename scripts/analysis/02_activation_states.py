#!/usr/bin/env python3
"""
Microglial Activation State Analysis for Alzheimer's Meta-Analysis
Identifies conserved activation states and gene signatures
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import decoupler as dc
from sklearn.metrics import adjusted_rand_score, silhouette_score
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=300, facecolor='white')

class ActivationStateAnalyzer:
    def __init__(self, processed_data_dir='data/processed', results_dir='results'):
        self.processed_data_dir = Path(processed_data_dir)
        self.results_dir = Path(results_dir)
        
        # Create results subdirectories
        (self.results_dir / 'figures' / 'activation_states').mkdir(parents=True, exist_ok=True)
        (self.results_dir / 'tables' / 'gene_signatures').mkdir(parents=True, exist_ok=True)
        
    def load_integrated_data(self, method='scvi'):
        """Load integrated microglial data"""
        print(f"Loading integrated data (method: {method})...")
        
        try:
            file_path = self.processed_data_dir / f'microglia_integrated_{method}.h5ad'
            adata = sc.read_h5ad(file_path)
            print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
            return adata
        except FileNotFoundError:
            print(f"Integrated data not found for {method}. Try alternative method.")
            return None
            
    def identify_activation_states(self, adata):
        """Identify microglial activation states through clustering"""
        print("\nIdentifying microglial activation states...")
        
        # Perform clustering at multiple resolutions
        resolutions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
        
        for res in resolutions:
            sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
            
        # Choose optimal resolution based on silhouette score
        best_resolution = None
        best_score = -1
        
        for res in resolutions:
            cluster_key = f'leiden_{res}'
            n_clusters = len(adata.obs[cluster_key].unique())
            
            if n_clusters > 1 and n_clusters < adata.n_obs // 10:
                if 'X_scVI' in adata.obsm:
                    embedding = adata.obsm['X_scVI']
                else:
                    embedding = adata.obsm['X_pca']
                    
                sil_score = silhouette_score(embedding, adata.obs[cluster_key])
                
                print(f"Resolution {res}: {n_clusters} clusters, silhouette: {sil_score:.3f}")
                
                if sil_score > best_score:
                    best_score = sil_score
                    best_resolution = res
                    
        if best_resolution is not None:
            print(f"\nSelected resolution {best_resolution} (silhouette: {best_score:.3f})")
            adata.obs['activation_state'] = adata.obs[f'leiden_{best_resolution}']
        else:
            print("Using default resolution 0.4")
            adata.obs['activation_state'] = adata.obs['leiden_0.4']
            
        return adata
        
    def analyze_differential_expression(self, adata):
        """Find differentially expressed genes between activation states"""
        print("\nAnalyzing differential gene expression...")
        
        # Calculate differential expression between activation states
        sc.tl.rank_genes_groups(
            adata, 
            groupby='activation_state',
            method='wilcoxon',
            key_added='rank_genes_activation',
            n_genes=100
        )
        
        # Extract results as DataFrame
        de_results = sc.get.rank_genes_groups_df(adata, key='rank_genes_activation', group=None)
        
        # Save results
        de_results.to_csv(
            self.results_dir / 'tables/gene_signatures/activation_state_markers.csv',
            index=False
        )
        
        # Also analyze by dataset to check consistency
        sc.tl.rank_genes_groups(
            adata,
            groupby='activation_state', 
            method='wilcoxon',
            key_added='rank_genes_by_dataset'
        )
        
        return de_results
        
    def perform_pathway_enrichment(self, adata, de_results):
        """Perform pathway enrichment analysis using decoupler"""
        print("\nPerforming pathway enrichment analysis...")
        
        try:
            # Load pathway databases
            msigdb = dc.get_resource('MSigDB')
            
            # Filter pathways by size
            msigdb = msigdb[(msigdb.groupby('source').source.transform('size') >= 15) & 
                           (msigdb.groupby('source').source.transform('size') <= 500)]
            
            # Prepare data for GSEA
            activation_states = adata.obs['activation_state'].unique()
            
            enrichment_results = {}
            
            for state in activation_states:
                print(f"Analyzing enrichment for state {state}...")
                
                # Get marker genes for this state
                state_markers = de_results[de_results['group'] == state].copy()
                
                if len(state_markers) > 0:
                    # Create gene ranking
                    gene_stats = state_markers.set_index('names')['scores']
                    gene_stats = gene_stats.sort_values(ascending=False)
                    
                    # Run GSEA
                    gsea_res = dc.run_gsea(
                        mat=gene_stats.to_frame().T,
                        net=msigdb,
                        source='source',
                        target='target',
                        verbose=False
                    )
                    
                    enrichment_results[state] = gsea_res
                    
            # Save enrichment results
            for state, results in enrichment_results.items():
                results.to_csv(
                    self.results_dir / f'tables/gene_signatures/enrichment_state_{state}.csv'
                )
                
            return enrichment_results
            
        except Exception as e:
            print(f"Pathway enrichment failed: {e}")
            return None
            
    def identify_conserved_signatures(self, adata):
        """Identify gene signatures conserved across datasets"""
        print("\nIdentifying conserved activation signatures...")
        
        # Get marker genes for each activation state in each dataset
        datasets = adata.obs['batch'].unique()
        states = adata.obs['activation_state'].unique()
        
        conserved_signatures = {}
        
        for state in states:
            print(f"Analyzing conservation for state {state}...")
            
            dataset_markers = {}
            
            # Get top markers for this state in each dataset
            for dataset in datasets:
                subset = adata[(adata.obs['batch'] == dataset) & 
                              (adata.obs['activation_state'] == state)]
                
                if subset.n_obs > 10:  # Minimum cells required
                    # Compare this state vs others within dataset
                    temp_adata = adata[adata.obs['batch'] == dataset].copy()
                    
                    # Create binary comparison (this state vs others)
                    temp_adata.obs['is_target_state'] = (
                        temp_adata.obs['activation_state'] == state
                    )
                    
                    if temp_adata.obs['is_target_state'].sum() > 5:
                        sc.tl.rank_genes_groups(
                            temp_adata,
                            groupby='is_target_state',
                            method='wilcoxon',
                            reference='rest'
                        )
                        
                        # Extract top genes
                        markers_df = sc.get.rank_genes_groups_df(
                            temp_adata, 
                            group=True,
                            key='rank_genes_groups'
                        )
                        
                        if len(markers_df) > 0:
                            top_markers = markers_df[
                                (markers_df['pvals_adj'] < 0.01) & 
                                (markers_df['logfoldchanges'] > 0.5)
                            ]['names'].head(50).tolist()
                            
                            dataset_markers[dataset] = set(top_markers)
                            
            # Find genes present in multiple datasets
            if len(dataset_markers) >= 2:
                all_genes = set.union(*dataset_markers.values()) if dataset_markers else set()
                conserved_genes = set.intersection(*dataset_markers.values()) if dataset_markers else set()
                
                # Calculate conservation score for each gene
                gene_conservation = {}
                for gene in all_genes:
                    presence_count = sum(1 for markers in dataset_markers.values() 
                                       if gene in markers)
                    conservation_score = presence_count / len(dataset_markers)
                    gene_conservation[gene] = conservation_score
                    
                # Filter for genes present in at least half of datasets
                min_presence = max(2, len(dataset_markers) // 2)
                conserved_genes = {
                    gene: score for gene, score in gene_conservation.items()
                    if score >= min_presence / len(dataset_markers)
                }
                
                conserved_signatures[state] = {
                    'genes': conserved_genes,
                    'n_datasets': len(dataset_markers),
                    'total_genes': len(all_genes)
                }
                
                print(f"  State {state}: {len(conserved_genes)} conserved genes "
                      f"across {len(dataset_markers)} datasets")
                      
        # Save conserved signatures
        conserved_df_list = []
        for state, info in conserved_signatures.items():
            for gene, score in info['genes'].items():
                conserved_df_list.append({
                    'activation_state': state,
                    'gene': gene,
                    'conservation_score': score,
                    'n_datasets': info['n_datasets']
                })
                
        if conserved_df_list:
            conserved_df = pd.DataFrame(conserved_df_list)
            conserved_df.to_csv(
                self.results_dir / 'tables/gene_signatures/conserved_signatures.csv',
                index=False
            )
            
        return conserved_signatures
        
    def validate_across_species(self, adata):
        """Validate conservation across human and mouse studies"""
        print("\nValidating activation states across species...")
        
        # Get activation states in human vs mouse
        human_data = adata[adata.obs['species'] == 'human']
        mouse_data = adata[adata.obs['species'] == 'mouse']
        
        validation_results = {}
        
        if human_data.n_obs > 0 and mouse_data.n_obs > 0:
            # Compare activation state distributions
            human_states = human_data.obs['activation_state'].value_counts()
            mouse_states = mouse_data.obs['activation_state'].value_counts()
            
            # Create comparison DataFrame
            comparison_df = pd.DataFrame({
                'human': human_states,
                'mouse': mouse_states
            }).fillna(0)
            
            # Calculate correlation
            if len(comparison_df) > 1:
                correlation = comparison_df['human'].corr(comparison_df['mouse'])
                validation_results['state_distribution_correlation'] = correlation
                
            # Save species comparison
            comparison_df.to_csv(
                self.results_dir / 'tables/species_activation_comparison.csv'
            )
            
            print(f"Species state correlation: {correlation:.3f}")
            
        return validation_results
        
    def generate_activation_plots(self, adata):
        """Generate comprehensive activation state visualization"""
        print("\nGenerating activation state plots...")
        
        # Create figure with multiple subplots
        fig, axes = plt.subplots(3, 3, figsize=(20, 18))
        fig.suptitle('Microglial Activation States Analysis', fontsize=20)
        
        # 1. UMAP colored by activation states
        sc.pl.umap(adata, color='activation_state', ax=axes[0,0], 
                  show=False, frameon=False, legend_fontsize=10)
        axes[0,0].set_title('Activation States (UMAP)')
        
        # 2. UMAP colored by dataset
        sc.pl.umap(adata, color='batch', ax=axes[0,1],
                  show=False, frameon=False, legend_fontsize=8)
        axes[0,1].set_title('Datasets')
        
        # 3. UMAP colored by species
        sc.pl.umap(adata, color='species', ax=axes[0,2],
                  show=False, frameon=False)
        axes[0,2].set_title('Species')
        
        # 4. Activation state composition by dataset
        comp_df = pd.crosstab(adata.obs['batch'], adata.obs['activation_state'])
        comp_df_norm = comp_df.div(comp_df.sum(axis=1), axis=0)
        
        sns.heatmap(comp_df_norm, annot=True, fmt='.2f', ax=axes[1,0], cmap='viridis')
        axes[1,0].set_title('State Composition by Dataset')
        axes[1,0].set_xlabel('Activation State')
        
        # 5. State composition by species
        species_comp = pd.crosstab(adata.obs['species'], adata.obs['activation_state'])
        species_comp_norm = species_comp.div(species_comp.sum(axis=1), axis=0)
        
        sns.heatmap(species_comp_norm, annot=True, fmt='.2f', ax=axes[1,1], cmap='plasma')
        axes[1,1].set_title('State Composition by Species')
        
        # 6. Cell counts per state
        state_counts = adata.obs['activation_state'].value_counts()
        state_counts.plot(kind='bar', ax=axes[1,2])
        axes[1,2].set_title('Cells per Activation State')
        axes[1,2].tick_params(axis='x', rotation=45)
        
        # 7. State distribution across conditions
        if 'condition' in adata.obs.columns:
            condition_comp = pd.crosstab(adata.obs['condition'], adata.obs['activation_state'])
            condition_comp_norm = condition_comp.div(condition_comp.sum(axis=1), axis=0)
            
            sns.heatmap(condition_comp_norm, annot=True, fmt='.2f', ax=axes[2,0], cmap='coolwarm')
            axes[2,0].set_title('State Composition by Condition')
            
        # 8. Expression of key microglial markers
        microglia_markers = ['CX3CR1', 'P2RY12', 'TMEM119', 'AIF1', 'TREM2']
        available_markers = [m for m in microglia_markers if m in adata.var_names]
        
        if available_markers:
            marker_expr = pd.DataFrame(index=adata.obs_names)
            for marker in available_markers:
                marker_expr[marker] = adata[:, marker].X.toarray().flatten()
                
            marker_expr['activation_state'] = adata.obs['activation_state']
            
            marker_means = marker_expr.groupby('activation_state')[available_markers].mean()
            
            sns.heatmap(marker_means.T, annot=True, fmt='.2f', ax=axes[2,1], cmap='Reds')
            axes[2,1].set_title('Microglial Marker Expression')
            
        # 9. Dataset representation per state
        dataset_state = pd.crosstab(adata.obs['activation_state'], adata.obs['batch'])
        dataset_state.plot(kind='bar', stacked=True, ax=axes[2,2])
        axes[2,2].set_title('Dataset Representation per State')
        axes[2,2].tick_params(axis='x', rotation=45)
        axes[2,2].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        
        plt.tight_layout()
        plt.savefig(
            self.results_dir / 'figures/activation_states/comprehensive_analysis.pdf',
            dpi=300, bbox_inches='tight'
        )
        plt.close()
        
    def create_summary_report(self, adata, conserved_signatures, validation_results):
        """Create comprehensive summary report"""
        print("\nCreating summary report...")
        
        report = f"""
# Microglial Activation States Meta-Analysis Report

## Dataset Summary
- Total cells analyzed: {adata.n_obs:,}
- Total genes: {adata.n_vars:,} 
- Number of datasets: {len(adata.obs['batch'].unique())}
- Species: {', '.join(adata.obs['species'].unique())}

## Activation States Identified
- Number of activation states: {len(adata.obs['activation_state'].unique())}
- States: {', '.join(sorted(adata.obs['activation_state'].unique()))}

### State Distribution:
"""
        
        # Add state distribution table
        state_dist = adata.obs.groupby('activation_state').agg({
            'batch': 'nunique',
            'species': lambda x: ', '.join(x.unique())
        }).rename(columns={'batch': 'n_datasets'})
        
        state_counts = adata.obs['activation_state'].value_counts()
        state_dist['n_cells'] = state_counts
        state_dist['percentage'] = (state_counts / adata.n_obs * 100).round(1)
        
        report += state_dist.to_string()
        
        # Add conserved signatures summary
        if conserved_signatures:
            report += f"""

## Conserved Gene Signatures
- Total conserved signatures identified: {len(conserved_signatures)}

### Signature Summary:
"""
            for state, info in conserved_signatures.items():
                n_genes = len(info['genes'])
                n_datasets = info['n_datasets']
                report += f"\n- State {state}: {n_genes} conserved genes across {n_datasets} datasets"
                
        # Add validation results
        if validation_results:
            report += f"\n\n## Cross-Species Validation\n"
            for metric, value in validation_results.items():
                report += f"- {metric}: {value:.3f}\n"
                
        # Save report
        with open(self.results_dir / 'reports/activation_states_summary.md', 'w') as f:
            f.write(report)
            
        print("Summary report saved!")

def main():
    """Main activation state analysis pipeline"""
    analyzer = ActivationStateAnalyzer()
    
    # Load integrated data
    adata = analyzer.load_integrated_data(method='scvi')
    
    if adata is None:
        # Try alternative integration method
        adata = analyzer.load_integrated_data(method='scanorama')
        
    if adata is None:
        print("No integrated data found. Run integration first.")
        return
        
    print(f"\n{'='*60}")
    print("MICROGLIAL ACTIVATION STATE ANALYSIS")
    print(f"{'='*60}")
    
    # 1. Identify activation states
    adata = analyzer.identify_activation_states(adata)
    
    # 2. Differential expression analysis
    de_results = analyzer.analyze_differential_expression(adata)
    
    # 3. Pathway enrichment
    enrichment_results = analyzer.perform_pathway_enrichment(adata, de_results)
    
    # 4. Identify conserved signatures
    conserved_signatures = analyzer.identify_conserved_signatures(adata)
    
    # 5. Cross-species validation
    validation_results = analyzer.validate_across_species(adata)
    
    # 6. Generate comprehensive plots
    analyzer.generate_activation_plots(adata)
    
    # 7. Create summary report
    analyzer.create_summary_report(adata, conserved_signatures, validation_results)
    
    # Save final annotated data
    adata.write(analyzer.processed_data_dir / 'microglia_activation_states_final.h5ad')
    
    print(f"\n{'='*60}")
    print("ACTIVATION STATE ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print("\nKey findings:")
    print(f"- {len(adata.obs['activation_state'].unique())} activation states identified")
    print(f"- {len(conserved_signatures)} conserved signatures found" if conserved_signatures else "- No conserved signatures analysis completed")
    print(f"- Analysis covers {len(adata.obs['batch'].unique())} datasets")
    print(f"- {adata.n_obs:,} total microglial cells analyzed")
    
    print("\nNext step: Run meta-analysis validation")

if __name__ == "__main__":
    main()
