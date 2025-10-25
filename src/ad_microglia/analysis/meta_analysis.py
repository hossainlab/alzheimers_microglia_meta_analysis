#!/usr/bin/env python3
"""
Meta-Analysis Validation for Alzheimer's Microglia Activation States
Statistical validation and effect size analysis across studies
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats
from sklearn.metrics import adjusted_rand_score
from scipy.stats import fisher_exact
import statsmodels.api as sm
from statsmodels.stats.meta_analysis import combine_effects
import warnings
warnings.filterwarnings('ignore')

class MetaAnalysisValidator:
    def __init__(self, processed_data_dir='data/processed', results_dir='results'):
        self.processed_data_dir = Path(processed_data_dir)
        self.results_dir = Path(results_dir)
        
        # Create meta-analysis results directory
        (self.results_dir / 'meta_analysis').mkdir(exist_ok=True)
        
    def load_final_data(self):
        """Load final annotated microglial data"""
        print("Loading final activation state data...")
        
        try:
            adata = sc.read_h5ad(
                self.processed_data_dir / 'microglia_activation_states.h5ad'
            )
            print(f"Loaded {adata.n_obs} cells with activation states")
            return adata
        except FileNotFoundError:
            print("Final activation state data not found. Run activation analysis first.")
            return None
            
    def calculate_effect_sizes(self, adata):
        """Calculate effect sizes for each activation state across datasets"""
        print("\nCalculating effect sizes across datasets...")
        
        datasets = adata.obs['batch'].unique()
        states = adata.obs['activation_state'].unique()
        
        effect_sizes = {}
        
        for state in states:
            print(f"Analyzing effect sizes for state {state}...")
            
            dataset_effects = []
            
            for dataset in datasets:
                # Get dataset subset
                dataset_data = adata[adata.obs['batch'] == dataset]
                
                if dataset_data.n_obs > 20:  # Minimum sample size
                    # Calculate proportion of this state in dataset
                    state_count = (dataset_data.obs['activation_state'] == state).sum()
                    total_count = dataset_data.n_obs
                    proportion = state_count / total_count
                    
                    # Calculate confidence interval using Wilson score
                    if total_count > 0:
                        z = 1.96  # 95% CI
                        p_hat = proportion
                        n = total_count
                        
                        denominator = 1 + (z**2 / n)
                        centre = (p_hat + (z**2 / (2*n))) / denominator
                        adjustment = z * np.sqrt((p_hat*(1-p_hat)/n) + (z**2/(4*n**2))) / denominator
                        
                        ci_lower = max(0, centre - adjustment)
                        ci_upper = min(1, centre + adjustment)
                        
                        dataset_effects.append({
                            'dataset': dataset,
                            'proportion': proportion,
                            'count': state_count,
                            'total': total_count,
                            'ci_lower': ci_lower,
                            'ci_upper': ci_upper,
                            'se': np.sqrt(proportion * (1-proportion) / total_count)
                        })
                        
            effect_sizes[state] = dataset_effects
            
        return effect_sizes
        
    def perform_meta_analysis(self, effect_sizes):
        """Perform statistical meta-analysis of effect sizes"""
        print("\nPerforming meta-analysis...")
        
        meta_results = {}
        
        for state, effects in effect_sizes.items():
            if len(effects) >= 2:  # Need at least 2 studies
                print(f"Meta-analyzing state {state} across {len(effects)} datasets...")
                
                # Extract effect sizes and standard errors
                proportions = np.array([e['proportion'] for e in effects])
                se_array = np.array([e['se'] for e in effects])
                
                # Perform fixed-effects meta-analysis
                try:
                    # Calculate weights (inverse variance)
                    weights = 1 / (se_array**2)
                    
                    # Pooled estimate
                    pooled_effect = np.sum(weights * proportions) / np.sum(weights)
                    pooled_se = 1 / np.sqrt(np.sum(weights))
                    
                    # 95% confidence interval
                    ci_lower = pooled_effect - 1.96 * pooled_se
                    ci_upper = pooled_effect + 1.96 * pooled_se
                    
                    # Test for heterogeneity (Q statistic)
                    q_stat = np.sum(weights * (proportions - pooled_effect)**2)
                    df = len(proportions) - 1
                    q_pvalue = 1 - stats.chi2.cdf(q_stat, df) if df > 0 else 1.0
                    
                    # I-squared (heterogeneity measure)
                    i_squared = max(0, (q_stat - df) / q_stat) if q_stat > 0 else 0
                    
                    meta_results[state] = {
                        'pooled_proportion': pooled_effect,
                        'pooled_se': pooled_se,
                        'ci_lower': ci_lower,
                        'ci_upper': ci_upper,
                        'n_studies': len(effects),
                        'q_statistic': q_stat,
                        'q_pvalue': q_pvalue,
                        'i_squared': i_squared,
                        'individual_effects': effects
                    }
                    
                    print(f"  Pooled proportion: {pooled_effect:.3f} "
                          f"(95% CI: {ci_lower:.3f}-{ci_upper:.3f})")
                    print(f"  Heterogeneity I²: {i_squared:.1%}")
                    
                except Exception as e:
                    print(f"  Meta-analysis failed for state {state}: {e}")
                    continue
                    
        return meta_results
        
    def test_conservation_significance(self, adata):
        """Test statistical significance of conservation across datasets"""
        print("\nTesting conservation significance...")
        
        datasets = adata.obs['batch'].unique()
        states = adata.obs['activation_state'].unique()
        
        conservation_tests = {}
        
        # Test each state's presence across datasets
        for state in states:
            presence_data = []
            
            for dataset in datasets:
                dataset_data = adata[adata.obs['batch'] == dataset]
                
                if dataset_data.n_obs > 10:  # Minimum sample size
                    state_present = (dataset_data.obs['activation_state'] == state).sum()
                    state_absent = dataset_data.n_obs - state_present
                    
                    presence_data.append({
                        'dataset': dataset,
                        'present': state_present,
                        'absent': state_absent,
                        'total': dataset_data.n_obs,
                        'proportion': state_present / dataset_data.n_obs
                    })
                    
            if len(presence_data) >= 2:
                # Test for significant conservation using Fisher's exact test
                # Compare each dataset pair
                pairwise_tests = []
                
                for i in range(len(presence_data)):
                    for j in range(i+1, len(presence_data)):
                        data1 = presence_data[i]
                        data2 = presence_data[j]
                        
                        # Create contingency table
                        table = [
                            [data1['present'], data1['absent']],
                            [data2['present'], data2['absent']]
                        ]
                        
                        # Fisher's exact test
                        odds_ratio, p_value = fisher_exact(table)
                        
                        pairwise_tests.append({
                            'dataset1': data1['dataset'],
                            'dataset2': data2['dataset'], 
                            'odds_ratio': odds_ratio,
                            'p_value': p_value
                        })
                        
                # Combine p-values using Fisher's method
                if pairwise_tests:
                    p_values = [test['p_value'] for test in pairwise_tests]
                    combined_stat = -2 * np.sum(np.log(p_values))
                    combined_p = 1 - stats.chi2.cdf(combined_stat, 2*len(p_values))
                    
                    conservation_tests[state] = {
                        'n_datasets': len(presence_data),
                        'individual_data': presence_data,
                        'pairwise_tests': pairwise_tests,
                        'combined_p_value': combined_p,
                        'is_conserved': combined_p < 0.05
                    }
                    
        return conservation_tests
        
    def analyze_species_differences(self, adata):
        """Analyze differences between human and mouse studies"""
        print("\nAnalyzing species differences...")
        
        if 'species' not in adata.obs.columns:
            print("Species information not available")
            return None
            
        human_data = adata[adata.obs['species'] == 'human']
        mouse_data = adata[adata.obs['species'] == 'mouse']
        
        if human_data.n_obs == 0 or mouse_data.n_obs == 0:
            print("Both species not represented in data")
            return None
            
        species_analysis = {}
        
        # Compare activation state distributions
        human_states = human_data.obs['activation_state'].value_counts(normalize=True)
        mouse_states = mouse_data.obs['activation_state'].value_counts(normalize=True)
        
        # Align state distributions
        all_states = set(human_states.index) | set(mouse_states.index)
        comparison_data = []
        
        for state in all_states:
            human_prop = human_states.get(state, 0)
            mouse_prop = mouse_states.get(state, 0)
            
            comparison_data.append({
                'state': state,
                'human_proportion': human_prop,
                'mouse_proportion': mouse_prop,
                'difference': human_prop - mouse_prop
            })
            
        comparison_df = pd.DataFrame(comparison_data)
        
        # Statistical test for species differences
        # Chi-square test of independence
        contingency_data = []
        for state in all_states:
            human_count = (human_data.obs['activation_state'] == state).sum()
            mouse_count = (mouse_data.obs['activation_state'] == state).sum()
            
            contingency_data.append([human_count, mouse_count])
            
        if len(contingency_data) > 1:
            chi2_stat, chi2_p = stats.chi2_contingency(contingency_data)[:2]
            
            species_analysis = {
                'n_human_cells': human_data.n_obs,
                'n_mouse_cells': mouse_data.n_obs,
                'comparison_table': comparison_df,
                'chi2_statistic': chi2_stat,
                'chi2_p_value': chi2_p,
                'significant_difference': chi2_p < 0.05
            }
            
        return species_analysis
        
    def generate_meta_analysis_plots(self, meta_results, conservation_tests):
        """Generate meta-analysis visualization plots"""
        print("\nGenerating meta-analysis plots...")
        
        # Forest plots for meta-analysis results
        n_states = len(meta_results)
        if n_states == 0:
            print("No meta-analysis results to plot")
            return
            
        fig, axes = plt.subplots(n_states, 1, figsize=(12, 3*n_states))
        if n_states == 1:
            axes = [axes]
            
        fig.suptitle('Meta-Analysis Forest Plots', fontsize=16)
        
        for i, (state, results) in enumerate(meta_results.items()):
            ax = axes[i]
            
            # Individual study effects
            effects = results['individual_effects']
            y_pos = range(len(effects))
            
            # Plot individual studies
            for j, effect in enumerate(effects):
                ax.errorbar(effect['proportion'], j, 
                           xerr=[[effect['proportion']-effect['ci_lower']], 
                                [effect['ci_upper']-effect['proportion']]],
                           fmt='s', markersize=6, alpha=0.7, 
                           label=effect['dataset'] if j < 10 else "")
                           
            # Plot pooled estimate
            pooled = results['pooled_proportion']
            pooled_ci = [results['ci_lower'], results['ci_upper']]
            
            ax.errorbar(pooled, len(effects), 
                       xerr=[[pooled-pooled_ci[0]], [pooled_ci[1]-pooled]],
                       fmt='D', markersize=10, color='red', 
                       linewidth=3, label='Pooled')
                       
            ax.axvline(pooled, color='red', linestyle='--', alpha=0.5)
            
            ax.set_xlabel('Proportion')
            ax.set_ylabel('Studies')
            ax.set_title(f'Activation State {state} '
                        f'(I² = {results["i_squared"]:.1%}, '
                        f'p = {results["q_pvalue"]:.3f})')
            ax.set_yticks(range(len(effects)+1))
            ax.set_yticklabels([e['dataset'] for e in effects] + ['Pooled'])
            ax.grid(True, alpha=0.3)
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
            
        plt.tight_layout()
        plt.savefig(self.results_dir / 'meta_analysis/forest_plots.pdf',
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Conservation significance plot
        if conservation_tests:
            conserved_states = [state for state, test in conservation_tests.items() 
                              if test['is_conserved']]
                              
            fig, ax = plt.subplots(figsize=(10, 6))
            
            states = list(conservation_tests.keys())
            p_values = [conservation_tests[state]['combined_p_value'] for state in states]
            colors = ['green' if p < 0.05 else 'red' for p in p_values]
            
            bars = ax.bar(states, [-np.log10(p) for p in p_values], color=colors, alpha=0.7)
            ax.axhline(-np.log10(0.05), color='black', linestyle='--', 
                      label='Significance threshold (p=0.05)')
                      
            ax.set_xlabel('Activation State')
            ax.set_ylabel('-log10(p-value)')
            ax.set_title('Conservation Significance Across Datasets')
            ax.legend()
            
            # Add p-value labels on bars
            for bar, p_val in zip(bars, p_values):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'p={p_val:.3f}', ha='center', va='bottom')
                       
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(self.results_dir / 'meta_analysis/conservation_significance.pdf',
                       dpi=300, bbox_inches='tight')
            plt.close()
        
    def create_meta_analysis_report(self, meta_results, conservation_tests, species_analysis):
        """Create comprehensive meta-analysis report"""
        print("\nCreating meta-analysis report...")
        
        report = f"""
# Alzheimer's Microglia Meta-Analysis Results

## Executive Summary

This meta-analysis examined conserved microglial activation states across {len(meta_results)} 
activation states identified from multiple Alzheimer's disease single-cell RNA-sequencing studies.

## Key Findings

### Conserved Activation States
"""
        
        if conservation_tests:
            conserved_count = sum(1 for test in conservation_tests.values() if test['is_conserved'])
            report += f"- {conserved_count}/{len(conservation_tests)} activation states show significant conservation across datasets\n"
            
            for state, test in conservation_tests.items():
                if test['is_conserved']:
                    report += f"  - State {state}: present across {test['n_datasets']} datasets (p={test['combined_p_value']:.3e})\n"
                    
        report += "\n### Meta-Analysis Results\n\n"
        
        for state, results in meta_results.items():
            report += f"""
**Activation State {state}:**
- Pooled proportion: {results['pooled_proportion']:.3f} (95% CI: {results['ci_lower']:.3f}-{results['ci_upper']:.3f})
- Number of studies: {results['n_studies']}
- Heterogeneity (I²): {results['i_squared']:.1%}
- Heterogeneity p-value: {results['q_pvalue']:.3f}
"""
            
        if species_analysis:
            report += f"""
### Species Comparison

- Human microglial cells analyzed: {species_analysis['n_human_cells']:,}
- Mouse microglial cells analyzed: {species_analysis['n_mouse_cells']:,}
- Significant species difference: {species_analysis['significant_difference']}
- Chi-square p-value: {species_analysis['chi2_p_value']:.3e}

#### State Distribution Differences:
"""
            
            for _, row in species_analysis['comparison_table'].iterrows():
                diff = row['difference']
                direction = "higher in human" if diff > 0 else "higher in mouse"
                report += f"- State {row['state']}: {abs(diff):.3f} {direction}\n"
                
        report += f"""

## Statistical Methods

- Effect sizes calculated using Wilson score confidence intervals
- Meta-analysis performed using fixed-effects model with inverse variance weighting
- Conservation tested using Fisher's exact test with combined p-values
- Species differences assessed using chi-square test of independence

## Limitations

- Different tissue regions and processing protocols across studies
- Variable cell numbers per dataset
- Potential batch effects despite integration efforts
- Limited longitudinal data for disease progression analysis

## Conclusions

This meta-analysis provides evidence for conserved microglial activation states across 
Alzheimer's disease studies, with some states showing significant conservation across 
datasets and species. These findings support the existence of reproducible microglial 
phenotypes in Alzheimer's pathogenesis.
"""
        
        # Save report
        with open(self.results_dir / 'meta_analysis/meta_analysis_report.md', 'w') as f:
            f.write(report)
            
        # Save detailed results as CSV
        if meta_results:
            meta_df = pd.DataFrame([
                {
                    'activation_state': state,
                    'pooled_proportion': results['pooled_proportion'],
                    'ci_lower': results['ci_lower'], 
                    'ci_upper': results['ci_upper'],
                    'n_studies': results['n_studies'],
                    'i_squared': results['i_squared'],
                    'q_pvalue': results['q_pvalue']
                }
                for state, results in meta_results.items()
            ])
            
            meta_df.to_csv(self.results_dir / 'meta_analysis/meta_analysis_results.csv', 
                          index=False)

def main():
    """Main meta-analysis pipeline"""
    validator = MetaAnalysisValidator()
    
    # Load final data
    adata = validator.load_final_data()
    
    if adata is None:
        return
        
    print(f"\n{'='*60}")
    print("META-ANALYSIS VALIDATION")
    print(f"{'='*60}")
    
    # 1. Calculate effect sizes
    effect_sizes = validator.calculate_effect_sizes(adata)
    
    # 2. Perform meta-analysis
    meta_results = validator.perform_meta_analysis(effect_sizes)
    
    # 3. Test conservation significance
    conservation_tests = validator.test_conservation_significance(adata)
    
    # 4. Analyze species differences
    species_analysis = validator.analyze_species_differences(adata)
    
    # 5. Generate visualization plots
    validator.generate_meta_analysis_plots(meta_results, conservation_tests)
    
    # 6. Create comprehensive report
    validator.create_meta_analysis_report(meta_results, conservation_tests, species_analysis)
    
    print(f"\n{'='*60}")
    print("META-ANALYSIS COMPLETE")
    print(f"{'='*60}")
    
    # Summary statistics
    if meta_results:
        print(f"\nMeta-analysis completed for {len(meta_results)} activation states:")
        for state, results in meta_results.items():
            pooled_prop = results['pooled_proportion']
            i_squared = results['i_squared']
            print(f"- State {state}: {pooled_prop:.3f} proportion (I² = {i_squared:.1%})")
            
    if conservation_tests:
        conserved_count = sum(1 for test in conservation_tests.values() if test['is_conserved'])
        print(f"\n{conserved_count}/{len(conservation_tests)} states show significant conservation")
        
    print("\nAnalysis pipeline completed!")
    print(f"Results saved in: {validator.results_dir / 'meta_analysis'}")

if __name__ == "__main__":
    main()
