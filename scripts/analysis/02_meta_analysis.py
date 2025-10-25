#!/usr/bin/env python3
"""
Meta-Analysis Script for Alzheimer's Microglia Meta-Analysis
Performs statistical meta-analysis to validate conserved activation states
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Try importing statistics packages
try:
    import statsmodels.api as sm
    from statsmodels.stats.meta_analysis import combine_effects
    STATSMODELS_AVAILABLE = True
except ImportError:
    print("Warning: statsmodels not available. Install with: pip install statsmodels")
    STATSMODELS_AVAILABLE = False

try:
    from scipy import stats
    from scipy.stats import chi2_contingency, fisher_exact
    SCIPY_AVAILABLE = True
except ImportError:
    print("Warning: scipy not available.")
    SCIPY_AVAILABLE = False

# Configure scanpy
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=300, facecolor='white')

class MetaAnalyzer:
    def __init__(self, processed_data_dir='data/processed', results_dir='results'):
        self.processed_data_dir = Path(processed_data_dir)
        self.results_dir = Path(results_dir)
        self.meta_results_dir = self.results_dir / 'meta_analysis'
        self.tables_dir = self.results_dir / 'tables'
        
        # Create directories
        self.meta_results_dir.mkdir(parents=True, exist_ok=True)
        
    def load_activation_data(self):
        """Load activation states data"""
        print("Loading activation states data...")
        
        activation_file = self.processed_data_dir / 'microglia_activation_states.h5ad'
        if not activation_file.exists():
            print("Activation states data not found. Run activation analysis first.")
            return None
            
        adata = sc.read_h5ad(activation_file)
        print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
        print(f"Datasets: {', '.join(adata.obs['batch'].unique())}")
        print(f"Activation states: {', '.join(sorted(adata.obs['activation_state'].unique()))}")
        
        return adata
    
    def calculate_effect_sizes(self, adata):
        """Calculate effect sizes for activation state proportions across datasets"""
        print("\nCalculating effect sizes...")
        
        results = []
        
        for state in adata.obs['activation_state'].unique():
            print(f"Analyzing activation state {state}...")
            
            # Calculate proportions in each dataset
            for dataset in adata.obs['batch'].unique():
                dataset_mask = adata.obs['batch'] == dataset
                dataset_cells = dataset_mask.sum()
                
                if dataset_cells < 20:  # Skip datasets with too few cells
                    continue
                
                state_cells = (dataset_mask & (adata.obs['activation_state'] == state)).sum()
                proportion = state_cells / dataset_cells
                
                # Calculate confidence intervals using Wilson score
                if dataset_cells > 0:
                    ci_lower, ci_upper = self.wilson_confidence_interval(state_cells, dataset_cells)
                else:
                    ci_lower, ci_upper = 0, 0
                
                # Calculate standard error
                se = np.sqrt((proportion * (1 - proportion)) / dataset_cells)
                
                results.append({
                    'activation_state': state,
                    'dataset': dataset,
                    'n_cells_total': dataset_cells,
                    'n_cells_state': state_cells,
                    'proportion': proportion,
                    'se': se,
                    'ci_lower': ci_lower,
                    'ci_upper': ci_upper
                })
        
        effect_sizes_df = pd.DataFrame(results)
        
        # Save results
        effect_sizes_df.to_csv(self.meta_results_dir / 'effect_sizes.csv', index=False)
        
        return effect_sizes_df
    
    def wilson_confidence_interval(self, successes, trials, confidence=0.95):
        """Calculate Wilson score confidence interval"""
        if trials == 0:
            return 0, 0
            
        z = stats.norm.ppf(1 - (1 - confidence) / 2)
        p = successes / trials
        
        denominator = 1 + z**2 / trials
        center = (p + z**2 / (2 * trials)) / denominator
        margin = z * np.sqrt((p * (1 - p) + z**2 / (4 * trials)) / trials) / denominator
        
        return max(0, center - margin), min(1, center + margin)
    
    def perform_fixed_effects_meta_analysis(self, effect_sizes_df):
        """Perform fixed-effects meta-analysis"""
        print("\nPerforming fixed-effects meta-analysis...")
        
        if not STATSMODELS_AVAILABLE:
            print("Statsmodels not available, using manual calculation")
            return self.manual_meta_analysis(effect_sizes_df)
        
        meta_results = []
        
        for state in effect_sizes_df['activation_state'].unique():
            state_data = effect_sizes_df[effect_sizes_df['activation_state'] == state]
            
            if len(state_data) < 2:
                continue
                
            # Remove rows with zero standard error
            state_data = state_data[state_data['se'] > 0]
            
            if len(state_data) < 2:
                continue
            
            try:
                # Fixed-effects meta-analysis
                effects = state_data['proportion'].values
                se = state_data['se'].values
                
                # Calculate weights (inverse variance)
                weights = 1 / (se**2)
                
                # Pooled estimate
                pooled_effect = np.sum(weights * effects) / np.sum(weights)
                pooled_se = np.sqrt(1 / np.sum(weights))
                
                # Confidence interval
                z = stats.norm.ppf(0.975)
                ci_lower = pooled_effect - z * pooled_se
                ci_upper = pooled_effect + z * pooled_se
                
                # Z-test for significance
                z_score = pooled_effect / pooled_se
                p_value = 2 * (1 - stats.norm.cdf(abs(z_score)))
                
                # Heterogeneity statistics (I²)
                Q = np.sum(weights * (effects - pooled_effect)**2)
                df = len(effects) - 1
                I2 = max(0, (Q - df) / Q) if Q > 0 else 0
                
                meta_results.append({
                    'activation_state': state,
                    'pooled_proportion': pooled_effect,
                    'se': pooled_se,
                    'ci_lower': ci_lower,
                    'ci_upper': ci_upper,
                    'z_score': z_score,
                    'p_value': p_value,
                    'Q_statistic': Q,
                    'I2_heterogeneity': I2,
                    'n_studies': len(state_data)
                })
                
            except Exception as e:
                print(f"Error in meta-analysis for state {state}: {e}")
                continue
        
        meta_results_df = pd.DataFrame(meta_results)
        
        # Adjust p-values for multiple testing
        if len(meta_results_df) > 1:
            from statsmodels.stats.multitest import multipletests
            _, meta_results_df['p_value_adj'], _, _ = multipletests(
                meta_results_df['p_value'], 
                method='fdr_bh'
            )
        
        # Save results
        meta_results_df.to_csv(self.meta_results_dir / 'meta_analysis_results.csv', index=False)
        
        return meta_results_df
    
    def manual_meta_analysis(self, effect_sizes_df):
        """Manual meta-analysis calculation when statsmodels not available"""
        print("Performing manual meta-analysis...")
        
        meta_results = []
        
        for state in effect_sizes_df['activation_state'].unique():
            state_data = effect_sizes_df[effect_sizes_df['activation_state'] == state]
            
            if len(state_data) < 2:
                continue
                
            # Remove rows with zero standard error
            state_data = state_data[state_data['se'] > 0]
            
            if len(state_data) < 2:
                continue
            
            effects = state_data['proportion'].values
            se = state_data['se'].values
            
            # Calculate weights (inverse variance)
            weights = 1 / (se**2)
            
            # Pooled estimate
            pooled_effect = np.sum(weights * effects) / np.sum(weights)
            pooled_se = np.sqrt(1 / np.sum(weights))
            
            # Confidence interval (95%)
            ci_lower = pooled_effect - 1.96 * pooled_se
            ci_upper = pooled_effect + 1.96 * pooled_se
            
            # Basic significance test
            z_score = pooled_effect / pooled_se
            p_value = 2 * (1 - stats.norm.cdf(abs(z_score))) if SCIPY_AVAILABLE else np.nan
            
            meta_results.append({
                'activation_state': state,
                'pooled_proportion': pooled_effect,
                'se': pooled_se,
                'ci_lower': ci_lower,
                'ci_upper': ci_upper,
                'z_score': z_score,
                'p_value': p_value,
                'n_studies': len(state_data)
            })
        
        return pd.DataFrame(meta_results)
    
    def test_conservation_significance(self, adata, effect_sizes_df):
        """Test statistical significance of conservation across datasets"""
        print("\nTesting conservation significance...")
        
        conservation_results = []
        
        for state in adata.obs['activation_state'].unique():
            print(f"Testing conservation for state {state}...")
            
            # Create contingency table
            datasets = adata.obs['batch'].unique()
            contingency_data = []
            
            for dataset in datasets:
                dataset_mask = adata.obs['batch'] == dataset
                state_cells = (dataset_mask & (adata.obs['activation_state'] == state)).sum()
                other_cells = dataset_mask.sum() - state_cells
                contingency_data.append([state_cells, other_cells])
            
            contingency_table = np.array(contingency_data)
            
            # Skip if any dataset has no cells of this state
            if np.any(contingency_table[:, 0] == 0):
                continue
            
            try:
                # Chi-square test for independence
                if SCIPY_AVAILABLE:
                    chi2, p_chi2, dof, expected = chi2_contingency(contingency_table)
                else:
                    chi2, p_chi2, dof = np.nan, np.nan, np.nan
                
                # Calculate coefficient of variation for proportions
                state_proportions = contingency_table[:, 0] / contingency_table.sum(axis=1)
                cv = stats.variation(state_proportions) if SCIPY_AVAILABLE else np.std(state_proportions) / np.mean(state_proportions)
                
                conservation_results.append({
                    'activation_state': state,
                    'chi2_statistic': chi2,
                    'chi2_p_value': p_chi2,
                    'degrees_freedom': dof,
                    'coefficient_variation': cv,
                    'mean_proportion': np.mean(state_proportions),
                    'std_proportion': np.std(state_proportions),
                    'n_datasets': len(datasets)
                })
                
            except Exception as e:
                print(f"Error testing conservation for state {state}: {e}")
                continue
        
        conservation_df = pd.DataFrame(conservation_results)
        
        # Adjust p-values
        if len(conservation_df) > 1 and STATSMODELS_AVAILABLE:
            from statsmodels.stats.multitest import multipletests
            _, conservation_df['chi2_p_value_adj'], _, _ = multipletests(
                conservation_df['chi2_p_value'].fillna(1), 
                method='fdr_bh'
            )
        
        # Save results
        conservation_df.to_csv(self.meta_results_dir / 'conservation_significance.csv', index=False)
        
        return conservation_df
    
    def create_forest_plots(self, effect_sizes_df, meta_results_df):
        """Create forest plots for meta-analysis results"""
        print("\nCreating forest plots...")
        
        # Determine number of states
        states = sorted(meta_results_df['activation_state'].unique())
        n_states = len(states)
        
        if n_states == 0:
            print("No states to plot")
            return
        
        # Create subplot grid
        fig, axes = plt.subplots(n_states, 1, figsize=(12, 4 * n_states))
        if n_states == 1:
            axes = [axes]
        
        for i, state in enumerate(states):
            ax = axes[i]
            
            # Get data for this state
            state_effects = effect_sizes_df[effect_sizes_df['activation_state'] == state]
            state_meta = meta_results_df[meta_results_df['activation_state'] == state].iloc[0]
            
            # Plot individual studies
            y_positions = range(len(state_effects))
            
            for j, (_, row) in enumerate(state_effects.iterrows()):
                # Plot point estimate and confidence interval
                ax.errorbar(row['proportion'], j, 
                           xerr=[[row['proportion'] - row['ci_lower']], 
                                [row['ci_upper'] - row['proportion']]], 
                           fmt='s', color='blue', alpha=0.7, 
                           markersize=6, capsize=3)
                
                # Add dataset label
                ax.text(-0.05, j, row['dataset'], ha='right', va='center', fontsize=10)
            
            # Add pooled estimate
            pooled_y = len(state_effects) + 0.5
            ax.errorbar(state_meta['pooled_proportion'], pooled_y,
                       xerr=[[state_meta['pooled_proportion'] - state_meta['ci_lower']],
                            [state_meta['ci_upper'] - state_meta['pooled_proportion']]], 
                       fmt='D', color='red', markersize=8, capsize=4, linewidth=2)
            ax.text(-0.05, pooled_y, 'Pooled', ha='right', va='center', 
                   fontsize=10, fontweight='bold')
            
            # Add vertical line at pooled estimate
            ax.axvline(state_meta['pooled_proportion'], color='red', 
                      linestyle='--', alpha=0.5)
            
            # Formatting
            ax.set_xlim(0, max(state_effects['ci_upper'].max(), 
                              state_meta['ci_upper']) * 1.1)
            ax.set_ylim(-0.5, len(state_effects) + 1)
            ax.set_xlabel('Proportion')
            ax.set_title(f'Activation State {state} - Forest Plot\n'
                        f'Pooled: {state_meta["pooled_proportion"]:.3f} '
                        f'(95% CI: {state_meta["ci_lower"]:.3f}-{state_meta["ci_upper"]:.3f})\n'
                        f'P-value: {state_meta["p_value"]:.2e}', fontsize=12)
            ax.grid(True, alpha=0.3, axis='x')
        
        plt.tight_layout()
        plt.savefig(self.meta_results_dir / 'forest_plots.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Forest plots saved")
    
    def create_conservation_plots(self, conservation_df):
        """Create conservation analysis plots"""
        print("\nCreating conservation plots...")
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot 1: Conservation significance
        ax1 = axes[0, 0]
        if 'chi2_p_value' in conservation_df.columns:
            # -log10 p-values
            log_p = -np.log10(conservation_df['chi2_p_value'].fillna(1))
            bars = ax1.bar(conservation_df['activation_state'], log_p, 
                          color=['red' if p < 0.05 else 'blue' 
                                for p in conservation_df['chi2_p_value'].fillna(1)])
            ax1.axhline(-np.log10(0.05), color='red', linestyle='--', alpha=0.7, 
                       label='p = 0.05')
            ax1.set_xlabel('Activation State')
            ax1.set_ylabel('-log10(p-value)')
            ax1.set_title('Conservation Significance Test')
            ax1.legend()
        
        # Plot 2: Coefficient of variation
        ax2 = axes[0, 1]
        if 'coefficient_variation' in conservation_df.columns:
            bars = ax2.bar(conservation_df['activation_state'], 
                          conservation_df['coefficient_variation'],
                          color='green', alpha=0.7)
            ax2.set_xlabel('Activation State')
            ax2.set_ylabel('Coefficient of Variation')
            ax2.set_title('Proportion Variability Across Datasets')
            
            # Add values on bars
            for bar, cv in zip(bars, conservation_df['coefficient_variation']):
                height = bar.get_height()
                ax2.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                        f'{cv:.3f}', ha='center', va='bottom', fontsize=10)
        
        # Plot 3: Mean proportions
        ax3 = axes[1, 0]
        if 'mean_proportion' in conservation_df.columns:
            bars = ax3.bar(conservation_df['activation_state'], 
                          conservation_df['mean_proportion'],
                          yerr=conservation_df['std_proportion'],
                          capsize=5, color='purple', alpha=0.7)
            ax3.set_xlabel('Activation State')
            ax3.set_ylabel('Mean Proportion')
            ax3.set_title('Mean Proportions Across Datasets')
        
        # Plot 4: Sample sizes
        ax4 = axes[1, 1]
        if 'n_datasets' in conservation_df.columns:
            bars = ax4.bar(conservation_df['activation_state'], 
                          conservation_df['n_datasets'],
                          color='orange', alpha=0.7)
            ax4.set_xlabel('Activation State')
            ax4.set_ylabel('Number of Datasets')
            ax4.set_title('Dataset Coverage')
            
            # Add values on bars
            for bar, n in zip(bars, conservation_df['n_datasets']):
                height = bar.get_height()
                ax4.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                        str(int(n)), ha='center', va='bottom', fontsize=10)
        
        plt.suptitle('Conservation Analysis Results', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(self.meta_results_dir / 'conservation_significance.pdf', 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        print("Conservation plots saved")
    
    def generate_meta_analysis_report(self, meta_results_df, conservation_df, effect_sizes_df):
        """Generate comprehensive meta-analysis report"""
        print("\nGenerating meta-analysis report...")
        
        # Calculate summary statistics
        n_states = len(meta_results_df)
        significant_states = (meta_results_df['p_value'] < 0.05).sum() if 'p_value' in meta_results_df.columns else 0
        
        if 'p_value_adj' in meta_results_df.columns:
            significant_adj = (meta_results_df['p_value_adj'] < 0.05).sum()
        else:
            significant_adj = 0
        
        # Most conserved states
        if 'coefficient_variation' in conservation_df.columns:
            most_conserved = conservation_df.nsmallest(3, 'coefficient_variation')
        else:
            most_conserved = pd.DataFrame()
        
        report = f"""# Meta-Analysis Report: Alzheimer's Microglia Activation States

## Summary
- **Total activation states analyzed**: {n_states}
- **States with significant pooled effects**: {significant_states}
- **States significant after FDR correction**: {significant_adj}
- **Datasets included**: {len(effect_sizes_df['dataset'].unique())}
- **Total cells analyzed**: {effect_sizes_df['n_cells_total'].sum():,}

## Meta-Analysis Results

### Significant Activation States
"""
        
        # Add significant states
        if 'p_value' in meta_results_df.columns:
            sig_states = meta_results_df[meta_results_df['p_value'] < 0.05].sort_values('p_value')
            
            for _, row in sig_states.iterrows():
                report += f"""
**State {row['activation_state']}**:
- Pooled proportion: {row['pooled_proportion']:.3f} (95% CI: {row['ci_lower']:.3f}-{row['ci_upper']:.3f})
- P-value: {row['p_value']:.2e}
- Studies included: {int(row['n_studies'])}
"""
                if 'I2_heterogeneity' in row:
                    report += f"- Heterogeneity (I²): {row['I2_heterogeneity']:.1f}%\n"
        
        # Conservation analysis
        report += f"""
## Conservation Analysis

### Most Conserved States (lowest variability):
"""
        
        if not most_conserved.empty and 'coefficient_variation' in most_conserved.columns:
            for _, row in most_conserved.iterrows():
                report += f"""
**State {row['activation_state']}**:
- Coefficient of variation: {row['coefficient_variation']:.3f}
- Mean proportion: {row['mean_proportion']:.3f} ± {row['std_proportion']:.3f}
"""
                if 'chi2_p_value' in row:
                    report += f"- Conservation p-value: {row['chi2_p_value']:.2e}\n"
        
        # Interpretation
        report += f"""
## Interpretation

### Statistical Significance
- States with p < 0.05 show consistent presence across datasets
- Low coefficient of variation indicates stable proportions
- Chi-square tests assess independence from dataset effects

### Conservation Criteria
Strong conservation evidence requires:
1. Significant pooled effect (p < 0.05)
2. Low coefficient of variation (< 0.5)
3. Present in multiple datasets (≥ 3)
4. Non-significant chi-square test (p > 0.05)

## Files Generated
- **Effect sizes**: `{self.meta_results_dir / 'effect_sizes.csv'}`
- **Meta-analysis results**: `{self.meta_results_dir / 'meta_analysis_results.csv'}`
- **Conservation tests**: `{self.meta_results_dir / 'conservation_significance.csv'}`
- **Forest plots**: `{self.meta_results_dir / 'forest_plots.pdf'}`
- **Conservation plots**: `{self.meta_results_dir / 'conservation_significance.pdf'}`

## Conclusions
"""
        
        # Add specific conclusions
        if significant_states > 0:
            report += f"\n- {significant_states} activation states show significant conservation across studies"
        
        if not most_conserved.empty:
            most_conserved_state = most_conserved.iloc[0]['activation_state']
            report += f"\n- State {most_conserved_state} shows the highest conservation"
        
        report += f"\n- Results support the existence of reproducible microglial activation states in Alzheimer's disease"
        
        # Save report
        report_path = self.meta_results_dir / 'meta_analysis_report.md'
        with open(report_path, 'w') as f:
            f.write(report)
        
        print(f"Meta-analysis report saved to {report_path}")
    
    def run_meta_analysis(self):
        """Run complete meta-analysis pipeline"""
        print("Starting meta-analysis...")
        
        # Load data
        adata = self.load_activation_data()
        if adata is None:
            return
        
        # Calculate effect sizes
        effect_sizes_df = self.calculate_effect_sizes(adata)
        
        # Perform meta-analysis
        meta_results_df = self.perform_fixed_effects_meta_analysis(effect_sizes_df)
        
        # Test conservation significance
        conservation_df = self.test_conservation_significance(adata, effect_sizes_df)
        
        # Create plots
        if not meta_results_df.empty:
            self.create_forest_plots(effect_sizes_df, meta_results_df)
        
        if not conservation_df.empty:
            self.create_conservation_plots(conservation_df)
        
        # Generate report
        self.generate_meta_analysis_report(meta_results_df, conservation_df, effect_sizes_df)
        
        print(f"\nMeta-analysis complete!")
        print(f"Results saved in: {self.meta_results_dir}")
        
        return meta_results_df, conservation_df

def main():
    analyzer = MetaAnalyzer()
    analyzer.run_meta_analysis()

if __name__ == "__main__":
    main()