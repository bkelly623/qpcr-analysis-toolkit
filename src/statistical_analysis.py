"""
Statistical Analysis Module for qPCR Data

This module provides comprehensive statistical testing for qPCR analysis results,
including t-tests, multiple comparison corrections, and confidence intervals.

Author: Your Name
Date: June 2025
"""

import pandas as pd
import numpy as np
from scipy import stats
from scipy.stats import ttest_ind, ttest_rel
import warnings


class qPCRStatistics:
    """
    Statistical analysis class for qPCR data.
    
    Provides methods for:
    - Student's t-tests (independent and paired)
    - Multiple comparison corrections
    - Confidence interval calculations
    - Effect size measurements
    - Power analysis
    """
    
    def __init__(self, alpha=0.05):
        """
        Initialize the statistics analyzer.
        
        Parameters:
        -----------
        alpha : float
            Significance level for statistical tests (default: 0.05)
        """
        self.alpha = alpha
        self.results = {}
        
    def perform_ttest(self, data, gene_column='Target_Gene', condition_column='Condition', 
                     value_column='Fold_Change', control_condition='Control', 
                     test_type='independent'):
        """
        Perform t-tests comparing treatment groups to control.
        
        Parameters:
        -----------
        data : pandas.DataFrame
            Data containing fold change values
        gene_column : str
            Column name containing gene identifiers
        condition_column : str
            Column name containing experimental conditions
        value_column : str
            Column name containing values to test (default: 'Fold_Change')
        control_condition : str
            Name of the control condition
        test_type : str
            Type of t-test ('independent' or 'paired')
            
        Returns:
        --------
        pandas.DataFrame
            Statistical test results for each gene and condition
        """
        results = []
        
        genes = data[gene_column].unique()
        conditions = data[condition_column].unique()
        treatment_conditions = [c for c in conditions if c != control_condition]
        
        for gene in genes:
            gene_data = data[data[gene_column] == gene]
            
            # Get control data for this gene
            control_data = gene_data[gene_data[condition_column] == control_condition][value_column]
            
            if len(control_data) == 0:
                warnings.warn(f"No control data found for gene {gene}")
                continue
                
            for treatment in treatment_conditions:
                treatment_data = gene_data[gene_data[condition_column] == treatment][value_column]
                
                if len(treatment_data) == 0:
                    continue
                    
                # Remove NaN values
                control_clean = control_data.dropna()
                treatment_clean = treatment_data.dropna()
                
                if len(control_clean) < 2 or len(treatment_clean) < 2:
                    warnings.warn(f"Insufficient data for {gene} {treatment} comparison")
                    continue
                
                # Perform appropriate t-test
                if test_type == 'independent':
                    t_stat, p_value = ttest_ind(treatment_clean, control_clean)
                elif test_type == 'paired':
                    if len(treatment_clean) != len(control_clean):
                        warnings.warn(f"Unequal sample sizes for paired test: {gene} {treatment}")
                        t_stat, p_value = ttest_ind(treatment_clean, control_clean)
                    else:
                        t_stat, p_value = ttest_rel(treatment_clean, control_clean)
                else:
                    raise ValueError("test_type must be 'independent' or 'paired'")
                
                # Calculate descriptive statistics
                control_mean = control_clean.mean()
                control_std = control_clean.std()
                treatment_mean = treatment_clean.mean()
                treatment_std = treatment_clean.std()
                
                # Calculate effect size (Cohen's d)
                pooled_std = np.sqrt(((len(control_clean) - 1) * control_std**2 + 
                                    (len(treatment_clean) - 1) * treatment_std**2) / 
                                   (len(control_clean) + len(treatment_clean) - 2))
                cohens_d = (treatment_mean - control_mean) / pooled_std if pooled_std > 0 else np.nan
                
                # Calculate 95% confidence interval for the difference
                se_diff = np.sqrt((control_std**2 / len(control_clean)) + 
                                (treatment_std**2 / len(treatment_clean)))
                df = len(control_clean) + len(treatment_clean) - 2
                t_critical = stats.t.ppf(1 - self.alpha/2, df)
                mean_diff = treatment_mean - control_mean
                ci_lower = mean_diff - t_critical * se_diff
                ci_upper = mean_diff + t_critical * se_diff
                
                results.append({
                    'Gene': gene,
                    'Comparison': f"{treatment} vs {control_condition}",
                    'Control_Mean': round(control_mean, 3),
                    'Control_Std': round(control_std, 3),
                    'Control_N': len(control_clean),
                    'Treatment_Mean': round(treatment_mean, 3),
                    'Treatment_Std': round(treatment_std, 3),
                    'Treatment_N': len(treatment_clean),
                    'Mean_Difference': round(mean_diff, 3),
                    'T_Statistic': round(t_stat, 3),
                    'P_Value': p_value,
                    'Significant': p_value < self.alpha,
                    'Cohens_D': round(cohens_d, 3),
                    'Effect_Size_Interpretation': self._interpret_effect_size(cohens_d),
                    'CI_Lower': round(ci_lower, 3),
                    'CI_Upper': round(ci_upper, 3),
                    'Degrees_Freedom': df
                })
        
        self.results['ttest'] = pd.DataFrame(results)
        return self.results['ttest']
    
    def apply_multiple_comparison_correction(self, method='bonferroni'):
        """
        Apply multiple comparison correction to p-values.
        
        Parameters:
        -----------
        method : str
            Correction method ('bonferroni', 'holm', 'fdr_bh')
            
        Returns:
        --------
        pandas.DataFrame
            Results with corrected p-values
        """
        if 'ttest' not in self.results:
            raise ValueError("No t-test results found. Run perform_ttest() first.")
        
        data = self.results['ttest'].copy()
        p_values = data['P_Value'].values
        
        if method == 'bonferroni':
            corrected_p = p_values * len(p_values)
            corrected_p = np.minimum(corrected_p, 1.0)  # Cap at 1.0
        elif method == 'holm':
            corrected_p = self._holm_correction(p_values)
        elif method == 'fdr_bh':
            corrected_p = self._fdr_correction(p_values)
        else:
            raise ValueError("method must be 'bonferroni', 'holm', or 'fdr_bh'")
        
        data[f'P_Value_{method.upper()}'] = corrected_p
        data[f'Significant_{method.upper()}'] = corrected_p < self.alpha
        
        self.results[f'ttest_{method}'] = data
        return data
    
    def _holm_correction(self, p_values):
        """Holm-Bonferroni correction."""
        sorted_indices = np.argsort(p_values)
        sorted_p = p_values[sorted_indices]
        n = len(p_values)
        
        corrected_p = np.zeros_like(p_values)
        for i, p in enumerate(sorted_p):
            corrected_p[sorted_indices[i]] = min(1.0, p * (n - i))
        
        # Ensure monotonicity
        for i in range(1, n):
            idx = sorted_indices[i]
            prev_idx = sorted_indices[i-1]
            corrected_p[idx] = max(corrected_p[idx], corrected_p[prev_idx])
        
        return corrected_p
    
    def _fdr_correction(self, p_values):
        """Benjamini-Hochberg FDR correction."""
        sorted_indices = np.argsort(p_values)
        sorted_p = p_values[sorted_indices]
        n = len(p_values)
        
        corrected_p = np.zeros_like(p_values)
        for i in range(n-1, -1, -1):
            corrected_p[sorted_indices[i]] = min(1.0, sorted_p[i] * n / (i + 1))
            if i < n - 1:
                corrected_p[sorted_indices[i]] = min(corrected_p[sorted_indices[i]], 
                                                   corrected_p[sorted_indices[i+1]])
        
        return corrected_p
    
    def _interpret_effect_size(self, cohens_d):
        """Interpret Cohen's d effect size."""
        if pd.isna(cohens_d):
            return "Unknown"
        
        abs_d = abs(cohens_d)
        if abs_d < 0.2:
            return "Negligible"
        elif abs_d < 0.5:
            return "Small"
        elif abs_d < 0.8:
            return "Medium"
        else:
            return "Large"
    
    def calculate_confidence_intervals(self, data, value_column='Fold_Change', 
                                     gene_column='Target_Gene', condition_column='Condition',
                                     confidence_level=0.95):
        """
        Calculate confidence intervals for fold change values.
        
        Parameters:
        -----------
        data : pandas.DataFrame
            Data containing fold change values
        value_column : str
            Column containing values for CI calculation
        gene_column : str
            Column containing gene identifiers
        condition_column : str
            Column containing experimental conditions
        confidence_level : float
            Confidence level (default: 0.95)
            
        Returns:
        --------
        pandas.DataFrame
            Confidence intervals for each gene-condition combination
        """
        alpha = 1 - confidence_level
        results = []
        
        groups = data.groupby([gene_column, condition_column])
        
        for (gene, condition), group in groups:
            values = group[value_column].dropna()
            
            if len(values) < 2:
                continue
                
            mean_val = values.mean()
            std_val = values.std()
            n = len(values)
            se = std_val / np.sqrt(n)
            
            # Calculate t-critical value
            df = n - 1
            t_critical = stats.t.ppf(1 - alpha/2, df)
            
            # Calculate confidence interval
            ci_lower = mean_val - t_critical * se
            ci_upper = mean_val + t_critical * se
            
            results.append({
                'Gene': gene,
                'Condition': condition,
                'Mean': round(mean_val, 3),
                'Std': round(std_val, 3),
                'N': n,
                'SE': round(se, 3),
                'CI_Lower': round(ci_lower, 3),
                'CI_Upper': round(ci_upper, 3),
                'CI_Width': round(ci_upper - ci_lower, 3),
                'Confidence_Level': confidence_level
            })
        
        return pd.DataFrame(results)
    
    def perform_comprehensive_analysis(self, data, control_condition='Control'):
        """
        Perform complete statistical analysis pipeline.
        
        Parameters:
        -----------
        data : pandas.DataFrame
            qPCR analysis results with fold change data
        control_condition : str
            Name of the control condition
            
        Returns:
        --------
        dict
            Complete statistical analysis results
        """
        print("📊 Running Comprehensive Statistical Analysis")
        print("=" * 50)
        
        # Perform t-tests
        ttest_results = self.perform_ttest(data, control_condition=control_condition)
        print(f"✓ Completed t-tests for {len(ttest_results)} comparisons")
        
        # Apply multiple comparison corrections
        bonferroni_results = self.apply_multiple_comparison_correction('bonferroni')
        fdr_results = self.apply_multiple_comparison_correction('fdr_bh')
        print("✓ Applied multiple comparison corrections")
        
        # Calculate confidence intervals
        ci_results = self.calculate_confidence_intervals(data)
        print(f"✓ Calculated confidence intervals for {len(ci_results)} groups")
        
        # Generate summary
        significant_uncorrected = (ttest_results['P_Value'] < self.alpha).sum()
        significant_bonferroni = (bonferroni_results['P_Value_BONFERRONI'] < self.alpha).sum()
        significant_fdr = (fdr_results['P_Value_FDR_BH'] < self.alpha).sum()
        
        print(f"\n📈 Statistical Summary:")
        print(f"  Significant (uncorrected): {significant_uncorrected}/{len(ttest_results)}")
        print(f"  Significant (Bonferroni): {significant_bonferroni}/{len(ttest_results)}")
        print(f"  Significant (FDR): {significant_fdr}/{len(ttest_results)}")
        
        return {
            'ttest_results': ttest_results,
            'bonferroni_corrected': bonferroni_results,
            'fdr_corrected': fdr_results,
            'confidence_intervals': ci_results,
            'summary': {
                'total_comparisons': len(ttest_results),
                'significant_uncorrected': significant_uncorrected,
                'significant_bonferroni': significant_bonferroni,
                'significant_fdr': significant_fdr,
                'alpha_level': self.alpha
            }
        }