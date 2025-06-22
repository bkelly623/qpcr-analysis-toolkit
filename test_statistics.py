#!/usr/bin/env python3
"""
Test statistical analysis functionality
"""

import sys
sys.path.append('src')

from qpcr_analyzer import qPCRAnalyzer
from statistical_analysis import qPCRStatistics

def main():
    print("🧬 Testing qPCR Statistical Analysis")
    print("=" * 50)
    
    # Step 1: Run qPCR analysis
    analyzer = qPCRAnalyzer(reference_gene='GAPDH', control_condition='Control')
    qpcr_results = analyzer.run_complete_analysis('data/sample_experiment.csv')
    
    print("\n📊 Running Statistical Tests...")
    print("-" * 30)
    
    # Step 2: Initialize statistical analyzer
    stats_analyzer = qPCRStatistics(alpha=0.05)
    
    # Step 3: Run comprehensive statistical analysis
    fold_change_data = qpcr_results['delta_delta_ct_data']
    stats_results = stats_analyzer.perform_comprehensive_analysis(fold_change_data)
    
    # Step 4: Display key results
    print("\n🔬 Statistical Test Results:")
    print(stats_results['ttest_results'][['Gene', 'Comparison', 'Treatment_Mean', 'P_Value', 'Significant', 'Cohens_D', 'Effect_Size_Interpretation']])
    
    print("\n📈 Multiple Comparison Corrections:")
    bonferroni = stats_results['bonferroni_corrected'][['Gene', 'Comparison', 'P_Value', 'P_Value_BONFERRONI', 'Significant_BONFERRONI']]
    print(bonferroni)
    
    print("\n📊 Confidence Intervals:")
    ci_results = stats_results['confidence_intervals']
    print(ci_results[['Gene', 'Condition', 'Mean', 'CI_Lower', 'CI_Upper']])
    
    print("\n✅ Statistical analysis complete!")
    
    return stats_results

if __name__ == "__main__":
    main()