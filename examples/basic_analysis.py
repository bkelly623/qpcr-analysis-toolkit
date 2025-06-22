#!/usr/bin/env python3
"""
Basic qPCR Analysis Example

This script demonstrates how to use the qPCR Analysis Toolkit for 
analyzing gene expression data from raw Ct values to publication-ready results.

Usage:
    python examples/basic_analysis.py

Requirements:
    - Virtual environment activated: source venv/bin/activate
    - Dependencies installed: pip install -r requirements.txt
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from qpcr_analyzer import qPCRAnalyzer
from statistical_analysis import qPCRStatistics
from visualization import qPCRVisualizer


def analyze_inflammatory_response():
    """
    Example analysis of inflammatory gene expression.
    
    This example demonstrates:
    1. Loading qPCR data from CSV
    2. Calculating ΔΔCt and fold changes
    3. Statistical significance testing
    4. Creating publication-ready plots
    """
    
    print("🧬 qPCR Analysis Toolkit - Basic Example")
    print("=" * 50)
    print("Analyzing inflammatory gene expression data...")
    
    # Step 1: Initialize the qPCR analyzer
    print("\n📋 Step 1: Setting up analysis parameters")
    analyzer = qPCRAnalyzer(
        reference_gene='GAPDH',        # Housekeeping gene for normalization
        control_condition='Control'     # Baseline condition for comparison
    )
    
    # Step 2: Load and process qPCR data
    print("\n🔬 Step 2: Loading and processing qPCR data")
    data_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'sample_experiment.csv')
    qpcr_results = analyzer.run_complete_analysis(data_file)
    
    # Step 3: Perform statistical analysis
    print("\n📊 Step 3: Statistical testing and validation")
    stats_analyzer = qPCRStatistics(alpha=0.05)
    stats_results = stats_analyzer.perform_comprehensive_analysis(
        qpcr_results['delta_delta_ct_data']
    )
    
    # Step 4: Create visualizations
    print("\n🎨 Step 4: Generating publication-ready plots")
    visualizer = qPCRVisualizer()
    figures = visualizer.create_comprehensive_report(qpcr_results, stats_results)
    
    # Step 5: Interpret results
    print("\n📋 Step 5: Results Summary")
    print("=" * 30)
    
    print("\n🔬 Experimental Design:")
    print(f"  • Reference Gene: {analyzer.reference_gene}")
    print(f"  • Control Condition: {analyzer.control_condition}")
    print(f"  • Total Samples: {len(qpcr_results['delta_delta_ct_data']['Sample_ID'].unique())}")
    print(f"  • Target Genes: {len(qpcr_results['delta_delta_ct_data']['Target_Gene'].unique())}")
    
    print("\n📈 Key Findings:")
    for _, row in stats_results['ttest_results'].iterrows():
        gene = row['Gene']
        fold_change = row['Treatment_Mean']
        p_value = row['P_Value']
        effect_size = row['Effect_Size_Interpretation']
        significant = "✓" if row['Significant'] else "✗"
        
        print(f"  {significant} {gene}: {fold_change:.1f}-fold change (p={p_value:.3f}, {effect_size} effect)")
    
    print("\n📊 Statistical Validation:")
    summary = stats_results['summary']
    print(f"  • Comparisons tested: {summary['total_comparisons']}")
    print(f"  • Significant (raw p-values): {summary['significant_uncorrected']}/{summary['total_comparisons']}")
    print(f"  • Significant (Bonferroni): {summary['significant_bonferroni']}/{summary['total_comparisons']}")
    print(f"  • Significant (FDR): {summary['significant_fdr']}/{summary['total_comparisons']}")
    
    print("\n🎯 Clinical Interpretation:")
    if summary['significant_fdr'] > 0:
        print("  • Statistically significant inflammatory response detected")
        print("  • Results are robust to multiple comparison correction")
        print("  • Large effect sizes suggest biologically meaningful changes")
    else:
        print("  • No significant changes detected after correction")
        print("  • May require larger sample sizes or different conditions")
    
    print("\n📁 Output Files:")
    print("  • Plots saved to: output/plots/")
    print("  • fold_change_bars.png - Main results visualization")
    print("  • ct_distributions.png - Data quality assessment")
    print("  • delta_ct_heatmap.png - Sample-wise expression patterns")
    print("  • qc_dashboard.png - Quality control metrics")
    
    print("\n✅ Analysis Complete!")
    print("📖 Tip: Open the generated plots to see publication-ready visualizations")
    
    return qpcr_results, stats_results, figures


def custom_analysis_example():
    """
    Example showing how to customize the analysis for different experimental designs.
    """
    
    print("\n" + "="*60)
    print("🔧 CUSTOM ANALYSIS EXAMPLE")
    print("="*60)
    
    # Example: Different reference gene and conditions
    custom_analyzer = qPCRAnalyzer(
        reference_gene='ACTB',          # Different housekeeping gene
        control_condition='Baseline'     # Different control name
    )
    
    print("\n📋 Custom Parameters:")
    print(f"  • Reference Gene: {custom_analyzer.reference_gene}")
    print(f"  • Control Condition: {custom_analyzer.control_condition}")
    
    # Example: Different statistical parameters
    custom_stats = qPCRStatistics(alpha=0.01)  # More stringent significance level
    
    print(f"  • Significance Level: {custom_stats.alpha}")
    
    print("\n💡 Customization Options:")
    print("  • Change reference genes (GAPDH, ACTB, 18S, etc.)")
    print("  • Modify significance thresholds (0.05, 0.01, 0.001)")
    print("  • Use different statistical corrections (Bonferroni, FDR, Holm)")
    print("  • Customize plot colors and styles")
    print("  • Add multiple treatment conditions")


if __name__ == "__main__":
    # Run the basic analysis
    qpcr_results, stats_results, figures = analyze_inflammatory_response()
    
    # Show customization options
    custom_analysis_example()
    
    print("\n🎓 Next Steps:")
    print("  1. Try modifying the sample data with your own experiments")
    print("  2. Explore different reference genes and conditions")
    print("  3. Adjust statistical parameters for your requirements")
    print("  4. Customize visualizations for your publications")
    print("  5. Use this toolkit for your biotech portfolio projects!")