#!/usr/bin/env python3
"""
Complete qPCR Analysis Toolkit Demonstration
Shows all features: Analysis + Statistics + Visualization + Reporting
"""

import sys
sys.path.append('src')

from qpcr_analyzer import qPCRAnalyzer
from statistical_analysis import qPCRStatistics
from visualization import qPCRVisualizer
from report_generator import qPCRReportGenerator

def main():
    print("🧬 COMPLETE qPCR ANALYSIS TOOLKIT DEMONSTRATION")
    print("=" * 70)
    print("🎯 This script demonstrates all toolkit capabilities:")
    print("   • qPCR data processing with ΔΔCt methodology")
    print("   • Statistical testing with multiple comparison corrections")
    print("   • Publication-ready visualizations")
    print("   • Automated report generation (HTML, Excel, JSON)")
    
    # Step 1: Core qPCR Analysis
    print("\n" + "="*50)
    print("📊 STEP 1: qPCR DATA ANALYSIS")
    print("="*50)
    
    analyzer = qPCRAnalyzer(reference_gene='GAPDH', control_condition='Control')
    qpcr_results = analyzer.run_complete_analysis('data/sample_experiment.csv')
    
    # Step 2: Statistical Analysis
    print("\n" + "="*50)
    print("📈 STEP 2: STATISTICAL VALIDATION")
    print("="*50)
    
    stats_analyzer = qPCRStatistics(alpha=0.05)
    stats_results = stats_analyzer.perform_comprehensive_analysis(
        qpcr_results['delta_delta_ct_data']
    )
    
    # Step 3: Visualization
    print("\n" + "="*50)
    print("🎨 STEP 3: PUBLICATION-READY PLOTS")
    print("="*50)
    
    visualizer = qPCRVisualizer()
    figures = visualizer.create_comprehensive_report(qpcr_results, stats_results)
    
    # Step 4: Report Generation
    print("\n" + "="*50)
    print("📋 STEP 4: AUTOMATED REPORTING")
    print("="*50)
    
    reporter = qPCRReportGenerator()
    report_paths = reporter.generate_complete_report_package(qpcr_results, stats_results)
    
    # Final Summary
    print("\n" + "="*70)
    print("🏆 TOOLKIT DEMONSTRATION COMPLETE!")
    print("="*70)
    
    print("\n📊 ANALYSIS SUMMARY:")
    print(f"  • Samples processed: {len(qpcr_results['raw_data']['Sample_ID'].unique())}")
    print(f"  • Genes analyzed: {len(qpcr_results['delta_delta_ct_data']['Target_Gene'].unique())}")
    print(f"  • Statistical tests: {len(stats_results['ttest_results'])}")
    print(f"  • Plots generated: {len(figures)}")
    print(f"  • Reports created: {len(report_paths)}")
    
    print("\n🔬 KEY BIOLOGICAL FINDINGS:")
    for _, row in stats_results['ttest_results'].iterrows():
        gene = row['Gene']
        fold_change = row['Treatment_Mean']
        p_value = row['P_Value']
        significant = "✓ SIGNIFICANT" if row['Significant'] else "✗ Not significant"
        print(f"  • {gene}: {fold_change:.1f}-fold upregulation (p={p_value:.3f}) {significant}")
    
    print("\n📁 OUTPUT FILES GENERATED:")
    print("  📈 Visualizations:")
    print("    • output/plots/fold_change_bars.png")
    print("    • output/plots/ct_distributions.png")
    print("    • output/plots/delta_ct_heatmap.png")
    print("    • output/plots/qc_dashboard.png")
    
    print("  📋 Reports:")
    for format_type, path in report_paths.items():
        print(f"    • {path}")
    
    print("\n🎯 PORTFOLIO IMPACT:")
    print("  ✅ Demonstrates molecular biology expertise")
    print("  ✅ Shows statistical analysis proficiency")
    print("  ✅ Creates publication-ready visualizations")
    print("  ✅ Generates professional reports")
    print("  ✅ Follows industry-standard methodologies")
    
    print("\n🚀 READY FOR:")
    print("  • Biotech job applications")
    print("  • Academic publications")
    print("  • Regulatory submissions")
    print("  • Client presentations")
    print("  • Portfolio showcase")
    
    return {
        'qpcr_results': qpcr_results,
        'stats_results': stats_results,
        'figures': figures,
        'reports': report_paths
    }

if __name__ == "__main__":
    results = main()