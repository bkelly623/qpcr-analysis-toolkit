#!/usr/bin/env python3
"""
Quick test of the qPCR analyzer functionality
"""

import sys
sys.path.append('src')

from qpcr_analyzer import qPCRAnalyzer

def main():
    print("🧬 Testing qPCR Analyzer")
    print("=" * 40)
    
    # Initialize analyzer
    analyzer = qPCRAnalyzer(reference_gene='GAPDH', control_condition='Control')
    
    # Run complete analysis
    try:
        results = analyzer.run_complete_analysis('data/sample_experiment.csv')
        
        print("\n📊 Sample Results:")
        print("-" * 30)
        
        # Show some key results
        print("\n🔬 Summary Statistics:")
        print(results['summary_statistics'])
        
        print("\n🧪 Sample Fold Changes:")
        sample_results = results['delta_delta_ct_data'][['Target_Gene', 'Condition', 'Fold_Change']].head(10)
        print(sample_results)
        
        print("\n✅ Core analyzer working correctly!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return False
    
    return True

if __name__ == "__main__":
    main()