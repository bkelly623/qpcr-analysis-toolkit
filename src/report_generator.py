"""
Automated Report Generator for qPCR Analysis

This module creates comprehensive analysis reports in multiple formats
including HTML summaries and Excel exports with embedded statistics.

Author: Your Name
Date: June 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
import json


class qPCRReportGenerator:
    """
    Generates comprehensive analysis reports from qPCR results.
    
    Features:
    - HTML summary reports with embedded plots
    - Excel exports with multiple sheets
    - JSON data exports for further analysis
    - Executive summary generation
    """
    
    def __init__(self, output_dir='output/reports'):
        """
        Initialize the report generator.
        
        Parameters:
        -----------
        output_dir : str
            Directory to save reports (default: 'output/reports')
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def generate_html_report(self, qpcr_results, stats_results, save_path=None):
        """
        Generate an HTML summary report with embedded results.
        
        Parameters:
        -----------
        qpcr_results : dict
            Results from qPCR analysis
        stats_results : dict
            Results from statistical analysis
        save_path : str, optional
            Path to save the HTML report
            
        Returns:
        --------
        str
            Path to the generated HTML file
        """
        
        # Generate HTML content
        html_content = self._create_html_template(qpcr_results, stats_results)
        
        # Save HTML file
        if save_path is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            save_path = self.output_dir / f"qpcr_analysis_report_{timestamp}.html"
        
        with open(save_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"✓ HTML report saved: {save_path}")
        return str(save_path)
    
    def generate_excel_report(self, qpcr_results, stats_results, save_path=None):
        """
        Generate Excel report with multiple sheets containing analysis results.
        
        Parameters:
        -----------
        qpcr_results : dict
            Results from qPCR analysis
        stats_results : dict
            Results from statistical analysis
        save_path : str, optional
            Path to save the Excel file
            
        Returns:
        --------
        str
            Path to the generated Excel file
        """
        
        if save_path is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            save_path = self.output_dir / f"qpcr_analysis_data_{timestamp}.xlsx"
        
        with pd.ExcelWriter(save_path, engine='openpyxl') as writer:
            # Sheet 1: Executive Summary
            summary_df = self._create_summary_dataframe(qpcr_results, stats_results)
            summary_df.to_excel(writer, sheet_name='Executive_Summary', index=False)
            
            # Sheet 2: Statistical Results
            stats_results['ttest_results'].to_excel(writer, sheet_name='Statistical_Tests', index=False)
            
            # Sheet 3: Fold Change Data
            qpcr_results['delta_delta_ct_data'].to_excel(writer, sheet_name='Fold_Changes', index=False)
            
            # Sheet 4: Raw Data
            qpcr_results['raw_data'].to_excel(writer, sheet_name='Raw_Data', index=False)
            
            # Sheet 5: Quality Control
            if 'confidence_intervals' in stats_results:
                stats_results['confidence_intervals'].to_excel(writer, sheet_name='Confidence_Intervals', index=False)
        
        print(f"✓ Excel report saved: {save_path}")
        return str(save_path)
    
    def generate_json_export(self, qpcr_results, stats_results, save_path=None):
        """
        Export analysis results to JSON format for further processing.
        
        Parameters:
        -----------
        qpcr_results : dict
            Results from qPCR analysis
        stats_results : dict
            Results from statistical analysis
        save_path : str, optional
            Path to save the JSON file
            
        Returns:
        --------
        str
            Path to the generated JSON file
        """
        
        # Convert dataframes to dictionaries for JSON serialization
        export_data = {
            'metadata': {
                'analysis_date': datetime.now().isoformat(),
                'toolkit_version': '1.0.0',
                'analysis_type': 'qPCR_ddCt_analysis'
            },
            'experimental_design': {
                'total_samples': len(qpcr_results['raw_data']['Sample_ID'].unique()),
                'target_genes': list(qpcr_results['delta_delta_ct_data']['Target_Gene'].unique()),
                'conditions': list(qpcr_results['raw_data']['Condition'].unique()),
                'biological_replicates': int(qpcr_results['raw_data'].groupby(['Condition', 'Gene']).size().max())
            },
            'statistical_summary': stats_results['summary'],
            'fold_changes': qpcr_results['delta_delta_ct_data'].to_dict('records'),
            'statistical_tests': stats_results['ttest_results'].to_dict('records')
        }
        
        if save_path is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            save_path = self.output_dir / f"qpcr_analysis_export_{timestamp}.json"
        
        with open(save_path, 'w', encoding='utf-8') as f:
            json.dump(export_data, f, indent=2, default=str)
        
        print(f"✓ JSON export saved: {save_path}")
        return str(save_path)
    
    def _create_html_template(self, qpcr_results, stats_results):
        """Create HTML template with embedded results."""
        
        # Get key statistics
        total_genes = len(qpcr_results['delta_delta_ct_data']['Target_Gene'].unique())
        total_samples = len(qpcr_results['raw_data']['Sample_ID'].unique())
        significant_genes = stats_results['summary']['significant_fdr']
        
        # Create statistical results table
        stats_table = stats_results['ttest_results'].to_html(classes='table table-striped', index=False)
        
        html_template = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>qPCR Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
        .header {{ background-color: #2E86AB; color: white; padding: 20px; border-radius: 8px; }}
        .summary {{ background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin: 20px 0; }}
        .section {{ margin: 30px 0; }}
        .highlight {{ background-color: #e7f3ff; padding: 15px; border-left: 4px solid #2E86AB; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        .significant {{ color: #28a745; font-weight: bold; }}
        .not-significant {{ color: #dc3545; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 qPCR Analysis Report</h1>
        <p>Generated on {datetime.now().strftime("%B %d, %Y at %I:%M %p")}</p>
    </div>
    
    <div class="summary">
        <h2>📊 Executive Summary</h2>
        <ul>
            <li><strong>Total Genes Analyzed:</strong> {total_genes}</li>
            <li><strong>Total Samples:</strong> {total_samples}</li>
            <li><strong>Statistically Significant Genes:</strong> {significant_genes}/{total_genes}</li>
            <li><strong>Analysis Method:</strong> ΔΔCt relative quantification</li>
            <li><strong>Statistical Correction:</strong> FDR (False Discovery Rate)</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>🔬 Key Findings</h2>
        <div class="highlight">
            <h3>Biological Interpretation:</h3>
            {self._generate_biological_interpretation(stats_results)}
        </div>
    </div>
    
    <div class="section">
        <h2>📈 Statistical Results</h2>
        {stats_table}
    </div>
    
    <div class="section">
        <h2>📁 Generated Files</h2>
        <ul>
            <li><strong>Fold Change Plot:</strong> output/plots/fold_change_bars.png</li>
            <li><strong>Ct Distributions:</strong> output/plots/ct_distributions.png</li>
            <li><strong>ΔCt Heatmap:</strong> output/plots/delta_ct_heatmap.png</li>
            <li><strong>QC Dashboard:</strong> output/plots/qc_dashboard.png</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>🔧 Methodology</h2>
        <p><strong>ΔΔCt Method:</strong> This analysis uses the industry-standard ΔΔCt method for relative gene expression quantification:</p>
        <ol>
            <li><strong>ΔCt = Target Ct - Reference Ct</strong> (normalizes against housekeeping gene)</li>
            <li><strong>ΔΔCt = Treatment ΔCt - Control ΔCt</strong> (compares to baseline)</li>
            <li><strong>Fold Change = 2^(-ΔΔCt)</strong> (converts to linear scale)</li>
        </ol>
        <p><strong>Statistical Testing:</strong> Student's t-tests with False Discovery Rate (FDR) correction for multiple comparisons.</p>
    </div>
    
    <footer style="margin-top: 50px; padding-top: 20px; border-top: 1px solid #ddd; color: #666;">
        <p>Generated by qPCR Analysis Toolkit v1.0.0 | For research use only</p>
    </footer>
</body>
</html>
        """
        
        return html_template
    
    def _generate_biological_interpretation(self, stats_results):
        """Generate biological interpretation text."""
        
        significant_results = stats_results['ttest_results'][stats_results['ttest_results']['Significant']]
        
        if len(significant_results) == 0:
            return "<p>No statistically significant changes detected in gene expression.</p>"
        
        interpretation = "<ul>"
        for _, row in significant_results.iterrows():
            gene = row['Gene']
            fold_change = row['Treatment_Mean']
            p_value = row['P_Value']
            
            direction = "upregulated" if fold_change > 1 else "downregulated"
            interpretation += f"<li><strong>{gene}</strong>: {fold_change:.1f}-fold {direction} (p={p_value:.3f})</li>"
        
        interpretation += "</ul>"
        
        return interpretation
    
    def _create_summary_dataframe(self, qpcr_results, stats_results):
        """Create executive summary dataframe."""
        
        summary_data = []
        
        # Add general information
        summary_data.append(['Analysis Date', datetime.now().strftime("%Y-%m-%d %H:%M:%S")])
        summary_data.append(['Total Samples', len(qpcr_results['raw_data']['Sample_ID'].unique())])
        summary_data.append(['Target Genes', len(qpcr_results['delta_delta_ct_data']['Target_Gene'].unique())])
        summary_data.append(['Statistical Method', 'Student\'s t-test with FDR correction'])
        summary_data.append(['Significance Level', '0.05'])
        
        # Add gene-specific results
        for _, row in stats_results['ttest_results'].iterrows():
            gene = row['Gene']
            fold_change = row['Treatment_Mean']
            p_value = row['P_Value']
            significant = 'Yes' if row['Significant'] else 'No'
            
            summary_data.append([f'{gene} Fold Change', f'{fold_change:.2f}'])
            summary_data.append([f'{gene} P-value', f'{p_value:.4f}'])
            summary_data.append([f'{gene} Significant', significant])
        
        return pd.DataFrame(summary_data, columns=['Parameter', 'Value'])
    
    def generate_complete_report_package(self, qpcr_results, stats_results):
        """
        Generate complete report package with all formats.
        
        Parameters:
        -----------
        qpcr_results : dict
            Results from qPCR analysis
        stats_results : dict
            Results from statistical analysis
            
        Returns:
        --------
        dict
            Paths to all generated report files
        """
        
        print("📋 Generating Complete Report Package")
        print("=" * 40)
        
        report_paths = {}
        
        # Generate HTML report
        report_paths['html'] = self.generate_html_report(qpcr_results, stats_results)
        
        # Generate Excel report
        report_paths['excel'] = self.generate_excel_report(qpcr_results, stats_results)
        
        # Generate JSON export
        report_paths['json'] = self.generate_json_export(qpcr_results, stats_results)
        
        print(f"\n✅ Complete report package generated!")
        print(f"📁 All reports saved to: {self.output_dir}")
        
        return report_paths