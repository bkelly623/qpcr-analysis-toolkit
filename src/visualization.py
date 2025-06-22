"""
Visualization Module for qPCR Data

This module creates publication-ready visualizations for qPCR analysis results,
including fold change plots, statistical annotations, and quality control charts.

Author: Your Name
Date: June 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings

# Set publication-quality defaults
plt.rcParams.update({
    'figure.figsize': (10, 6),
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

class qPCRVisualizer:
    """
    Visualization class for qPCR analysis results.
    
    Creates publication-ready plots including:
    - Fold change bar charts with error bars
    - Statistical significance annotations
    - Ct value distributions
    - Heatmaps for multi-gene comparisons
    - Quality control plots
    """
    
    def __init__(self, output_dir='output/plots'):
        """
        Initialize the visualizer.
        
        Parameters:
        -----------
        output_dir : str
            Directory to save plots (default: 'output/plots')
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Color palette for conditions
        self.colors = {
            'Control': '#2E86AB',
            'Treatment': '#A23B72',
            'Reference': '#F18F01'
        }
        
    def plot_fold_change_bars(self, fold_change_data, stats_data=None, 
                             save_path=None, show_plot=True):
        """
        Create bar chart of fold changes with error bars and statistical annotations.
        
        Parameters:
        -----------
        fold_change_data : pandas.DataFrame
            Data containing fold change values
        stats_data : pandas.DataFrame, optional
            Statistical test results for significance annotations
        save_path : str, optional
            Path to save the plot
        show_plot : bool
            Whether to display the plot
            
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure object
        """
        # Calculate summary statistics
        summary = fold_change_data.groupby(['Target_Gene', 'Condition']).agg({
            'Fold_Change': ['mean', 'std', 'sem']
        }).reset_index()
        
        # Flatten column names
        summary.columns = ['Gene', 'Condition', 'Mean_FC', 'Std_FC', 'SEM_FC']
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        genes = summary['Gene'].unique()
        conditions = summary['Condition'].unique()
        
        x = np.arange(len(genes))
        width = 0.35
        
        # Plot bars for each condition
        for i, condition in enumerate(conditions):
            condition_data = summary[summary['Condition'] == condition]
            
            # Ensure genes are in the same order
            means = []
            errors = []
            for gene in genes:
                gene_data = condition_data[condition_data['Gene'] == gene]
                if len(gene_data) > 0:
                    means.append(gene_data['Mean_FC'].iloc[0])
                    errors.append(gene_data['SEM_FC'].iloc[0])
                else:
                    means.append(0)
                    errors.append(0)
            
            bars = ax.bar(x + i*width, means, width, 
                         yerr=errors, capsize=5,
                         label=condition,
                         color=self.colors.get(condition, f'C{i}'),
                         alpha=0.8,
                         edgecolor='black', linewidth=1)
        
        # Add statistical significance annotations
        if stats_data is not None:
            self._add_significance_annotations(ax, stats_data, genes, x, width)
        
        # Customize the plot
        ax.set_xlabel('Gene', fontweight='bold')
        ax.set_ylabel('Fold Change (2^-ΔΔCt)', fontweight='bold')
        ax.set_title('Gene Expression Fold Changes\n(Relative to Control)', 
                    fontweight='bold', pad=20)
        ax.set_xticks(x + width/2)
        ax.set_xticklabels(genes)
        ax.legend(title='Condition', frameon=True, shadow=True)
        
        # Add horizontal line at y=1 (no change)
        ax.axhline(y=1, color='gray', linestyle='--', alpha=0.7, linewidth=1)
        
        # Set y-axis to log scale if fold changes are large
        max_fc = summary['Mean_FC'].max()
        if max_fc > 10:
            ax.set_yscale('log')
            ax.set_ylabel('Fold Change (2^-ΔΔCt) [Log Scale]', fontweight='bold')
        
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save plot
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved fold change plot: {save_path}")
        else:
            default_path = self.output_dir / 'fold_change_bars.png'
            plt.savefig(default_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved fold change plot: {default_path}")
        
        if show_plot:
            plt.show()
        
        return fig
    
    def plot_ct_distributions(self, raw_data, save_path=None, show_plot=True):
        """
        Create box plots showing Ct value distributions by gene and condition.
        
        Parameters:
        -----------
        raw_data : pandas.DataFrame
            Raw qPCR data with Ct values
        save_path : str, optional
            Path to save the plot
        show_plot : bool
            Whether to display the plot
            
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure object
        """
        fig, ax = plt.subplots(figsize=(14, 8))
        
        # Create box plot
        sns.boxplot(data=raw_data, x='Gene', y='Ct_Value', hue='Condition',
                   palette=[self.colors.get(c, f'C{i}') for i, c in enumerate(raw_data['Condition'].unique())],
                   ax=ax)
        
        # Customize the plot
        ax.set_xlabel('Gene', fontweight='bold')
        ax.set_ylabel('Ct Value (PCR Cycles)', fontweight='bold')
        ax.set_title('Raw Ct Value Distributions\n(Lower Ct = Higher Expression)', 
                    fontweight='bold', pad=20)
        
        # Invert y-axis (lower Ct = higher expression)
        ax.invert_yaxis()
        
        # Add grid
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save plot
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved Ct distribution plot: {save_path}")
        else:
            default_path = self.output_dir / 'ct_distributions.png'
            plt.savefig(default_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved Ct distribution plot: {default_path}")
        
        if show_plot:
            plt.show()
        
        return fig
    
    def plot_delta_ct_heatmap(self, delta_ct_data, save_path=None, show_plot=True):
        """
        Create heatmap of ΔCt values across samples and genes.
        
        Parameters:
        -----------
        delta_ct_data : pandas.DataFrame
            Data containing ΔCt values
        save_path : str, optional
            Path to save the plot
        show_plot : bool
            Whether to display the plot
            
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure object
        """
        # Pivot data for heatmap
        heatmap_data = delta_ct_data.pivot_table(
            index=['Sample_ID', 'Condition'], 
            columns='Target_Gene', 
            values='Delta_Ct',
            fill_value=np.nan
        )
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create heatmap
        sns.heatmap(heatmap_data, 
                   annot=True, 
                   fmt='.2f',
                   cmap='RdYlBu_r',
                   center=0,
                   cbar_kws={'label': 'ΔCt Value'},
                   ax=ax)
        
        # Customize the plot
        ax.set_title('ΔCt Values Heatmap\n(Target - Reference Gene)', 
                    fontweight='bold', pad=20)
        ax.set_xlabel('Target Gene', fontweight='bold')
        ax.set_ylabel('Sample (Condition)', fontweight='bold')
        
        plt.tight_layout()
        
        # Save plot
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved ΔCt heatmap: {save_path}")
        else:
            default_path = self.output_dir / 'delta_ct_heatmap.png'
            plt.savefig(default_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved ΔCt heatmap: {default_path}")
        
        if show_plot:
            plt.show()
        
        return fig
    
    def plot_qc_dashboard(self, raw_data, processed_data=None, save_path=None, show_plot=True):
        """
        Create quality control dashboard with multiple diagnostic plots.
        
        Parameters:
        -----------
        raw_data : pandas.DataFrame
            Raw qPCR data
        processed_data : pandas.DataFrame, optional
            Processed data with technical replicate means
        save_path : str, optional
            Path to save the plot
        show_plot : bool
            Whether to display the plot
            
        Returns:
        --------
        matplotlib.figure.Figure
            The created figure object
        """
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('qPCR Quality Control Dashboard', fontsize=18, fontweight='bold')
        
        # Plot 1: Ct value ranges by gene
        ax1 = axes[0, 0]
        raw_data.boxplot(column='Ct_Value', by='Gene', ax=ax1)
        ax1.set_title('Ct Value Ranges by Gene')
        ax1.set_xlabel('Gene')
        ax1.set_ylabel('Ct Value')
        ax1.grid(True, alpha=0.3)
        
        # Plot 2: Technical replicate correlation
        ax2 = axes[0, 1]
        if 'Technical_Rep' in raw_data.columns:
            tech_rep_data = raw_data.pivot_table(
                index=['Sample_ID', 'Gene', 'Condition'],
                columns='Technical_Rep',
                values='Ct_Value'
            ).reset_index()
            
            if len(tech_rep_data.columns) >= 3:  # Has at least 2 technical reps
                tech_rep_data.plot.scatter(x=1, y=2, ax=ax2, alpha=0.7)
                ax2.set_xlabel('Technical Rep 1 (Ct)')
                ax2.set_ylabel('Technical Rep 2 (Ct)')
                ax2.set_title('Technical Replicate Correlation')
                
                # Add correlation line
                x_min, x_max = ax2.get_xlim()
                ax2.plot([x_min, x_max], [x_min, x_max], 'r--', alpha=0.8)
        
        # Plot 3: Sample-wise Ct distributions
        ax3 = axes[1, 0]
        sample_means = raw_data.groupby('Sample_ID')['Ct_Value'].mean()
        ax3.hist(sample_means, bins=15, alpha=0.7, edgecolor='black')
        ax3.set_xlabel('Mean Ct Value per Sample')
        ax3.set_ylabel('Frequency')
        ax3.set_title('Sample Quality Distribution')
        ax3.grid(True, alpha=0.3)
        
        # Plot 4: Failed reactions summary
        ax4 = axes[1, 1]
        failed_reactions = raw_data[
            (raw_data['Ct_Value'].isna()) | 
            (raw_data['Ct_Value'] <= 0) | 
            (raw_data['Ct_Value'] > 40)
        ]
        
        if len(failed_reactions) > 0:
            failure_counts = failed_reactions.groupby('Gene').size()
            failure_counts.plot(kind='bar', ax=ax4, color='red', alpha=0.7)
            ax4.set_title('Failed Reactions by Gene')
            ax4.set_xlabel('Gene')
            ax4.set_ylabel('Number of Failed Reactions')
        else:
            ax4.text(0.5, 0.5, 'No Failed Reactions\n✓ Excellent Data Quality', 
                    transform=ax4.transAxes, ha='center', va='center',
                    fontsize=14, fontweight='bold', color='green')
            ax4.set_title('Failed Reactions Summary')
        
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved QC dashboard: {save_path}")
        else:
            default_path = self.output_dir / 'qc_dashboard.png'
            plt.savefig(default_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved QC dashboard: {default_path}")
        
        if show_plot:
            plt.show()
        
        return fig
    
    def _add_significance_annotations(self, ax, stats_data, genes, x_positions, width):
        """Add statistical significance annotations to bar plot."""
        for i, gene in enumerate(genes):
            gene_stats = stats_data[stats_data['Gene'] == gene]
            
            if len(gene_stats) > 0:
                p_value = gene_stats['P_Value'].iloc[0]
                
                # Determine significance level
                if p_value < 0.001:
                    sig_text = '***'
                elif p_value < 0.01:
                    sig_text = '**'
                elif p_value < 0.05:
                    sig_text = '*'
                else:
                    sig_text = 'ns'
                
                # Get the height for annotation
                bars_at_position = [bar for bar in ax.patches if abs(bar.get_x() - (x_positions[i] + width/2)) < width]
                if bars_at_position:
                    max_height = max([bar.get_height() + bar.get_yerr() if hasattr(bar, 'get_yerr') else bar.get_height() for bar in bars_at_position])
                    
                    # Add annotation
                    ax.annotate(sig_text, 
                              xy=(x_positions[i] + width/2, max_height),
                              xytext=(0, 5), textcoords='offset points',
                              ha='center', va='bottom',
                              fontweight='bold', fontsize=14)
    
    def create_comprehensive_report(self, qpcr_results, stats_results):
        """
        Create a comprehensive visualization report with all plots.
        
        Parameters:
        -----------
        qpcr_results : dict
            Results from qPCR analysis
        stats_results : dict
            Results from statistical analysis
            
        Returns:
        --------
        dict
            Dictionary of created figure objects
        """
        print("🎨 Creating Comprehensive Visualization Report")
        print("=" * 50)
        
        figures = {}
        
        # 1. Fold change bar plot with statistics
        figures['fold_change'] = self.plot_fold_change_bars(
            qpcr_results['delta_delta_ct_data'],
            stats_results['ttest_results'],
            show_plot=False
        )
        
        # 2. Ct distribution plots
        figures['ct_distributions'] = self.plot_ct_distributions(
            qpcr_results['raw_data'],
            show_plot=False
        )
        
        # 3. ΔCt heatmap
        figures['delta_ct_heatmap'] = self.plot_delta_ct_heatmap(
            qpcr_results['delta_ct_data'],
            show_plot=False
        )
        
        # 4. Quality control dashboard
        figures['qc_dashboard'] = self.plot_qc_dashboard(
            qpcr_results['raw_data'],
            qpcr_results['processed_data'],
            show_plot=False
        )
        
        print(f"✅ Created {len(figures)} publication-ready plots")
        print(f"📁 All plots saved to: {self.output_dir}")
        
        return figures