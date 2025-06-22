"""
qPCR Data Analysis Module

This module provides comprehensive analysis of qPCR (quantitative PCR) data,
including ΔCt and ΔΔCt calculations, fold change analysis, and quality control.

Author: Your Name
Date: June 2025
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings


class qPCRAnalyzer:
    """
    A comprehensive qPCR data analysis class that handles:
    - Data loading and validation
    - ΔCt calculation (target - reference)
    - ΔΔCt calculation (treated - control)  
    - Fold change calculation (2^(-ΔΔCt))
    - Quality control checks
    """
    
    def __init__(self, reference_gene='GAPDH', control_condition='Control'):
        """
        Initialize the qPCR analyzer.
        
        Parameters:
        -----------
        reference_gene : str
            Name of the reference/housekeeping gene (default: 'GAPDH')
        control_condition : str
            Name of the control condition (default: 'Control')
        """
        self.reference_gene = reference_gene
        self.control_condition = control_condition
        self.raw_data = None
        self.processed_data = None
        self.delta_ct_data = None
        self.delta_delta_ct_data = None
        self.fold_change_data = None
        
    def load_data(self, filepath):
        """
        Load qPCR data from CSV file.
        
        Parameters:
        -----------
        filepath : str or Path
            Path to the CSV file containing qPCR data
            
        Returns:
        --------
        pandas.DataFrame
            Loaded qPCR data
        """
        try:
            self.raw_data = pd.read_csv(filepath)
            print(f"✓ Successfully loaded {len(self.raw_data)} rows of data")
            return self.raw_data
        except FileNotFoundError:
            raise FileNotFoundError(f"Could not find file: {filepath}")
        except Exception as e:
            raise Exception(f"Error loading data: {str(e)}")
    
    def validate_data(self):
        """
        Validate the loaded qPCR data for required columns and data quality.
        
        Returns:
        --------
        dict
            Validation results and warnings
        """
        if self.raw_data is None:
            raise ValueError("No data loaded. Please run load_data() first.")
        
        required_columns = ['Sample_ID', 'Gene', 'Condition', 'Ct_Value']
        missing_columns = [col for col in required_columns if col not in self.raw_data.columns]
        
        validation_results = {
            'valid': len(missing_columns) == 0,
            'missing_columns': missing_columns,
            'warnings': []
        }
        
        if not validation_results['valid']:
            raise ValueError(f"Missing required columns: {missing_columns}")
        
        # Check for missing reference gene
        if self.reference_gene not in self.raw_data['Gene'].values:
            validation_results['warnings'].append(f"Reference gene '{self.reference_gene}' not found in data")
        
        # Check for missing control condition
        if self.control_condition not in self.raw_data['Condition'].values:
            validation_results['warnings'].append(f"Control condition '{self.control_condition}' not found in data")
        
        # Check for invalid Ct values
        invalid_ct = self.raw_data['Ct_Value'].isna() | (self.raw_data['Ct_Value'] <= 0) | (self.raw_data['Ct_Value'] > 40)
        if invalid_ct.any():
            validation_results['warnings'].append(f"Found {invalid_ct.sum()} invalid Ct values (NaN, ≤0, or >40)")
        
        # Print validation summary
        if validation_results['valid']:
            print("✓ Data validation passed")
            if validation_results['warnings']:
                print("⚠ Warnings:")
                for warning in validation_results['warnings']:
                    print(f"  - {warning}")
        
        return validation_results
    
    def calculate_technical_replicates_mean(self):
        """
        Calculate mean Ct values across technical replicates.
        
        Returns:
        --------
        pandas.DataFrame
            Data with technical replicates averaged
        """
        if self.raw_data is None:
            raise ValueError("No data loaded. Please run load_data() first.")
        
        # Group by everything except technical replicate and Ct_Value
        grouping_columns = ['Sample_ID', 'Gene', 'Condition', 'Biological_Rep']
        
        # Calculate mean and standard deviation of technical replicates
        self.processed_data = self.raw_data.groupby(grouping_columns).agg({
            'Ct_Value': ['mean', 'std', 'count']
        }).reset_index()
        
        # Flatten column names
        self.processed_data.columns = grouping_columns + ['Ct_Mean', 'Ct_Std', 'Tech_Rep_Count']
        
        print(f"✓ Calculated technical replicate means for {len(self.processed_data)} samples")
        return self.processed_data
    
    def calculate_delta_ct(self):
        """
        Calculate ΔCt values (Target Ct - Reference Ct) for each sample.
        
        Returns:
        --------
        pandas.DataFrame
            Data with ΔCt values calculated
        """
        if self.processed_data is None:
            self.calculate_technical_replicates_mean()
        
        # Pivot to get reference and target genes in separate columns
        pivot_data = self.processed_data.pivot_table(
            index=['Sample_ID', 'Condition', 'Biological_Rep'],
            columns='Gene',
            values='Ct_Mean',
            fill_value=np.nan
        ).reset_index()
        
        # Check if reference gene exists
        if self.reference_gene not in pivot_data.columns:
            raise ValueError(f"Reference gene '{self.reference_gene}' not found in data")
        
        # Calculate ΔCt for each target gene
        target_genes = [col for col in pivot_data.columns if col not in ['Sample_ID', 'Condition', 'Biological_Rep', self.reference_gene]]
        
        delta_ct_results = []
        
        for _, row in pivot_data.iterrows():
            reference_ct = row[self.reference_gene]
            
            for target_gene in target_genes:
                target_ct = row[target_gene]
                
                if pd.notna(reference_ct) and pd.notna(target_ct):
                    delta_ct = target_ct - reference_ct
                    
                    delta_ct_results.append({
                        'Sample_ID': row['Sample_ID'],
                        'Condition': row['Condition'],
                        'Biological_Rep': row['Biological_Rep'],
                        'Target_Gene': target_gene,
                        'Reference_Gene': self.reference_gene,
                        'Target_Ct': target_ct,
                        'Reference_Ct': reference_ct,
                        'Delta_Ct': delta_ct
                    })
        
        self.delta_ct_data = pd.DataFrame(delta_ct_results)
        print(f"✓ Calculated ΔCt values for {len(self.delta_ct_data)} measurements")
        return self.delta_ct_data
    
    def calculate_delta_delta_ct(self):
        """
        Calculate ΔΔCt values (Treated ΔCt - Control ΔCt) and fold changes.
        
        Returns:
        --------
        pandas.DataFrame
            Data with ΔΔCt and fold change values
        """
        if self.delta_ct_data is None:
            self.calculate_delta_ct()
        
        # Calculate mean ΔCt for control condition
        control_delta_ct = self.delta_ct_data[
            self.delta_ct_data['Condition'] == self.control_condition
        ].groupby('Target_Gene')['Delta_Ct'].mean().to_dict()
        
        # Calculate ΔΔCt and fold change
        results = []
        
        for target_gene in self.delta_ct_data['Target_Gene'].unique():
            gene_data = self.delta_ct_data[self.delta_ct_data['Target_Gene'] == target_gene]
            control_mean_delta_ct = control_delta_ct.get(target_gene, np.nan)
            
            for _, row in gene_data.iterrows():
                delta_delta_ct = row['Delta_Ct'] - control_mean_delta_ct
                fold_change = 2 ** (-delta_delta_ct) if pd.notna(delta_delta_ct) else np.nan
                
                results.append({
                    'Sample_ID': row['Sample_ID'],
                    'Condition': row['Condition'],
                    'Biological_Rep': row['Biological_Rep'],
                    'Target_Gene': row['Target_Gene'],
                    'Delta_Ct': row['Delta_Ct'],
                    'Control_Mean_Delta_Ct': control_mean_delta_ct,
                    'Delta_Delta_Ct': delta_delta_ct,
                    'Fold_Change': fold_change
                })
        
        self.delta_delta_ct_data = pd.DataFrame(results)
        print(f"✓ Calculated ΔΔCt and fold change values for {len(self.delta_delta_ct_data)} measurements")
        return self.delta_delta_ct_data
    
    def get_summary_statistics(self):
        """
        Generate summary statistics for the analysis.
        
        Returns:
        --------
        dict
            Summary statistics including means, standard deviations, and fold changes
        """
        if self.delta_delta_ct_data is None:
            self.calculate_delta_delta_ct()
        
        summary = self.delta_delta_ct_data.groupby(['Target_Gene', 'Condition']).agg({
            'Fold_Change': ['mean', 'std', 'count'],
            'Delta_Ct': ['mean', 'std'],
            'Delta_Delta_Ct': ['mean', 'std']
        }).round(3)
        
        print("✓ Generated summary statistics")
        return summary
    
    def run_complete_analysis(self, filepath):
        """
        Run the complete qPCR analysis pipeline.
        
        Parameters:
        -----------
        filepath : str or Path
            Path to the CSV file containing qPCR data
            
        Returns:
        --------
        dict
            Complete analysis results
        """
        print("🧬 Starting qPCR Analysis Pipeline")
        print("=" * 50)
        
        # Load and validate data
        self.load_data(filepath)
        self.validate_data()
        
        # Run analysis steps
        self.calculate_technical_replicates_mean()
        self.calculate_delta_ct()
        self.calculate_delta_delta_ct()
        
        # Generate summary
        summary = self.get_summary_statistics()
        
        print("\n✅ Analysis Complete!")
        print(f"📊 Analyzed {len(self.delta_delta_ct_data['Target_Gene'].unique())} target genes")
        print(f"🧪 Processed {len(self.delta_delta_ct_data['Sample_ID'].unique())} samples")
        
        return {
            'raw_data': self.raw_data,
            'processed_data': self.processed_data,
            'delta_ct_data': self.delta_ct_data,
            'delta_delta_ct_data': self.delta_delta_ct_data,
            'summary_statistics': summary
        }