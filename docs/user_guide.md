# qPCR Analysis Toolkit - User Guide

## Table of Contents
1. [Getting Started](#getting-started)
2. [Data Format Requirements](#data-format-requirements)
3. [Basic Analysis Workflow](#basic-analysis-workflow)
4. [Statistical Analysis](#statistical-analysis)
5. [Visualization Options](#visualization-options)
6. [Report Generation](#report-generation)
7. [Troubleshooting](#troubleshooting)
8. [Advanced Features](#advanced-features)

## Getting Started

### Prerequisites
- Python 3.8 or higher
- Virtual environment (recommended)
- Basic understanding of qPCR methodology

### Installation
```bash
# Clone the repository
git clone <your-repo-url>
cd qpcr-analysis-toolkit

# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Quick Start
```bash
# Run the example analysis
python examples/basic_analysis.py

# Or test the complete pipeline
python test_complete_pipeline.py
```

## Data Format Requirements

### CSV File Structure
Your qPCR data should be formatted as a CSV file with the following required columns:

| Column | Description | Example |
|--------|-------------|---------|
| `Sample_ID` | Unique sample identifier | S001, S002, etc. |
| `Gene` | Gene name (target or reference) | GAPDH, IL6, TNF |
| `Condition` | Experimental condition | Control, Treatment |
| `Ct_Value` | PCR cycle threshold value | 18.45, 25.23, etc. |

### Optional Columns
| Column | Description | Purpose |
|--------|-------------|---------|
| `Biological_Rep` | Biological replicate number | Statistical analysis |
| `Technical_Rep` | Technical replicate number | Quality control |
| `Well_Position` | PCR plate well position | Traceability |

### Example Data Format
```csv
Sample_ID,Gene,Condition,Biological_Rep,Technical_Rep,Ct_Value,Well_Position
S001,GAPDH,Control,1,1,18.45,A1
S001,GAPDH,Control,1,2,18.52,A2
S001,IL6,Control,1,1,28.34,A3
S001,IL6,Control,1,2,28.41,A4
```

## Basic Analysis Workflow

### 1. Initialize the Analyzer
```python
from src.qpcr_analyzer import qPCRAnalyzer

analyzer = qPCRAnalyzer(
    reference_gene='GAPDH',        # Your housekeeping gene
    control_condition='Control'     # Your baseline condition
)
```

### 2. Load and Process Data
```python
# Run complete analysis pipeline
results = analyzer.run_complete_analysis('path/to/your/data.csv')

# Access specific results
raw_data = results['raw_data']
fold_changes = results['delta_delta_ct_data']
summary_stats = results['summary_statistics']
```

### 3. Quality Control Checks
The analyzer automatically performs:
- ✅ Data validation (required columns, valid Ct values)
- ✅ Technical replicate averaging
- ✅ Missing data handling
- ✅ Outlier detection warnings

## Statistical Analysis

### Running Statistical Tests
```python
from src.statistical_analysis import qPCRStatistics

# Initialize with significance level
stats_analyzer = qPCRStatistics(alpha=0.05)

# Perform comprehensive statistical analysis
stats_results = stats_analyzer.perform_comprehensive_analysis(
    fold_change_data=results['delta_delta_ct_data']
)
```

### Available Statistical Methods

#### T-Tests
- **Independent samples**: Compares treatment vs control groups
- **Paired samples**: For matched experimental designs
- Automatically calculates effect sizes (Cohen's d)

#### Multiple Comparison Corrections
- **Bonferroni**: Conservative correction (multiplies p-values)
- **FDR (Benjamini-Hochberg)**: Controls false discovery rate
- **Holm**: Step-down method balancing power and control

#### Effect Size Interpretation
| Cohen's d | Interpretation |
|-----------|---------------|
| < 0.2 | Negligible |
| 0.2 - 0.5 | Small |
| 0.5 - 0.8 | Medium |
| > 0.8 | Large |

## Visualization Options

### Creating Publication-Ready Plots
```python
from src.visualization import qPCRVisualizer

visualizer = qPCRVisualizer()

# Generate all plots
figures = visualizer.create_comprehensive_report(qpcr_results, stats_results)

# Or create individual plots
fold_change_plot = visualizer.plot_fold_change_bars(
    fold_change_data=results['delta_delta_ct_data'],
    stats_data=stats_results['ttest_results']
)
```

### Available Plot Types

#### 1. Fold Change Bar Charts
- Bar plots with error bars (SEM)
- Statistical significance annotations (*, **, ***)
- Customizable colors per condition
- Log scale for large fold changes

#### 2. Ct Distribution Plots
- Box plots showing data quality
- Grouped by gene and condition
- Outlier detection visual
- Inverted y-axis (lower Ct = higher expression)

#### 3. ΔCt Heatmaps
- Sample-wise expression patterns
- Color-coded by ΔCt values
- Hierarchical clustering options
- Missing data handling

#### 4. Quality Control Dashboard
- Technical replicate correlations
- Failed reaction summaries
- Sample quality distributions
- Data completeness metrics

## Report Generation

### Automated Reports
```python
from src.report_generator import qPCRReportGenerator

reporter = qPCRReportGenerator()

# Generate complete report package
report_paths = reporter.generate_complete_report_package(
    qpcr_results, stats_results
)
```

### Report Formats

#### HTML Reports
- Executive summary with key findings
- Embedded statistical tables
- Methodology explanations
- Professional formatting

#### Excel Exports
- Multiple sheets: Summary, Statistics, Fold Changes, Raw Data
- Ready for further analysis
- Formatted tables with headers

#### JSON Exports
- Machine-readable format
- Complete analysis metadata
- API integration ready

## Troubleshooting

### Common Issues and Solutions

#### "FileNotFoundError: Could not find file"
```python
# Solution: Check file path and format
import os
print(os.path.exists('your_file.csv'))  # Should return True
```

#### "Missing required columns"
```python
# Solution: Verify your CSV has these columns:
required_columns = ['Sample_ID', 'Gene', 'Condition', 'Ct_Value']
```

#### "Reference gene not found"
```python
# Solution: Check gene names in your data
print(data['Gene'].unique())  # List all genes in your dataset
```

#### "No significant results"
Possible causes:
- Small effect sizes (biological variation)
- Insufficient sample sizes
- High variability in technical replicates
- Need for different reference genes

#### "Technical replicate issues"
```python
# Check technical replicate quality
tech_rep_summary = data.groupby(['Sample_ID', 'Gene'])['Ct_Value'].agg(['mean', 'std'])
print(tech_rep_summary[tech_rep_summary['std'] > 0.5])  # High variability
```

### Data Quality Guidelines

#### Acceptable Ct Value Ranges
- **Reference genes**: 15-25 (stable, moderate expression)
- **Target genes**: 15-35 (depending on expression level)
- **Failed reactions**: Ct > 40 or no amplification

#### Technical Replicate Guidelines
- **Standard deviation**: < 0.5 Ct between replicates
- **Coefficient of variation**: < 5% for fold change calculations
- **Minimum replicates**: 2 technical + 3 biological replicates

## Advanced Features

### Custom Reference Genes
```python
# Multiple reference genes (advanced)
analyzer = qPCRAnalyzer(
    reference_gene=['GAPDH', 'ACTB'],  # Geometric mean normalization
    control_condition='Baseline'
)
```

### Batch Processing
```python
# Process multiple experiments
experiments = ['exp1.csv', 'exp2.csv', 'exp3.csv']

results = {}
for exp in experiments:
    results[exp] = analyzer.run_complete_analysis(exp)
```

### Custom Statistical Parameters
```python
# Adjust statistical parameters
stats_analyzer = qPCRStatistics(
    alpha=0.01,                    # More stringent significance level
    test_type='paired'             # For matched samples
)
```

### Publication-Ready Customization
```python
# Customize plot appearance
visualizer = qPCRVisualizer()
visualizer.colors = {
    'Control': '#1f77b4',
    'Treatment': '#ff7f0e',