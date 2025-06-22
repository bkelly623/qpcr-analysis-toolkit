# qPCR Analysis Toolkit

A comprehensive Python toolkit for analyzing quantitative PCR (qPCR) data using industry-standard methodologies. This project implements the ΔΔCt method for relative gene expression analysis with proper statistical testing and publication-ready visualizations.

## 🧬 Project Overview

This toolkit processes raw qPCR Ct (cycle threshold) values and performs:
- ΔCt normalization against reference genes
- ΔΔCt calculation for relative quantification  
- Fold change analysis using the 2^(-ΔΔCt) method
- Statistical significance testing
- Quality control and data validation
- Professional scientific visualizations

## 🎯 Use Cases

- **Drug Discovery**: Analyze gene expression changes in response to treatments
- **Clinical Diagnostics**: Quantify viral loads or biomarker expression
- **Research Studies**: Compare gene expression across experimental conditions
- **Quality Control**: Validate qPCR data integrity and experimental design

## 📊 Example Results

The toolkit successfully analyzes inflammatory response data showing:
- **IL6**: 14.9-fold upregulation (p < 0.001)
- **TNF**: 7.0-fold upregulation (p < 0.001)

Perfect for demonstrating the analysis of cytokine responses in biomedical research.

## 🚀 Quick Start

### Prerequisites
- Python 3.8+
- Virtual environment (recommended)

### Installation

```bash
# Clone the repository
git clone <your-repo-url>
cd qpcr-analysis-toolkit

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

```python
from src.qpcr_analyzer import qPCRAnalyzer

# Initialize analyzer
analyzer = qPCRAnalyzer(reference_gene='GAPDH', control_condition='Control')

# Run complete analysis
results = analyzer.run_complete_analysis('data/sample_experiment.csv')

# View summary statistics
print(results['summary_statistics'])
```

## 📁 Project Structure

```
qpcr-analysis-toolkit/
├── README.md                 # Project documentation
├── requirements.txt          # Python dependencies
├── .gitignore               # Git ignore rules
├── data/                    # Sample datasets
│   ├── sample_experiment.csv    # Example qPCR data
│   └── metadata.csv            # Gene information
├── src/                     # Source code modules
│   ├── qpcr_analyzer.py        # Core analysis engine
│   ├── statistical_analysis.py # Statistical testing
│   ├── visualization.py        # Plotting functions
│   └── utils.py                # Utility functions
├── examples/                # Usage examples
│   └── basic_analysis.py       # Getting started guide
├── output/                  # Generated results
│   ├── plots/                  # Visualization outputs
│   └── reports/                # Analysis reports
└── tests/                   # Unit tests
    └── test_qpcr_analyzer.py   # Core functionality tests
```

## 🔬 Scientific Methodology

### The ΔΔCt Method

This toolkit implements the gold-standard ΔΔCt method for relative quantification:

1. **ΔCt Calculation**: `Target_Ct - Reference_Ct`
   - Normalizes target gene expression against housekeeping genes
   - Accounts for sample-to-sample variation in RNA/DNA quality

2. **ΔΔCt Calculation**: `Treated_ΔCt - Control_ΔCt`  
   - Compares experimental conditions to control baseline
   - Enables relative quantification between groups

3. **Fold Change**: `2^(-ΔΔCt)`
   - Converts log₂ scale to linear fold change
   - Interpretable as "X-fold increase/decrease"

### Quality Control Features

- **Technical Replicate Averaging**: Reduces PCR-specific noise
- **Data Validation**: Checks for invalid Ct values (NaN, ≤0, >40)
- **Missing Data Handling**: Robust analysis with incomplete datasets
- **Outlier Detection**: Statistical identification of aberrant values

## 📈 Features

### Current Implementation
- ✅ Core ΔΔCt analysis pipeline
- ✅ Technical replicate processing
- ✅ Data validation and quality control
- ✅ Summary statistics generation
- ✅ Comprehensive error handling

### In Development
- 🔄 Statistical significance testing (t-tests, multiple comparisons)
- 🔄 Publication-ready visualizations (bar charts, heatmaps, box plots)
- 🔄 Automated report generation (PDF/Excel outputs)
- 🔄 Batch processing for multiple experiments

## 🧪 Sample Data Format

The toolkit expects CSV data with the following structure:

```csv
Sample_ID,Gene,Condition,Biological_Rep,Technical_Rep,Ct_Value,Well_Position
S001,GAPDH,Control,1,1,18.45,A1
S001,IL6,Control,1,1,28.34,A3
S002,GAPDH,Treatment,1,1,18.78,D1
S002,IL6,Treatment,1,1,25.23,D3
```

### Required Columns
- `Sample_ID`: Unique sample identifier
- `Gene`: Target or reference gene name  
- `Condition`: Experimental condition (e.g., Control, Treatment)
- `Ct_Value`: PCR cycle threshold value
- `Biological_Rep`: Biological replicate number
- `Technical_Rep`: Technical replicate number (optional)

## 🔧 Development Setup

### Running Tests
```bash
# Run the basic functionality test
python test_analyzer.py

# Run full test suite (when implemented)
pytest tests/
```

### Virtual Environment Management
```bash
# Activate environment (required for each session)
source venv/bin/activate

# Deactivate when done
deactivate
```

## 📚 Dependencies

- **pandas**: Data manipulation and analysis
- **numpy**: Numerical computing and array operations
- **scipy**: Statistical functions and hypothesis testing
- **matplotlib**: Publication-quality plotting
- **seaborn**: Statistical data visualization
- **openpyxl**: Excel file input/output
- **pytest**: Unit testing framework

## 🎓 Learning Objectives

This project demonstrates proficiency in:

### Molecular Biology
- qPCR principles and data interpretation
- Gene expression analysis methodologies
- Experimental design for biological studies
- Reference gene selection and validation

### Data Science
- Statistical analysis of biological data
- Hypothesis testing and multiple comparisons
- Data visualization and scientific plotting
- Quality control and data validation

### Software Engineering  
- Object-oriented Python programming
- Modular code architecture and design patterns
- Error handling and input validation
- Unit testing and test-driven development
- Version control with Git
- Professional documentation

## 🤝 Contributing

This is a portfolio project demonstrating biotech data analysis capabilities. For suggestions or improvements, please open an issue or submit a pull request.

## 📄 License

This project is available under the MIT License. See LICENSE file for details.

## 👨‍💻 Author

**Your Name**  
*Biotech Data Scientist*  
- GitHub: [@yourusername](https://github.com/yourusername)
- LinkedIn: [Your LinkedIn](https://linkedin.com/in/yourprofile)

---

*Built with modern Python tools for the biotech industry. Demonstrates real-world molecular biology data analysis capabilities.*