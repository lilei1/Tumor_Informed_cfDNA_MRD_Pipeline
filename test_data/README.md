# Test Data for Tumor-Informed cfDNA MRD Pipeline

This directory contains test data and workflows for validating the Tumor-Informed cfDNA MRD Pipeline.

## Directory Structure

```
test_data/
├── wes/                    # WES test data (Step 1)
│   ├── tumor/             # Tumor WES FASTQ files
│   └── normal/            # Normal WBC WES FASTQ files
├── plasma/                 # Plasma cfDNA test data (Step 2)
│   ├── patients/          # Patient plasma samples (multiple timepoints)
│   └── healthy/           # Healthy donor plasma samples (Step 2.5)
├── refs/                   # Reference genome
├── resources/              # Resource files
│   ├── exome.interval_list # Exome capture intervals
│   ├── gnomad.af-only.vcf.gz # gnomAD germline resource
│   ├── tss.bed            # TSS regions for fragmentomics
│   ├── gc_hg38.wig        # GC content for ichorCNA
│   ├── map_hg38.wig       # Mappability for ichorCNA
│   └── mock_truthset.bed  # Mock truth set for error model testing
├── pon/                    # Panel of Normals
├── nextflow.config         # Test-specific configuration
├── test_pipeline.nf        # Test workflow for Step 1
├── test_step2.nf          # Test workflow for Step 2
├── test_step2_5.nf        # Test workflow for Step 2.5
├── test_step3.nf          # Test workflow for Step 3
├── simple_step2_test.nf   # Simple Step 2 validation
├── test_step2_5_comprehensive.nf # Comprehensive Step 2.5 test
├── test_step2_5_simple.nf # Simple Step 2.5 validation
├── run_test.sh             # Test execution script
└── README.md               # This file
```

## Test Data Description

### Step 1: WES Preprocessing and Somatic Variant Calling

**Input Data:**
- `wes/tumor/test_patient_T0_R{1,2}.fastq.gz`: Mock tumor WES data
- `wes/normal/test_patient_T0_R{1,2}.fastq.gz`: Mock normal WBC WES data

**Expected Outputs:**
- `results/wes/`: Preprocessed BAM files, QC metrics
- `results/somatic/`: Somatic variant calls, truth set BED

**Test Command:**
```bash
cd test_data
nextflow run test_pipeline.nf
```

### Step 2: Plasma cfDNA UMI-based WGS Processing

**Input Data:**
- `plasma/patients/test_patient_T{0,1}_R{1,2}.fastq.gz`: Mock patient plasma samples
- Multiple timepoints (T0, T1) for longitudinal analysis

**Expected Outputs:**
- `results/plasma/consensus/`: UMI consensus BAM files
- `results/plasma/variants/`: Variant evidence at truth set loci
- `results/plasma/features/`: Fragmentomics, end-motifs, TSS coverage
- `results/plasma/cnv/`: CNV analysis results

**Test Command:**
```bash
cd test_data
nextflow run test_step2.nf
```

### Step 2.5: Background Error Model from Healthy Donors

**Input Data:**
- `plasma/healthy/healthy_donor_{1,2,3}_R{1,2}.fastq.gz`: Mock healthy donor plasma
- 3 donors with diverse UMI patterns (ATCG, CGAT, TAGC)
- Each donor has 3 UMI families for robust error modeling

**Expected Outputs:**
- `results/error_model/`: Context-specific error rates (JSON/TSV)
- `results/error_model/validation/`: Error model validation results

**Test Commands:**
```bash
cd test_data
# Basic validation
nextflow run test_step2_5.nf

# Enhanced validation
nextflow run test_step2_5_simple.nf

# Comprehensive testing
nextflow run test_step2_5_comprehensive.nf
```

### Step 3: Feature Integration with Error-Aware Scoring

**Input Data:**
- Integrated features from Steps 1, 2, and 2.5
- Background error model for context-aware variant scoring
- Multi-modal features: variants, fragmentomics, CNV, methylation

**Expected Outputs:**
- `results/features/integrated/`: Unified feature matrices
- `results/features/mrd_scores/`: MRD probability scores
- `results/features/longitudinal/`: Longitudinal MRD trends
- `results/reports/`: Comprehensive MRD analysis reports

**Test Command:**
```bash
cd test_data
nextflow run test_step3.nf
```

## Mock Data Details

### Enhanced Healthy Donor Data (Step 2.5)
- **3 Healthy Donors**: Each with unique UMI patterns
- **9 Sample Files**: 3 donors × 3 UMI families each
- **Diverse Sequences**: ATCG, CGAT, TAGC patterns for error model diversity
- **Realistic Structure**: Proper paired-end FASTQ format with quality scores

### Enhanced Truth Set (Step 2.5)
- **12 Regions**: Across 3 chromosomes (chr1, chr2, chr3)
- **Varied Contexts**: Different sequence contexts for error model training
- **Realistic Coordinates**: Simulated genomic positions

### Enhanced Reference Genome (Step 2.5)
- **3 Chromosomes**: chr1, chr2, chr3
- **780 bp Total**: Sufficient for testing genomic operations
- **Realistic Sequences**: TGCATGCATGC... patterns

## Running Tests

### Prerequisites
- Nextflow 20.0 or later
- Docker (for full pipeline execution)
- Python 3.7+ with scientific packages

### Test Execution
```bash
# Navigate to test directory
cd test_data

# Run individual step tests
nextflow run test_pipeline.nf      # Step 1
nextflow run test_step2.nf         # Step 2
nextflow run test_step2_5.nf       # Step 2.5
nextflow run test_step3.nf         # Step 3

# Run comprehensive tests
nextflow run test_step2_5_simple.nf
nextflow run test_step2_5_comprehensive.nf

# Run all tests with script
bash run_test.sh
```

### Expected Results
Each test should complete successfully and generate:
- Validation reports in `results/validate/`
- Mock data files for testing
- Process logs and execution traces

## Validation Status

- ✅ **Step 1**: WES preprocessing and somatic variant calling
- ✅ **Step 2**: Plasma cfDNA UMI-based processing
- ✅ **Step 2.5**: Background error model from healthy donors
- ✅ **Step 3**: Feature integration with error-aware scoring
- 🔄 **Step 4**: MRD classification (ready to implement)
- 🔄 **Step 5**: Validation and reporting (ready to implement)

## Notes

- All test workflows use mock data to validate pipeline logic
- External tools (samtools, bwa, etc.) are not required for basic testing
- Full pipeline execution requires Docker container with all bioinformatics tools
- Test data is designed to be biologically plausible while remaining lightweight
- Enhanced toy data provides robust foundation for testing complex workflows

## Troubleshooting

### Common Issues
1. **Missing files**: Ensure all mock data files are present
2. **Channel errors**: Check Nextflow syntax and channel definitions
3. **Process failures**: Review process input/output specifications
4. **Memory issues**: Adjust resource requirements in nextflow.config

### Debug Mode
Run tests with verbose output:
```bash
nextflow run test_step3.nf -v
```

### Log Files
Check `.nextflow.log` for detailed execution information and error messages.
