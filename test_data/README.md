# Test Data for Tumor-Informed cfDNA MRD Pipeline

This directory contains test data and validation scripts for the Tumor-Informed cfDNA MRD Pipeline.

## Directory Structure

```
test_data/
├── wes/                    # WES test data (Step 1)
│   ├── tumor/             # Tumor WES FASTQ files
│   └── normal/            # Normal WBC WES FASTQ files
├── plasma/                 # Plasma cfDNA test data (Step 2)
│   ├── patients/          # Patient plasma samples (multiple timepoints)
│   └── healthy/           # Healthy donor plasma samples
├── refs/                   # Reference genome
├── resources/              # Resource files
│   ├── exome.interval_list # Exome capture intervals
│   ├── gnomad.af-only.vcf.gz # gnomAD germline resource
│   ├── tss.bed            # TSS regions for fragmentomics
│   ├── gc_hg38.wig        # GC content for ichorCNA
│   └── map_hg38.wig       # Mappability for ichorCNA
├── pon/                    # Panel of Normals
├── nextflow.config         # Test-specific configuration
├── test_pipeline.nf        # Test workflow for Step 1
├── test_step2.nf          # Test workflow for Step 2
├── run_test.sh             # Test execution script
└── README.md               # This file
```

## Test Data Description

### Step 1: WES Preprocessing and Somatic Variant Calling

- **Tumor WES**: `wes/tumor/test_patient_T0_R{1,2}.fastq.gz`
- **Normal WES**: `wes/normal/test_patient_T0_R{1,2}.fastq.gz`
- **Reference**: `refs/GRCh38.fa` (minimal mock genome)
- **Resources**: Exome intervals, gnomAD, PoN

### Step 2: Plasma cfDNA UMI-based WGS Processing

- **Patient Plasma**: `plasma/patients/test_patient_T{0,1}_R{1,2}.fastq.gz`
  - T0: Baseline timepoint
  - T1: Follow-up timepoint
  - UMI-tagged reads for consensus calling
- **Healthy Plasma**: `plasma/healthy/` (for control analysis)
- **Additional Resources**: TSS regions, GC content, mappability

## Running Tests

### Test Step 1 (WES Processing)

```bash
cd test_data
./run_test.sh
```

This will:
1. Validate input data structure
2. Test channel creation and data flow
3. Generate test outputs

### Test Step 2 (Plasma Processing)

```bash
cd test_data
nextflow run test_step2.nf -profile docker
```

This will:
1. Validate plasma data inputs
2. Test UMI extraction workflow
3. Validate reference file availability

## Expected Outputs

### Step 1 Outputs
- `results/wes/`: Preprocessed BAM files and QC metrics
- `results/variants/`: Somatic variants and truth set
- `results/pipeline_report.html`: Execution report

### Step 2 Outputs
- `results/plasma/umi/`: UMI-extracted FASTQ files
- `results/plasma/consensus/`: UMI consensus BAM files
- `results/plasma/variants/`: Variant calls at truth set loci
- `results/plasma/features/`: Fragmentomics and other features
- `results/plasma/cnv/`: CNV analysis results

## Test Limitations

- **Mock Data**: All FASTQ files contain minimal, non-biological sequences
- **Small Size**: Files are intentionally small for quick testing
- **Tool Execution**: Full bioinformatics tools are not executed in test mode
- **Docker Required**: Full pipeline testing requires Docker container

## Validation

The test workflows validate:
1. **Data Structure**: Correct file organization and naming
2. **Channel Logic**: Proper data flow between processes
3. **Parameter Passing**: Configuration and reference file handling
4. **Workflow Logic**: Process dependencies and execution order

## Next Steps

After successful testing:
1. Build Docker container with all tools
2. Run full pipeline with real data
3. Validate outputs against expected results
4. Scale to production AWS environment
