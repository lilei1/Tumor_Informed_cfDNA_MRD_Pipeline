# Tumor-Informed cfDNA MRD Pipeline

A comprehensive pipeline for detecting Minimal Residual Disease (MRD) in colorectal cancer patients using tumor-informed circulating cell-free DNA analysis.

## Project Overview

This pipeline integrates multiple data modalities to detect ultra-low frequency variants and predict MRD status:
- **90 Stage II/III CRC patients** with matched tumor/WBC exome data
- **50 healthy controls** for baseline establishment
- **Longitudinal plasma cfDNA samples** (t0-t4) with UMI barcoding
- **Multi-modal analysis**: variants, methylation, and fragmentomics

## Objectives

1. **Develop a statistically principled caller** for ultra-low VAF variants with context-aware error suppression and UMI consensus
2. **Classifier**: Integrate variant, methylation, and fragmentomics features into a single MRD probability with calibrated thresholds
3. **Validation**: Benchmark LoD, sensitivity/specificity, and longitudinal precision using spike-ins, reference standards, and retrospective cohorts
4. **Platform**: Ship a Nextflow + Docker pipeline running on AWS Batch/EC2/S3, with automated QC, audit-ready logs, and one-click reports

## Pipeline Architecture

```
├── workflows/           # Nextflow workflow definitions
├── modules/            # Individual processing modules
├── scripts/            # Utility scripts and analysis tools
├── configs/            # Configuration files
├── containers/         # Docker container definitions
├── tests/              # Test data and validation scripts
├── docs/               # Documentation and protocols
└── results/            # Output directory structure
```

## Implementation Phases

### Phase 1: Core Infrastructure (Weeks 1-4)
- [ ] UMI-aware preprocessing pipeline
- [ ] Tumor-informed variant calling
- [ ] Context-aware error suppression
- [ ] Quality control modules

### Phase 2: Feature Extraction (Weeks 5-8)
- [ ] Variant feature extraction
- [ ] Methylation analysis pipeline
- [ ] Fragmentomics feature extraction
- [ ] Feature integration framework

### Phase 3: MRD Classification (Weeks 9-12)
- [ ] Machine learning classifier development
- [ ] Feature selection and optimization
- [ ] Threshold calibration
- [ ] Model validation

### Phase 4: Validation & Benchmarking (Weeks 13-16)
- [ ] Spike-in experiments
- [ ] Reference standard validation
- [ ] Retrospective cohort analysis
- [ ] Performance benchmarking

### Phase 5: Production Platform (Weeks 17-20)
- [ ] Nextflow workflow implementation
- [ ] Docker containerization
- [ ] AWS infrastructure setup
- [ ] Automated QC and reporting

## Getting Started

1. Install Nextflow: `curl -s https://get.nextflow.io | bash`
2. Install Docker: Follow [Docker installation guide](https://docs.docker.com/get-docker/)
3. Configure AWS credentials
4. Run the pipeline: `nextflow run main.nf`

## Citation

If you use this pipeline in your research, please cite:
[Citation information to be added]

## License

[License information to be added]
