# Test Data for Tumor-Informed cfDNA MRD Pipeline

This directory contains minimal test data to validate **Step 1** of the Tumor-Informed cfDNA MRD Pipeline: **Tumor-Normal WES → Somatic Truth Set**.

## Test Data Structure

```
test_data/
├── wes/
│   ├── tumor/
│   │   ├── test_patient_T0_R1.fastq.gz    # Test tumor R1 FASTQ
│   │   └── test_patient_T0_R2.fastq.gz    # Test tumor R2 FASTQ
│   └── normal/
│       ├── test_patient_T0_R1.fastq.gz    # Test normal R1 FASTQ
│       └── test_patient_T0_R2.fastq.gz    # Test normal R2 FASTQ
├── refs/
│   └── GRCh38.fa                          # Minimal reference genome
├── resources/
│   ├── exome.interval_list                # Exome capture intervals
│   └── gnomad.af-only.vcf.gz             # Minimal gnomAD resource
├── pon/
│   └── pon.vcf.gz                         # Minimal Panel of Normals
├── nextflow.config                         # Test configuration
├── test_pipeline.nf                        # Test pipeline script
├── run_test.sh                             # Test execution script
└── README.md                               # This file
```

## Test Data Description

### FASTQ Files
- **Minimal paired-end reads** (3 reads per file)
- **64bp sequences** with high quality scores
- **Naming convention**: `test_patient_T0_R{1,2}.fastq.gz`
- **Format**: Standard FASTQ with quality scores

### Reference Files
- **GRCh38.fa**: Minimal genome with 2 chromosomes (chr1, chr2)
- **exome.interval_list**: 2 target regions for exome capture
- **gnomad.af-only.vcf.gz**: 2 common variants for filtering
- **pon.vcf.gz**: 2 variants for Panel of Normals

## Running the Test

### Prerequisites
1. **Nextflow** (≥22.10.0)
2. **Docker** (running)
3. **Bash shell**

### Quick Test
```bash
# Navigate to test data directory
cd test_data

# Run the test
./run_test.sh
```

### Manual Test
```bash
# Navigate to test data directory
cd test_data

# Run with Nextflow
nextflow run test_pipeline.nf -profile docker
```

## Expected Outputs

After successful execution, you should see:

```
results/
├── wes/
│   ├── test_patient_T0.bam                # Aligned tumor BAM
│   ├── test_patient_T0.bai                # Tumor BAM index
│   ├── test_patient_T0.dedup.bam          # Deduplicated tumor BAM
│   ├── test_patient_T0.dedup.bai          # Deduplicated tumor BAM index
│   ├── normal_T0.bam                      # Aligned normal BAM
│   ├── normal_T0.bai                      # Normal BAM index
│   ├── normal_T0.dedup.bam                # Deduplicated normal BAM
│   ├── normal_T0.dedup.bai                # Deduplicated normal BAM index
│   └── qc/                                # Quality control reports
│       ├── tumor_fastp.html               # Tumor FastQC report
│       ├── tumor_fastp.json               # Tumor FastQC data
│       ├── normal_fastp.html              # Normal FastQC report
│       ├── normal_fastp.json              # Normal FastQC data
│       ├── tumor.dup.txt                  # Tumor duplicate metrics
│       ├── normal.dup.txt                 # Normal duplicate metrics
│       ├── tumor.hs.txt                   # Tumor HS metrics
│       └── normal.hs.txt                  # Normal HS metrics
├── variants/
│   ├── test_patient.truthset.vcf.gz       # Somatic truth set VCF
│   ├── test_patient.truthset.vcf.gz.tbi   # Truth set VCF index
│   ├── test_patient.truthset.bed          # Truth set BED file
│   └── tumor_vs_normal.filtered.vcf.gz    # Filtered somatic variants
├── pipeline_report.html                    # Pipeline execution report
├── timeline_report.html                    # Timeline report
└── trace.txt                              # Execution trace
```

## Test Validation

### What to Check
1. **Pipeline completion** - No errors in execution
2. **Output files** - All expected files generated
3. **File sizes** - BAM files should be >0 bytes
4. **Quality metrics** - Check QC reports for reasonable values

### Common Issues
1. **Docker not running** - Start Docker daemon
2. **Memory issues** - Increase Docker memory allocation
3. **File permissions** - Ensure read/write access to directories

## Test Data Limitations

⚠️ **Important Notes**:
- **Minimal data**: Only 3 reads per file (not realistic for production)
- **Small genome**: Only 2 chromosomes with 64bp each
- **Limited variants**: Only 2 variants in gnomAD and PoN
- **No real mutations**: Test data doesn't contain actual somatic variants

This test data is designed to **validate the pipeline workflow** and **test tool integration**, not for **biological analysis**.

## Next Steps

After successful testing:
1. **Replace with real data** - Use your actual tumor-normal WES pairs
2. **Update references** - Use full GRCh38 genome and real resources
3. **Scale resources** - Adjust CPU/memory for production runs
4. **Implement Step 2** - Plasma cfDNA feature extraction

## Troubleshooting

### Pipeline Fails
- Check Docker is running: `docker info`
- Verify Nextflow version: `nextflow -version`
- Check file permissions and paths
- Review error logs in `work/` directory

### Missing Outputs
- Check `results/` directory structure
- Verify input file paths in configuration
- Review process logs for errors
- Check resource allocation (CPU/memory)

### Performance Issues
- Increase Docker memory allocation
- Adjust process resource requirements in config
- Use AWS Batch for cloud execution
- Optimize input file sizes and formats
