# Step 1 Implementation: Tumor-Normal WES → Somatic Truth Set

## Overview

This document details the implementation of Step 1 of the Tumor-Informed cfDNA MRD Pipeline, which processes tumor-normal WES pairs to create a high-confidence somatic "truth set" for each patient.

## Workflow Summary

```
Tumor/Normal FASTQ → FastQC → fastp → BWA-MEM2 → MarkDuplicates → HS Metrics → Mutect2 → Truth Set
```

## Detailed Implementation

### 1.1 FastQC → fastp (trim)

**Process**: `WES_PREPROCESSING`

**Inputs**:
- `tumor_T0_R{1,2}.fastq.gz` (tumor tissue WES)
- `normalWBC_T0_R{1,2}.fastq.gz` (matched normal WES)

**Outputs**:
- `wes/qc/*_fastp.html` - FastQC reports
- `wes/qc/*_fastp.json` - FastQC JSON data
- `wes/trim/*_R{1,2}.fq.gz` - Trimmed FASTQ files

**Commands**:
```bash
# FastQC
fastqc -o wes/qc tumor_T0_R1.fastq.gz tumor_T0_R2.fastq.gz
fastqc -o wes/qc normalWBC_T0_R1.fastq.gz normalWBC_T0_R2.fastq.gz

# fastp trimming
fastp -i tumor_T0_R1.fastq.gz -I tumor_T0_R2.fastq.gz \
      -o wes/trim/tumor_R1.fq.gz -O wes/trim/tumor_R2.fq.gz \
      -h wes/qc/tumor_fastp.html -j wes/qc/tumor_fastp.json

fastp -i normalWBC_T0_R1.fastq.gz -I normalWBC_T0_R2.fastq.gz \
      -o wes/trim/normal_R1.fq.gz -O wes/trim/normal_R2.fq.gz \
      -h wes/qc/normal_fastp.html -j wes/qc/normal_fastp.json
```

### 1.2 Align (BWA-MEM2) → sort/index

**Process**: `WES_PREPROCESSING`

**Inputs**: Trimmed FASTQ files from step 1.1

**Outputs**:
- `wes/tumor_T0.bam` - Aligned tumor BAM
- `wes/normal_T0.bam` - Aligned normal BAM
- `wes/*.bai` - BAM index files

**Commands**:
```bash
# BWA-MEM2 alignment
bwa-mem2 mem -t 16 refs/GRCh38.fa wes/trim/tumor_R1.fq.gz wes/trim/tumor_R2.fq.gz \
  | samtools sort -@8 -o wes/tumor_T0.bam
samtools index wes/tumor_T0.bam

bwa-mem2 mem -t 16 refs/GRCh38.fa wes/trim/normal_R1.fq.gz wes/trim/normal_R2.fq.gz \
  | samtools sort -@8 -o wes/normal_T0.bam
samtools index wes/normal_T0.bam
```

### 1.3 Mark dups, metrics, target-restricted depth

**Process**: `WES_PREPROCESSING`

**Inputs**: Aligned BAM files from step 1.2

**Outputs**:
- `wes/tumor_T0.dedup.bam` - Deduplicated tumor BAM
- `wes/normal_T0.dedup.bam` - Deduplicated normal BAM
- `wes/qc/*.dup.txt` - Duplicate metrics
- `wes/qc/*.hs.txt` - Hybrid selection metrics

**Commands**:
```bash
# Mark duplicates
gatk MarkDuplicatesSpark -I wes/tumor_T0.bam -O wes/tumor_T0.dedup.bam -M wes/qc/tumor.dup.txt
gatk MarkDuplicatesSpark -I wes/normal_T0.bam -O wes/normal_T0.dedup.bam -M wes/qc/normal.dup.txt
samtools index wes/tumor_T0.dedup.bam; samtools index wes/normal_T0.dedup.bam

# Collect HS metrics
gatk CollectHsMetrics I=wes/tumor_T0.dedup.bam O=wes/qc/tumor.hs.txt R=refs/GRCh38.fa \
   BAIT_INTERVALS=resources/exome.interval_list TARGET_INTERVALS=resources/exome.interval_list
gatk CollectHsMetrics I=wes/normal_T0.dedup.bam O=wes/qc/normal.hs.txt R=refs/GRCh38.fa \
   BAIT_INTERVALS=resources/exome.interval_list TARGET_INTERVALS=resources/exome.interval_list
```

### 1.4 Mutect2 (tumor–normal mode) with PoN + gnomAD

**Process**: `SOMATIC_VARIANT_CALLING`

**Inputs**: Deduplicated BAM files from step 1.3

**Outputs**:
- `wes/tumor_vs_normal.filtered.vcf.gz` - Filtered somatic variants
- `wes/tumor.contamination.table` - Contamination estimates

**Commands**:
```bash
# Pileups for contamination/orientation
gatk GetPileupSummaries -I wes/tumor_T0.dedup.bam -V resources/gnomad.af-only.vcf.gz \
  -L resources/exome.interval_list -O wes/tumor.pileups.table
gatk GetPileupSummaries -I wes/normal_T0.dedup.bam -V resources/gnomad.af-only.vcf.gz \
  -L resources/exome.interval_list -O wes/normal.pileups.table
gatk CalculateContamination -I wes/tumor.pileups.table -matched wes/normal.pileups.table \
  -O wes/tumor.contamination.table

# Mutect2 calling
gatk Mutect2 \
  -R refs/GRCh38.fa \
  -I wes/tumor_T0.dedup.bam -I wes/normal_T0.dedup.bam \
  -normal normal_sample_id \
  --panel-of-normals pon/pon.vcf.gz \
  --germline-resource resources/gnomad.af-only.vcf.gz \
  -L resources/exome.interval_list \
  --f1r2-tar-gz wes/tumor.f1r2.tar.gz \
  -O wes/tumor_vs_normal.unfiltered.vcf.gz

# Learn read orientation model
gatk LearnReadOrientationModel -I wes/tumor.f1r2.tar.gz -O wes/tumor.read-orientation-model.tar.gz

# Filter variants
gatk FilterMutectCalls \
  -R refs/GRCh38.fa \
  -V wes/tumor_vs_normal.unfiltered.vcf.gz \
  --contamination-table wes/tumor.contamination.table \
  --ob-priors wes/tumor.read-orientation-model.tar.gz \
  -O wes/tumor_vs_normal.filtered.vcf.gz
```

### 1.5 Truth set derivation (filters → BED of loci)

**Process**: `SOMATIC_VARIANT_CALLING`

**Inputs**: Filtered VCF from step 1.4

**Outputs**:
- `wes/patientX.truthset.vcf.gz` - High-confidence somatic variants
- `wes/patientX.truthset.bed` - BED file of truth set loci

**Commands**:
```bash
# Keep PASS, exonic/splice, min tumor depth/VAF, remove germline/CHIP
bcftools view -f PASS wes/tumor_vs_normal.filtered.vcf.gz \
  | bcftools filter -i 'FORMAT[AD][0:1]>=5 && FORMAT[DP]>=50 && INFO[TLOD]>=6' \
  | bcftools annotate -x INFO,^INFO/CSQ \
  -Oz -o wes/patientX.truthset.vcf.gz
tabix -p vcf wes/patientX.truthset.vcf.gz

# BED of loci for fast scanning in plasma
bcftools query -f'%CHROM\t%POS0\t%END\n' wes/patientX.truthset.vcf.gz | sort -k1,1 -k2,2n \
  > wes/patientX.truthset.bed
```

## Key Outputs

### Per-Patient Somatic Truth Set
- **VCF**: `patientX.truthset.vcf.gz` - High-confidence somatic variants
- **BED**: `patientX.truthset.bed` - Genomic coordinates for plasma scanning

### Quality Metrics
- **FastQC reports**: HTML and JSON quality reports
- **Duplicate metrics**: Duplicate rate and statistics
- **HS metrics**: Hybrid selection and coverage statistics
- **Contamination estimates**: Sample contamination levels

## Configuration Parameters

```groovy
params {
    // Quality thresholds
    minDepth = 50        // Minimum depth for truth set
    minTlod = 6          // Minimum tumor LOD score
    minAd = 5            // Minimum alternate allele depth
    
    // File paths
    genome = "$baseDir/refs/GRCh38.fa"
    exomeBed = "$baseDir/resources/exome.interval_list"
    gnomadVcf = "$baseDir/resources/gnomad.af-only.vcf.gz"
    ponVcf = "$baseDir/pon/pon.vcf.gz"
}
```

## Resource Requirements

### WES_PREPROCESSING
- **CPUs**: 16
- **Memory**: 64 GB
- **Time**: 6 hours

### SOMATIC_VARIANT_CALLING
- **CPUs**: 16
- **Memory**: 64 GB
- **Time**: 8 hours

## Testing

Run the test workflow:
```bash
nextflow run test_step1.nf -profile docker
```

## Next Steps

After completing Step 1:
1. **Validate truth sets** for all 90 patients
2. **Implement Step 2**: Plasma cfDNA feature extraction
3. **Implement Step 3**: MRD classification using truth set loci

## Troubleshooting

### Common Issues
1. **Memory errors**: Increase Java heap size in GATK commands
2. **File not found**: Check reference file paths in configuration
3. **Permission errors**: Ensure proper file permissions for output directories

### Quality Checks
1. **Coverage**: Verify target coverage >100x for tumor samples
2. **Variant quality**: Check TLOD scores and depth distributions
3. **Contamination**: Monitor contamination estimates (<5% recommended)
