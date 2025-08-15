#!/usr/bin/env nextflow

/*
 * Preprocessing Module
 * Handles WES preprocessing for tumor-normal pairs to create somatic truth set
 */

process WES_PREPROCESSING {
    tag "${sample}_${type}"
    
    publishDir "${params.outdir}/wes", mode: 'copy'
    
    cpus 16
    memory '64 GB'
    time '6h'
    
    input:
    tuple val(sample), val(type), path(fastq1), path(fastq2)
    path ref_genome
    path exome_intervals
    
    output:
    tuple val(sample), val(type), path("*.bam"), path("*.bai"), emit: bam
    tuple val(sample), val(type), path("*.dedup.bam"), path("*.dedup.bai"), emit: dedup_bam
    tuple val(sample), val(type), path("*.hs.txt"), emit: hs_metrics
    tuple val(sample), val(type), path("*.dup.txt"), emit: dup_metrics
    path "*.log", emit: logs
    
    script:
    def sample_id = type == 'tumor' ? 'tumor_sample_id' : 'normal_sample_id'
    
    """
    #!/bin/bash
    set -e
    
    # Create output directories
    mkdir -p wes/qc wes/trim
    
    echo "Processing ${type} sample: ${sample}"
    
    # 1.1 FastQC and fastp trimming
    echo "Running FastQC..."
    fastqc -o wes/qc ${fastq1} ${fastq2}
    
    echo "Running fastp trimming..."
    fastp -i ${fastq1} -I ${fastq2} \
          -o wes/trim/${type}_R1.fq.gz -O wes/trim/${type}_R2.fq.gz \
          -h wes/qc/${type}_fastp.html -j wes/qc/${type}_fastp.json \
          --qualified_quality_phred 20 \
          --length_required 50 \
          --detect_adapter_for_pe \
          --thread 8
    
    # 1.2 BWA-MEM2 alignment and sorting
    echo "Running BWA-MEM2 alignment..."
    bwa-mem2 mem -t 16 ${ref_genome} \
        wes/trim/${type}_R1.fq.gz wes/trim/${type}_R2.fq.gz \
        | samtools sort -@8 -o wes/${type}_T0.bam
    
    echo "Indexing BAM..."
    samtools index wes/${type}_T0.bam
    
    # 1.3 Mark duplicates and collect metrics
    echo "Marking duplicates..."
    gatk MarkDuplicatesSpark \
        -I wes/${type}_T0.bam \
        -O wes/${type}_T0.dedup.bam \
        -M wes/qc/${type}.dup.txt \
        --java-options "-Xmx32g"
    
    echo "Indexing deduplicated BAM..."
    samtools index wes/${type}_T0.dedup.bam
    
    echo "Collecting HS metrics..."
    gatk CollectHsMetrics \
        I=wes/${type}_T0.dedup.bam \
        O=wes/qc/${type}.hs.txt \
        R=${ref_genome} \
        BAIT_INTERVALS=${exome_intervals} \
        TARGET_INTERVALS=${exome_intervals} \
        --java-options "-Xmx32g"
    
    # Clean up intermediate files
    rm -f wes/trim/${type}_R*.fq.gz
    
    echo "${type} preprocessing completed successfully"
    """
}

process SOMATIC_VARIANT_CALLING {
    tag "${sample}"
    
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    cpus 16
    memory '64 GB'
    time '8h'
    
    input:
    tuple val(sample), path(tumor_dedup_bam), path(tumor_dedup_bai)
    tuple val(sample), path(normal_dedup_bam), path(normal_dedup_bai)
    path ref_genome
    path exome_intervals
    path gnomad_vcf
    path pon_vcf
    
    output:
    tuple val(sample), path("*.truthset.vcf.gz"), path("*.truthset.vcf.gz.tbi"), emit: truth_set
    tuple val(sample), path("*.truthset.bed"), emit: truth_set_bed
    tuple val(sample), path("*.filtered.vcf.gz"), path("*.filtered.vcf.gz.tbi"), emit: filtered_variants
    path "*.log", emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "Running somatic variant calling for ${sample}"
    
    # 1.4 Mutect2 with PoN + gnomAD
    echo "Getting pileup summaries for contamination..."
    gatk GetPileupSummaries \
        -I ${tumor_dedup_bam} \
        -V ${gnomad_vcf} \
        -L ${exome_intervals} \
        -O wes/tumor.pileups.table \
        --java-options "-Xmx32g"
    
    gatk GetPileupSummaries \
        -I ${normal_dedup_bam} \
        -V ${gnomad_vcf} \
        -L ${exome_intervals} \
        -O wes/normal.pileups.table \
        --java-options "-Xmx32g"
    
    echo "Calculating contamination..."
    gatk CalculateContamination \
        -I wes/tumor.pileups.table \
        -matched wes/normal.pileups.table \
        -O wes/tumor.contamination.table \
        --java-options "-Xmx32g"
    
    echo "Running Mutect2..."
    gatk Mutect2 \
        -R ${ref_genome} \
        -I ${tumor_dedup_bam} \
        -I ${normal_dedup_bam} \
        -normal normal_sample_id \
        --panel-of-normals ${pon_vcf} \
        --germline-resource ${gnomad_vcf} \
        -L ${exome_intervals} \
        --f1r2-tar-gz wes/tumor.f1r2.tar.gz \
        -O wes/tumor_vs_normal.unfiltered.vcf.gz \
        --java-options "-Xmx32g"
    
    echo "Learning read orientation model..."
    gatk LearnReadOrientationModel \
        -I wes/tumor.f1r2.tar.gz \
        -O wes/tumor.read-orientation-model.tar.gz \
        --java-options "-Xmx32g"
    
    echo "Filtering Mutect2 calls..."
    gatk FilterMutectCalls \
        -R ${ref_genome} \
        -V wes/tumor_vs_normal.unfiltered.vcf.gz \
        --contamination-table wes/tumor.contamination.table \
        --ob-priors wes/tumor.read-orientation-model.tar.gz \
        -O wes/tumor_vs_normal.filtered.vcf.gz \
        --java-options "-Xmx32g"
    
    # 1.5 Truth set derivation
    echo "Creating truth set..."
    bcftools view -f PASS wes/tumor_vs_normal.filtered.vcf.gz \
        | bcftools filter -i 'FORMAT[AD][0:1]>=5 && FORMAT[DP]>=50 && INFO[TLOD]>=6' \
        | bcftools annotate -x INFO,^INFO/CSQ \
        -Oz -o wes/${sample}.truthset.vcf.gz
    
    echo "Indexing truth set VCF..."
    tabix -p vcf wes/${sample}.truthset.vcf.gz
    
    echo "Creating truth set BED..."
    bcftools query -f'%CHROM\\t%POS0\\t%END\\n' wes/${sample}.truthset.vcf.gz \
        | sort -k1,1 -k2,2n > wes/${sample}.truthset.bed
    
    # Clean up intermediate files
    rm -f wes/*.pileups.table wes/*.contamination.table wes/*.f1r2.tar.gz wes/*.read-orientation-model.tar.gz
    
    echo "Somatic variant calling completed for ${sample}"
    """
}
