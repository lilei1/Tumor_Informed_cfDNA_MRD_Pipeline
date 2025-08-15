#!/usr/bin/env nextflow

/*
 * Simplified Test Workflow for Tumor-Informed cfDNA MRD Pipeline
 * Tests workflow structure and channel logic without full tool execution
 */

// Test data channels
test_tumor = Channel
    .fromPath("wes/tumor/*_R{1,2}.fastq.gz")
    .map { file -> 
        def sample = file.name.replaceAll(/.*?(\w+)_T0_R[12]\.fastq\.gz/, '$1')
        [sample, file]
    }
    .groupTuple()
    .map { sample, files -> 
        def r1 = files.find { it.name.contains('_R1.fastq.gz') }
        def r2 = files.find { it.name.contains('_R2.fastq.gz') }
        [sample, r1, r2]
    }

test_normal = Channel
    .fromPath("wes/normal/*_R{1,2}.fastq.gz")
    .map { file -> 
        def sample = file.name.replaceAll(/.*?(\w+)_T0_R[12]\.fastq\.gz/, '$1')
        [sample, file]
    }
    .groupTuple()
    .map { sample, files -> 
        def r1 = files.find { it.name.contains('_R1.fastq.gz') }
        def r2 = files.find { it.name.contains('_R2.fastq.gz') }
        [sample, r1, r2]
    }

// Reference files
ref_genome = Channel.fromPath(params.genome)
exome_intervals = Channel.fromPath(params.exomeBed)

// Mock processes for testing workflow structure
process MOCK_WES_PREPROCESSING {
    tag "${sample}_${type}"
    
    publishDir "${params.outdir}/wes", mode: 'copy'
    
    input:
    tuple val(sample), val(type), path(fastq1), path(fastq2)
    path ref_genome
    path exome_intervals
    
    output:
    tuple val(sample), val(type), path("*.mock.bam"), path("*.mock.bai"), emit: bam
    tuple val(sample), val(type), path("*.mock.dedup.bam"), path("*.mock.bai"), emit: dedup_bam
    tuple val(sample), val(type), path("*.mock.hs.txt"), emit: hs_metrics
    tuple val(sample), val(type), path("*.mock.dup.txt"), emit: dup_metrics
    
    script:
    """
    #!/bin/bash
    echo "Mock WES preprocessing for ${sample} ${type}"
    
    # Create mock output files
    echo "Mock BAM content" > ${sample}_${type}_T0.mock.bam
    echo "Mock BAM index" > ${sample}_${type}_T0.mock.bai
    echo "Mock dedup BAM" > ${sample}_${type}_T0.mock.dedup.bam
    echo "Mock dedup BAM index" > ${sample}_${type}_T0.mock.bai
    
    # Create mock metrics
    echo "Mock HS metrics for ${sample} ${type}" > ${sample}.mock.hs.txt
    echo "Mock duplicate metrics for ${sample} ${type}" > ${sample}.mock.dup.txt
    
    echo "Mock WES preprocessing completed for ${sample} ${type}"
    """
}

process MOCK_SOMATIC_VARIANT_CALLING {
    tag "${sample}"
    
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    tuple val(sample), path(tumor_dedup_bam), path(tumor_dedup_bai)
    tuple val(sample), path(normal_dedup_bam), path(normal_dedup_bai)
    path ref_genome
    path exome_intervals
    path gnomad_vcf
    path pon_vcf
    
    output:
    tuple val(sample), path("*.mock.truthset.vcf.gz"), path("*.mock.truthset.vcf.gz.tbi"), emit: truth_set
    tuple val(sample), path("*.mock.truthset.bed"), emit: truth_set_bed
    tuple val(sample), path("*.mock.filtered.vcf.gz"), path("*.mock.filtered.vcf.gz.tbi"), emit: filtered_variants
    
    script:
    """
    #!/bin/bash
    echo "Mock somatic variant calling for ${sample}"
    
    # Create mock output files
    echo "Mock truth set VCF" > ${sample}.mock.truthset.vcf.gz
    echo "Mock truth set VCF index" > ${sample}.mock.truthset.vcf.gz.tbi
    echo "Mock truth set BED" > ${sample}.mock.truthset.bed
    echo "Mock filtered variants VCF" > ${sample}.mock.filtered.vcf.gz
    echo "Mock filtered variants VCF index" > ${sample}.mock.filtered.vcf.gz.tbi
    
    echo "Mock somatic variant calling completed for ${sample}"
    """
}

workflow {
    // Test WES preprocessing
    MOCK_WES_PREPROCESSING(
        test_tumor.map { sample, r1, r2 -> [sample, 'tumor', r1, r2] },
        ref_genome,
        exome_intervals
    )
    
    MOCK_WES_PREPROCESSING(
        test_normal.map { sample, r1, r2 -> [sample, 'normal', r1, r2] },
        ref_genome,
        exome_intervals
    )
    
    // Test somatic variant calling
    MOCK_SOMATIC_VARIANT_CALLING(
        MOCK_WES_PREPROCESSING.out.dedup_bam.filter { sample, type, bam, bai -> type == 'tumor' },
        MOCK_WES_PREPROCESSING.out.dedup_bam.filter { sample, type, bam, bai -> type == 'normal' },
        ref_genome,
        exome_intervals,
        Channel.fromPath(params.gnomadVcf),
        Channel.fromPath(params.ponVcf)
    )
}

// Output channels - define after workflow
workflow.out.truth_set = MOCK_SOMATIC_VARIANT_CALLING.out.truth_set
workflow.out.truth_set_bed = MOCK_SOMATIC_VARIANT_CALLING.out.truth_set_bed
workflow.out.filtered_variants = MOCK_SOMATIC_VARIANT_CALLING.out.filtered_variants
workflow.out.wes_bams = MOCK_WES_PREPROCESSING.out.dedup_bam
workflow.out.wes_metrics = MOCK_WES_PREPROCESSING.out.hs_metrics
