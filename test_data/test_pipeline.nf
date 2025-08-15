#!/usr/bin/env nextflow

/*
 * Test Pipeline for Tumor-Informed cfDNA MRD Pipeline
 * Tests Step 1: Tumor-Normal WES → Somatic Truth Set
 */

include { WES_PREPROCESSING } from '../modules/preprocessing'
include { SOMATIC_VARIANT_CALLING } from '../modules/preprocessing'

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
gnomad_vcf = Channel.fromPath(params.gnomadVcf)
pon_vcf = Channel.fromPath(params.ponVcf)

workflow {
    // Test WES preprocessing
    WES_PREPROCESSING(
        test_tumor.map { sample, r1, r2 -> [sample, 'tumor', r1, r2] },
        ref_genome,
        exome_intervals
    )
    
    WES_PREPROCESSING(
        test_normal.map { sample, r1, r2 -> [sample, 'normal', r1, r2] },
        ref_genome,
        exome_intervals
    )
    
    // Test somatic variant calling
    SOMATIC_VARIANT_CALLING(
        WES_PREPROCESSING.out.dedup_bam.filter { sample, type, bam, bai -> type == 'tumor' },
        WES_PREPROCESSING.out.dedup_bam.filter { sample, type, bam, bai -> type == 'normal' },
        ref_genome,
        exome_intervals,
        gnomad_vcf,
        pon_vcf
    )
}

// Output channels
workflow.out.truth_set = SOMATIC_VARIANT_CALLING.out.truth_set
workflow.out.truth_set_bed = SOMATIC_VARIANT_CALLING.out.truth_set_bed
workflow.out.filtered_variants = SOMATIC_VARIANT_CALLING.out.filtered_variants
workflow.out.wes_bams = WES_PREPROCESSING.out.dedup_bam
workflow.out.wes_metrics = WES_PREPROCESSING.out.hs_metrics
