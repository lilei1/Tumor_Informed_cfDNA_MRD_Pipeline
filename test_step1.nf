#!/usr/bin/env nextflow

/*
 * Test script for Step 1: Tumor-Normal WES → Somatic Truth Set
 */

include { WES_PREPROCESSING } from './modules/preprocessing'
include { SOMATIC_VARIANT_CALLING } from './modules/preprocessing'

// Test data channels
test_tumor = Channel.of(['test_patient', 'test_tumor_R1.fastq.gz', 'test_tumor_R2.fastq.gz'])
test_normal = Channel.of(['test_patient', 'test_normal_R1.fastq.gz', 'test_normal_R2.fastq.gz'])

// Reference files
ref_genome = Channel.fromPath(params.genome)
exome_intervals = Channel.fromPath(params.exomeBed)

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
        Channel.fromPath(params.gnomadVcf),
        Channel.fromPath(params.ponVcf)
    )
}

workflow.out.truth_set = SOMATIC_VARIANT_CALLING.out.truth_set
workflow.out.truth_set_bed = SOMATIC_VARIANT_CALLING.out.truth_set_bed
