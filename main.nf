#!/usr/bin/env nextflow

/*
 * Tumor-Informed cfDNA MRD Pipeline
 * Main workflow for detecting Minimal Residual Disease in CRC patients
 */

// Import modules
include { WES_PREPROCESSING } from './modules/preprocessing'
include { SOMATIC_VARIANT_CALLING } from './modules/preprocessing'
include { FEATURE_EXTRACTION } from './modules/feature_extraction'
include { MRD_CLASSIFICATION } from './modules/mrd_classification'
include { VALIDATION } from './modules/validation'
include { REPORTING } from './modules/reporting'

// Input channels for WES samples
Channel
    .fromPath(params.tumorDir + "/*_R{1,2}.fastq.gz")
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
    .set { tumor_ch }

Channel
    .fromPath(params.normalDir + "/*_R{1,2}.fastq.gz")
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
    .set { normal_ch }

// Input channels for plasma samples (for later steps)
Channel
    .fromPath(params.plasmaDir + "/*_T*_R{1,2}.fastq.gz")
    .map { file -> 
        def sample = file.name.replaceAll(/.*?(\w+)_T(\d+)_R[12]\.fastq\.gz/, '$1\t$2')
        def parts = sample.split('\t')
        [parts[0], parts[1].toInteger(), file]
    }
    .groupTuple()
    .map { sample, timepoints, files -> 
        def timepointFiles = [:]
        timepoints.eachWithIndex { tp, i ->
            def r1 = files.find { it.name.contains("_T${tp}_R1.fastq.gz") }
            def r2 = files.find { it.name.contains("_T${tp}_R2.fastq.gz") }
            timepointFiles[tp] = [r1, r2]
        }
        [sample, timepointFiles]
    }
    .set { plasma_ch }

Channel
    .fromPath(params.healthyDir + "/*_R{1,2}.fastq.gz")
    .map { file -> 
        def sample = file.name.replaceAll(/.*?(\w+)_R[12]\.fastq\.gz/, '$1')
        [sample, file]
    }
    .groupTuple()
    .map { sample, files -> 
        def r1 = files.find { it.name.contains('_R1.fastq.gz') }
        def r2 = files.find { it.name.contains('_R2.fastq.gz') }
        [sample, r1, r2]
    }
    .set { healthy_ch }

// Reference files
ref_genome = Channel.fromPath(params.genome)
exome_intervals = Channel.fromPath(params.exomeBed)
gnomad_vcf = Channel.fromPath(params.gnomadVcf)
pon_vcf = Channel.fromPath(params.ponVcf)
dmr_bed = Channel.fromPath(params.dmrBed)
tss_bed = Channel.fromPath(params.tssBed)

// Main workflow
workflow {
    
    // Step 1: WES Preprocessing and Somatic Variant Calling
    // Process tumor samples
    WES_PREPROCESSING(
        tumor_ch.map { sample, r1, r2 -> [sample, 'tumor', r1, r2] },
        ref_genome,
        exome_intervals
    )
    
    // Process normal samples
    WES_PREPROCESSING(
        normal_ch.map { sample, r1, r2 -> [sample, 'normal', r1, r2] },
        ref_genome,
        exome_intervals
    )
    
    // Somatic variant calling (tumor vs normal)
    SOMATIC_VARIANT_CALLING(
        WES_PREPROCESSING.out.dedup_bam.filter { sample, type, bam, bai -> type == 'tumor' },
        WES_PREPROCESSING.out.dedup_bam.filter { sample, type, bam, bai -> type == 'normal' },
        ref_genome,
        exome_intervals,
        gnomad_vcf,
        pon_vcf
    )
    
    // Step 2: Feature Extraction (for plasma samples - to be implemented)
    // FEATURE_EXTRACTION(
    //     plasma_ch,
    //     SOMATIC_VARIANT_CALLING.out.truth_set,
    //     dmr_bed,
    //     tss_bed,
    //     ref_genome
    // )
    
    // Step 3: MRD Classification (to be implemented)
    // MRD_CLASSIFICATION(
    //     FEATURE_EXTRACTION.out.features,
    //     SOMATIC_VARIANT_CALLING.out.truth_set,
    //     plasma_ch
    // )
    
    // Step 4: Validation (to be implemented)
    // VALIDATION(
    //     MRD_CLASSIFICATION.out.mrd_scores,
    //     FEATURE_EXTRACTION.out.features,
    //     params.outdir
    // )
    
    // Step 5: Reporting (to be implemented)
    // REPORTING(
    //     SOMATIC_VARIANT_CALLING.out.truth_set,
    //     FEATURE_EXTRACTION.out.features,
    //     MRD_CLASSIFICATION.out.mrd_scores,
    //     VALIDATION.out.validation_metrics,
    //     params.outdir
    // )
}

// Output channels for Step 1
workflow.out.truth_set = SOMATIC_VARIANT_CALLING.out.truth_set
workflow.out.truth_set_bed = SOMATIC_VARIANT_CALLING.out.truth_set_bed
workflow.out.filtered_variants = SOMATIC_VARIANT_CALLING.out.filtered_variants
workflow.out.wes_bams = WES_PREPROCESSING.out.dedup_bam
workflow.out.wes_metrics = WES_PREPROCESSING.out.hs_metrics
