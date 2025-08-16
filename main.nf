#!/usr/bin/env nextflow

// Tumor-Informed cfDNA MRD Pipeline
// Main workflow orchestrator

// Include modules
include { WES_PREPROCESSING; SOMATIC_VARIANT_CALLING } from './modules/preprocessing'
include { UMI_EXTRACTION; PLASMA_ALIGNMENT; UMI_CONSENSUS; PLASMA_QC; TRUTHSET_VARIANT_CALLING; FRAGMENTOMICS_FEATURES; METHYLATION_FEATURES; CNV_ANALYSIS } from './modules/plasma_processing'
include { BUILD_ERROR_MODEL; VALIDATE_ERROR_MODEL } from './modules/error_model'
include { FEATURE_EXTRACTION } from './modules/feature_extraction'
include { MRD_CLASSIFICATION } from './modules/mrd_classification'
include { VALIDATION } from './modules/validation'
include { REPORTING } from './modules/reporting'

// Input channels
tumor_ch = Channel
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

normal_ch = Channel
    .fromPath(params.normalDir + "/*_R{1,2}.fastq.gz")
    .map { file -> 
        def sample = file.name.replaceAll(/.*?(\w+)_T0_R{1,2}\.fastq\.gz/, '$1')
        [sample, file]
    }
    .groupTuple()
    .map { sample, files -> 
        def r1 = files.find { it.name.contains('_R1.fastq.gz') }
        def r2 = files.find { it.name.contains('_R2.fastq.gz') }
        [sample, r1, r2]
    }

// Plasma samples (patients) - multiple timepoints
plasma_patients_ch = Channel
    .fromPath(params.plasmaDir + "/*_T{0,1,2,3,4}_R{1,2}.fastq.gz")
    .map { file -> 
        def sample = file.name.replaceAll(/.*?(\w+)_T(\d+)_R[12]\.fastq\.gz/, '$1')
        def timepoint = file.name.replaceAll(/.*?(\w+)_T(\d+)_R[12]\.fastq\.gz/, 'T$2')
        [sample, timepoint, file]
    }
    .groupTuple()
    .map { sample, timepoint, files -> 
        def r1 = files.find { it.name.contains('_R1.fastq.gz') }
        def r2 = files.find { it.name.contains('_R2.fastq.gz') }
        [sample, timepoint, r1, r2]
    }

// Healthy plasma samples
healthy_plasma_ch = Channel
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

// Reference files
ref_genome = Channel.fromPath(params.genome)
exome_intervals = Channel.fromPath(params.exomeBed)
gnomad_vcf = Channel.fromPath(params.gnomadVcf)
pon_vcf = Channel.fromPath(params.ponVcf)
dmr_bed = Channel.fromPath(params.dmrBed)
tss_bed = Channel.fromPath(params.tssBed)

// Additional resources for Step 2
gc_wig = Channel.fromPath(params.gcWig)
map_wig = Channel.fromPath(params.mapWig)

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
    
    // Step 2: Plasma cfDNA Processing
    // UMI extraction for plasma samples
    UMI_EXTRACTION(plasma_patients_ch)
    
    // Align UMI-tagged reads
    PLASMA_ALIGNMENT(
        UMI_EXTRACTION.out.umi_fastq,
        ref_genome
    )
    
    // Generate UMI consensus
    UMI_CONSENSUS(PLASMA_ALIGNMENT.out.raw_bam)
    
    // Quality control for consensus BAMs
    PLASMA_QC(
        UMI_CONSENSUS.out.consensus_bam,
        ref_genome
    )
    
    // Variant calling at truth set loci
    TRUTHSET_VARIANT_CALLING(
        UMI_CONSENSUS.out.consensus_bam,
        ref_genome,
        SOMATIC_VARIANT_CALLING.out.truth_set_bed,
        gnomad_vcf,
        pon_vcf
    )
    
    // Extract fragmentomics features
    FRAGMENTOMICS_FEATURES(
        UMI_CONSENSUS.out.consensus_bam,
        ref_genome,
        tss_bed
    )
    
    // Optional: Methylation features (if EM-seq data available)
    // METHYLATION_FEATURES(...)
    
    // CNV analysis using ichorCNA
    CNV_ANALYSIS(
        PLASMA_QC.out.coverage_wig,
        gc_wig,
        map_wig
    )
    
    // Step 2.5: Background Error Model from Healthy Donors
    // Process healthy plasma samples for error model
    UMI_EXTRACTION(healthy_plasma_ch.map { sample, r1, r2 -> [sample, 'healthy', r1, r2] })
    
    PLASMA_ALIGNMENT(
        UMI_EXTRACTION.out.umi_fastq.filter { sample, type, r1, r2 -> type == 'healthy' },
        ref_genome
    )
    
    UMI_CONSENSUS(PLASMA_ALIGNMENT.out.raw_bam)
    
    // Build background error model
    BUILD_ERROR_MODEL(
        UMI_CONSENSUS.out.consensus_bam.filter { sample, type, bam, bai -> type == 'healthy' },
        SOMATIC_VARIANT_CALLING.out.truth_set_bed,
        ref_genome
    )
    
    // Validate error model
    VALIDATE_ERROR_MODEL(
        BUILD_ERROR_MODEL.out.error_model_json,
        BUILD_ERROR_MODEL.out.error_model_tsv
    )
    
    // Step 3: Feature Extraction (to be implemented)
    // FEATURE_EXTRACTION(...)
    
    // Step 4: MRD Classification (to be implemented)
    // MRD_CLASSIFICATION(...)
    
    // Step 5: Validation (to be implemented)
    // VALIDATION(...)
    
    // Step 6: Reporting (to be implemented)
    // REPORTING(...)
}

// Output channels for Step 1
workflow.out.truth_set = SOMATIC_VARIANT_CALLING.out.truth_set
workflow.out.truth_set_bed = SOMATIC_VARIANT_CALLING.out.truth_set_bed
workflow.out.filtered_variants = SOMATIC_VARIANT_CALLING.out.filtered_variants
workflow.out.wes_bams = WES_PREPROCESSING.out.dedup_bam
workflow.out.wes_metrics = WES_PREPROCESSING.out.hs_metrics

// Output channels for Step 2
workflow.out.plasma_consensus = UMI_CONSENSUS.out.consensus_bam
workflow.out.plasma_variants = TRUTHSET_VARIANT_CALLING.out.filtered_variants
workflow.out.plasma_variant_counts = TRUTHSET_VARIANT_CALLING.out.variant_counts
workflow.out.plasma_fragmentomics = FRAGMENTOMICS_FEATURES.out.fragmentomics
workflow.out.plasma_endmotifs = FRAGMENTOMICS_FEATURES.out.endmotifs
workflow.out.plasma_tss_coverage = FRAGMENTOMICS_FEATURES.out.tss_coverage
workflow.out.plasma_cnv = CNV_ANALYSIS.out.tumor_fraction
workflow.out.plasma_qc = PLASMA_QC.out.insert_metrics

// Output channels for Step 2.5
workflow.out.error_model = BUILD_ERROR_MODEL.out.error_model_json
workflow.out.error_model_tsv = BUILD_ERROR_MODEL.out.error_model_tsv
workflow.out.error_model_summary = BUILD_ERROR_MODEL.out.summary
workflow.out.error_model_validation = VALIDATE_ERROR_MODEL.out.validation
