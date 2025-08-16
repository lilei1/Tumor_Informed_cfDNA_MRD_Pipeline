#!/usr/bin/env nextflow

// Tumor-Informed cfDNA MRD Pipeline
// Main workflow orchestrator

// Include modules
include { WES_PREPROCESSING; SOMATIC_VARIANT_CALLING } from './modules/preprocessing'
include { UMI_EXTRACTION; PLASMA_ALIGNMENT; UMI_CONSENSUS; PLASMA_QC; TRUTHSET_VARIANT_CALLING; FRAGMENTOMICS_FEATURES; METHYLATION_FEATURES; CNV_ANALYSIS } from './modules/plasma_processing'
include { BUILD_ERROR_MODEL; VALIDATE_ERROR_MODEL } from './modules/error_model'
include { INTEGRATE_FEATURES; CALCULATE_LONGITUDINAL_MRD; GENERATE_MRD_REPORT } from './modules/feature_integration'

// Input channels
tumor_ch = Channel.fromPath(params.tumorDir).map { file ->
    def sample = file.name.replaceAll(/.*?(\w+)_T0_R[12]\.fastq\.gz/, '$1')
    [sample, file]
}.groupTuple()

normal_ch = Channel.fromPath(params.normalDir).map { file ->
    def sample = file.name.replaceAll(/.*?(\w+)_T0_R[12]\.fastq\.gz/, '$1')
    [sample, file]
}.groupTuple()

plasma_patients_ch = Channel.fromPath(params.plasmaDir).map { file ->
    def sample = file.name.replaceAll(/.*?(\w+)_T(\d+)_R[12]\.fastq\.gz/, '$1\t$2')
    def parts = sample.split('\t')
    [parts[0], parts[1], file]
}.groupTuple().map { sample, timepoint, files ->
    def r1 = files.find { it.name.contains('_R1.fastq.gz') }
    def r2 = files.find { it.name.contains('_R2.fastq.gz') }
    [sample, timepoint, r1, r2]
}

healthy_plasma_ch = Channel.fromPath(params.healthyDir).map { file ->
    def sample = file.name.replaceAll(/.*?(\w+)_R[12]\.fastq\.gz/, '$1')
    [sample, file]
}.groupTuple().map { sample, files ->
    def r1 = files.find { it.name.contains('_R1.fastq.gz') }
    def r2 = files.find { it.name.contains('_R2.fastq.gz') }
    [sample, r1, r2]
}

// Reference and resource channels
ref_genome = Channel.fromPath(params.genome)
exome_bed = Channel.fromPath(params.exomeBed)
gnomad_vcf = Channel.fromPath(params.gnomadVcf)
pon_vcf = Channel.fromPath(params.ponVcf)
dmr_bed = Channel.fromPath(params.dmrBed)
tss_bed = Channel.fromPath(params.tssBed)
gc_wig = Channel.fromPath(params.gcWig)
map_wig = Channel.fromPath(params.mapWig)

workflow {
    // Step 1: WES Preprocessing and Somatic Variant Calling
    WES_PREPROCESSING(tumor_ch, normal_ch, ref_genome, exome_bed, gnomad_vcf, pon_vcf)
    SOMATIC_VARIANT_CALLING(WES_PREPROCESSING.out.dedup_bam, ref_genome, exome_bed, gnomad_vcf, pon_vcf)
    
    // Step 2: Plasma cfDNA Processing
    UMI_EXTRACTION(plasma_patients_ch)
    PLASMA_ALIGNMENT(UMI_EXTRACTION.out.umi_fastq, ref_genome)
    UMI_CONSENSUS(PLASMA_ALIGNMENT.out.raw_bam)
    PLASMA_QC(UMI_CONSENSUS.out.consensus_bam, ref_genome)
    TRUTHSET_VARIANT_CALLING(UMI_CONSENSUS.out.consensus_bam, ref_genome, SOMATIC_VARIANT_CALLING.out.truth_set_bed, gnomad_vcf, pon_vcf)
    FRAGMENTOMICS_FEATURES(UMI_CONSENSUS.out.consensus_bam, ref_genome, tss_bed)
    CNV_ANALYSIS(PLASMA_QC.out.coverage_wig, gc_wig, map_wig)
    
    // Step 2.5: Background Error Model from Healthy Donors
    UMI_EXTRACTION(healthy_plasma_ch.map { sample, r1, r2 -> [sample, 'healthy', r1, r2] })
    PLASMA_ALIGNMENT(UMI_EXTRACTION.out.umi_fastq.filter { sample, type, r1, r2 -> type == 'healthy' }, ref_genome)
    UMI_CONSENSUS(PLASMA_ALIGNMENT.out.raw_bam)
    BUILD_ERROR_MODEL(UMI_CONSENSUS.out.consensus_bam.filter { sample, type, bam, bai -> type == 'healthy' }, SOMATIC_VARIANT_CALLING.out.truth_set_bed, ref_genome)
    VALIDATE_ERROR_MODEL(BUILD_ERROR_MODEL.out.error_model_json, BUILD_ERROR_MODEL.out.error_model_tsv)
    
    // Step 3: Feature Integration with Error-Aware Scoring
    // Prepare inputs for feature integration
    feature_inputs = TRUTHSET_VARIANT_CALLING.out.filtered_variants
        .join(FRAGMENTOMICS_FEATURES.out.fragmentomics)
        .join(FRAGMENTOMICS_FEATURES.out.endmotifs)
        .join(FRAGMENTOMICS_FEATURES.out.tss_coverage)
        .join(CNV_ANALYSIS.out.tumor_fraction)
        .join(BUILD_ERROR_MODEL.out.error_model_json)
        .map { sample, timepoint, variant, fragmentomics, endmotifs, tss_coverage, cnv_results, error_model ->
            [sample, timepoint, variant, fragmentomics, endmotifs, tss_coverage, cnv_results, error_model]
        }
    
    // Integrate features with error-aware scoring
    INTEGRATE_FEATURES(
        feature_inputs,
        SOMATIC_VARIANT_CALLING.out.truth_set_bed,
        ref_genome
    )
    
    // Calculate longitudinal MRD trends
    CALCULATE_LONGITUDINAL_MRD(
        INTEGRATE_FEATURES.out.mrd_scores
            .groupTuple()
            .map { sample, scores -> [sample, scores.flatten()] }
    )
    
    // Generate comprehensive MRD reports
    GENERATE_MRD_REPORT(
        INTEGRATE_FEATURES.out.integrated_features.collect(),
        INTEGRATE_FEATURES.out.mrd_scores.collect(),
        CALCULATE_LONGITUDINAL_MRD.out.longitudinal_mrd.collect(),
        INTEGRATE_FEATURES.out.feature_importance.collect()
    )
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

// Output channels for Step 3
workflow.out.integrated_features = INTEGRATE_FEATURES.out.integrated_features
workflow.out.mrd_scores = INTEGRATE_FEATURES.out.mrd_scores
workflow.out.feature_importance = INTEGRATE_FEATURES.out.feature_importance
workflow.out.longitudinal_mrd = CALCULATE_LONGITUDINAL_MRD.out.longitudinal_mrd
workflow.out.mrd_trends = CALCULATE_LONGITUDINAL_MRD.out.mrd_trend
workflow.out.mrd_reports = GENERATE_MRD_REPORT.out.mrd_report
workflow.out.mrd_reports_pdf = GENERATE_MRD_REPORT.out.mrd_report_pdf
