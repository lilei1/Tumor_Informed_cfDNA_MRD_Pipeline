#!/usr/bin/env nextflow

/*
 * Comprehensive Test Workflow for Full Pipeline including Step 5
 * Tests all steps: WES preprocessing, plasma processing, feature integration, 
 * MRD classification, and QC Gates & MultiQC
 */

// Include all modules
include { WES_PREPROCESSING; SOMATIC_VARIANT_CALLING } from '../modules/preprocessing'
include { UMI_EXTRACTION; PLASMA_ALIGNMENT; UMI_CONSENSUS; PLASMA_QC; TRUTHSET_VARIANT_CALLING; FRAGMENTOMICS_FEATURES; METHYLATION_FEATURES; CNV_ANALYSIS } from '../modules/plasma_processing'
include { BUILD_ERROR_MODEL; VALIDATE_ERROR_MODEL } from '../modules/error_model'
include { INTEGRATE_FEATURES; CALCULATE_LONGITUDINAL_MRD; GENERATE_MRD_REPORT } from '../modules/feature_integration_simple'
include { TRAIN_MRD_CLASSIFIER; CALIBRATE_THRESHOLDS; CLASSIFY_MRD_SAMPLES; GENERATE_CLINICAL_REPORT } from '../modules/mrd_classification_simple'
include { WES_QC_GATES; PLASMA_CONSENSUS_QC_GATES; PER_RUN_QC_GATES; MULTIQC_INTEGRATION } from '../modules/qc_gates_simple'

// Create comprehensive mock data
process CREATE_COMPREHENSIVE_MOCK_DATA {
    tag "comprehensive_mock_data"
    
    output:
    path "mock_tumor_fastq_R1.fastq.gz", emit: tumor_r1
    path "mock_tumor_fastq_R2.fastq.gz", emit: tumor_r2
    path "mock_normal_fastq_R1.fastq.gz", emit: normal_r1
    path "mock_normal_fastq_R2.fastq.gz", emit: normal_r2
    path "mock_plasma_fastq_R1.fastq.gz", emit: plasma_r1
    path "mock_plasma_fastq_R2.fastq.gz", emit: plasma_r2
    path "mock_healthy_fastq_R1.fastq.gz", emit: healthy_r1
    path "mock_healthy_fastq_R2.fastq.gz", emit: healthy_r2
    path "mock_ref_genome.fa", emit: ref_genome
    path "mock_exome_bed.bed", emit: exome_bed
    path "mock_gnomad.vcf.gz", emit: gnomad_vcf
    path "mock_pon.vcf.gz", emit: pon_vcf
    path "mock_dmr.bed", emit: dmr_bed
    path "mock_tss.bed", emit: tss_bed
    path "mock_gc.wig", emit: gc_wig
    path "mock_map.wig", emit: map_wig
    
    script:
    """
    #!/bin/bash
    set -e
    echo "Creating comprehensive mock data for full pipeline testing..."
    
    # Create mock FASTQ files
    echo "Mock tumor R1" | gzip > mock_tumor_fastq_R1.fastq.gz
    echo "Mock tumor R2" | gzip > mock_tumor_fastq_R2.fastq.gz
    echo "Mock normal R1" | gzip > mock_normal_fastq_R1.fastq.gz
    echo "Mock normal R2" | gzip > mock_normal_fastq_R2.fastq.gz
    echo "Mock plasma R1" | gzip > mock_plasma_fastq_R1.fastq.gz
    echo "Mock plasma R2" | gzip > mock_plasma_fastq_R2.fastq.gz
    echo "Mock healthy R1" | gzip > mock_healthy_fastq_R1.fastq.gz
    echo "Mock healthy R2" | gzip > mock_healthy_fastq_R2.fastq.gz
    
    # Create mock reference files
    echo ">chr1" > mock_ref_genome.fa
    echo "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT" >> mock_ref_genome.fa
    
    echo "chr1\t1000\t2000" > mock_exome_bed.bed
    echo "chr1\t3000\t4000" >> mock_exome_bed.bed
    
    echo "##fileformat=VCFv4.2" > mock_gnomad.vcf
    echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> mock_gnomad.vcf
    echo "chr1\t1500\t.\tA\tT\t100\tPASS\tAF=0.1" >> mock_gnomad.vcf
    
    echo "##fileformat=VCFv4.2" > mock_pon.vcf
    echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> mock_pon.vcf
    echo "chr1\t1500\t.\tA\tT\t100\tPASS\tAF=0.05" >> mock_pon.vcf
    
    echo "chr1\t1000\t2000" > mock_dmr.bed
    echo "chr1\t3000\t4000" >> mock_dmr.bed
    
    echo "chr1\t1000\t2000" > mock_tss.bed
    echo "chr1\t3000\t4000" >> mock_tss.bed
    
    echo "fixedStep chrom=chr1 start=1000 step=100" > mock_gc.wig
    echo "0.5" >> mock_gc.wig
    echo "0.6" >> mock_gc.wig
    
    echo "fixedStep chrom=chr1 start=1000 step=100" > mock_map.wig
    echo "1.0" >> mock_map.wig
    echo "1.0" >> mock_map.wig
    
    echo "Comprehensive mock data creation completed"
    """
}

workflow {
    // Create comprehensive mock data
    CREATE_COMPREHENSIVE_MOCK_DATA()
    
    // Step 1: WES Preprocessing and Somatic Variant Calling
    WES_PREPROCESSING(
        Channel.of(['patient_01', CREATE_COMPREHENSIVE_MOCK_DATA.out.tumor_r1, CREATE_COMPREHENSIVE_MOCK_DATA.out.tumor_r2]),
        Channel.of(['patient_01', CREATE_COMPREHENSIVE_MOCK_DATA.out.normal_r1, CREATE_COMPREHENSIVE_MOCK_DATA.out.normal_r2]),
        CREATE_COMPREHENSIVE_MOCK_DATA.out.ref_genome,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.exome_bed
    )
    
    SOMATIC_VARIANT_CALLING(
        WES_PREPROCESSING.out.dedup_bam,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.ref_genome,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.exome_bed,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.gnomad_vcf,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.pon_vcf
    )
    
    // Step 2: Plasma cfDNA Processing
    UMI_EXTRACTION(
        Channel.of(['patient_01', 'T0', CREATE_COMPREHENSIVE_MOCK_DATA.out.plasma_r1, CREATE_COMPREHENSIVE_MOCK_DATA.out.plasma_r2])
    )
    
    PLASMA_ALIGNMENT(
        UMI_EXTRACTION.out.umi_fastq,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.ref_genome
    )
    
    UMI_CONSENSUS(
        PLASMA_ALIGNMENT.out.raw_bam
    )
    
    PLASMA_QC(
        UMI_CONSENSUS.out.consensus_bam,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.ref_genome
    )
    
    TRUTHSET_VARIANT_CALLING(
        UMI_CONSENSUS.out.consensus_bam,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.ref_genome,
        SOMATIC_VARIANT_CALLING.out.truth_set_bed,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.gnomad_vcf,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.pon_vcf
    )
    
    FRAGMENTOMICS_FEATURES(
        UMI_CONSENSUS.out.consensus_bam,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.ref_genome,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.tss_bed
    )
    
    CNV_ANALYSIS(
        PLASMA_QC.out.coverage_wig,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.gc_wig,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.map_wig
    )
    
    // Step 2.5: Background Error Model from Healthy Donors
    UMI_EXTRACTION(
        Channel.of(['healthy_01', 'healthy', CREATE_COMPREHENSIVE_MOCK_DATA.out.healthy_r1, CREATE_COMPREHENSIVE_MOCK_DATA.out.healthy_r2])
    )
    
    PLASMA_ALIGNMENT(
        UMI_EXTRACTION.out.umi_fastq.filter { sample, type, r1, r2 -> type == 'healthy' },
        CREATE_COMPREHENSIVE_MOCK_DATA.out.ref_genome
    )
    
    UMI_CONSENSUS(
        PLASMA_ALIGNMENT.out.raw_bam
    )
    
    BUILD_ERROR_MODEL(
        UMI_CONSENSUS.out.consensus_bam.filter { sample, type, bam, bai -> type == 'healthy' },
        SOMATIC_VARIANT_CALLING.out.truth_set_bed,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.ref_genome
    )
    
    VALIDATE_ERROR_MODEL(
        BUILD_ERROR_MODEL.out.error_model_json,
        BUILD_ERROR_MODEL.out.error_model_tsv
    )
    
    // Step 3: Feature Integration with Error-Aware Scoring
    INTEGRATE_FEATURES(
        TRUTHSET_VARIANT_CALLING.out.filtered_variants
            .join(FRAGMENTOMICS_FEATURES.out.fragmentomics)
            .join(FRAGMENTOMICS_FEATURES.out.endmotifs)
            .join(FRAGMENTOMICS_FEATURES.out.tss_coverage)
            .join(CNV_ANALYSIS.out.tumor_fraction)
            .join(BUILD_ERROR_MODEL.out.error_model_json)
            .map { sample, timepoint, variant, fragmentomics, endmotifs, tss_coverage, cnv_results, error_model ->
                [sample, timepoint, variant, fragmentomics, endmotifs, tss_coverage, cnv_results, error_model]
            },
        SOMATIC_VARIANT_CALLING.out.truth_set_bed,
        CREATE_COMPREHENSIVE_MOCK_DATA.out.ref_genome
    )
    
    CALCULATE_LONGITUDINAL_MRD(
        INTEGRATE_FEATURES.out.mrd_scores
            .groupTuple()
            .map { sample, scores -> [sample, scores.flatten()] }
    )
    
    GENERATE_MRD_REPORT(
        INTEGRATE_FEATURES.out.integrated_features.collect(),
        INTEGRATE_FEATURES.out.mrd_scores.collect(),
        CALCULATE_LONGITUDINAL_MRD.out.longitudinal_mrd.collect(),
        INTEGRATE_FEATURES.out.feature_importance.collect()
    )
    
    // Step 4: MRD Classification and Threshold Calibration
    TRAIN_MRD_CLASSIFIER(
        INTEGRATE_FEATURES.out.integrated_features.collect(),
        INTEGRATE_FEATURES.out.mrd_scores.collect(),
        INTEGRATE_FEATURES.out.integrated_features.collect(),
        INTEGRATE_FEATURES.out.mrd_scores.collect()
    )
    
    CALIBRATE_THRESHOLDS(
        TRAIN_MRD_CLASSIFIER.out.trained_model,
        INTEGRATE_FEATURES.out.integrated_features.collect(),
        INTEGRATE_FEATURES.out.mrd_scores.collect()
    )
    
    CLASSIFY_MRD_SAMPLES(
        TRAIN_MRD_CLASSIFIER.out.trained_model,
        CALIBRATE_THRESHOLDS.out.calibrated_thresholds,
        INTEGRATE_FEATURES.out.integrated_features.collect()
    )
    
    GENERATE_CLINICAL_REPORT(
        CLASSIFY_MRD_SAMPLES.out.mrd_classification,
        CLASSIFY_MRD_SAMPLES.out.confidence_scores,
        CLASSIFY_MRD_SAMPLES.out.risk_assessment,
        TRAIN_MRD_CLASSIFIER.out.feature_importance.collect()
    )
    
    // Step 5: QC Gates & MultiQC
    WES_QC_GATES(
        WES_PREPROCESSING.out.hs_metrics.filter { sample, type, metrics -> type == 'tumor' },
        WES_PREPROCESSING.out.hs_metrics.filter { sample, type, metrics -> type == 'normal' },
        WES_PREPROCESSING.out.dup_metrics.filter { sample, type, metrics -> type == 'tumor' },
        WES_PREPROCESSING.out.dup_metrics.filter { sample, type, metrics -> type == 'normal' },
        WES_PREPROCESSING.out.depth_stats.filter { sample, type, stats -> type == 'tumor' },
        WES_PREPROCESSING.out.depth_stats.filter { sample, type, stats -> type == 'normal' }
    )
    
    PLASMA_CONSENSUS_QC_GATES(
        UMI_CONSENSUS.out.consensus_bam,
        UMI_CONSENSUS.out.umi_stats,
        PLASMA_QC.out.coverage_wig,
        PLASMA_QC.out.insert_metrics
    )
    
    PER_RUN_QC_GATES(
        SOMATIC_VARIANT_CALLING.out.contamination_table,
        SOMATIC_VARIANT_CALLING.out.sex_check_results,
        SOMATIC_VARIANT_CALLING.out.pon_fp_results
    )
    
    MULTIQC_INTEGRATION(
        WES_QC_GATES.out.qc_summary.collect(),
        PLASMA_CONSENSUS_QC_GATES.out.qc_summary.collect(),
        PER_RUN_QC_GATES.out.qc_summary.collect()
    )
}
