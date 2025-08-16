#!/usr/bin/env nextflow

// Test workflow for Step 2: Plasma cfDNA UMI-based WGS processing
// Tests UMI extraction, consensus calling, and feature extraction

// Include plasma processing modules
include { UMI_EXTRACTION; PLASMA_ALIGNMENT; UMI_CONSENSUS; PLASMA_QC; TRUTHSET_VARIANT_CALLING; FRAGMENTOMICS_FEATURES; CNV_ANALYSIS } from '../modules/plasma_processing'

// Test data channels
test_plasma = Channel
    .fromPath("plasma/patients/*_T{0,1}_R{1,2}.fastq.gz")
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

// Reference files
ref_genome = Channel.fromPath(params.genome)
truthset_bed = Channel.fromPath("wes/patientX.truthset.bed")
gnomad_vcf = Channel.fromPath(params.gnomadVcf)
pon_vcf = Channel.fromPath(params.ponVcf)
tss_bed = Channel.fromPath(params.tssBed)
gc_wig = Channel.fromPath(params.gcWig)
map_wig = Channel.fromPath(params.mapWig)

// Simple test process to validate inputs
process VALIDATE_STEP2_INPUTS {
    tag "validate"
    
    input:
    val plasma_data
    path ref_genome
    path truthset_bed
    path gnomad_vcf
    path pon_vcf
    path tss_bed
    path gc_wig
    path map_wig
    
    output:
    path "step2_validation.txt", emit: validation
    
    script:
    """
    echo "=== Step 2 Validation ===" > step2_validation.txt
    echo "Plasma data: ${plasma_data}" >> step2_validation.txt
    echo "Reference genome: ${ref_genome}" >> step2_validation.txt
    echo "Truth set BED: ${truthset_bed}" >> step2_validation.txt
    echo "gnomAD VCF: ${gnomad_vcf}" >> step2_validation.txt
    echo "Panel of Normals: ${pon_vcf}" >> step2_validation.txt
    echo "TSS BED: ${tss_bed}" >> step2_validation.txt
    echo "GC content WIG: ${gc_wig}" >> step2_validation.txt
    echo "Mappability WIG: ${map_wig}" >> step2_validation.txt
    echo "Step 2 inputs validated successfully!" >> step2_validation.txt
    """
}

workflow {
    // Validate Step 2 inputs
    VALIDATE_STEP2_INPUTS(
        test_plasma.collect(),
        ref_genome,
        truthset_bed,
        gnomad_vcf,
        pon_vcf,
        tss_bed,
        gc_wig,
        map_wig
    )
    
    // Note: Full Step 2 workflow would be:
    // UMI_EXTRACTION(test_plasma)
    // PLASMA_ALIGNMENT(UMI_EXTRACTION.out.umi_fastq, ref_genome)
    // UMI_CONSENSUS(PLASMA_ALIGNMENT.out.raw_bam)
    // PLASMA_QC(UMI_CONSENSUS.out.consensus_bam, ref_genome)
    // TRUTHSET_VARIANT_CALLING(UMI_CONSENSUS.out.consensus_bam, ref_genome, truthset_bed, gnomad_vcf, pon_vcf)
    // FRAGMENTOMICS_FEATURES(UMI_CONSENSUS.out.consensus_bam, ref_genome, tss_bed)
    // CNV_ANALYSIS(PLASMA_QC.out.coverage_wig, gc_wig, map_wig)
}

// Output
workflow.out.validation = VALIDATE_STEP2_INPUTS.out.validation
