#!/usr/bin/env nextflow

/*
 * Variant Calling Module
 * Tumor-informed somatic variant detection with ultra-low VAF sensitivity
 */

process VARIANT_CALLING {
    tag "${sample}"
    
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    cpus 16
    memory '64 GB'
    time '8h'
    
    input:
    tuple val(sample), path(tumor_bam), path(tumor_bai)
    tuple val(sample), path(normal_bam), path(normal_bai)
    path ref_genome
    path exome_bed
    path gnomad_vcf
    path pon_vcf
    
    output:
    tuple val(sample), path("*.truthset.vcf.gz"), path("*.truthset.vcf.gz.tbi"), emit: truth_set
    tuple val(sample), path("*.somatic.vcf.gz"), path("*.somatic.vcf.gz.tbi"), emit: somatic_variants
    path "*.log", emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    
    # Create tumor-normal pair VCF
    gatk Mutect2 \
        --reference ${ref_genome} \
        --tumor ${tumor_bam} \
        --normal ${normal_bam} \
        --intervals ${exome_bed} \
        --germline-resource ${gnomad_vcf} \
        --panel-of-normals ${pon_vcf} \
        --output ${sample}.somatic.vcf.gz \
        --f1r2-tar-gz ${sample}.f1r2.tar.gz \
        --bam-output ${sample}.tumor-normal.bam \
        --java-options "-Xmx32g"
    
    # Learn read orientation model
    gatk LearnReadOrientationModel \
        --input ${sample}.f1r2.tar.gz \
        --output ${sample}.read-orientation-model.tar.gz
    
    # Calculate contamination
    gatk CalculateContamination \
        --input ${sample}.f1r2.tar.gz \
        --tumor-segmentation ${sample}.segments.table \
        --output ${sample}.contamination.table
    
    # Filter variants
    gatk FilterMutectCalls \
        --reference ${ref_genome} \
        --variant ${sample}.somatic.vcf.gz \
        --contamination-table ${sample}.contamination.table \
        --orientation-bias-artifact-priors ${sample}.read-orientation-model.tar.gz \
        --output ${sample}.filtered.vcf.gz
    
    # Create high-confidence truth set
    gatk SelectVariants \
        --reference ${ref_genome} \
        --variant ${sample}.filtered.vcf.gz \
        --output ${sample}.truthset.vcf.gz \
        --select-type-to-include SNP \
        --select-type-to-include INDEL \
        --restrict-alleles-to BIALLELIC \
        --exclude-filtered true \
        --exclude-non-variants true
    
    # Index VCF
    tabix -p vcf ${sample}.truthset.vcf.gz
    
    # Generate summary statistics
    gatk VariantsToTable \
        --reference ${ref_genome} \
        --variant ${sample}.truthset.vcf.gz \
        --output ${sample}.truthset.table \
        -F CHROM -F POS -F REF -F ALT -F QUAL -F FILTER -F AF -F DP
    
    # Create BED file of truth set loci
    gatk VariantsToTable \
        --reference ${ref_genome} \
        --variant ${sample}.truthset.vcf.gz \
        --output ${sample}.truthset.bed \
        -F CHROM -F POS -F POS -F REF -F ALT \
        --split-multi-allelic
    
    # Clean up intermediate files
    rm -f *.f1r2.tar.gz *.read-orientation-model.tar.gz *.contamination.table *.segments.table
    """
}
