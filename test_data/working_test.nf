#!/usr/bin/env nextflow

/*
 * Working Test Workflow for Tumor-Informed cfDNA MRD Pipeline
 * Tests basic workflow structure and channel logic
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

// Mock process for tumor samples
process MOCK_TUMOR_PROCESS {
    tag "${sample}_tumor"
    
    publishDir "${params.outdir}/mock", mode: 'copy'
    
    input:
    tuple val(sample), path(fastq1), path(fastq2)
    path ref_genome
    path exome_intervals
    
    output:
    tuple val(sample), path("*.txt"), emit: results
    
    script:
    """
    echo "Processing tumor sample ${sample}" > ${sample}_tumor_result.txt
    echo "Input files: ${fastq1}, ${fastq2}" >> ${sample}_tumor_result.txt
    echo "Reference: ${ref_genome}" >> ${sample}_tumor_result.txt
    echo "Intervals: ${exome_intervals}" >> ${sample}_tumor_result.txt
    echo "Mock tumor processing completed successfully!" >> ${sample}_tumor_result.txt
    """
}

// Mock process for normal samples
process MOCK_NORMAL_PROCESS {
    tag "${sample}_normal"
    
    publishDir "${params.outdir}/mock", mode: 'copy'
    
    input:
    tuple val(sample), path(fastq1), path(fastq2)
    path ref_genome
    path exome_intervals
    
    output:
    tuple val(sample), path("*.txt"), emit: results
    
    script:
    """
    echo "Processing normal sample ${sample}" > ${sample}_normal_result.txt
    echo "Input files: ${fastq1}, ${fastq2}" >> ${sample}_normal_result.txt
    echo "Reference: ${ref_genome}" >> ${sample}_normal_result.txt
    echo "Intervals: ${exome_intervals}" >> ${sample}_normal_result.txt
    echo "Mock normal processing completed successfully!" >> ${sample}_normal_result.txt
    """
}

workflow {
    // Process tumor samples
    tumor_results = MOCK_TUMOR_PROCESS(
        test_tumor.map { sample, r1, r2 -> [sample, r1, r2] },
        ref_genome,
        exome_intervals
    )
    
    // Process normal samples
    normal_results = MOCK_NORMAL_PROCESS(
        test_normal.map { sample, r1, r2 -> [sample, r1, r2] },
        ref_genome,
        exome_intervals
    )
    
    // Combine results
    all_results = tumor_results.mix(normal_results)
    
    // Output channels
    workflow.out.results = all_results
}
