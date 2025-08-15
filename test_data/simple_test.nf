#!/usr/bin/env nextflow

/*
 * Simple Test Workflow for Tumor-Informed cfDNA MRD Pipeline
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

// Simple mock process
process MOCK_PROCESS {
    tag "${sample}_${type}"
    
    publishDir "${params.outdir}/mock", mode: 'copy'
    
    input:
    tuple val(sample), val(type), path(fastq1), path(fastq2)
    path ref_genome
    path exome_intervals
    
    output:
    tuple val(sample), val(type), path("*.txt"), emit: results
    
    script:
    """
    echo "Processing ${sample} ${type}" > ${sample}_${type}_result.txt
    echo "Input files: ${fastq1}, ${fastq2}" >> ${sample}_${type}_result.txt
    echo "Reference: ${ref_genome}" >> ${sample}_${type}_result.txt
    echo "Intervals: ${exome_intervals}" >> ${sample}_${type}_result.txt
    echo "Mock processing completed successfully!" >> ${sample}_${type}_result.txt
    """
}

workflow {
    // Process tumor samples
    MOCK_PROCESS(
        test_tumor.map { sample, r1, r2 -> [sample, 'tumor', r1, r2] },
        ref_genome,
        exome_intervals
    )
    
    // Process normal samples
    MOCK_PROCESS(
        test_normal.map { sample, r1, r2 -> [sample, 'normal', r1, r2] },
        ref_genome,
        exome_intervals
    )
}

// Output channels
workflow.out.results = MOCK_PROCESS.out.results
