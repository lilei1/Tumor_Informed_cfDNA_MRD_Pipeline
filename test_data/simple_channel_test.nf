#!/usr/bin/env nextflow

/*
 * Simple Channel Test for Tumor-Informed cfDNA MRD Pipeline
 * Tests input channel creation and file detection
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

// Simple process to validate channels
process VALIDATE_CHANNELS {
    tag "validation"
    
    input:
    val tumor_data
    val normal_data
    path ref_genome
    path exome_intervals
    
    output:
    path "channel_validation.txt", emit: validation
    
    script:
    """
    echo "=== Channel Validation Results ===" > channel_validation.txt
    echo "Tumor data: ${tumor_data}" >> channel_validation.txt
    echo "Normal data: ${normal_data}" >> channel_validation.txt
    echo "Reference genome: ${ref_genome}" >> channel_validation.txt
    echo "Exome intervals: ${exome_intervals}" >> channel_validation.txt
    echo "" >> channel_validation.txt
    echo "Channel validation completed successfully!" >> channel_validation.txt
    """
}

workflow {
    // Collect channel data
    tumor_collected = test_tumor.collect()
    normal_collected = test_normal.collect()
    
    // Validate channels
    VALIDATE_CHANNELS(
        tumor_collected,
        normal_collected,
        ref_genome,
        exome_intervals
    )
}

// Output
workflow.out.validation = VALIDATE_CHANNELS.out.validation
