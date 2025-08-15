#!/usr/bin/env nextflow

/*
 * Minimal Test for Tumor-Informed cfDNA MRD Pipeline
 * Tests basic channel creation and file detection
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

// Simple process
process PRINT_CHANNELS {
    tag "print"
    
    input:
    val tumor_data
    val normal_data
    path ref_genome
    path exome_intervals
    
    output:
    path "channel_info.txt", emit: info
    
    script:
    """
    echo "=== Channel Information ===" > channel_info.txt
    echo "Tumor data: ${tumor_data}" >> channel_info.txt
    echo "Normal data: ${normal_data}" >> channel_info.txt
    echo "Reference genome: ${ref_genome}" >> channel_info.txt
    echo "Exome intervals: ${exome_intervals}" >> channel_info.txt
    echo "Test completed successfully!" >> channel_info.txt
    """
}

workflow {
    // Collect and print channel information
    tumor_collected = test_tumor.collect()
    normal_collected = test_normal.collect()
    
    PRINT_CHANNELS(
        tumor_collected,
        normal_collected,
        ref_genome,
        exome_intervals
    )
}
