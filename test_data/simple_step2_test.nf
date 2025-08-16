#!/usr/bin/env nextflow

// Simple test workflow for Step 2: Plasma cfDNA processing
// Tests input validation and channel creation

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
tss_bed = Channel.fromPath(params.tssBed)
gc_wig = Channel.fromPath(params.gcWig)
map_wig = Channel.fromPath(params.mapWig)

// Simple validation process
process VALIDATE_STEP2 {
    tag "validate"
    
    input:
    val plasma_data
    path ref_genome
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
    echo "TSS BED: ${tss_bed}" >> step2_validation.txt
    echo "GC content WIG: ${gc_wig}" >> step2_validation.txt
    echo "Mappability WIG: ${map_wig}" >> step2_validation.txt
    echo "Step 2 inputs validated successfully!" >> step2_validation.txt
    """
}

workflow {
    // Run the validation
    VALIDATE_STEP2(
        test_plasma.collect(),
        ref_genome,
        tss_bed,
        gc_wig,
        map_wig
    )
}
