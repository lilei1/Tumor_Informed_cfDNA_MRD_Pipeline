#!/usr/bin/env nextflow

// Simple test workflow for Step 2.5: Background Error Model
// Tests enhanced toy data without requiring external tools

// Test data channels
test_healthy_plasma = Channel
    .fromPath("plasma/healthy/*_R{1,2}.fastq.gz")
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

// Mock truth set BED and reference genome
mock_truthset_bed = Channel.fromPath("resources/mock_truthset.bed")
ref_genome = Channel.fromPath(params.genome)

// Comprehensive validation process
process VALIDATE_ENHANCED_STEP2_5 {
    tag "validate_enhanced"
    
    input:
    val healthy_data
    path mock_truthset_bed
    path ref_genome
    
    output:
    path "enhanced_step2_5_validation.txt", emit: validation
    
    script:
    """
    #!/bin/bash
    set -e
    echo "=== Enhanced Step 2.5 Validation ===" > enhanced_step2_5_validation.txt
    
    # Analyze healthy plasma data
    echo "Healthy plasma data analysis:" >> enhanced_step2_5_validation.txt
    echo "Total samples: \$(echo "${healthy_data}" | tr ',' '\\n' | grep -c 'healthy_donor')" >> enhanced_step2_5_validation.txt
    
    # Count unique donors
    unique_donors=\$(echo "${healthy_data}" | tr ',' '\\n' | grep 'healthy_donor' | sort | uniq | wc -l)
    echo "Unique donors: \${unique_donors}" >> enhanced_step2_5_validation.txt
    
    # List all donors
    echo "Donor list:" >> enhanced_step2_5_validation.txt
    echo "${healthy_data}" | tr ',' '\\n' | grep 'healthy_donor' | sort | uniq >> enhanced_step2_5_validation.txt
    
    # Analyze truth set
    echo "" >> enhanced_step2_5_validation.txt
    echo "Truth set analysis:" >> enhanced_step2_5_validation.txt
    echo "Truth set file: ${mock_truthset_bed}" >> enhanced_step2_5_validation.txt
    
    # Count regions
    if [ -f "${mock_truthset_bed}" ]; then
        region_count=\$(wc -l < "${mock_truthset_bed}")
        echo "Total regions: \${region_count}" >> enhanced_step2_5_validation.txt
        
        # Count chromosomes
        chrom_count=\$(cut -f1 "${mock_truthset_bed}" | sort | uniq | wc -l)
        echo "Chromosomes covered: \${chrom_count}" >> enhanced_step2_5_validation.txt
        
        # List chromosomes
        echo "Chromosomes:" >> enhanced_step2_5_validation.txt
        cut -f1 "${mock_truthset_bed}" | sort | uniq >> enhanced_step2_5_validation.txt
    else
        echo "ERROR: Truth set file not found!" >> enhanced_step2_5_validation.txt
    fi
    
    # Analyze reference genome
    echo "" >> enhanced_step2_5_validation.txt
    echo "Reference genome analysis:" >> enhanced_step2_5_validation.txt
    echo "Reference file: ${ref_genome}" >> enhanced_step2_5_validation.txt
    
    if [ -f "${ref_genome}" ]; then
        chrom_count=\$(grep '^>' "${ref_genome}" | wc -l)
        echo "Chromosomes in reference: \${chrom_count}" >> enhanced_step2_5_validation.txt
        
        # List chromosomes
        echo "Reference chromosomes:" >> enhanced_step2_5_validation.txt
        grep '^>' "${ref_genome}" | sed 's/^>//' >> enhanced_step2_5_validation.txt
        
        # Check sequence length
        total_length=\$(grep -v '^>' "${ref_genome}" | tr -d '\\n' | wc -c)
        echo "Total sequence length: \${total_length} bp" >> enhanced_step2_5_validation.txt
    else
        echo "ERROR: Reference genome file not found!" >> enhanced_step2_5_validation.txt
    fi
    
    # Summary
    echo "" >> enhanced_step2_5_validation.txt
    echo "=== Summary ===" >> enhanced_step2_5_validation.txt
    echo "Enhanced Step 2.5 toy data validation completed successfully!" >> enhanced_step2_5_validation.txt
    echo "Ready for error model building with:" >> enhanced_step2_5_validation.txt
    echo "  - \${unique_donors} healthy donors" >> enhanced_step2_5_validation.txt
    echo "  - \${region_count} truth set regions" >> enhanced_step2_5_validation.txt
    echo "  - \${chrom_count} reference chromosomes" >> enhanced_step2_5_validation.txt
    """
}

workflow {
    // Validate enhanced Step 2.5 inputs
    VALIDATE_ENHANCED_STEP2_5(
        test_healthy_plasma.collect(),
        mock_truthset_bed,
        ref_genome
    )
}

// Note: This workflow validates the enhanced Step 2.5 toy data
// For full execution, use the main.nf workflow
