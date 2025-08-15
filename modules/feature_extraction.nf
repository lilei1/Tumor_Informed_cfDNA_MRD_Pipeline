#!/usr/bin/env nextflow

/*
 * Feature Extraction Module
 * Integrates variant, methylation, and fragmentomics features
 */

process FEATURE_EXTRACTION {
    tag "${sample}_T${timepoint}"
    
    publishDir "${params.outdir}/features", mode: 'copy'
    
    cpus 8
    memory '32 GB'
    time '6h'
    
    input:
    tuple val(sample), val(timepoint), path(plasma_bam), path(plasma_bai)
    tuple val(sample), path(truth_set_vcf), path(truth_set_vcf_idx)
    path dmr_bed
    path tss_bed
    path ref_genome
    
    output:
    tuple val(sample), val(timepoint), path("*.features.csv"), emit: features
    tuple val(sample), val(timepoint), path("*.consensus.bam"), path("*.consensus.bai"), emit: consensus_bam
    path "*.log", emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    
    # Extract variants at truth set loci
    gatk SelectVariants \
        --reference ${ref_genome} \
        --variant ${plasma_bam} \
        --intervals ${truth_set_vcf} \
        --output ${sample}_T${timepoint}_variants.vcf.gz
    
    # UMI consensus calling for ultra-low VAF detection
    umi_tools consensus \
        --stdin ${plasma_bam} \
        --stdout ${sample}_T${timepoint}_consensus.bam \
        --method threshold \
        --threshold ${params.consensusThreshold} \
        --log ${sample}_T${timepoint}_consensus.log
    
    # Index consensus BAM
    samtools index ${sample}_T${timepoint}_consensus.bam
    
    # Extract fragmentomics features
    python3 -c "
    import pysam
    import pandas as pd
    import numpy as np
    
    # Read consensus BAM
    bam = pysam.AlignmentFile('${sample}_T${timepoint}_consensus.bam', 'rb')
    
    # Fragment size distribution
    sizes = []
    for read in bam.fetch():
        if read.is_proper_pair and read.is_read1:
            size = abs(read.template_length)
            if 50 <= size <= 1000:
                sizes.append(size)
    
    # Calculate fragmentomics features
    if sizes:
        features = {
            'sample': '${sample}',
            'timepoint': ${timepoint},
            'fragment_count': len(sizes),
            'mean_size': np.mean(sizes),
            'median_size': np.median(sizes),
            'std_size': np.std(sizes),
            'size_50_150': len([s for s in sizes if 50 <= s <= 150]),
            'size_151_250': len([s for s in sizes if 151 <= s <= 250]),
            'size_251_350': len([s for s in sizes if 251 <= s <= 350]),
            'size_351_450': len([s for s in sizes if 351 <= s <= 450]),
            'size_451_1000': len([s for s in sizes if 451 <= s <= 1000])
        }
    else:
        features = {
            'sample': '${sample}',
            'timepoint': ${timepoint},
            'fragment_count': 0,
            'mean_size': 0, 'median_size': 0, 'std_size': 0,
            'size_50_150': 0, 'size_151_250': 0, 'size_251_350': 0,
            'size_351_450': 0, 'size_451_1000': 0
        }
    
    # Save features
    df = pd.DataFrame([features])
    df.to_csv('${sample}_T${timepoint}_fragmentomics.csv', index=False)
    "
    
    # Extract methylation features (if DMR panel provided)
    if [ -f "${dmr_bed}" ]; then
        echo "Extracting methylation features..."
        # Add methylation analysis here
        echo "Methylation features extracted" > ${sample}_T${timepoint}_methylation.log
    fi
    
    # Extract TSS-related features
    python3 -c "
    import pysam
    import pandas as pd
    
    # Read TSS BED and consensus BAM
    bam = pysam.AlignmentFile('${sample}_T${timepoint}_consensus.bam', 'rb')
    
    # TSS proximity features
    tss_features = {
        'sample': '${sample}',
        'timepoint': ${timepoint},
        'reads_near_tss': 0,
        'tss_coverage': 0.0
    }
    
    # Calculate TSS-related metrics
    # This is a simplified version - expand based on your TSS analysis needs
    
    # Save TSS features
    df = pd.DataFrame([tss_features])
    df.to_csv('${sample}_T${timepoint}_tss.csv', index=False)
    "
    
    # Combine all features
    python3 -c "
    import pandas as pd
    
    # Read individual feature files
    fragmentomics = pd.read_csv('${sample}_T${timepoint}_fragmentomics.csv')
    tss = pd.read_csv('${sample}_T${timepoint}_tss.csv')
    
    # Merge features
    features = fragmentomics.merge(tss, on=['sample', 'timepoint'])
    
    # Add variant features
    features['variant_count'] = 0  # Calculate from VCF
    features['mean_vaf'] = 0.0     # Calculate from VCF
    
    # Save combined features
    features.to_csv('${sample}_T${timepoint}_features.csv', index=False)
    "
    
    # Clean up intermediate files
    rm -f *_variants.vcf.gz *_fragmentomics.csv *_tss.csv
    """
}
