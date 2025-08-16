#!/usr/bin/env nextflow

// Plasma cfDNA UMI-based WGS processing module
// Implements Step 2: UMI extraction → consensus → variant calling → features

process UMI_EXTRACTION {
    tag "${sample}_${timepoint}"
    publishDir "${params.outdir}/plasma/umi", mode: 'copy'
    cpus 8
    memory '32 GB'
    time '4h'
    
    input:
    tuple val(sample), val(timepoint), path(fastq1), path(fastq2)
    
    output:
    tuple val(sample), val(timepoint), path("*.umi.fq.gz"), path("*.umi.fq.gz"), emit: umi_fastq
    tuple val(sample), val(timepoint), path("*.log"), emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    mkdir -p plasma/qc
    echo "Extracting UMIs for ${sample} at ${timepoint}"
    
    # UMI extraction using UMI-tools
    # Adjust bc-pattern based on your UMI kit design
    umi_tools extract \\
        --bc-pattern=NNNNNNNN \\
        -I ${fastq1} --read2-in=${fastq2} \\
        -S plasma/${timepoint}_R1.umi.fq.gz --read2-out=plasma/${timepoint}_R2.umi.fq.gz \\
        --log=plasma/qc/${timepoint}.umi_extract.log
    
    echo "UMI extraction completed for ${sample}_${timepoint}"
    """
}

process PLASMA_ALIGNMENT {
    tag "${sample}_${timepoint}"
    publishDir "${params.outdir}/plasma/alignment", mode: 'copy'
    cpus 16
    memory '64 GB'
    time '6h'
    
    input:
    tuple val(sample), val(timepoint), path(umi_fastq1), path(umi_fastq2)
    path ref_genome
    
    output:
    tuple val(sample), val(timepoint), path("*.raw.bam"), path("*.raw.bai"), emit: raw_bam
    tuple val(sample), val(timepoint), path("*.log"), emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    echo "Aligning plasma reads for ${sample} at ${timepoint}"
    
    # Align raw UMI-tagged reads
    bwa-mem2 mem -t 16 ${ref_genome} \\
        ${umi_fastq1} ${umi_fastq2} \\
        | samtools sort -@8 -o plasma/${timepoint}.raw.bam
    
    samtools index plasma/${timepoint}.raw.bam
    
    echo "Alignment completed for ${sample}_${timepoint}"
    """
}

process UMI_CONSENSUS {
    tag "${sample}_${timepoint}"
    publishDir "${params.outdir}/plasma/consensus", mode: 'copy'
    cpus 16
    memory '64 GB'
    time '8h'
    
    input:
    tuple val(sample), val(timepoint), path(raw_bam), path(raw_bai)
    
    output:
    tuple val(sample), val(timepoint), path("*.consensus.bam"), path("*.consensus.bai"), emit: consensus_bam
    tuple val(sample), val(timepoint), path("*.grouped.bam"), emit: grouped_bam
    tuple val(sample), val(timepoint), path("*.log"), emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    echo "Generating UMI consensus for ${sample} at ${timepoint}"
    
    # Group reads by UMI (adjacency method)
    fgbio GroupReadsByUmi \\
        -i ${raw_bam} -o plasma/${timepoint}.grouped.bam \\
        -s adjacency -e 1 -f PRIMER_SEQUENCE_TAG=RX
    
    # Call consensus per UMI family
    fgbio CallMolecularConsensusReads \\
        -i plasma/${timepoint}.grouped.bam -o plasma/${timepoint}.consensus.unmapped.bam \\
        --min-reads 2 --error-rate-post-umi 45 --min-input-base-quality 20
    
    # Filter consensus reads
    fgbio FilterConsensusReads \\
        -i plasma/${timepoint}.consensus.unmapped.bam -o plasma/${timepoint}.consensus.filtered.bam \\
        --min-base-quality 30 --max-no-call-fraction 0.2
    
    # Re-align consensus reads
    samtools sort -@8 -o plasma/${timepoint}.consensus.bam plasma/${timepoint}.consensus.filtered.bam
    samtools index plasma/${timepoint}.consensus.bam
    
    echo "UMI consensus completed for ${sample}_${timepoint}"
    """
}

process PLASMA_QC {
    tag "${sample}_${timepoint}"
    publishDir "${params.outdir}/plasma/qc", mode: 'copy'
    cpus 8
    memory '32 GB'
    time '4h'
    
    input:
    tuple val(sample), val(timepoint), path(consensus_bam), path(consensus_bai)
    path ref_genome
    
    output:
    tuple val(sample), val(timepoint), path("*.insert.txt"), path("*.insert.pdf"), emit: insert_metrics
    tuple val(sample), val(timepoint), path("*.per-base.wig"), emit: coverage_wig
    tuple val(sample), val(timepoint), path("*.log"), emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    mkdir -p plasma/qc
    echo "Running QC for ${sample} at ${timepoint}"
    
    # Insert size metrics
    gatk CollectInsertSizeMetrics \\
        I=${consensus_bam} O=plasma/qc/${timepoint}.insert.txt H=plasma/qc/${timepoint}.insert.pdf \\
        R=${ref_genome} --java-options "-Xmx16g"
    
    # Coverage analysis with mosdepth
    mosdepth -t 8 plasma/qc/${timepoint} ${consensus_bam}
    
    echo "QC completed for ${sample}_${timepoint}"
    """
}

process TRUTHSET_VARIANT_CALLING {
    tag "${sample}_${timepoint}"
    publishDir "${params.outdir}/plasma/variants", mode: 'copy'
    cpus 16
    memory '64 GB'
    time '6h'
    
    input:
    tuple val(sample), val(timepoint), path(consensus_bam), path(consensus_bai)
    path ref_genome
    path truthset_bed
    path gnomad_vcf
    path pon_vcf
    
    output:
    tuple val(sample), val(timepoint), path("*.truthset.unfiltered.vcf.gz"), emit: unfiltered_variants
    tuple val(sample), val(timepoint), path("*.truthset.filtered.vcf.gz"), emit: filtered_variants
    tuple val(sample), val(timepoint), path("*.variant_counts.tsv"), emit: variant_counts
    tuple val(sample), val(timepoint), path("*.log"), emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    echo "Calling variants at truth set loci for ${sample} at ${timepoint}"
    
    # Option A: Mutect2 tumor-only restricted to loci
    gatk Mutect2 \\
        -R ${ref_genome} \\
        -I ${consensus_bam} \\
        --germline-resource ${gnomad_vcf} \\
        --panel-of-normals ${pon_vcf} \\
        -L ${truthset_bed} \\
        -O plasma/${timepoint}.truthset.unfiltered.vcf.gz \\
        --java-options "-Xmx32g"
    
    gatk FilterMutectCalls \\
        -R ${ref_genome} \\
        -V plasma/${timepoint}.truthset.unfiltered.vcf.gz \\
        -O plasma/${timepoint}.truthset.filtered.vcf.gz \\
        --java-options "-Xmx32g"
    
    # Option B: Exact-count pileup (UMI-aware) with bcftools
    bcftools mpileup -f ${ref_genome} -R ${truthset_bed} \\
        -a FORMAT/AD,FORMAT/DP ${consensus_bam} \\
        | bcftools call -mv -Ob -o plasma/${timepoint}.truthset.bcf
    
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%DP\\t%AD]\\n' plasma/${timepoint}.truthset.bcf \\
        > features/${sample}_${timepoint}.variant_counts.tsv
    
    echo "Variant calling completed for ${sample}_${timepoint}"
    """
}

process FRAGMENTOMICS_FEATURES {
    tag "${sample}_${timepoint}"
    publishDir "${params.outdir}/plasma/features", mode: 'copy'
    cpus 8
    memory '32 GB'
    time '4h'
    
    input:
    tuple val(sample), val(timepoint), path(consensus_bam), path(consensus_bai)
    path ref_genome
    path tss_bed
    
    output:
    tuple val(sample), val(timepoint), path("*.fragmentomics.tsv"), emit: fragmentomics
    tuple val(sample), val(timepoint), path("*.endmotifs.txt"), emit: endmotifs
    tuple val(sample), val(timepoint), path("*.tss_cov.tsv"), emit: tss_coverage
    tuple val(sample), val(timepoint), path("*.log"), emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    mkdir -p features
    echo "Extracting fragmentomics features for ${sample} at ${timepoint}"
    
    # Length distribution (paired-end)
    samtools view -f 3 ${consensus_bam} \\
        | python3 - << 'PY'
import sys, numpy as np
L=[int(x.strip().split()[8]) for x in sys.stdin if len(x.strip().split()) > 8]
short=sum(1 for x in L if x<150); long=sum(1 for x in L if 150<=x<=220)
print("short/long_ratio\t{:.3f}".format(short/max(1,long)))
PY
        > features/${sample}_${timepoint}.fragmentomics.tsv
    
    # End-motifs (5' 4-mers)
    bedtools bamtobed -i ${consensus_bam} \\
        | python3 - << 'PY'
import sys
for line in sys.stdin:
    fields = line.strip().split()
    if len(fields) >= 6:
        if fields[5] == "+":
            print(f"{fields[0]}\t{fields[1]}\t{int(fields[1])+1}\t+")
        else:
            print(f"{fields[0]}\t{int(fields[2])-1}\t{fields[2]}\t-")
PY
        | bedtools getfasta -fi ${ref_genome} -bed - -s \\
        | paste - - | cut -f2 \\
        | python3 - << 'PY'
import sys
for line in sys.stdin:
    seq = line.strip()
    if len(seq) >= 4:
        print(seq[:4].upper())
PY
        | sort | uniq -c | sort -nr \\
        > features/${sample}_${timepoint}.endmotifs.txt
    
    # TSS enrichment
    bedtools coverage -a ${tss_bed} -b ${consensus_bam} \\
        | python3 - << 'PY'
import sys
cov = {}
for line in sys.stdin:
    fields = line.strip().split()
    if len(fields) >= 7:
        key = fields[3]
        value = float(fields[6])
        cov[key] = cov.get(key, 0) + value
for k, v in cov.items():
    print(f"{k}\t{v}")
PY
        > features/${sample}_${timepoint}.tss_cov.tsv
    
    echo "Fragmentomics features extracted for ${sample}_${timepoint}"
    """
}

process METHYLATION_FEATURES {
    tag "${sample}_${timepoint}"
    publishDir "${params.outdir}/plasma/methylation", mode: 'copy'
    cpus 16
    memory '64 GB'
    time '8h'
    
    input:
    tuple val(sample), val(timepoint), path(em_fastq1), path(em_fastq2)
    path ref_genome
    path dmr_bed
    
    output:
    tuple val(sample), val(timepoint), path("*.dmr_scores.tsv"), emit: dmr_scores
    tuple val(sample), val(timepoint), path("*.bedGraph"), emit: methylation_bedgraph
    tuple val(sample), val(timepoint), path("*.log"), emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    mkdir -p plasma/bs features
    echo "Processing methylation for ${sample} at ${timepoint}"
    
    # Align with Bismark
    bismark -1 ${em_fastq1} -2 ${em_fastq2} -o plasma/bs/ ${ref_genome}
    
    # Deduplicate
    deduplicate_bismark --bam plasma/bs/${timepoint}_EM_bismark_bt2_pe.bam
    
    # Extract methylation
    bismark_methylation_extractor --gzip --bedGraph --CX_context --cytosine_report \\
        --genome_folder ${ref_genome} plasma/bs/${timepoint}_EM_bismark_bt2_pe.deduplicated.bam
    
    # Score DMRs against cancer DMR panel
    python3 scripts/score_dmrs.py \\
        --dmr-bed ${dmr_bed} \\
        --meth-bed plasma/bs/${timepoint}_EM.bedGraph \\
        --out features/${sample}_${timepoint}.dmr_scores.tsv
    
    echo "Methylation processing completed for ${sample}_${timepoint}"
    """
}

process CNV_ANALYSIS {
    tag "${sample}_${timepoint}"
    publishDir "${params.outdir}/plasma/cnv", mode: 'copy'
    cpus 8
    memory '32 GB'
    time '6h'
    
    input:
    tuple val(sample), val(timepoint), path(coverage_wig)
    path gc_wig
    path map_wig
    
    output:
    tuple val(sample), val(timepoint), path("*.tumor_fraction.txt"), emit: tumor_fraction
    tuple val(sample), val(timepoint), path("*.cnv_segments.txt"), emit: cnv_segments
    tuple val(sample), val(timepoint), path("*.log"), emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    mkdir -p features/ichorCNA_${timepoint}
    echo "Running CNV analysis for ${sample} at ${timepoint}"
    
    # Run ichorCNA for tumor fraction and CNV segments
    Rscript ichorCNA/scripts/runIchorCNA.R \\
        --id ${sample}_${timepoint} \\
        --WIG ${coverage_wig} \\
        --gcWig ${gc_wig} \\
        --mapWig ${map_wig} \\
        --normal 0.95 --ploidy "2" \\
        --outDir features/ichorCNA_${timepoint}/
    
    # Parse tumor fraction and CNV segments
    echo "CNV analysis completed for ${sample}_${timepoint}"
    """
}
