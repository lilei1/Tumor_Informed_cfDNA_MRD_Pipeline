#!/usr/bin/env nextflow

// Background Error Model Module
// Step 2.5: Learn context-specific error rates from healthy donor plasma

process BUILD_ERROR_MODEL {
    tag "error_model"
    publishDir "${params.outdir}/error_model", mode: 'copy'
    cpus 16
    memory '64 GB'
    time '12h'
    
    input:
    path healthy_bams
    path truthset_bed
    path ref_genome
    
    output:
    path "error_model_tricontext.json", emit: error_model_json
    path "error_model_tricontext.tsv", emit: error_model_tsv
    path "error_model_summary.txt", emit: summary
    path "*.log", emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    mkdir -p error_model
    echo "Building background error model from healthy donor plasma"
    
    # Create BAM list file
    echo "Creating BAM list for error model..."
    for bam in *.bam; do
        if [[ -f "\$bam" ]]; then
            echo "\$bam" >> error_model/healthy_bams.txt
        fi
    done
    
    # Build error model using Python script
    python3 << 'PY'
import json
import pandas as pd
import pysam
from collections import defaultdict
import sys

def get_trinucleotide_context(ref_genome, chrom, pos):
    """Get trinucleotide context at given position"""
    try:
        # Get 3bp context centered on position
        start = max(0, pos - 1)
        end = pos + 2
        context = ref_genome.fetch(chrom, start, end)
        return context.upper()
    except:
        return "NNN"

def build_error_model(healthy_bams, truthset_bed, ref_genome_path):
    """Build context-specific error model from healthy plasma"""
    
    # Initialize error rate storage
    error_rates = defaultdict(lambda: defaultdict(int))
    total_counts = defaultdict(lambda: defaultdict(int))
    
    # Load reference genome
    ref_genome = pysam.FastaFile(ref_genome_path)
    
    # Load truth set regions
    truth_regions = []
    with open(truthset_bed, 'r') as f:
        for line in f:
            if line.strip():
                fields = line.strip().split()
                if len(fields) >= 3:
                    truth_regions.append((fields[0], int(fields[1]), int(fields[2])))
    
    print(f"Processing {len(truth_regions)} truth set regions")
    
    # Process each healthy BAM
    for bam_file in healthy_bams:
        print(f"Processing {bam_file}")
        bam = pysam.AlignmentFile(bam_file, 'rb')
        
        for chrom, start, end in truth_regions:
            try:
                # Get pileup at this region
                for pileup_column in bam.pileup(chrom, start, end, truncate=True):
                    pos = pileup_column.pos
                    
                    # Get trinucleotide context
                    context = get_trinucleotide_context(ref_genome, chrom, pos)
                    if context == "NNN":
                        continue
                    
                    # Count reads and errors
                    for read in pileup_column.pileups:
                        if read.is_del or read.is_refskip:
                            continue
                        
                        # Get reference and alternate bases
                        ref_base = ref_genome.fetch(chrom, pos, pos + 1).upper()
                        alt_base = read.alignment.query_sequence[read.query_position]
                        
                        # Skip if same as reference
                        if alt_base == ref_base:
                            continue
                        
                        # Get strand and position info
                        strand = "+" if not read.alignment.is_reverse else "-"
                        read_pos = read.query_position
                        read_len = read.alignment.query_length
                        
                        # Normalize read position to 0-1 scale
                        norm_pos = read_pos / read_len if read_len > 0 else 0.5
                        
                        # Create context key
                        context_key = f"{context}_{strand}_{norm_pos:.2f}"
                        
                        # Count errors
                        error_rates[context_key][f"{ref_base}->{alt_base}"] += 1
                        total_counts[context_key]["total"] += 1
                        
            except Exception as e:
                print(f"Error processing region {chrom}:{start}-{end}: {e}")
                continue
        
        bam.close()
    
    ref_genome.close()
    
    # Calculate error rates
    results = []
    for context_key, errors in error_rates.items():
        total = total_counts[context_key]["total"]
        if total > 0:
            for mutation, count in errors.items():
                error_rate = count / total
                context, strand, pos = context_key.rsplit("_", 2)
                
                results.append({
                    "context": context,
                    "strand": strand,
                    "read_position": float(pos),
                    "mutation": mutation,
                    "error_count": count,
                    "total_count": total,
                    "error_rate": error_rate
                })
    
    return results

def main():
    # Get input files
    import glob
    healthy_bams = glob.glob("*.bam")
    truthset_bed = "truthset.bed"
    ref_genome = "ref.fa"
    
    print(f"Found {len(healthy_bams)} healthy BAM files")
    
    # Build error model
    error_model = build_error_model(healthy_bams, truthset_bed, ref_genome)
    
    # Save as JSON
    with open("error_model_tricontext.json", "w") as f:
        json.dump(error_model, f, indent=2)
    
    # Save as TSV
    df = pd.DataFrame(error_model)
    df.to_csv("error_model_tricontext.tsv", sep="\\t", index=False)
    
    # Create summary
    with open("error_model_summary.txt", "w") as f:
        f.write("=== Background Error Model Summary ===\\n")
        f.write(f"Total contexts analyzed: {len(set(item['context'] for item in error_model))}\\n")
        f.write(f"Total error rate entries: {len(error_model)}\\n")
        
        # Summary by context
        context_summary = defaultdict(list)
        for item in error_model:
            context_summary[item['context']].append(item['error_rate'])
        
        f.write("\\nError rates by context:\\n")
        for context, rates in context_summary.items():
            avg_rate = sum(rates) / len(rates)
            f.write(f"{context}: {avg_rate:.6f} (n={len(rates)})\\n")
    
    print("Error model built successfully!")

if __name__ == "__main__":
    main()
PY
    
    echo "Background error model completed"
    """
}

process VALIDATE_ERROR_MODEL {
    tag "validate"
    publishDir "${params.outdir}/error_model", mode: 'copy'
    cpus 4
    memory '16 GB'
    time '2h'
    
    input:
    path error_model_json
    path error_model_tsv
    
    output:
    path "error_model_validation.txt", emit: validation
    
    script:
    """
    #!/bin/bash
    set -e
    echo "Validating background error model"
    
    # Validate JSON structure
    python3 -c "
import json
import pandas as pd

# Load error model
with open('${error_model_json}', 'r') as f:
    model = json.load(f)

# Check structure
required_fields = ['context', 'strand', 'read_position', 'mutation', 'error_rate']
for item in model:
    for field in required_fields:
        if field not in item:
            print(f'Missing field: {field}')
            exit(1)

# Load TSV
df = pd.read_csv('${error_model_tsv}', sep='\\t')

# Validate data
print(f'JSON entries: {len(model)}')
print(f'TSV entries: {len(df)}')
print(f'Contexts: {df.context.nunique()}')
print(f'Average error rate: {df.error_rate.mean():.6f}')
print(f'Min error rate: {df.error_rate.min():.6f}')
print(f'Max error rate: {df.error_rate.max():.6f}')

print('Error model validation passed!')
"
    
    # Create validation report
    echo "=== Error Model Validation ===" > error_model_validation.txt
    echo "JSON file: ${error_model_json}" >> error_model_validation.txt
    echo "TSV file: ${error_model_tsv}" >> error_model_validation.txt
    echo "Validation completed successfully!" >> error_model_validation.txt
    """
}
