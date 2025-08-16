#!/usr/bin/env python3
"""
Enhanced Mock Data Generator for Tumor-Informed cfDNA MRD Pipeline
Generates biologically plausible mock data for comprehensive pipeline testing
"""

import os
import random
import json
import gzip
from pathlib import Path
import numpy as np
from datetime import datetime

class EnhancedMockDataGenerator:
    def __init__(self, output_dir="test_data"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Biological parameters
        self.chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        self.bases = ["A", "C", "G", "T"]
        self.quality_scores = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
        
        # UMI patterns for different kits
        self.umi_patterns = {
            "UMI8": "NNNNNNNN",
            "UMI10": "NNNNNNNNNN",
            "UMI12": "NNNNNNNNNNNN"
        }
        
        # Fragment length distributions (cfDNA typical)
        self.fragment_lengths = {
            "short": (100, 150),    # Apoptotic fragments
            "medium": (150, 220),   # Normal cfDNA
            "long": (220, 400)      # Necrotic fragments
        }
        
        # Error rates by context (simulated)
        self.context_error_rates = {
            "CpG": 0.02,      # Higher error in CpG contexts
            "AT": 0.005,      # Lower error in AT contexts
            "GC": 0.008,      # Medium error in GC contexts
            "default": 0.01   # Default error rate
        }
    
    def generate_reference_genome(self, size_mb=100):
        """Generate a more realistic reference genome"""
        print(f"Generating {size_mb}MB reference genome...")
        
        ref_dir = self.output_dir / "refs"
        ref_dir.mkdir(exist_ok=True)
        
        with open(ref_dir / "GRCh38.fa", "w") as f:
            for chrom in self.chromosomes[:5]:  # Limit to 5 chromosomes for testing
                # Generate chromosome sequence
                seq_length = random.randint(100000, 500000)  # 100-500kb per chromosome
                sequence = "".join(random.choices(self.bases, k=seq_length))
                
                f.write(f">{chrom}\n")
                # Write in 80-character lines
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + "\n")
        
        print(f"Reference genome generated: {ref_dir / 'GRCh38.fa'}")
    
    def generate_wes_data(self, num_patients=3):
        """Generate WES data for multiple patients"""
        print(f"Generating WES data for {num_patients} patients...")
        
        wes_dir = self.output_dir / "wes"
        wes_dir.mkdir(exist_ok=True)
        
        for patient_id in range(1, num_patients + 1):
            # Tumor samples
            self._generate_fastq(
                wes_dir / "tumor" / f"patient_{patient_id}_T0_R1.fastq.gz",
                wes_dir / "tumor" / f"patient_{patient_id}_T0_R2.fastq.gz",
                num_reads=10000,
                sample_type="tumor",
                patient_id=patient_id
            )
            
            # Normal samples
            self._generate_fastq(
                wes_dir / "normal" / f"patient_{patient_id}_T0_R1.fastq.gz",
                wes_dir / "normal" / f"patient_{patient_id}_T0_R2.fastq.gz",
                num_reads=10000,
                sample_type="normal",
                patient_id=patient_id
            )
        
        print(f"WES data generated for {num_patients} patients")
    
    def generate_plasma_data(self, num_patients=3, timepoints=3):
        """Generate plasma cfDNA data with UMI tags"""
        print(f"Generating plasma data for {num_patients} patients at {timepoints} timepoints...")
        
        plasma_dir = self.output_dir / "plasma"
        plasma_dir.mkdir(exist_ok=True)
        
        for patient_id in range(1, num_patients + 1):
            for timepoint in range(timepoints):
                # Patient plasma samples
                self._generate_fastq(
                    plasma_dir / "patients" / f"patient_{patient_id}_T{timepoint}_R1.fastq.gz",
                    plasma_dir / "patients" / f"patient_{patient_id}_T{timepoint}_R2.fastq.gz",
                    num_reads=50000,  # More reads for plasma
                    sample_type="plasma",
                    patient_id=patient_id,
                    timepoint=timepoint,
                    umi_tagged=True
                )
        
        print(f"Plasma data generated for {num_patients} patients at {timepoints} timepoints")
    
    def generate_healthy_donor_data(self, num_donors=5):
        """Generate healthy donor plasma data for error model building"""
        print(f"Generating healthy donor data for {num_donors} donors...")
        
        healthy_dir = self.output_dir / "plasma" / "healthy"
        healthy_dir.mkdir(exist_ok=True)
        
        for donor_id in range(1, num_donors + 1):
            self._generate_fastq(
                healthy_dir / f"healthy_donor_{donor_id}_R1.fastq.gz",
                healthy_dir / f"healthy_donor_{donor_id}_R2.fastq.gz",
                num_reads=30000,
                sample_type="healthy",
                donor_id=donor_id,
                umi_tagged=True
            )
        
        print(f"Healthy donor data generated for {num_donors} donors")
    
    def generate_truth_set(self, num_variants=100):
        """Generate a realistic truth set based on simulated variants"""
        print(f"Generating truth set with {num_variants} variants...")
        
        resources_dir = self.output_dir / "resources"
        resources_dir.mkdir(exist_ok=True)
        
        # Generate variants across chromosomes
        variants = []
        for i in range(num_variants):
            chrom = random.choice(self.chromosomes[:5])
            pos = random.randint(1000, 100000)
            end = pos + random.randint(1, 10)
            strand = random.choice(["+", "-"])
            variant_id = f"variant_{i+1}"
            
            variants.append([chrom, pos, end, strand, variant_id])
        
        # Sort by chromosome and position
        variants.sort(key=lambda x: (x[0], x[1]))
        
        # Write BED file
        with open(resources_dir / "mock_truthset.bed", "w") as f:
            for variant in variants:
                f.write("\t".join(map(str, variant)) + "\n")
        
        # Generate VCF file
        self._generate_vcf(resources_dir / "mock_truthset.vcf.gz", variants)
        
        print(f"Truth set generated: {resources_dir / 'mock_truthset.bed'}")
    
    def generate_resource_files(self):
        """Generate additional resource files"""
        print("Generating resource files...")
        
        resources_dir = self.output_dir / "resources"
        resources_dir.mkdir(exist_ok=True)
        
        # Exome intervals
        self._generate_exome_intervals(resources_dir / "exome.interval_list")
        
        # TSS regions
        self._generate_tss_regions(resources_dir / "tss.bed")
        
        # GC content and mappability (simplified)
        self._generate_gc_content(resources_dir / "gc_hg38.wig")
        self._generate_mappability(resources_dir / "map_hg38.wig")
        
        # Mock gnomAD resource
        self._generate_gnomad_resource(resources_dir / "gnomad.af-only.vcf.gz")
        
        print("Resource files generated")
    
    def _generate_fastq(self, r1_path, r2_path, num_reads, sample_type, **kwargs):
        """Generate paired FASTQ files"""
        r1_path.parent.mkdir(exist_ok=True)
        r2_path.parent.mkdir(exist_ok=True)
        
        with gzip.open(r1_path, 'wt') as r1, gzip.open(r2_path, 'wt') as r2:
            for read_id in range(1, num_reads + 1):
                # Generate read sequence
                read_length = random.randint(100, 150)
                sequence = "".join(random.choices(self.bases, k=read_length))
                
                # Generate quality scores
                qualities = "".join(random.choices(self.quality_scores, k=read_length))
                
                # Read header
                if kwargs.get('umi_tagged'):
                    umi = "".join(random.choices(self.bases, k=8))
                    header = f"@{sample_type}_{kwargs.get('patient_id', kwargs.get('donor_id', 'unknown'))}_read_{read_id}:UMI:{umi}"
                else:
                    header = f"@{sample_type}_{kwargs.get('patient_id', 'unknown')}_read_{read_id}"
                
                # Write R1
                r1.write(f"{header}\n{sequence}\n+\n{qualities}\n")
                
                # Generate complementary R2 sequence
                complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
                r2_sequence = "".join(complement.get(base, base) for base in sequence[::-1])
                r2_qualities = qualities[::-1]
                
                r2.write(f"{header}\n{r2_sequence}\n+\n{r2_qualities}\n")
    
    def _generate_vcf(self, vcf_path, variants):
        """Generate a mock VCF file"""
        with gzip.open(vcf_path, 'wt') as f:
            # VCF header
            f.write("##fileformat=VCFv4.2\n")
            f.write("##source=MockDataGenerator\n")
            f.write("##reference=GRCh38.fa\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            
            for chrom, pos, end, strand, variant_id in variants:
                ref = random.choice(self.bases)
                alt = random.choice([b for b in self.bases if b != ref])
                qual = random.randint(30, 100)
                
                f.write(f"{chrom}\t{pos}\t{variant_id}\t{ref}\t{alt}\t{qual}\tPASS\t.\tGT:AD:DP\t0/1:10,5:15\n")
    
    def _generate_exome_intervals(self, interval_path):
        """Generate exome capture intervals"""
        with open(interval_path, 'w') as f:
            f.write("##interval_list\n")
            f.write("##source=MockDataGenerator\n")
            f.write("##reference=GRCh38.fa\n")
            f.write("##feature_type=exon\n")
            f.write("##feature_description=Mock exome capture intervals\n")
            
            for chrom in self.chromosomes[:5]:
                for exon_id in range(1, 21):  # 20 exons per chromosome
                    start = exon_id * 5000
                    end = start + random.randint(100, 500)
                    f.write(f"{chrom}\t{start}\t{end}\t+\texon_{chrom}_{exon_id}\n")
    
    def _generate_tss_regions(self, tss_path):
        """Generate TSS regions for fragmentomics"""
        with open(tss_path, 'w') as f:
            for chrom in self.chromosomes[:5]:
                for tss_id in range(1, 11):  # 10 TSS per chromosome
                    center = tss_id * 10000
                    start = center - 1000
                    end = center + 1000
                    f.write(f"{chrom}\t{start}\t{end}\tTSS_{chrom}_{tss_id}\n")
    
    def _generate_gc_content(self, gc_path):
        """Generate GC content file"""
        with open(gc_path, 'w') as f:
            for chrom in self.chromosomes[:5]:
                for pos in range(0, 100000, 1000):
                    gc_content = random.uniform(0.3, 0.7)
                    f.write(f"fixedStep chrom={chrom} start={pos} step=1000\n{gc_content:.3f}\n")
    
    def _generate_mappability(self, map_path):
        """Generate mappability file"""
        with open(map_path, 'w') as f:
            for chrom in self.chromosomes[:5]:
                for pos in range(0, 100000, 1000):
                    mappability = random.uniform(0.8, 1.0)
                    f.write(f"fixedStep chrom={chrom} start={pos} step=1000\n{mappability:.3f}\n")
    
    def _generate_gnomad_resource(self, gnomad_path):
        """Generate mock gnomAD resource"""
        with gzip.open(gnomad_path, 'wt') as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##source=MockDataGenerator\n")
            f.write("##reference=GRCh38.fa\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            for chrom in self.chromosomes[:5]:
                for pos in range(1000, 100000, 2000):
                    ref = random.choice(self.bases)
                    alt = random.choice([b for b in self.bases if b != ref])
                    af = random.uniform(0.001, 0.5)
                    
                    f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\tAF={af:.6f}\n")
    
    def generate_panel_of_normals(self):
        """Generate a mock Panel of Normals"""
        print("Generating Panel of Normals...")
        
        pon_dir = self.output_dir / "pon"
        pon_dir.mkdir(exist_ok=True)
        
        with gzip.open(pon_dir / "pon.vcf.gz", 'wt') as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##source=MockDataGenerator\n")
            f.write("##reference=GRCh38.fa\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            
            # Generate common variants found in normals
            for chrom in self.chromosomes[:3]:
                for pos in range(1000, 50000, 5000):
                    ref = random.choice(self.bases)
                    alt = random.choice([b for b in self.bases if b != ref])
                    af = random.uniform(0.01, 0.3)  # Common variants
                    
                    f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\tAF={af:.6f}\n")
        
        print(f"Panel of Normals generated: {pon_dir / 'pon.vcf.gz'}")
    
    def generate_all_data(self):
        """Generate all mock data"""
        print("=== Enhanced Mock Data Generation ===")
        print(f"Output directory: {self.output_dir}")
        
        # Generate all data types
        self.generate_reference_genome()
        self.generate_wes_data()
        self.generate_plasma_data()
        self.generate_healthy_donor_data()
        self.generate_truth_set()
        self.generate_resource_files()
        self.generate_panel_of_normals()
        
        # Generate configuration file
        self._generate_nextflow_config()
        
        print("\n=== Mock Data Generation Complete ===")
        print("All mock data has been generated and is ready for pipeline testing!")
        print(f"Total size: {self._calculate_directory_size(self.output_dir):.1f} MB")
    
    def _generate_nextflow_config(self):
        """Generate a test-specific Nextflow configuration"""
        config_path = self.output_dir / "nextflow.config"
        
        config_content = f"""// Test-specific Nextflow configuration
// Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

// Import base configuration
includeConfig '../nextflow.config'

// Override parameters for testing
params {{
    // Use test data directory
    baseDir = "{self.output_dir.absolute()}"
    
    // Test-specific parameters
    maxCpus = 4
    maxMemory = '8 GB'
    
    // Reduced resource requirements for testing
    cpus = 2
    memory = '4 GB'
    time = '1h'
}}

// Test-specific profiles
profiles {{
    test {{
        // Use local executor for testing
        executor = 'local'
        
        // Enable test mode
        params.testMode = true
    }}
    
    docker {{
        // Docker configuration for testing
        docker.enabled = true
        docker.runOptions = '-u $(id -u):$(id -g)'
    }}
}}
"""
        
        with open(config_path, 'w') as f:
            f.write(config_content)
        
        print(f"Test configuration generated: {config_path}")
    
    def _calculate_directory_size(self, directory):
        """Calculate directory size in MB"""
        total_size = 0
        for dirpath, dirnames, filenames in os.walk(directory):
            for filename in filenames:
                filepath = os.path.join(dirpath, filename)
                if os.path.exists(filepath):
                    total_size += os.path.getsize(filepath)
        return total_size / (1024 * 1024)  # Convert to MB

def main():
    """Main function to generate enhanced mock data"""
    generator = EnhancedMockDataGenerator()
    generator.generate_all_data()

if __name__ == "__main__":
    main()
