#!/bin/bash

# Test script for Tumor-Informed cfDNA MRD Pipeline
# Tests Step 1: Tumor-Normal WES → Somatic Truth Set

set -e

echo "=== Testing Tumor-Informed cfDNA MRD Pipeline ==="
echo "Step 1: Tumor-Normal WES → Somatic Truth Set"
echo ""

# Check if Nextflow is installed
if ! command -v nextflow &> /dev/null; then
    echo "Error: Nextflow is not installed. Please install Nextflow first."
    echo "Install with: curl -s https://get.nextflow.io | bash"
    exit 1
fi

# Check if Docker is running
if ! docker info &> /dev/null; then
    echo "Error: Docker is not running. Please start Docker first."
    exit 1
fi

echo "✓ Nextflow version: $(nextflow -version)"
echo "✓ Docker is running"
echo ""

# Create results directory
mkdir -p results

echo "=== Running Test Pipeline ==="
echo "This will test the complete Step 1 workflow with minimal test data."
echo ""

# Run the test pipeline
echo "Executing: nextflow run test_pipeline.nf -profile docker"
echo ""

nextflow run test_pipeline.nf -profile docker

echo ""
echo "=== Test Pipeline Completed ==="
echo "Check the 'results' directory for outputs:"
echo "  - wes/: Preprocessed BAM files and QC metrics"
echo "  - variants/: Somatic variants and truth set"
echo "  - pipeline_report.html: Pipeline execution report"
echo ""

# List results
if [ -d "results" ]; then
    echo "Generated files:"
    find results -type f -name "*.bam" -o -name "*.vcf.gz" -o -name "*.bed" -o -name "*.html" | sort
fi

echo ""
echo "Test completed successfully!"
