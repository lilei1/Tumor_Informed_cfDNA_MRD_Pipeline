#!/usr/bin/env nextflow

/*
 * Test Workflow for Step 6: One-Click Reporting Dashboard
 * Tests the comprehensive reporting system with interactive visualizations
 */

// Include reporting dashboard module (simplified version)
include { GENERATE_PATIENT_DASHBOARD; GENERATE_BATCH_REPORTS; GENERATE_CLINICAL_RMARKDOWN } from '../modules/reporting_dashboard_simple'

// Create mock data for testing
process CREATE_MOCK_REPORTING_DATA {
    tag "mock_reporting_data"
    
    output:
    path "mock_patient_id.txt", emit: patient_id
    path "mock_mrd_scores.tsv", emit: mrd_scores
    path "mock_truth_set_loci.tsv", emit: truth_set_loci
    path "mock_fragmentomics.tsv", emit: fragmentomics
    path "mock_dmr_scores.tsv", emit: dmr_scores
    path "mock_cnv_results.tsv", emit: cnv_results
    path "mock_qc_results.tsv", emit: qc_results
    path "mock_clinical_data.tsv", emit: clinical_data
    path "mock_patient_list.txt", emit: patient_list
    path "mock_all_mrd_scores.tsv", emit: all_mrd_scores
    path "mock_all_truth_set.tsv", emit: all_truth_set
    path "mock_all_fragmentomics.tsv", emit: all_fragmentomics
    path "mock_all_cnv_results.tsv", emit: all_cnv_results
    path "mock_all_qc_results.tsv", emit: all_qc_results
    path "mock_clinical_metadata.tsv", emit: clinical_metadata
    path "mock_rmarkdown_template.Rmd", emit: rmarkdown_template
    
    script:
    """
    #!/bin/bash
    set -e
    echo "Creating mock reporting data for Step 6 testing..."
    
    # Create mock patient ID
    echo "patient_01" > mock_patient_id.txt
    
    # Create mock MRD scores
    cat > mock_mrd_scores.tsv << 'EOF'
    timepoint	mrd_probability
    T0	0.05
    T1	0.15
    T2	0.35
    T3	0.65
    EOF
    
    # Create mock truth set loci
    cat > mock_truth_set_loci.tsv << 'EOF'
    locus	depth	umi_families	vaf	p_value	status
    chr1:1000	150	5	0.02	1e-6	PASS
    chr1:2000	200	8	0.03	1e-8	PASS
    chr1:3000	180	6	0.025	1e-7	PASS
    EOF
    
    # Create mock fragmentomics
    cat > mock_fragmentomics.tsv << 'EOF'
    fragment_length	gc_content	coverage	end_motif	motif_count	tss_position	tss_coverage
    166	0.5	50	AA	25	-1000	40
    180	0.6	60	AT	30	-800	45
    150	0.4	40	AG	20	-600	35
    EOF
    
    # Create mock DMR scores
    cat > mock_dmr_scores.tsv << 'EOF'
    region	methylation_score	status
    chr1:1000-2000	0.8	HYPERMETHYLATED
    chr1:3000-4000	0.2	HYPOMETHYLATED
    EOF
    
    # Create mock CNV results
    cat > mock_cnv_results.tsv << 'EOF'
    chromosome	copy_number	timepoint	tumor_fraction
    chr1	2.1	T0	0.05
    chr2	1.9	T1	0.06
    chr3	2.0	T2	0.07
    EOF
    
    # Create mock QC results
    cat > mock_qc_results.tsv << 'EOF'
    metric	status	value
    WES_QC	PASS	85.0
    Plasma_QC	PASS	90.0
    Per_Run_QC	PASS	95.0
    EOF
    
    # Create mock clinical data
    cat > mock_clinical_data.tsv << 'EOF'
    patient_id	diagnosis_date	stage	treatment
    patient_01	2024-01-15	IIIA	FOLFOX + Surgery
    EOF
    
    # Create mock patient list
    cat > mock_patient_list.txt << 'EOF'
    patient_01
    patient_02
    patient_03
    EOF
    
    # Create mock all MRD scores
    cat > mock_all_mrd_scores.tsv << 'EOF'
    patient_id	timepoint	mrd_probability
    patient_01	T0	0.05
    patient_01	T1	0.15
    patient_02	T0	0.10
    patient_02	T1	0.25
    patient_03	T0	0.08
    patient_03	T1	0.20
    EOF
    
    # Create mock all truth set
    cat > mock_all_truth_set.tsv << 'EOF'
    patient_id	locus	depth	umi_families	vaf	p_value	status
    patient_01	chr1:1000	150	5	0.02	1e-6	PASS
    patient_02	chr1:2000	200	8	0.03	1e-8	PASS
    patient_03	chr1:3000	180	6	0.025	1e-7	PASS
    EOF
    
    # Create mock all fragmentomics
    cat > mock_all_fragmentomics.tsv << 'EOF'
    patient_id	fragment_length	gc_content	coverage	end_motif	motif_count
    patient_01	166	0.5	50	AA	25
    patient_02	180	0.6	60	AT	30
    patient_03	150	0.4	40	AG	20
    EOF
    
    # Create mock all CNV results
    cat > mock_all_cnv_results.tsv << 'EOF'
    patient_id	chromosome	copy_number	timepoint	tumor_fraction
    patient_01	chr1	2.1	T0	0.05
    patient_02	chr2	1.9	T1	0.06
    patient_03	chr3	2.0	T2	0.07
    EOF
    
    # Create mock all QC results
    cat > mock_all_qc_results.tsv << 'EOF'
    patient_id	metric	status	value
    patient_01	WES_QC	PASS	85.0
    patient_02	Plasma_QC	PASS	90.0
    patient_03	Per_Run_QC	PASS	95.0
    EOF
    
    # Create mock clinical metadata
    cat > mock_clinical_metadata.tsv << 'EOF'
    patient_id	diagnosis_date	stage	treatment
    patient_01	2024-01-15	IIIA	FOLFOX + Surgery
    patient_02	2024-02-01	IIB	FOLFOX + Surgery
    patient_03	2024-02-15	IIIA	FOLFOX + Surgery
    EOF
    
    # Create mock R Markdown template
    cat > mock_rmarkdown_template.Rmd << 'EOF'
    ---
    title: "Clinical MRD Report"
    author: "Tumor-Informed cfDNA MRD Pipeline"
    date: "`r Sys.Date()`"
    output:
      html_document:
        toc: true
        toc_float: true
        theme: cosmo
        highlight: tango
    ---
    
    # Clinical MRD Assessment
    
    This is a test R Markdown template for clinical reporting.
    EOF
    
    echo "Mock reporting data creation completed"
    """
}

workflow {
    // Create mock reporting data
    CREATE_MOCK_REPORTING_DATA()
    
    // Test individual patient dashboard generation
    GENERATE_PATIENT_DASHBOARD(
        CREATE_MOCK_REPORTING_DATA.out.patient_id,
        CREATE_MOCK_REPORTING_DATA.out.mrd_scores,
        CREATE_MOCK_REPORTING_DATA.out.truth_set_loci,
        CREATE_MOCK_REPORTING_DATA.out.fragmentomics,
        CREATE_MOCK_REPORTING_DATA.out.dmr_scores,
        CREATE_MOCK_REPORTING_DATA.out.cnv_results,
        CREATE_MOCK_REPORTING_DATA.out.qc_results,
        CREATE_MOCK_REPORTING_DATA.out.clinical_data
    )
    
    // Test batch report generation
    GENERATE_BATCH_REPORTS(
        CREATE_MOCK_REPORTING_DATA.out.patient_list,
        CREATE_MOCK_REPORTING_DATA.out.all_mrd_scores,
        CREATE_MOCK_REPORTING_DATA.out.all_truth_set,
        CREATE_MOCK_REPORTING_DATA.out.all_fragmentomics,
        CREATE_MOCK_REPORTING_DATA.out.all_cnv_results,
        CREATE_MOCK_REPORTING_DATA.out.all_qc_results,
        CREATE_MOCK_REPORTING_DATA.out.clinical_metadata
    )
    
    // Test clinical R Markdown report generation
    GENERATE_CLINICAL_RMARKDOWN(
        CREATE_MOCK_REPORTING_DATA.out.all_mrd_scores,
        CREATE_MOCK_REPORTING_DATA.out.rmarkdown_template
    )
}
