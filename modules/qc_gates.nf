#!/usr/bin/env nextflow

/*
 * QC Gates & MultiQC Module
 * Implements quality control gates for clinical-grade MRD detection
 */

// WES Quality Control Gates
process WES_QC_GATES {
    tag "wes_qc_gates"
    publishDir "${params.outdir}/qc/wes", mode: 'copy'
    cpus 4
    memory '8 GB'
    time '2h'
    
    input:
    path tumor_metrics      // GATK CollectHsMetrics output
    path normal_metrics     // GATK CollectHsMetrics output
    path tumor_dup_metrics // MarkDuplicates metrics
    path normal_dup_metrics // MarkDuplicates metrics
    path tumor_depth       // Depth statistics
    path normal_depth      // Depth statistics
    
    output:
    path "wes_qc_summary.tsv", emit: qc_summary
    path "wes_qc_report.html", emit: qc_report
    path "wes_qc_failures.tsv", emit: qc_failures
    path "wes_qc_passed.tsv", emit: qc_passed
    
    script:
    """
    #!/bin/bash
    set -e
    echo "Running WES QC Gates..."
    
    # Parse GATK metrics
    python3 << 'PY'
    import pandas as pd
    import numpy as np
    import sys
    
    def parse_gatk_metrics(metrics_file):
        """Parse GATK CollectHsMetrics output"""
        with open(metrics_file, 'r') as f:
            lines = f.readlines()
        
        # Find metrics section
        for i, line in enumerate(lines):
            if line.startswith('## METRICS CLASS'):
                break
        
        # Parse metrics
        metrics = {}
        for line in lines[i+1:]:
            if line.startswith('#'):
                continue
            if line.strip() == '':
                break
            parts = line.strip().split('\\t')
            if len(parts) >= 2:
                try:
                    metrics[parts[0]] = float(parts[1])
                except:
                    continue
        
        return metrics
    
    def parse_dup_metrics(dup_file):
        """Parse MarkDuplicates metrics"""
        with open(dup_file, 'r') as f:
            lines = f.readlines()
        
        metrics = {}
        for line in lines:
            if line.startswith('LIBRARY'):
                continue
            if line.startswith('UNPAIRED_READ_DUPLICATES'):
                parts = line.strip().split('\\t')
                if len(parts) >= 2:
                    metrics['duplication_rate'] = float(parts[1])
                break
        
        return metrics
    
    def check_wes_qc_gates(tumor_metrics, normal_metrics, tumor_dup, normal_dup, tumor_depth, normal_depth):
        """Check WES QC gates against clinical thresholds"""
        
        # QC thresholds
        thresholds = {
            'on_target_percent': 80.0,      # Minimum on-target %
            'fold_80': 2.0,                 # Maximum fold-80
            'tumor_depth': 120,             # Minimum median on-target depth (tumor)
            'normal_depth': 80,             # Minimum median on-target depth (normal)
            'bases_30x': 30.0,              # Minimum % bases ≥30x
            'dup_rate': 20.0                # Maximum duplication rate %
        }
        
        # Parse depth statistics
        tumor_depth_stats = parse_depth_stats(tumor_depth)
        normal_depth_stats = parse_depth_stats(normal_depth)
        
        # QC checks
        qc_results = {
            'sample': ['tumor', 'normal'],
            'on_target_percent': [],
            'fold_80': [],
            'median_depth': [],
            'bases_30x': [],
            'dup_rate': [],
            'qc_status': []
        }
        
        # Tumor QC
        tumor_on_target = tumor_metrics.get('PCT_TARGET_BASES_30X', 0) * 100
        tumor_fold80 = tumor_metrics.get('FOLD_80_BASE_PENALTY', 999)
        tumor_median_depth = tumor_depth_stats.get('median_depth', 0)
        tumor_bases30x = tumor_metrics.get('PCT_TARGET_BASES_30X', 0) * 100
        tumor_dup_rate = tumor_dup.get('duplication_rate', 0) * 100
        
        # Normal QC
        normal_on_target = normal_metrics.get('PCT_TARGET_BASES_30X', 0) * 100
        normal_fold80 = normal_metrics.get('FOLD_80_BASE_PENALTY', 999)
        normal_median_depth = normal_depth_stats.get('median_depth', 0)
        normal_bases30x = normal_metrics.get('PCT_TARGET_BASES_30X', 0) * 100
        normal_dup_rate = normal_dup.get('duplication_rate', 0) * 100
        
        # Check tumor QC
        tumor_qc_passed = (
            tumor_on_target >= thresholds['on_target_percent'] and
            tumor_fold80 <= thresholds['fold_80'] and
            tumor_median_depth >= thresholds['tumor_depth'] and
            tumor_bases30x >= thresholds['bases_30x'] and
            tumor_dup_rate <= thresholds['dup_rate']
        )
        
        # Check normal QC
        normal_qc_passed = (
            normal_on_target >= thresholds['on_target_percent'] and
            normal_fold80 <= thresholds['fold_80'] and
            normal_median_depth >= thresholds['normal_depth'] and
            normal_bases30x >= thresholds['bases_30x'] and
            normal_dup_rate <= thresholds['dup_rate']
        )
        
        # Compile results
        qc_results['on_target_percent'].extend([tumor_on_target, normal_on_target])
        qc_results['fold_80'].extend([tumor_fold80, normal_fold80])
        qc_results['median_depth'].extend([tumor_median_depth, normal_median_depth])
        qc_results['bases_30x'].extend([tumor_bases30x, normal_bases30x])
        qc_results['dup_rate'].extend([tumor_dup_rate, normal_dup_rate])
        qc_results['qc_status'].extend(['PASS' if tumor_qc_passed else 'FAIL', 'PASS' if normal_qc_passed else 'FAIL'])
        
        return qc_results, tumor_qc_passed and normal_qc_passed
    
    def parse_depth_stats(depth_file):
        """Parse depth statistics file"""
        stats = {}
        try:
            with open(depth_file, 'r') as f:
                for line in f:
                    if line.startswith('median_depth'):
                        parts = line.strip().split('\\t')
                        if len(parts) >= 2:
                            stats['median_depth'] = float(parts[1])
                        break
        except:
            stats['median_depth'] = 0
        return stats
    
    # Main execution
    try:
        tumor_metrics = parse_gatk_metrics('${tumor_metrics}')
        normal_metrics = parse_gatk_metrics('${normal_metrics}')
        tumor_dup = parse_dup_metrics('${tumor_dup_metrics}')
        normal_dup = parse_dup_metrics('${normal_dup_metrics}')
        
        qc_results, overall_passed = check_wes_qc_gates(
            tumor_metrics, normal_metrics, tumor_dup, normal_dup,
            '${tumor_depth}', '${normal_depth}'
        )
        
        # Save QC summary
        qc_df = pd.DataFrame(qc_results)
        qc_df.to_csv('wes_qc_summary.tsv', sep='\\t', index=False)
        
        # Save QC failures
        failures = qc_df[qc_df['qc_status'] == 'FAIL']
        if not failures.empty:
            failures.to_csv('wes_qc_failures.tsv', sep='\\t', index=False)
        else:
            pd.DataFrame({'sample': [], 'qc_status': []}).to_csv('wes_qc_failures.tsv', sep='\\t', index=False)
        
        # Save QC passed
        passed = qc_df[qc_df['qc_status'] == 'PASS']
        passed.to_csv('wes_qc_passed.tsv', sep='\\t', index=False)
        
        # Generate HTML report
        html_report = f'''
        <html>
        <head><title>WES QC Report</title></head>
        <body>
        <h1>WES Quality Control Report</h1>
        <h2>Overall Status: {'PASS' if overall_passed else 'FAIL'}</h2>
        <table border="1">
        <tr><th>Sample</th><th>On-Target %</th><th>Fold-80</th><th>Median Depth</th><th>Bases ≥30x %</th><th>Dup Rate %</th><th>Status</th></tr>
        '''
        
        for i in range(len(qc_results['sample'])):
            html_report += f'''
        <tr>
        <td>{qc_results['sample'][i]}</td>
        <td>{qc_results['on_target_percent'][i]:.1f}%</td>
        <td>{qc_results['fold_80'][i]:.2f}</td>
        <td>{qc_results['median_depth'][i]:.0f}x</td>
        <td>{qc_results['bases_30x'][i]:.1f}%</td>
        <td>{qc_results['dup_rate'][i]:.1f}%</td>
        <td style="color: {'green' if qc_results['qc_status'][i] == 'PASS' else 'red'}">{qc_results['qc_status'][i]}</td>
        </tr>
        '''
        
        html_report += '''
        </table>
        <h3>QC Thresholds:</h3>
        <ul>
        <li>On-target %: ≥80%</li>
        <li>Fold-80: ≤2.0</li>
        <li>Tumor depth: ≥120x</li>
        <li>Normal depth: ≥80x</li>
        <li>Bases ≥30x: ≥30%</li>
        <li>Duplication rate: ≤20%</li>
        </ul>
        </body>
        </html>
        '''
        
        with open('wes_qc_report.html', 'w') as f:
            f.write(html_report)
        
        print(f"WES QC completed. Overall status: {'PASS' if overall_passed else 'FAIL'}")
        
    except Exception as e:
        print(f"Error in WES QC: {e}")
        sys.exit(1)
    PY
    
    echo "WES QC Gates completed"
    """
}

// Plasma Consensus Quality Control Gates
process PLASMA_CONSENSUS_QC_GATES {
    tag "plasma_consensus_qc"
    publishDir "${params.outdir}/qc/plasma", mode: 'copy'
    cpus 4
    memory '8 GB'
    time '2h'
    
    input:
    path consensus_bam      // Consensus BAM file
    path umi_stats         // UMI family statistics
    path depth_stats       // Depth statistics
    path insert_size_stats // Insert size statistics
    
    output:
    path "plasma_qc_summary.tsv", emit: qc_summary
    path "plasma_qc_report.html", emit: qc_report
    path "plasma_qc_failures.tsv", emit: qc_failures
    path "plasma_qc_passed.tsv", emit: qc_passed
    
    script:
    """
    #!/bin/bash
    set -e
    echo "Running Plasma Consensus QC Gates..."
    
    # Parse UMI and depth statistics
    python3 << 'PY'
    import pandas as pd
    import numpy as np
    import sys
    
    def parse_umi_stats(umi_file):
        """Parse UMI family statistics"""
        stats = {}
        try:
            with open(umi_file, 'r') as f:
                for line in f:
                    if line.startswith('median_family_size'):
                        parts = line.strip().split('\\t')
                        if len(parts) >= 2:
                            stats['median_family_size'] = float(parts[1])
                        break
                    elif line.startswith('duplex_percent'):
                        parts = line.strip().split('\\t')
                        if len(parts) >= 2:
                            stats['duplex_percent'] = float(parts[1])
                        break
        except:
            stats['median_family_size'] = 0
            stats['duplex_percent'] = 0
        return stats
    
    def parse_depth_stats(depth_file):
        """Parse depth statistics"""
        stats = {}
        try:
            with open(depth_file, 'r') as f:
                for line in f:
                    if line.startswith('raw_depth'):
                        parts = line.strip().split('\\t')
                        if len(parts) >= 2:
                            stats['raw_depth'] = float(parts[1])
                        break
                    elif line.startswith('consensus_depth'):
                        parts = line.strip().split('\\t')
                        if len(parts) >= 2:
                            stats['consensus_depth'] = float(parts[1])
                        break
        except:
            stats['raw_depth'] = 0
            stats['consensus_depth'] = 0
        return stats
    
    def parse_insert_size_stats(insert_file):
        """Parse insert size statistics"""
        stats = {}
        try:
            with open(insert_file, 'r') as f:
                for line in f:
                    if line.startswith('insert_size_mode'):
                        parts = line.strip().split('\\t')
                        if len(parts) >= 2:
                            stats['insert_size_mode'] = float(parts[1])
                        break
        except:
            stats['insert_size_mode'] = 0
        return stats
    
    def check_plasma_qc_gates(umi_stats, depth_stats, insert_stats):
        """Check plasma consensus QC gates against clinical thresholds"""
        
        # QC thresholds
        thresholds = {
            'raw_depth_min': 400,           # Minimum raw depth
            'raw_depth_max': 500,           # Maximum raw depth (optimal range)
            'consensus_depth': 60,          # Minimum consensus depth
            'median_umi_family': 3,         # Minimum median UMI family size
            'duplex_percent': 50,           # Minimum duplex percentage (if duplex)
            'insert_size_mode': 166,        # Expected insert size mode (±20)
        }
        
        # QC checks
        raw_depth = depth_stats.get('raw_depth', 0)
        consensus_depth = depth_stats.get('consensus_depth', 0)
        median_family = umi_stats.get('median_family_size', 0)
        duplex_pct = umi_stats.get('duplex_percent', 0)
        insert_mode = insert_stats.get('insert_size_mode', 0)
        
        # Check QC gates
        raw_depth_ok = thresholds['raw_depth_min'] <= raw_depth <= thresholds['raw_depth_max']
        consensus_depth_ok = consensus_depth >= thresholds['consensus_depth']
        umi_family_ok = median_family >= thresholds['median_umi_family']
        duplex_ok = duplex_pct >= thresholds['duplex_percent'] if duplex_pct > 0 else True
        insert_size_ok = abs(insert_mode - thresholds['insert_size_mode']) <= 20
        
        overall_passed = all([raw_depth_ok, consensus_depth_ok, umi_family_ok, duplex_ok, insert_size_ok])
        
        # Compile results
        qc_results = {
            'metric': ['raw_depth', 'consensus_depth', 'median_umi_family', 'duplex_percent', 'insert_size_mode'],
            'value': [raw_depth, consensus_depth, median_family, duplex_pct, insert_mode],
            'threshold': [f"{thresholds['raw_depth_min']}-{thresholds['raw_depth_max']}", 
                         f"≥{thresholds['consensus_depth']}", 
                         f"≥{thresholds['median_umi_family']}", 
                         f"≥{thresholds['duplex_percent']}%", 
                         f"{thresholds['insert_size_mode']}±20"],
            'status': ['PASS' if raw_depth_ok else 'FAIL',
                      'PASS' if consensus_depth_ok else 'FAIL',
                      'PASS' if umi_family_ok else 'FAIL',
                      'PASS' if duplex_ok else 'FAIL',
                      'PASS' if insert_size_ok else 'FAIL']
        }
        
        return qc_results, overall_passed
    
    # Main execution
    try:
        umi_stats = parse_umi_stats('${umi_stats}')
        depth_stats = parse_depth_stats('${depth_stats}')
        insert_stats = parse_insert_size_stats('${insert_size_stats}')
        
        qc_results, overall_passed = check_plasma_qc_gates(umi_stats, depth_stats, insert_stats)
        
        # Save QC summary
        qc_df = pd.DataFrame(qc_results)
        qc_df.to_csv('plasma_qc_summary.tsv', sep='\\t', index=False)
        
        # Save QC failures
        failures = qc_df[qc_df['status'] == 'FAIL']
        if not failures.empty:
            failures.to_csv('plasma_qc_failures.tsv', sep='\\t', index=False)
        else:
            pd.DataFrame({'metric': [], 'status': []}).to_csv('plasma_qc_failures.tsv', sep='\\t', index=False)
        
        # Save QC passed
        passed = qc_df[qc_df['status'] == 'PASS']
        passed.to_csv('plasma_qc_passed.tsv', sep='\\t', index=False)
        
        # Generate HTML report
        html_report = f'''
        <html>
        <head><title>Plasma Consensus QC Report</title></head>
        <body>
        <h1>Plasma Consensus Quality Control Report</h1>
        <h2>Overall Status: {'PASS' if overall_passed else 'FAIL'}</h2>
        <table border="1">
        <tr><th>Metric</th><th>Value</th><th>Threshold</th><th>Status</th></tr>
        '''
        
        for i in range(len(qc_results['metric'])):
            html_report += f'''
        <tr>
        <td>{qc_results['metric'][i]}</td>
        <td>{qc_results['value'][i]:.1f}</td>
        <td>{qc_results['threshold'][i]}</td>
        <td style="color: {'green' if qc_results['status'][i] == 'PASS' else 'red'}">{qc_results['status'][i]}</td>
        </tr>
        '''
        
        html_report += '''
        </table>
        <h3>QC Thresholds:</h3>
        <ul>
        <li>Raw depth: 400-500x (optimal range)</li>
        <li>Consensus depth: ≥60x</li>
        <li>Median UMI family: ≥3</li>
        <li>Duplex %: ≥50% (if duplex)</li>
        <li>Insert size mode: 166±20 bp</li>
        </ul>
        </body>
        </html>
        '''
        
        with open('plasma_qc_report.html', 'w') as f:
            f.write(html_report)
        
        print(f"Plasma consensus QC completed. Overall status: {'PASS' if overall_passed else 'FAIL'}")
        
    except Exception as e:
        print(f"Error in plasma consensus QC: {e}")
        sys.exit(1)
    PY
    
    echo "Plasma Consensus QC Gates completed"
    """
}

// Per-Run Quality Control Gates
process PER_RUN_QC_GATES {
    tag "per_run_qc"
    publishDir "${params.outdir}/qc/per_run", mode: 'copy'
    cpus 4
    memory '8 GB'
    time '2h'
    
    input:
    path contamination_table  // GATK CalculateContamination output
    path sex_check_results    // Sex determination results
    path pon_fp_results       // Panel of Normals false positive results
    
    output:
    path "per_run_qc_summary.tsv", emit: qc_summary
    path "per_run_qc_report.html", emit: qc_report
    path "per_run_qc_failures.tsv", emit: qc_failures
    path "per_run_qc_passed.tsv", emit: qc_passed
    
    script:
    """
    #!/bin/bash
    set -e
    echo "Running Per-Run QC Gates..."
    
    # Parse per-run QC metrics
    python3 << 'PY'
    import pandas as pd
    import numpy as np
    import sys
    
    def parse_contamination(contam_file):
        """Parse contamination table"""
        contam = 0.0
        try:
            with open(contam_file, 'r') as f:
                for line in f:
                    if line.startswith('contamination'):
                        parts = line.strip().split('\\t')
                        if len(parts) >= 2:
                            contam = float(parts[1])
                        break
        except:
            contam = 0.0
        return contam
    
    def parse_sex_check(sex_file):
        """Parse sex determination results"""
        sex_result = 'UNKNOWN'
        try:
            with open(sex_file, 'r') as f:
                for line in f:
                    if line.startswith('predicted_sex'):
                        parts = line.strip().split('\\t')
                        if len(parts) >= 2:
                            sex_result = parts[1]
                        break
        except:
            sex_result = 'UNKNOWN'
        return sex_result
    
    def parse_pon_fp(pon_file):
        """Parse Panel of Normals false positive results"""
        fp_rate = 0.0
        try:
            with open(pon_file, 'r') as f:
                for line in f:
                    if line.startswith('false_positive_rate'):
                        parts = line.strip().split('\\t')
                        if len(parts) >= 2:
                            fp_rate = float(parts[1])
                        break
        except:
            fp_rate = 0.0
        return fp_rate
    
    def check_per_run_qc_gates(contamination, sex_result, fp_rate):
        """Check per-run QC gates against clinical thresholds"""
        
        # QC thresholds
        thresholds = {
            'contamination': 0.02,         # Maximum contamination (2%)
            'sex_check': 'CONSISTENT',     # Sex determination should be consistent
            'pon_fp_rate': 0.05            # Maximum PoN false positive rate (5%)
        }
        
        # QC checks
        contam_ok = contamination <= thresholds['contamination']
        sex_ok = sex_result != 'UNKNOWN' and sex_result != 'INCONSISTENT'
        fp_ok = fp_rate <= thresholds['pon_fp_rate']
        
        overall_passed = all([contam_ok, sex_ok, fp_ok])
        
        # Compile results
        qc_results = {
            'metric': ['contamination', 'sex_check', 'pon_fp_rate'],
            'value': [f"{contamination:.3f}", sex_result, f"{fp_rate:.3f}"],
            'threshold': [f"≤{thresholds['contamination']:.3f}", 
                         thresholds['sex_check'], 
                         f"≤{thresholds['pon_fp_rate']:.3f}"],
            'status': ['PASS' if contam_ok else 'FAIL',
                      'PASS' if sex_ok else 'FAIL',
                      'PASS' if fp_ok else 'FAIL']
        }
        
        return qc_results, overall_passed
    
    # Main execution
    try:
        contamination = parse_contamination('${contamination_table}')
        sex_result = parse_sex_check('${sex_check_results}')
        fp_rate = parse_pon_fp('${pon_fp_results}')
        
        qc_results, overall_passed = check_per_run_qc_gates(contamination, sex_result, fp_rate)
        
        # Save QC summary
        qc_df = pd.DataFrame(qc_results)
        qc_df.to_csv('per_run_qc_summary.tsv', sep='\\t', index=False)
        
        # Save QC failures
        failures = qc_df[qc_df['status'] == 'FAIL']
        if not failures.empty:
            failures.to_csv('per_run_qc_failures.tsv', sep='\\t', index=False)
        else:
            pd.DataFrame({'metric': [], 'status': []}).to_csv('per_run_qc_failures.tsv', sep='\\t', index=False)
        
        # Save QC passed
        passed = qc_df[qc_df['status'] == 'PASS']
        passed.to_csv('per_run_qc_passed.tsv', sep='\\t', index=False)
        
        # Generate HTML report
        html_report = f'''
        <html>
        <head><title>Per-Run QC Report</title></head>
        <body>
        <h1>Per-Run Quality Control Report</h1>
        <h2>Overall Status: {'PASS' if overall_passed else 'FAIL'}</h2>
        <table border="1">
        <tr><th>Metric</th><th>Value</th><th>Threshold</th><th>Status</th></tr>
        '''
        
        for i in range(len(qc_results['metric'])):
            html_report += f'''
        <tr>
        <td>{qc_results['metric'][i]}</td>
        <td>{qc_results['value'][i]}</td>
        <td>{qc_results['threshold'][i]}</td>
        <td style="color: {'green' if qc_results['status'][i] == 'PASS' else 'red'}">{qc_results['status'][i]}</td>
        </tr>
        '''
        
        html_report += '''
        </table>
        <h3>QC Thresholds:</h3>
        <ul>
        <li>Contamination: ≤2%</li>
        <li>Sex check: Consistent determination</li>
        <li>PoN FP rate: ≤5%</li>
        </ul>
        </body>
        </html>
        '''
        
        with open('per_run_qc_report.html', 'w') as f:
            f.write(html_report)
        
        print(f"Per-run QC completed. Overall status: {'PASS' if overall_passed else 'FAIL'}")
        
    except Exception as e:
        print(f"Error in per-run QC: {e}")
        sys.exit(1)
    PY
    
    echo "Per-Run QC Gates completed"
    """
}

// MultiQC Integration
process MULTIQC_INTEGRATION {
    tag "multiqc_integration"
    publishDir "${params.outdir}/reports/multiqc", mode: 'copy'
    cpus 2
    memory '4 GB'
    time '1h'
    
    input:
    path wes_qc_dir      // WES QC directory
    path plasma_qc_dir   // Plasma QC directory
    path per_run_qc_dir // Per-run QC directory
    
    output:
    path "multiqc_report.html", emit: multiqc_report
    path "multiqc_data/", emit: multiqc_data
    path "multiqc.log", emit: multiqc_log
    
    script:
    """
    #!/bin/bash
    set -e
    echo "Running MultiQC integration..."
    
    # Create MultiQC configuration
    cat > multiqc_config.yaml << 'EOF'
    report_comment: "Tumor-Informed cfDNA MRD Pipeline QC Report"
    report_section_order:
      - "WES_QC_Gates"
      - "Plasma_Consensus_QC"
      - "Per_Run_QC"
      - "General_Stats"
    
    # Custom section headers
    section_comments:
      WES_QC_Gates: "WES Quality Control Gates - Clinical Thresholds"
      Plasma_Consensus_QC: "Plasma Consensus Quality Control - UMI-based Metrics"
      Per_Run_QC: "Per-Run Quality Control - Contamination & Validation"
    
    # File patterns
    file_patterns:
      wes_qc: "wes_qc_summary.tsv"
      plasma_qc: "plasma_qc_summary.tsv"
      per_run_qc: "per_run_qc_summary.tsv"
    
    # Custom plots
    custom_plots:
      qc_gates_summary:
        plot_type: "table"
        data_labels:
          - "WES QC"
          - "Plasma QC"
          - "Per-Run QC"
    EOF
    
    # Run MultiQC
    multiqc \
        --config multiqc_config.yaml \
        --outdir . \
        --filename multiqc_report.html \
        --force \
        ${wes_qc_dir} \
        ${plasma_qc_dir} \
        ${per_run_qc_dir}
    
    # Generate custom QC summary
    python3 << 'PY'
    import pandas as pd
    import os
    import glob
    
    def compile_qc_summary():
        """Compile comprehensive QC summary from all modules"""
        
        qc_summary = {
            'module': [],
            'total_samples': [],
            'passed': [],
            'failed': [],
            'pass_rate': [],
            'critical_failures': []
        }
        
        # WES QC summary
        wes_qc_files = glob.glob('${wes_qc_dir}/*_qc_summary.tsv')
        if wes_qc_files:
            wes_df = pd.read_csv(wes_qc_files[0], sep='\\t')
            total = len(wes_df)
            passed = len(wes_df[wes_df['qc_status'] == 'PASS'])
            failed = total - passed
            pass_rate = (passed / total * 100) if total > 0 else 0
            
            qc_summary['module'].append('WES_QC')
            qc_summary['total_samples'].append(total)
            qc_summary['passed'].append(passed)
            qc_summary['failed'].append(failed)
            qc_summary['pass_rate'].append(f"{pass_rate:.1f}%")
            qc_summary['critical_failures'].append('None' if failed == 0 else f"{failed} samples")
        
        # Plasma QC summary
        plasma_qc_files = glob.glob('${plasma_qc_dir}/*_qc_summary.tsv')
        if plasma_qc_files:
            plasma_df = pd.read_csv(plasma_qc_files[0], sep='\\t')
            total = len(plasma_df)
            passed = len(plasma_df[plasma_df['status'] == 'PASS'])
            failed = total - passed
            pass_rate = (passed / total * 100) if total > 0 else 0
            
            qc_summary['module'].append('Plasma_QC')
            qc_summary['total_samples'].append(total)
            qc_summary['passed'].append(passed)
            qc_summary['failed'].append(failed)
            qc_summary['pass_rate'].append(f"{pass_rate:.1f}%")
            qc_summary['critical_failures'].append('None' if failed == 0 else f"{failed} metrics")
        
        # Per-run QC summary
        per_run_qc_files = glob.glob('${per_run_qc_dir}/*_qc_summary.tsv')
        if per_run_qc_files:
            per_run_df = pd.read_csv(per_run_qc_files[0], sep='\\t')
            total = len(per_run_df)
            passed = len(per_run_df[per_run_df['status'] == 'PASS'])
            failed = total - passed
            pass_rate = (passed / total * 100) if total > 0 else 0
            
            qc_summary['module'].append('Per_Run_QC')
            qc_summary['total_samples'].append(total)
            qc_summary['passed'].append(passed)
            qc_summary['failed'].append(failed)
            qc_summary['pass_rate'].append(f"{pass_rate:.1f}%")
            qc_summary['critical_failures'].append('None' if failed == 0 else f"{failed} metrics")
        
        # Save summary
        summary_df = pd.DataFrame(qc_summary)
        summary_df.to_csv('qc_summary_compiled.tsv', sep='\\t', index=False)
        
        print("QC summary compilation completed")
    
    compile_qc_summary()
    PY
    
    echo "MultiQC integration completed"
    """
}
