#!/usr/bin/env nextflow

/*
 * Reporting Module
 * Generates QC dashboards, audit-ready logs, and one-click reports
 */

process REPORTING {
    tag "reporting"
    
    publishDir "${params.outdir}/reports", mode: 'copy'
    
    cpus 4
    memory '16 GB'
    time '2h'
    
    input:
    tuple val(sample), path(truth_set_vcf), path(truth_set_vcf_idx)
    tuple val(sample), val(timepoint), path(features_csv)
    tuple val(sample), val(timepoint), path(mrd_scores_csv)
    tuple val("validation"), path(validation_metrics_csv)
    val outdir
    
    output:
    tuple val("reports"), path("*.html"), path("*.pdf"), emit: reports
    path "*.log", emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    
    # Create comprehensive reporting script
    cat > generate_reports.py << 'EOF'
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from datetime import datetime
    import os
    import glob
    
    def load_all_data():
        """Load all pipeline outputs"""
        data = {}
        
        # Load MRD scores
        mrd_files = glob.glob('*_mrd_scores.csv')
        if mrd_files:
            mrd_data = []
            for file in mrd_files:
                df = pd.read_csv(file)
                mrd_data.append(df)
            data['mrd_scores'] = pd.concat(mrd_data, ignore_index=True) if mrd_data else pd.DataFrame()
        
        # Load features
        feature_files = glob.glob('*_features.csv')
        if feature_files:
            feature_data = []
            for file in feature_files:
                df = pd.read_csv(file)
                feature_data.append(df)
            data['features'] = pd.concat(feature_data, ignore_index=True) if feature_data else pd.DataFrame()
        
        # Load validation metrics
        validation_files = glob.glob('*validation_metrics.csv')
        if validation_files:
            data['validation'] = pd.read_csv(validation_files[0])
        
        return data
    
    def generate_qc_dashboard(data):
        """Generate QC dashboard with key metrics"""
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('MRD Pipeline QC Dashboard', fontsize=16, fontweight='bold')
        
        # Sample count by timepoint
        if not data['mrd_scores'].empty:
            timepoint_counts = data['mrd_scores'].groupby('timepoint').size()
            axes[0, 0].bar(timepoint_counts.index, timepoint_counts.values, color='skyblue')
            axes[0, 0].set_xlabel('Timepoint')
            axes[0, 0].set_ylabel('Sample Count')
            axes[0, 0].set_title('Samples per Timepoint')
            axes[0, 0].grid(True, alpha=0.3)
        
        # MRD score distribution
        if not data['mrd_scores'].empty:
            axes[0, 1].hist(data['mrd_scores']['mrd_score'], bins=20, alpha=0.7, color='lightgreen', edgecolor='black')
            axes[0, 1].set_xlabel('MRD Score')
            axes[0, 1].set_ylabel('Frequency')
            axes[0, 1].set_title('MRD Score Distribution')
            axes[0, 1].grid(True, alpha=0.3)
        
        # Fragment size distribution
        if not data['features'].empty and 'fragment_count' in data['features'].columns:
            fragment_sizes = ['size_50_150', 'size_151_250', 'size_251_350', 'size_351_450', 'size_451_1000']
            size_labels = ['50-150', '151-250', '251-350', '351-450', '451-1000']
            
            if all(col in data['features'].columns for col in fragment_sizes):
                size_counts = [data['features'][col].sum() for col in fragment_sizes]
                axes[0, 2].pie(size_counts, labels=size_labels, autopct='%1.1f%%', startangle=90)
                axes[0, 2].set_title('Fragment Size Distribution')
        
        # Longitudinal MRD trends
        if not data['mrd_scores'].empty:
            sample_groups = data['mrd_scores'].groupby('sample')
            for sample, group in sample_groups:
                if len(group) > 1:
                    axes[1, 0].plot(group['timepoint'], group['mrd_score'], 'o-', label=sample, alpha=0.7)
            axes[1, 0].set_xlabel('Timepoint')
            axes[1, 0].set_ylabel('MRD Score')
            axes[1, 0].set_title('Longitudinal MRD Trends')
            axes[1, 0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            axes[1, 0].grid(True, alpha=0.3)
        
        # Validation metrics
        if 'validation' in data and not data['validation'].empty:
            metrics = data['validation'].iloc[0]
            metric_names = ['sensitivity', 'specificity', 'precision', 'recall']
            metric_values = [metrics.get(name, 0) for name in metric_names]
            
            axes[1, 1].bar(metric_names, metric_values, color=['red', 'blue', 'green', 'orange'])
            axes[1, 1].set_ylabel('Score')
            axes[1, 1].set_title('Performance Metrics')
            axes[1, 1].set_ylim(0, 1)
            axes[1, 1].grid(True, alpha=0.3)
        
        # Sample quality metrics
        if not data['features'].empty:
            quality_metrics = ['fragment_count', 'mean_size', 'std_size']
            if all(col in data['features'].columns for col in quality_metrics):
                quality_data = data['features'][quality_metrics].describe()
                axes[1, 2].text(0.1, 0.9, quality_data.to_string(), transform=axes[1, 2].transAxes, 
                               fontsize=10, verticalalignment='top', fontfamily='monospace')
                axes[1, 2].set_title('Quality Metrics Summary')
                axes[1, 2].axis('off')
        
        plt.tight_layout()
        plt.savefig('qc_dashboard.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        return 'qc_dashboard.png'
    
    def generate_html_report(data, dashboard_path):
        """Generate HTML report"""
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>MRD Pipeline Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
                .section {{ margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }}
                .metric {{ display: inline-block; margin: 10px; padding: 10px; background-color: #f9f9f9; border-radius: 3px; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
                .dashboard {{ text-align: center; margin: 20px 0; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Tumor-Informed cfDNA MRD Pipeline Report</h1>
                <p>Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                <p>Pipeline Version: 1.0.0</p>
            </div>
            
            <div class="section">
                <h2>Pipeline Summary</h2>
                <div class="metric">
                    <strong>Total Samples:</strong> {len(data['mrd_scores']['sample'].unique()) if not data['mrd_scores'].empty else 0}
                </div>
                <div class="metric">
                    <strong>Timepoints:</strong> {len(data['mrd_scores']['timepoint'].unique()) if not data['mrd_scores'].empty else 0}
                </div>
                <div class="metric">
                    <strong>MRD Positive:</strong> {len(data['mrd_scores'][data['mrd_scores']['mrd_status'] == 'Positive']) if not data['mrd_scores'].empty else 0}
                </div>
            </div>
            
            <div class="section">
                <h2>QC Dashboard</h2>
                <div class="dashboard">
                    <img src="{dashboard_path}" alt="QC Dashboard" style="max-width: 100%; height: auto;">
                </div>
            </div>
        """
        
        # Add validation metrics
        if 'validation' in data and not data['validation'].empty:
            metrics = data['validation'].iloc[0]
            html_content += f"""
            <div class="section">
                <h2>Validation Metrics</h2>
                <table>
                    <tr><th>Metric</th><th>Value</th></tr>
                    <tr><td>ROC AUC</td><td>{metrics.get('roc_auc', 'N/A'):.3f}</td></tr>
                    <tr><td>Sensitivity</td><td>{metrics.get('sensitivity', 'N/A'):.3f}</td></tr>
                    <tr><td>Specificity</td><td>{metrics.get('specificity', 'N/A'):.3f}</td></tr>
                    <tr><td>Precision</td><td>{metrics.get('precision', 'N/A'):.3f}</td></tr>
                    <tr><td>LOD (95%)</td><td>{metrics.get('lod_95', 'N/A'):.6f}</td></tr>
                </table>
            </div>
            """
        
        # Add sample-level results
        if not data['mrd_scores'].empty:
            html_content += """
            <div class="section">
                <h2>Sample Results</h2>
                <table>
                    <tr><th>Sample</th><th>Timepoint</th><th>MRD Score</th><th>Status</th><th>Confidence</th></tr>
            """
            
            for _, row in data['mrd_scores'].iterrows():
                html_content += f"""
                    <tr>
                        <td>{row['sample']}</td>
                        <td>T{row['timepoint']}</td>
                        <td>{row['mrd_score']:.3f}</td>
                        <td>{row['mrd_status']}</td>
                        <td>{row['confidence']}</td>
                    </tr>
                """
            
            html_content += """
                </table>
            </div>
            """
        
        html_content += """
            <div class="section">
                <h2>Pipeline Information</h2>
                <p><strong>Objective:</strong> Detect Minimal Residual Disease in CRC patients using tumor-informed cfDNA analysis</p>
                <p><strong>Features:</strong> Variants, methylation, and fragmentomics integration</p>
                <p><strong>Technology:</strong> UMI-tagged WGS with ultra-low VAF detection</p>
            </div>
        </body>
        </html>
        """
        
        with open('pipeline_report.html', 'w') as f:
            f.write(html_content)
        
        return 'pipeline_report.html'
    
    def main():
        print("Loading pipeline data...")
        data = load_all_data()
        
        print("Generating QC dashboard...")
        dashboard_path = generate_qc_dashboard(data)
        
        print("Generating HTML report...")
        html_path = generate_html_report(data, dashboard_path)
        
        print("Generating audit log...")
        audit_log = f"""
        MRD Pipeline Audit Log
        =======================
        Timestamp: {datetime.now().isoformat()}
        Pipeline Version: 1.0.0
        
        Data Summary:
        - MRD Scores: {len(data['mrd_scores']) if not data['mrd_scores'].empty else 0} records
        - Features: {len(data['features']) if not data['features'].empty else 0} records
        - Validation Metrics: {len(data['validation']) if 'validation' in data and not data['validation'].empty else 0} records
        
        Output Files:
        - QC Dashboard: {dashboard_path}
        - HTML Report: {html_path}
        
        Pipeline completed successfully.
        """
        
        with open('audit_log.txt', 'w') as f:
            f.write(audit_log)
        
        print("Report generation completed successfully!")
    
    if __name__ == "__main__":
        main()
    EOF
    
    # Run report generation
    python3 generate_reports.py
    
    # Convert HTML to PDF (if wkhtmltopdf is available)
    if command -v wkhtmltopdf >/dev/null 2>&1; then
        wkhtmltopdf --page-size A4 --orientation Landscape pipeline_report.html pipeline_report.pdf
    else
        echo "wkhtmltopdf not available, skipping PDF generation"
    fi
    
    # Clean up
    rm -f generate_reports.py
    """
}
