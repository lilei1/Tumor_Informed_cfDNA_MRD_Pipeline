#!/usr/bin/env nextflow

/*
 * Validation Module
 * Benchmarks LoD, sensitivity/specificity, and longitudinal precision
 */

process VALIDATION {
    tag "validation"
    
    publishDir "${params.outdir}/validation", mode: 'copy'
    
    cpus 8
    memory '32 GB'
    time '4h'
    
    input:
    tuple val(sample), val(timepoint), path(mrd_scores_csv)
    tuple val(sample), val(timepoint), path(features_csv)
    val outdir
    
    output:
    tuple val("validation"), path("*.validation_metrics.csv"), emit: validation_metrics
    path "*.log", emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    
    # Create validation script
    cat > validation_analysis.py << 'EOF'
    import pandas as pd
    import numpy as np
    from sklearn.metrics import roc_auc_score, precision_recall_curve, confusion_matrix
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os
    
    def calculate_lod_metrics(features_df):
        """Calculate Limit of Detection metrics"""
        # This is a simplified LOD calculation
        # In practice, you would use spike-in data with known concentrations
        
        # Calculate background noise from healthy samples
        healthy_features = features_df[features_df['sample'].str.contains('healthy', case=False)]
        
        if len(healthy_features) > 0:
            background_noise = healthy_features['variant_count'].mean()
            lod_95 = background_noise + (1.96 * healthy_features['variant_count'].std())
        else:
            lod_95 = 0.001  # Default LOD
        
        return {
            'background_noise': background_noise,
            'lod_95': lod_95,
            'lod_99': lod_95 * 1.5
        }
    
    def calculate_sensitivity_specificity(mrd_scores_df, truth_labels=None):
        """Calculate sensitivity and specificity"""
        if truth_labels is None:
            # Use timepoint T0 as baseline (assuming T0 should be positive for patients)
            # This is a simplified approach - in practice use known clinical outcomes
            
            # Separate patient and healthy samples
            patient_scores = mrd_scores_df[~mrd_scores_df['sample'].str.contains('healthy', case=False)]
            healthy_scores = mrd_scores_df[mrd_scores_df['sample'].str.contains('healthy', case=False)]
            
            if len(patient_scores) > 0 and len(healthy_scores) > 0:
                # Use T0 timepoint for patients as positive, healthy as negative
                patient_t0 = patient_scores[patient_scores['timepoint'] == 0]
                healthy_t0 = healthy_scores[healthy_scores['timepoint'] == 0]
                
                if len(patient_t0) > 0 and len(healthy_t0) > 0:
                    # Calculate ROC AUC
                    y_true = [1] * len(patient_t0) + [0] * len(healthy_t0)
                    y_scores = list(patient_t0['mrd_score']) + list(healthy_t0['mrd_score'])
                    
                    try:
                        auc = roc_auc_score(y_true, y_scores)
                    except:
                        auc = 0.5
                    
                    # Calculate optimal threshold
                    precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
                    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-8)
                    optimal_threshold = thresholds[np.argmax(f1_scores[:-1])]
                    
                    # Calculate sensitivity and specificity at optimal threshold
                    y_pred = [1 if score >= optimal_threshold else 0 for score in y_scores]
                    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
                    
                    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
                    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
                    
                    return {
                        'auc': auc,
                        'optimal_threshold': optimal_threshold,
                        'sensitivity': sensitivity,
                        'specificity': specificity,
                        'precision': tp / (tp + fp) if (tp + fp) > 0 else 0,
                        'recall': sensitivity
                    }
        
        return {
            'auc': 0.5,
            'optimal_threshold': 0.5,
            'sensitivity': 0.0,
            'specificity': 0.0,
            'precision': 0.0,
            'recall': 0.0
        }
    
    def calculate_longitudinal_precision(mrd_scores_df):
        """Calculate longitudinal precision metrics"""
        # Group by sample and calculate consistency across timepoints
        sample_groups = mrd_scores_df.groupby('sample')
        
        longitudinal_metrics = []
        for sample, group in sample_groups:
            if len(group) > 1:  # Multiple timepoints
                timepoints = sorted(group['timepoint'].unique())
                scores = [group[group['timepoint'] == tp]['mrd_score'].iloc[0] for tp in timepoints]
                
                # Calculate consistency metrics
                score_variance = np.var(scores)
                score_range = max(scores) - min(scores)
                trend_consistency = 1.0 if len(scores) <= 2 else np.corrcoef(timepoints, scores)[0, 1]
                
                longitudinal_metrics.append({
                    'sample': sample,
                    'timepoints': len(timepoints),
                    'score_variance': score_variance,
                    'score_range': score_range,
                    'trend_consistency': trend_consistency if not np.isnan(trend_consistency) else 0.0
                })
        
        if longitudinal_metrics:
            df = pd.DataFrame(longitudinal_metrics)
            return {
                'mean_variance': df['score_variance'].mean(),
                'mean_range': df['score_range'].mean(),
                'mean_consistency': df['trend_consistency'].mean(),
                'samples_with_multiple_timepoints': len(df)
            }
        
        return {
            'mean_variance': 0.0,
            'mean_range': 0.0,
            'mean_consistency': 0.0,
            'samples_with_multiple_timepoints': 0
        }
    
    def main():
        # Load data
        mrd_files = [f for f in os.listdir('.') if f.endswith('_mrd_scores.csv')]
        feature_files = [f for f in os.listdir('.') if f.endswith('_features.csv')]
        
        if not mrd_files or not feature_files:
            print("No MRD scores or features found")
            return
        
        # Load and combine data
        mrd_data = []
        for file in mrd_files:
            df = pd.read_csv(file)
            mrd_data.append(df)
        
        if mrd_data:
            mrd_df = pd.concat(mrd_data, ignore_index=True)
        else:
            print("No MRD data loaded")
            return
        
        # Load features
        feature_data = []
        for file in feature_files:
            df = pd.read_csv(file)
            feature_data.append(df)
        
        if feature_data:
            features_df = pd.concat(feature_data, ignore_index=True)
        else:
            features_df = pd.DataFrame()
        
        # Calculate validation metrics
        lod_metrics = calculate_lod_metrics(features_df)
        performance_metrics = calculate_sensitivity_specificity(mrd_df)
        longitudinal_metrics = calculate_longitudinal_precision(mrd_df)
        
        # Combine all metrics
        validation_results = {
            'metric_type': 'validation_summary',
            'lod_background_noise': lod_metrics['background_noise'],
            'lod_95': lod_metrics['lod_95'],
            'lod_99': lod_metrics['lod_99'],
            'roc_auc': performance_metrics['auc'],
            'optimal_threshold': performance_metrics['optimal_threshold'],
            'sensitivity': performance_metrics['sensitivity'],
            'specificity': performance_metrics['specificity'],
            'precision': performance_metrics['precision'],
            'recall': performance_metrics['recall'],
            'longitudinal_variance': longitudinal_metrics['mean_variance'],
            'longitudinal_range': longitudinal_metrics['mean_range'],
            'longitudinal_consistency': longitudinal_metrics['mean_consistency'],
            'samples_multiple_timepoints': longitudinal_metrics['samples_with_multiple_timepoints']
        }
        
        # Save results
        results_df = pd.DataFrame([validation_results])
        results_df.to_csv('validation_metrics.csv', index=False)
        
        # Generate validation plots
        if len(mrd_df) > 0:
            # ROC curve
            plt.figure(figsize=(12, 8))
            
            plt.subplot(2, 2, 1)
            # Sample distribution by MRD score
            plt.hist(mrd_df['mrd_score'], bins=20, alpha=0.7, edgecolor='black')
            plt.xlabel('MRD Score')
            plt.ylabel('Frequency')
            plt.title('Distribution of MRD Scores')
            
            plt.subplot(2, 2, 2)
            # Longitudinal consistency
            sample_groups = mrd_df.groupby('sample')
            for sample, group in sample_groups:
                if len(group) > 1:
                    plt.plot(group['timepoint'], group['mrd_score'], 'o-', label=sample, alpha=0.7)
            plt.xlabel('Timepoint')
            plt.ylabel('MRD Score')
            plt.title('Longitudinal MRD Scores')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            plt.tight_layout()
            plt.savefig('validation_summary.png', dpi=300, bbox_inches='tight')
            plt.close()
        
        print("Validation analysis completed")
    
    if __name__ == "__main__":
        main()
    EOF
    
    # Run validation analysis
    python3 validation_analysis.py
    
    # Clean up
    rm -f validation_analysis.py
    """
}
