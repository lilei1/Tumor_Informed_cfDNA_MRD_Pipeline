#!/usr/bin/env nextflow

/*
 * MRD Classification Module
 * Integrates features to predict MRD probability with calibrated thresholds
 */

process MRD_CLASSIFICATION {
    tag "${sample}"
    
    publishDir "${params.outdir}/classification", mode: 'copy'
    
    cpus 4
    memory '16 GB'
    time '2h'
    
    input:
    tuple val(sample), val(timepoint), path(features_csv)
    tuple val(sample), path(truth_set_vcf), path(truth_set_vcf_idx)
    tuple val(sample), val(timepoint), path(plasma_bam), path(plasma_bai)
    
    output:
    tuple val(sample), val(timepoint), path("*.mrd_score.csv"), emit: mrd_scores
    tuple val(sample), val(timepoint), path("*.mrd_probability.csv"), emit: mrd_probabilities
    path "*.log", emit: logs
    
    script:
    """
    #!/bin/bash
    set -e
    
    # Create Python script for MRD classification
    cat > mrd_classifier.py << 'EOF'
    import pandas as pd
    import numpy as np
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_auc_score, precision_recall_curve
    import joblib
    import os
    
    def load_features(feature_files):
        """Load and combine all feature files for a sample"""
        all_features = []
        for file in feature_files:
            if os.path.exists(file):
                df = pd.read_csv(file)
                all_features.append(df)
        
        if not all_features:
            return None
        
        return pd.concat(all_features, ignore_index=True)
    
    def extract_variant_features(vcf_file, bam_file):
        """Extract variant-related features from VCF and BAM"""
        # This is a simplified version - implement based on your needs
        features = {
            'variant_count': 0,
            'mean_vaf': 0.0,
            'max_vaf': 0.0,
            'variant_burden': 0.0
        }
        return features
    
    def calculate_mrd_score(features, model=None):
        """Calculate MRD score using trained model or heuristics"""
        if model is None:
            # Simple heuristic-based scoring
            score = 0.0
            
            # Fragmentomics contribution
            if features['fragment_count'] > 0:
                size_ratio = features['size_50_150'] / features['fragment_count']
                score += size_ratio * 0.3
            
            # Variant contribution
            if features['variant_count'] > 0:
                score += min(features['variant_count'] * 0.1, 0.4)
            
            # VAF contribution
            score += min(features['mean_vaf'] * 2.0, 0.3)
            
            return min(score, 1.0)
        else:
            # Use trained model
            X = features.drop(['sample', 'timepoint'], axis=1)
            return model.predict_proba(X)[0][1]
    
    def main():
        # Load features
        feature_files = [f for f in os.listdir('.') if f.endswith('_features.csv')]
        features = load_features(feature_files)
        
        if features is None:
            print("No feature files found")
            return
        
        # Extract sample and timepoint info
        sample = features['sample'].iloc[0]
        timepoints = sorted(features['timepoint'].unique())
        
        # Calculate MRD scores for each timepoint
        results = []
        for tp in timepoints:
            tp_features = features[features['timepoint'] == tp].iloc[0]
            
            # Calculate MRD score
            mrd_score = calculate_mrd_score(tp_features)
            
            # Determine MRD status based on threshold
            mrd_status = 'Positive' if mrd_score > 0.5 else 'Negative'
            confidence = 'High' if abs(mrd_score - 0.5) > 0.2 else 'Low'
            
            result = {
                'sample': sample,
                'timepoint': tp,
                'mrd_score': mrd_score,
                'mrd_status': mrd_status,
                'confidence': confidence,
                'fragment_count': tp_features['fragment_count'],
                'variant_count': tp_features.get('variant_count', 0),
                'mean_vaf': tp_features.get('mean_vaf', 0.0)
            }
            results.append(result)
        
        # Save results
        results_df = pd.DataFrame(results)
        results_df.to_csv(f'{sample}_mrd_scores.csv', index=False)
        
        # Create probability matrix
        prob_matrix = results_df.pivot(index='sample', columns='timepoint', values='mrd_score')
        prob_matrix.to_csv(f'{sample}_mrd_probability.csv')
        
        print(f"MRD classification completed for {sample}")
    
    if __name__ == "__main__":
        main()
    EOF
    
    # Run MRD classification
    python3 mrd_classifier.py
    
    # Generate longitudinal plots
    python3 -c "
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # Load MRD scores
    scores = pd.read_csv('${sample}_mrd_scores.csv')
    
    # Create longitudinal plot
    plt.figure(figsize=(10, 6))
    plt.plot(scores['timepoint'], scores['mrd_score'], 'o-', linewidth=2, markersize=8)
    plt.axhline(y=0.5, color='r', linestyle='--', alpha=0.7, label='MRD Threshold')
    plt.fill_between(scores['timepoint'], 0.5, 1.0, alpha=0.2, color='red', label='MRD Positive')
    plt.fill_between(scores['timepoint'], 0.0, 0.5, alpha=0.2, color='green', label='MRD Negative')
    
    plt.xlabel('Timepoint')
    plt.ylabel('MRD Score')
    plt.title(f'Longitudinal MRD Scores - {sample}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('${sample}_longitudinal_mrd.png', dpi=300, bbox_inches='tight')
    plt.close()
    "
    
    # Clean up
    rm -f mrd_classifier.py
    """
}
