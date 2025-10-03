# -*- coding: utf-8 -*-
"""
Linear Regression Model for Protein Expression Prediction from RNA Data
=====================================================================

This script test/infer the baseline for the STP Challenge.

Authors: Zijie Fang @ WangLAB in HKUST
Date: 2025-10-03
"""

import scanpy as sc
import numpy as np
import pandas as pd
import pickle
import os

from scipy.sparse import issparse

docker = False # set to True if running in docker
if not docker:
    INPUT_DIR = "./input"
    OUTPUT_DIR = "./output"
else:
    INPUT_DIR = "/app/input"
    OUTPUT_DIR = "/app/output"


def process_validation_data(val_data_path, significant_genes):
    """
    Process validation RNA data using the same pipeline as training data.
    
    Args:
        val_data_path: Path to validation data file
        significant_genes: List of significant genes selected during training
        
    Returns:
        tuple: (preprocessed_data_df, position_df)
    """
    print(f"📂 Loading validation data from: {val_data_path}")
    val_rna = sc.read_h5ad(val_data_path)
    print(f"   Validation data shape: {val_rna.shape[0]} spots × {val_rna.shape[1]} genes")
    
    # Apply the same preprocessing pipeline
    # 1. CPM normalization
    counts = val_rna.X
    if issparse(counts):
        total_counts = counts.sum(axis=1).A1
        cpm = counts.multiply(1e4 / (total_counts[:, None] + 1e-8))
    else:
        total_counts = counts.sum(axis=1)
        cpm = counts / (total_counts[:, None] + 1e-8) * 1e4
    
    # 2. log1p transformation
    if issparse(cpm):
        cpm_log1p = cpm.copy()
        cpm_log1p.data = np.log1p(cpm_log1p.data)
    else:
        cpm_log1p = np.log1p(cpm)
    
    # 3. Extract significant gene expression
    gene_indices = [val_rna.var_names.get_loc(g) for g in significant_genes if g in val_rna.var_names]
    
    if issparse(cpm_log1p):
        cpm_log1p_dense = cpm_log1p.tocsr()[:, gene_indices].toarray()
        val_expr_df = pd.DataFrame(
            cpm_log1p_dense,
            index=val_rna.obs_names,
            columns=[g for g in significant_genes if g in val_rna.var_names]
        )
    else:
        val_expr_df = pd.DataFrame(
            cpm_log1p[:, gene_indices],
            index=val_rna.obs_names,
            columns=[g for g in significant_genes if g in val_rna.var_names]
        )
    
    # Extract position information
    POSITION_COLS = ['pxl_row_in_fullres', 'pxl_col_in_fullres']
    position_df = val_rna.obs[POSITION_COLS].copy()
    
    print(f"✅ Validation data preprocessing completed: {val_expr_df.shape}")
    return val_expr_df, position_df


def predict_protein_expression(models, val_expr_df, target_proteins):
    """
    Predict protein expression for validation data using trained models.
    
    Args:
        models: Dictionary of trained models
        val_expr_df: Validation gene expression DataFrame
        target_proteins: List of target protein names
        
    Returns:
        pd.DataFrame: Predicted protein expression DataFrame
    """
    print("🔮 Predicting protein expression for validation data...")
    
    # Initialize prediction dataframe
    val_pred_protein = pd.DataFrame(
        index=val_expr_df.index,
        columns=target_proteins
    )
    
    # Make predictions for each protein
    X_val = val_expr_df.values
    for protein in target_proteins:
        if protein in models:
            print(f"   Predicting {protein}...")
            val_pred_protein.loc[:, protein] = models[protein].predict(X_val)
        else:
            print(f"   Warning: Model for {protein} not found")
            val_pred_protein.loc[:, protein] = 0.0
    
    print(f"✅ Protein predictions completed: {val_pred_protein.shape}")
    return val_pred_protein


def inverse_normalize_predictions(predictions_df, scaler_mean, scaler_var):
    """
    Apply inverse z-score normalization and log1p transformation to predictions.
    
    Args:
        predictions_df: DataFrame with normalized predictions
        scaler_mean: Mean values from StandardScaler
        scaler_var: Variance values from StandardScaler
        
    Returns:
        pd.DataFrame: Denormalized predictions
    """
    print("   Applying inverse z-score normalization...")
    
    # Initialize denormalized dataframe
    denorm_df = predictions_df.copy()
    
    # Apply inverse z-score normalization
    for i, col in enumerate(predictions_df.columns):
        if i < len(scaler_mean):
            # denorm_df[col] = (predictions_df[col].values * np.sqrt(scaler_var[i])) + scaler_mean[i]
            denorm_df[col] = (predictions_df[col].values * scaler_var[i]) + scaler_mean[i]
        else:
            print(f"   Warning: No standardization parameters for {col}")
            
    print("   Applying inverse log1p transformation...")
    
    # Apply inverse log1p transformation (exponm1)
    denorm_df = denorm_df.applymap(lambda x: np.expm1(x))
    
    return denorm_df


def create_submission_file(denorm_predictions, position_data, output_filename=None):
    """
    Create final submission file with position and protein predictions.
    
    Args:
        denorm_predictions: DataFrame with denormalized protein predictions
        position_data: DataFrame with position information
        output_filename: Optional custom filename for output
        
    Returns:
        str: Path to saved submission file
    """
    # Combine position and protein data
    submission_df = pd.concat([position_data, denorm_predictions], axis=1)
    
    # Generate filename with timestamp if not provided
    if output_filename is None:
        output_filename = f"submission.csv"
    
    # Set barcode as index name
    submission_df.index.name = 'barcode'
    
    # Save to CSV
    submission_df.to_csv(os.path.join(OUTPUT_DIR, output_filename), index=True)
    
    print(f"✅ Submission file saved: {output_filename}")
    print(f"   Shape: {submission_df.shape}")
    print(f"   Columns: {list(submission_df.columns)}")
    
    return output_filename


def load_model_results(model_path="checkpoints/model_results.pkl"):
    """
    Load trained model results and extract configuration parameters.
    
    Args:
        model_path: Path to model results pickle file
        
    Returns:
        tuple: (models, protein_mean_, protein_var_, FORCE_KEEP_GENES, significant_genes)
    """
    print(f"📂 Loading model results from: {model_path}")
    model_results = pickle.load(open(model_path, "rb"))
    
    protein_mean_ = model_results['protein_mean_']
    protein_var_ = model_results['protein_var_']
    FORCE_KEEP_GENES = model_results['target_proteins']
    significant_genes = model_results['significant_genes']
    
    print(f"✅ Model loaded: {len(model_results['models'])} protein models")
    print(f"   Target proteins: {len(FORCE_KEEP_GENES)}")
    print(f"   Significant genes: {len(significant_genes)}")
    
    return model_results['models'], protein_mean_, protein_var_, FORCE_KEEP_GENES, significant_genes


def main():
    """Main testing pipeline."""
    # Get validation data path
    val_rna_path = os.path.join(INPUT_DIR, "valid_rna.h5ad")
    
    print("🔍 Loading and preprocessing validation data...")
    
    # Load model results and extract parameters
    models, protein_mean_, protein_var_, FORCE_KEEP_GENES, significant_genes = load_model_results()
    
    # Process validation data
    val_expr_df, position_df = process_validation_data(val_rna_path, significant_genes)
    
    # Make predictions using trained models
    val_pred_protein = predict_protein_expression(
        models, 
        val_expr_df, 
        FORCE_KEEP_GENES
    )
    
    # %%
    # =============================================================================
    # Inverse Normalization and Final Output Generation
    # =============================================================================

    print("🔄 Applying inverse normalization to predictions...")
    
    # Apply inverse normalization
    print("📊 Applying inverse transformations...")
    val_pred_protein_denorm = inverse_normalize_predictions(
        val_pred_protein, 
        protein_mean_, 
        protein_var_
    )

    # Create final submission file
    submission_file = create_submission_file(val_pred_protein_denorm, position_df)

    print("\n🎯 Analysis completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()