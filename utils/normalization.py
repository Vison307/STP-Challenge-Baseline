# -*- coding: utf-8 -*-
"""
Normalization Module
===================

This module handles expression data normalization and transformation.

Authors: Zijie Fang @ WangLAB in HKUST
Date: 2025-10-03
"""

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from sklearn.preprocessing import StandardScaler


def cpm_normalization(counts_matrix):
    """
    Perform CPM (Counts Per Million) normalization.
    
    Args:
        counts_matrix: Expression count matrix (sparse or dense)
        
    Returns:
        CPM normalized matrix
    """
    if issparse(counts_matrix):
        total_counts = counts_matrix.sum(axis=1).A1
        cpm = counts_matrix.multiply(1e4 / (total_counts[:, None] + 1e-8))
    else:
        total_counts = counts_matrix.sum(axis=1)
        cpm = counts_matrix / (total_counts[:, None] + 1e-8) * 1e4
    
    return cpm


def log1p_transform(matrix):
    """
    Apply log1p transformation to expression matrix.
    
    Args:
        matrix: Expression matrix (sparse or dense)
        
    Returns:
        log1p transformed matrix
    """
    if issparse(matrix):
        log1p_matrix = matrix.copy()
        log1p_matrix.data = np.log1p(log1p_matrix.data)
    else:
        log1p_matrix = np.log1p(matrix)
    
    return log1p_matrix


def normalize_expression_data(data):
    """
    Normalize expression data using CPM + log1p.
    
    Args:
        data: AnnData object with filtered expression data
        
    Returns:
        tuple: (data object with normalized layers, cpm matrix, cpm_log1p matrix)
    """
    # Extract count matrix
    counts = data.X
    
    # 1. CPM normalization
    cpm = cpm_normalization(counts)
    
    # 2. log1p transformation
    cpm_log1p = log1p_transform(cpm)
    
    # Store normalized matrices in layers
    data.layers["cpm"] = cpm
    data.layers["cpm_log1p"] = cpm_log1p
    
    return data, cpm, cpm_log1p


def preprocess_protein_data(pro_adata, target_proteins=None, apply_zscore=True):
    """
    Preprocess protein expression data with log1p transformation and optional z-score.
    
    Args:
        pro_adata: AnnData object with protein data
        target_proteins: List of target proteins to select (if None, use all)
        apply_zscore: Whether to apply z-score standardization
        
    Returns:
        tuple: (processed DataFrame, scaler object for inverse transform)
    """
    # 1. Convert sparse matrix to dense if necessary
    if issparse(pro_adata.X):
        protein_expr = pro_adata.X.toarray()
    else:
        protein_expr = pro_adata.X
    
    # 2. Apply log1p transformation (common preprocessing for proteomics)
    protein_expr_log1p = np.log1p(protein_expr)
    
    # 3. Create DataFrame
    protein_expr_df = pd.DataFrame(
        protein_expr_log1p,
        index=pro_adata.obs_names,
        columns=pro_adata.var_names
    )
    
    # 4. Apply z-score standardization if requested
    scaler = None
    if apply_zscore:
        scaler = StandardScaler()
        protein_expr_zscore = scaler.fit_transform(protein_expr_log1p)
        protein_expr_df = pd.DataFrame(
            protein_expr_zscore,
            index=pro_adata.obs_names,
            columns=pro_adata.var_names
        )
        
        print("📊 Z-score standardization applied")
        print(f"   Mean values shape: {scaler.mean_.shape}")
        print(f"   Variance values shape: {scaler.var_.shape}")
    
    # 5. Select target proteins if specified
    if target_proteins:
        # Filter to only include proteins that exist in the data
        available_proteins = [col for col in target_proteins if col in protein_expr_df.columns]
        missing_proteins = [col for col in target_proteins if col not in protein_expr_df.columns]
        
        if missing_proteins:
            print(f"⚠️  Warning: Missing proteins in data: {missing_proteins}")
        
        protein_expr_df = protein_expr_df[available_proteins]
        print(f"✅ Selected {len(available_proteins)} target proteins")
    
    # 6. Summary of processed data
    print("📊 Summary of preprocessed protein data:")
    print(f"   Shape: {protein_expr_df.shape}")
    
    return protein_expr_df, scaler
