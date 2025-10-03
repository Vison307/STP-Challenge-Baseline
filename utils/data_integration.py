# -*- coding: utf-8 -*-
"""
Data Integration Module
======================

This module handles RNA and protein data integration and alignment.

Authors: Zijie Fang @ WangLAB in HKUST
Date: 2025-10-03
"""

import pandas as pd
import anndata as ad


def load_protein_data(data_path):
    """
    Load protein expression data from H5AD file.
    
    Args:
        data_path: Path to protein data file
        
    Returns:
        AnnData: AnnData object containing protein data
    """
    print(f"📂 Loading protein data from: {data_path}")
    pro_adata = ad.read_h5ad(data_path)
    
    print(f"   Protein data shape: {pro_adata.shape[0]} spots × {pro_adata.shape[1]} proteins")
    print("   Observation metadata columns:", list(pro_adata.obs.columns))
    print("   Variable metadata columns:", list(pro_adata.var.columns))
    
    return pro_adata


def create_protein_dataframe(pro_adata, protein_expr_df):
    """
    Create comprehensive protein dataframe with position and expression data.
    
    Args:
        pro_adata: AnnData object with protein data
        protein_expr_df: DataFrame with protein expression data
        
    Returns:
        pd.DataFrame: Combined dataframe with position and expression data
    """
    # Extract spot position information
    spot_pos_df = pro_adata.obs.copy()
    
    # Combine position and expression data
    protein_df = pd.concat([spot_pos_df, protein_expr_df], axis=1)
    
    # Ensure deterministic ordering
    protein_df = protein_df.sort_index()
    
    return protein_df


def align_rna_protein_data(rna_df, protein_df):
    """
    Align RNA and protein datasets by common spot indices.
    
    Args:
        rna_df: DataFrame with RNA expression data (including position columns)
        protein_df: DataFrame with protein expression data
        
    Returns:
        tuple: (aligned_rna_df, aligned_protein_df)
    """
    # Find common spot indices
    intersection_index = sorted(list(set(rna_df.index) & set(protein_df.index)))
    print(f"   Common spots: {len(intersection_index)} out of {len(rna_df)} RNA and {len(protein_df)} protein spots")
    
    # Align datasets
    aligned_rna_df = rna_df.loc[intersection_index].sort_index()
    aligned_protein_df = protein_df.loc[intersection_index].sort_index()
    
    # Remove position columns from RNA data (only needed for alignment)
    POSITION_COLUMNS = ['array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    aligned_rna_df = aligned_rna_df.drop(columns=POSITION_COLUMNS)
    
    return aligned_rna_df, aligned_protein_df
