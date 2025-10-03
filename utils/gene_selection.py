# -*- coding: utf-8 -*-
"""
Gene Selection Module
====================

This module handles gene selection based on statistical criteria.

Authors: Zijie Fang @ WangLAB in HKUST
Date: 2025-10-03
"""

import numpy as np
import pandas as pd
from scipy.sparse import issparse


def calculate_gene_statistics(expression_matrix, gene_names):
    """
    Calculate mean and variance statistics for genes.
    
    Args:
        expression_matrix: Gene expression matrix (sparse or dense)
        gene_names: List of gene names
        
    Returns:
        pd.DataFrame: DataFrame with gene statistics
    """
    if issparse(expression_matrix):
        # For sparse matrices, convert to dense for variance calculation
        # Note: This might be memory intensive for large datasets
        gene_means = np.array(expression_matrix.mean(axis=0)).ravel()
        gene_vars = np.array(expression_matrix.toarray().var(axis=0)).ravel()
    else:
        gene_means = np.array(expression_matrix.mean(axis=0)).ravel()
        gene_vars = np.array(expression_matrix.var(axis=0)).ravel()
    
    return pd.DataFrame({
        'gene': gene_names,
        'mean': gene_means,
        'var': gene_vars
    })


def select_significant_genes(gene_stats, top_n=2000, force_keep_genes=None):
    """
    Select genes based on high mean and variance, plus force-keep list.
    
    Args:
        gene_stats: DataFrame with gene statistics
        top_n: Number of top genes to select by mean and variance
        force_keep_genes: List of genes to force keep
        
    Returns:
        list: List of selected significant genes
    """
    # Select genes with both high mean and high variance
    top_mean_genes = set(gene_stats.sort_values('mean', ascending=False).head(top_n)['gene'])
    top_var_genes = set(gene_stats.sort_values('var', ascending=False).head(top_n)['gene'])
    significant_genes = list(top_mean_genes & top_var_genes)
    
    # Add force-keep genes if provided
    if force_keep_genes:
        force_keep_in_stats = set(gene_stats['gene']).intersection(force_keep_genes)
        significant_genes = list(set(significant_genes).union(force_keep_in_stats))
    
    return significant_genes


def extract_gene_expression(data, significant_genes, expression_layer="cpm_log1p"):
    """
    Extract expression data for significant genes.
    
    Args:
        data: AnnData object with normalized expression layers
        significant_genes: List of significant gene names
        expression_layer: Name of expression layer to use
        
    Returns:
        pd.DataFrame: DataFrame with expression data for significant genes
    """
    # Extract spot position information
    POSITION_COLS = ['array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres']
    spot_pos = data.obs[POSITION_COLS].copy()
    spot_pos.index.name = 'spot_id'
    
    # Get expression matrix
    expr_matrix = data.layers[expression_layer]
    
    # Extract significant gene expression
    gene_indices = [data.var_names.get_loc(g) for g in significant_genes]
    
    if issparse(expr_matrix):
        expr_dense = expr_matrix.tocsr()[:, gene_indices].toarray()
    else:
        expr_dense = expr_matrix[:, gene_indices]
    
    # Create expression DataFrame
    expr_df = pd.DataFrame(
        expr_dense,
        index=data.obs_names,
        columns=significant_genes
    )
    
    # Combine position and expression data
    combined_df = pd.concat([spot_pos, expr_df], axis=1)
    
    # Ensure deterministic ordering
    combined_df = combined_df.sort_index()
    
    return combined_df
