# -*- coding: utf-8 -*-
"""
Quality Control Module
=====================

This module handles quality control metrics calculation and data filtering.

Authors: Zijie Fang @ WangLAB in HKUST
Date: 2025-10-03
"""

import numpy as np
from scipy.sparse import issparse


def calculate_qc_metrics(data):
    """
    Calculate quality control metrics for spots/cells.
    
    Args:
        data: AnnData object
        
    Returns:
        AnnData: AnnData object with QC metrics added to obs
    """
    # 1. Calculate total UMI count (n_counts) and detected genes (n_genes)
    if issparse(data.X):
        n_counts = data.X.sum(axis=1).A1
        n_genes = (data.X > 0).sum(axis=1).A1
    else:
        n_counts = data.X.sum(axis=1)
        n_genes = (data.X > 0).sum(axis=1)
    
    data.obs['n_counts'] = n_counts
    data.obs['n_genes'] = n_genes
    
    # 2. Calculate mitochondrial gene percentage
    mito_genes = [gene for gene in data.var_names if gene.upper().startswith('MT-')]
    if len(mito_genes) > 0:
        if issparse(data.X):
            mito_counts = data[:, mito_genes].X.sum(axis=1).A1
        else:
            mito_counts = data[:, mito_genes].X.sum(axis=1)
        data.obs['pct_mito'] = mito_counts / (n_counts + 1e-8)  # Avoid division by zero
    else:
        data.obs['pct_mito'] = 0.0  # No mitochondrial genes
    
    return data


def filter_spots(data, min_counts=500, max_counts=25000, min_genes=200, max_pct_mito=0.2):
    """
    Filter spots based on quality control criteria.
    
    Args:
        data: AnnData object
        min_counts: Minimum UMI count threshold
        max_counts: Maximum UMI count threshold  
        min_genes: Minimum number of detected genes
        max_pct_mito: Maximum mitochondrial gene percentage
        
    Returns:
        AnnData: Filtered AnnData object
    """
    qc_mask = (
        (data.obs['n_counts'] >= min_counts) &
        (data.obs['n_counts'] <= max_counts) &
        (data.obs['n_genes'] >= min_genes) &
        (data.obs['pct_mito'] <= max_pct_mito) &
        (data.obs['in_tissue'] == 1)
    )
    
    return data[qc_mask].copy()
