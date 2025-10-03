# -*- coding: utf-8 -*-
"""
Utils package for protein expression prediction
==============================================

This package contains utility modules for the STP Challenge project.

Authors: Zijie Fang @ WangLAB in HKUST
Date: 2025-10-03
"""

# Make modules easily importable from the package
from .common import seed_everything, get_dataset_paths
from .gene_selection import (
    calculate_gene_statistics, 
    select_significant_genes, 
    extract_gene_expression
)
from .quality_control import calculate_qc_metrics, filter_spots
from .normalization import normalize_expression_data, preprocess_protein_data
from .data_integration import (
    load_protein_data, 
    create_protein_dataframe, 
    align_rna_protein_data
)
from .model_training import train_protein_models

__all__ = [
    # Utils
    'seed_everything',
    'get_dataset_paths',
    
    # Gene selection
    'calculate_gene_statistics',
    'select_significant_genes', 
    'extract_gene_expression',
    
    # Quality control
    'calculate_qc_metrics',
    'filter_spots',
    
    # Normalization
    'normalize_expression_data',
    'preprocess_protein_data',
    
    # Data integration
    'load_protein_data',
    'create_protein_dataframe',
    'align_rna_protein_data',
    
    # Model training
    'train_protein_models',
]
