# -*- coding: utf-8 -*-
"""
Common utilities for protein expression prediction
================================================

This module contains utility functions used across the project.

Authors: Zijie Fang @ WangLAB in HKUST
Date: 2025-10-03
"""

import os
import random
import numpy as np
import scanpy as sc


def seed_everything(seed=42):
    """
    Set random seeds for reproducible results across multiple libraries.
    
    Args:
        seed (int): Random seed value
    """
    # Python and NumPy seeds
    np.random.seed(seed)
    random.seed(seed)
    
    # Set environment variable
    os.environ['PYTHONHASHSEED'] = str(seed)
    
    # Scanpy (if used)
    if hasattr(sc._settings, 'random_state'):
        sc._settings.random_state = seed
    
    print(f"🔬 Deterministic execution enabled with seed: {seed}")


def get_dataset_paths(dataset_dir="./input"):
    """
    Get standardized dataset paths.
    
    Args:
        dataset_dir: Base directory containing dataset files
        
    Returns:
        tuple: (train_rna_path, train_pro_path)
    """
    train_rna_path = os.path.join(dataset_dir, "train_rna.h5ad")
    train_pro_path = os.path.join(dataset_dir, "train_pro.h5ad")
    return train_rna_path, train_pro_path
