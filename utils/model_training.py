# -*- coding: utf-8 -*-
"""
Model Training Module
=====================

This module handles protein prediction model training and evaluation.

Authors: Zijie Fang @ WangLAB in HKUST
Date: 2025-10-03
"""

import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from tqdm import tqdm


def train_protein_models(rna_data, protein_data, target_proteins, test_size=0.2, random_state=42):
    """
    Train individual linear regression models for each target protein.
    
    Args:
        rna_data: DataFrame with RNA expression data
        protein_data: DataFrame with protein expression data  
        target_proteins: List of target protein names
        test_size: Fraction of data to use for testing
        random_state: Random seed for reproducibility
        
    Returns:
        dict: Dictionary containing trained models and evaluation metrics
    """
    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(
        rna_data.values, 
        protein_data[target_proteins].values, 
        test_size=test_size, 
        random_state=random_state
    )
    
    print(f"   Training set: {X_train.shape[0]} spots × {X_train.shape[1]} genes")
    print(f"   Testing set: {X_test.shape[0]} spots × {len(target_proteins)} proteins")
    
    # Train individual models for each protein
    protein_models = {}
    predictions = []
    
    print("   Training models...")
    for i, protein in tqdm(enumerate(target_proteins)):
        reg = LinearRegression()
        reg.fit(X_train, y_train[:, i])
        y_pred_i = reg.predict(X_test)
        predictions.append(y_pred_i)
        protein_models[protein] = reg
    
    predictions = np.array(predictions).T  # Shape: (n_samples, n_proteins)
    
    # Evaluate models using Spearman correlation
    print("📊 Model evaluation (Spearman correlations):")
    correlations = {}
    for i, protein in enumerate(target_proteins):
        scc, _ = spearmanr(y_test[:, i], predictions[:, i])
        correlations[protein] = scc
        print(f"   {protein:12s}: {scc:.4f}")
    
    # Calculate average correlation
    mean_scc = np.mean(list(correlations.values()))
    print(f"   {'Average':12s}: {mean_scc:.4f}")
    
    return {
        'models': protein_models,
        'correlations': correlations,
        'mean_correlation': mean_scc,
        'X_test': X_test,
        'y_test': y_test,
        'predictions': predictions,
        'target_proteins': target_proteins
    }
