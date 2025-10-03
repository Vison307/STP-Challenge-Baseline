# -*- coding: utf-8 -*-
"""
Linear Regression Model for Protein Expression Prediction from RNA Data
=====================================================================

This script implements a baseline for the STP Challenge.

Authors: Zijie Fang @ WangLAB in HKUST
Date: 2025-10-03
"""

import os
import pickle
import scanpy as sc

# Configuration constants
RANDOM_SEED = 42
DATASET_DIR = "./input"

# FORCE_KEEP_GENES - genes that must be included (proteins we want to predict)
FORCE_KEEP_GENES = [
    "synd", "FOXP3", "CD16", "CD31", "CXCL13", "Ki67", "OLIG2", "CXCR5", 
    "HLA-A", "PD-L1", "PSD95", "CD20", "CD68", "CD44", "SMA", "MSH6", 
    "CD23", "GFAP", "SYNA", "Podoplanin", "Vimentin", "CD47", "CD74", 
    "SIRP", "Granzyme B", "IDH1", "MPO", "CD45", "CD21", "FIBR", "C-KIT", 
    "CD3e", "TOX", "PD-1", "PDGFR", "CD4", "MAP2", "CD8", "MGMT", "CD38", 
    "HLA-DR", "CD14", "ICOS", "Granzyme K"
]

# QC_THRESHOLDS - quality control filtering parameters
QC_THRESHOLDS = {
    'min_counts': 500,
    'max_counts': 25000,
    'min_genes': 200,
    'max_pct_mito': 0.2
}

# Import our refactored modules
from utils.common import seed_everything, get_dataset_paths
from utils.gene_selection import (
    calculate_gene_statistics, select_significant_genes, extract_gene_expression
)
from utils.quality_control import calculate_qc_metrics, filter_spots
from utils.normalization import normalize_expression_data, preprocess_protein_data
from utils.data_integration import (
    load_protein_data, create_protein_dataframe, align_rna_protein_data
)
from utils.model_training import train_protein_models


def main():
    """Main training pipeline."""
    
    # Set random seed for reproducibility
    seed_everything(RANDOM_SEED)
    
    # %%
    # =============================================================================
    # Data Loading and Configuration
    # =============================================================================
    
    print("📂 Loading training RNA dataset...")
    train_rna_path, train_pro_path = get_dataset_paths(DATASET_DIR)
    train_rna = sc.read_h5ad(train_rna_path)
    print(f"   Loaded RNA data: {train_rna.shape[0]} spots, {train_rna.shape[1]} genes")
    
    # %%
    # =============================================================================
    # Gene Selection: High Mean and Variance Genes
    # =============================================================================
    
    print("🧬 Selecting significant genes based on mean and variance...")
    
    # Calculate gene statistics
    gene_stats = calculate_gene_statistics(train_rna.X, train_rna.var_names)
    
    # Select significant genes
    significant_genes = select_significant_genes(
        gene_stats, 
        top_n=2000, 
        force_keep_genes=FORCE_KEEP_GENES
    )
    
    print(f"✅ Selected {len(significant_genes)} significant genes")
    
    # %%
    # =============================================================================
    # Quality Control and Data Filtering
    # =============================================================================
    
    print("🔍 Performing quality control filtering...")
    
    # Calculate QC metrics
    train_rna = calculate_qc_metrics(train_rna)
    
    # Apply filtering
    train_rna_qc = filter_spots(train_rna, **QC_THRESHOLDS)
    
    print(f"✅ After QC filtering: {train_rna_qc.n_obs} spots remained (originally {train_rna.n_obs})")
    
    # %%
    # =============================================================================
    # Data Normalization: CPM + log1p Transformation
    # ================================================================================
    
    print("🔄 Performing CPM normalization and log1p transformation...")
    
    # Apply normalization to QC-filtered data
    train_rna_qc, cpm_matrix, cpm_log1p_matrix = normalize_expression_data(train_rna_qc)
    
    print("✅ Normalization completed - CPM and log1p matrices stored in data.layers")
    
    # %%
    # =============================================================================
    # Extract Significant Gene Expression Data
    # =============================================================================
    
    print("🧬 Extracting significant gene expression data...")
    
    # Extract significant gene expression
    rna_df = extract_gene_expression(train_rna_qc, significant_genes)
    print(f"✅ RNA expression data extracted: {rna_df.shape[0]} spots × {len(significant_genes)} genes")
    
    # %%
    # =============================================================================
    # Protein Data Loading and Preprocessing
    # =============================================================================
    
    print("🧬 Loading and preprocessing protein expression data...")
    
    # Load protein data
    pro_adata = load_protein_data(train_pro_path)
    
    # Preprocess protein data with z-score standardization
    protein_expr_df, protein_scaler = preprocess_protein_data(
        pro_adata, 
        target_proteins=FORCE_KEEP_GENES,
        apply_zscore=True
    )
    
    # Store scaler parameters for later inverse transformation
    protein_mean_ = protein_scaler.mean_
    protein_var_ = protein_scaler.var_
    
    # %%
    # =============================================================================
    # Data Integration: Combine RNA and Protein Data
    # =============================================================================
    
    print("🔗 Integrating RNA and protein datasets...")
    
    # Create protein dataframe
    protein_df = create_protein_dataframe(pro_adata, protein_expr_df)
    
    # Align RNA and protein datasets
    rna_df_aligned, protein_df_aligned = align_rna_protein_data(rna_df, protein_df)
    
    print(f"✅ Data integration completed:")
    print(f"   RNA data shape: {rna_df_aligned.shape}")
    print(f"   Protein data shape: {protein_df_aligned.shape}")
    print(f"   RNA features: {list(rna_df_aligned.columns[:5])}...")
    print(f"   Protein features: {list(protein_df_aligned.columns[:5])}...")
    
    # %%
    # =============================================================================
    # Model Training: Linear Regression for Each Protein
    # =============================================================================
    
    print("🤖 Training linear regression models for protein prediction...")
    
    # Train models using aligned data
    model_results = train_protein_models(
        rna_df_aligned, 
        protein_df_aligned, 
        target_proteins=FORCE_KEEP_GENES,
        test_size=0.2, 
        random_state=RANDOM_SEED
    )
    
    print(f"✅ Model training completed - Average Spearman correlation: {model_results['mean_correlation']:.4f}")
    
    # Save results
    if not os.path.exists("checkpoints"):
        os.makedirs("checkpoints")
    
    model_results['protein_mean_'] = protein_mean_
    model_results['protein_var_'] = protein_var_
    model_results['significant_genes'] = significant_genes
    
    pickle.dump(model_results, open("checkpoints/model_results.pkl", "wb"))


if __name__ == "__main__":
    main()