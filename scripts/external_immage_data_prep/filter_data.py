import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
import mygene

## Filter cells and remove multivisit donors

# --------------------------- #
#        CONFIGURATION        #
# --------------------------- #

# QC thresholds
MIN_COUNTS = 1000
MAX_COUNTS = 20000
MIN_FEATURES = 300
MAX_FEATURES = 6000
MAX_PCT_MT = 10

dataset = snakemake.params.dataset

# --------------------------- #
#       HELPER FUNCTIONS      #
# --------------------------- #
def filter_cells(adata):
    """Filter cells based on QC metrics."""
    print("Filtering low-quality cells...")
    before = adata.n_obs

    keep = (
        (adata.obs['nCount_RNA'] > MIN_COUNTS) &
        (adata.obs['nCount_RNA'] < MAX_COUNTS) &
        (adata.obs['nFeature_RNA'] > MIN_FEATURES) &
        (adata.obs['nFeature_RNA'] < MAX_FEATURES) &
        (adata.obs['percent_mito'] < MAX_PCT_MT)
    )

    #adata = adata[keep]
    adata = adata[keep, :].to_memory()
    print(f"Filtered cells: {before - adata.n_obs} removed, {adata.n_obs} retained.")
    return adata


# --------------------------- #
#            MAIN             #
# --------------------------- #
if __name__ == "__main__":
    print("Starting preprocessing pipeline...")

    # Load AnnData
    adata = sc.read_h5ad(snakemake.input[0], backed='r')
    print(f"Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Filter by QC
    adata = filter_cells(adata)

    # Save processed outputs
    print("Saving processed AnnData objects...")
    adata.write(snakemake.output[0])

    print("Saving metadata (obs) table...")
    adata.obs.to_csv(snakemake.output[1])

    print("Preprocessing completed")

