import scanpy as sc
import episcanpy as epi
import numpy as np
import pandas as pd
import scipy.sparse as sp

# -----------------------
# Load in backed mode
# -----------------------
adata = sc.read_h5ad(snakemake.input[0])

print(adata)
print(adata.obs.columns)
print(adata.var.columns)

# -----------------------
# Keep only necessary metadata
# -----------------------
adata.obs = adata.obs[
    [
        "sample",
        "sex",
        "age",
        "menstruating"
    ]
]

# -----------------------
# QC metrics (SAFE)
# -----------------------

# nCount_peaks
adata.obs["nCount_peaks"] = np.array(adata.X.sum(axis=1)).ravel()
adata.obs["nFeature_peaks"] = np.diff(adata.X.indptr)

# % reads in peaks not available in data (because not using cellranger)
adata.obs["pct_reads_in_peaks"] = np.nan

# -----------------------
# Episcanpy QC
# (run once, cache results)
# -----------------------

epi.pp.tss_enrichment(adata)
epi.pp.nucleosome_signal(adata)
epi.pp.blacklist_ratio(adata)

# -----------------------
# Save QC metrics ONLY
# -----------------------
adata.obs.to_csv(snakemake.output[0])

