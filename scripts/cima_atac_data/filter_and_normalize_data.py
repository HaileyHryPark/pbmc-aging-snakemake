import scanpy as sc
import episcanpy as epi
import pandas as pd
import numpy as np
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD

# -----------------------
# Load QC table
# -----------------------
qc = pd.read_csv(snakemake.input[1], index_col=0)

keep_cells = qc.query(
    "nCount_peaks > 9000 & nCount_peaks < 100000 & "
    "blacklist_ratio < 0.01 & "
    "nucleosome_signal < 4 & "
    "tss_enrichment > 4"
).index

print(f"Keeping {len(keep_cells)} cells")

# -----------------------
# Reload & subset (now in memory)
# -----------------------
adata = sc.read_h5ad(snakemake.input[0])
adata = adata[keep_cells].copy()

# -----------------------
# Feature filtering (CRITICAL)
# -----------------------
epi.pp.filter_features(adata, min_cells=50)

# -----------------------
# TF-IDF normalization
# -----------------------
epi.pp.tfidf(adata)

# -----------------------
# LSI (Signac-style)
# -----------------------
svd = TruncatedSVD(n_components=50, random_state=0)
adata.obsm["X_lsi"] = svd.fit_transform(adata.X)

# -----------------------
# Save
# -----------------------
adata.write(snakemake.output[0])

