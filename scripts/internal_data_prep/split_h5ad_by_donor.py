import scanpy as sc
import numpy as np
import os
from pathlib import Path


# -----------------------------
# Load and Prepare
# -----------------------------
adata = sc.read_h5ad(snakemake.input[0])
dataset = snakemake.params.dataset
if dataset == "onek1k":
    split_size = 100
else:
    split_size = 30

adata.obs = adata.obs[['donor_id', 'disease', 'sex', 'age', 'self_reported_ethnicity', 'nCount_RNA', 'nFeature_RNA', 'percent_mito']].copy()

donor_list = adata.obs["donor_id"].unique()

n_donors = len(donor_list)
n_splits = int(np.ceil(n_donors / split_size))

log_lines = [f"Splitting {n_donors} donors into {n_splits} splits split_size per group)\n"]

# -----------------------------
# Split and Save
# -----------------------------
for i in range(n_splits):
    start = i * split_size
    end = min((i + 1) * split_size, n_donors)
    donors_in_group = donor_list[start:end]

    adata_subset = adata[adata.obs["donor_id"].isin(donors_in_group)].copy()
    out_path = f"data/internal_data_prep/{dataset}_filtered_split{i+1:02d}.h5ad"
    adata_subset.write(out_path)

    log_lines.append(f"Group {i+1}: donors {start+1} to {end} ({len(donors_in_group)} donors)\n")
    log_lines.append(f"   Saved: {out_path} with {adata_subset.n_obs} cells\n")
    log_lines.append(f"   Donors: {', '.join(map(str, donors_in_group))}\n")

# -----------------------------
# Save Log File
# -----------------------------
with open(snakemake.output[0], "w") as f:
    f.writelines(log_lines)

