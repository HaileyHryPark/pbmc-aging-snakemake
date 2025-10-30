import scanpy as sc
import numpy as np
import os
from pathlib import Path


# -----------------------------
# Load and Prepare
# -----------------------------
adata = sc.read_h5ad(snakemake.input[0])
dataset = snakemake.params.dataset
split_size = 30

adata.obs = adata.obs[['donor_id', 'sample_id', 'disease', 'sex', 'age', 'self_reported_ethnicity', 'nCount_RNA', 'nFeature_RNA', 'percent_mito']].copy()

sample_list = adata.obs["sample_id"].unique()

n_samples = len(sample_list)
n_splits = int(np.ceil(n_samples / split_size))

log_lines = [f"Splitting {n_samples} samples into {n_splits} splits split_size per group)\n"]

# -----------------------------
# Split and Save
# -----------------------------
for i in range(n_splits):
    start = i * split_size
    end = min((i + 1) * split_size, n_samples)
    samples_in_group = sample_list[start:end]

    adata_subset = adata[adata.obs["sample_id"].isin(samples_in_group)].copy()
    out_path = f"data/external_dis_data_prep/{dataset}_filtered_split{i+1:02d}.h5ad"
    adata_subset.write(out_path)

    log_lines.append(f"Group {i+1}: samples {start+1} to {end} ({len(samples_in_group)} samples)\n")
    log_lines.append(f"   Saved: {out_path} with {adata_subset.n_obs} cells\n")
    log_lines.append(f"   Donors: {', '.join(map(str, samples_in_group))}\n")

# -----------------------------
# Save Log File
# -----------------------------
with open(snakemake.output[0], "w") as f:
    f.writelines(log_lines)

