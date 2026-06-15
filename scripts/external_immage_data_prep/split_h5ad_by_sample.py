import scanpy as sc
import numpy as np
import gc

# -----------------------------
# CONFIG
# -----------------------------
dataset = snakemake.params.dataset
split_size = 10

# -----------------------------
# LOAD DATA
# -----------------------------
print("Loading AnnData...")
adata = sc.read_h5ad(snakemake.input[0])

print(adata)

# -----------------------------
# KEEP ONLY REQUIRED OBS COLUMNS
# (without full dataframe copy)
# -----------------------------
keep_cols = [
    'donor_id',
    'sample_id',
    'disease',
    'sex',
    'age',
    'self_reported_ethnicity',
    'nCount_RNA',
    'nFeature_RNA',
    'percent_mito'
]

existing_cols = [x for x in keep_cols if x in adata.obs.columns]

# replace obs with smaller dataframe
adata.obs = adata.obs.loc[:, existing_cols]

# -----------------------------
# OPTIONAL MEMORY REDUCTION
# -----------------------------
# remove raw if unnecessary
if adata.raw is not None:
    print("Removing raw matrix to save memory...")
    adata.raw = None

# remove layers if unnecessary
if len(adata.layers.keys()) > 0:
    print("Removing layers to save memory...")
    adata.layers.clear()

gc.collect()

# -----------------------------
# SPLIT SAMPLES
# -----------------------------
sample_list = adata.obs["sample_id"].unique()

n_samples = len(sample_list)
n_splits = int(np.ceil(n_samples / split_size))

log_lines = [
    f"Splitting {n_samples} samples into {n_splits} groups\n"
]

# -----------------------------
# SPLIT + SAVE
# -----------------------------
for i in range(n_splits):

    print(f"\nProcessing split {i+1}/{n_splits}")

    start = i * split_size
    end = min((i + 1) * split_size, n_samples)

    samples_in_group = sample_list[start:end]

    # boolean mask
    mask = adata.obs["sample_id"].isin(samples_in_group).values

    print(f"Creating view...")

    # IMPORTANT:
    # no .copy() yet
    adata_subset = adata[mask]

    out_path = (
        f"data/external_immage_data_prep/"
        f"{dataset}_filtered_split{i+1:02d}.h5ad"
    )

    print(f"Writing {out_path}")

    # write directly
    adata_subset.write(out_path)

    log_lines.append(
        f"Group {i+1}: samples {start+1} to {end}\n"
    )

    log_lines.append(
        f"Saved: {out_path} "
        f"with {adata_subset.n_obs} cells\n"
    )

    # cleanup
    del adata_subset
    gc.collect()

# -----------------------------
# SAVE LOG
# -----------------------------
with open(snakemake.output[0], "w") as f:
    f.writelines(log_lines)

print("Done.")
