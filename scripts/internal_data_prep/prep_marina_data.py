import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse
import mygene

# Load AnnData
adata = sc.read_h5ad(snakemake.input[0])
print(f"Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} genes")

meta = pd.read_csv(snakemake.input[1], index_col=0)
meta = meta[["Donor_id","Age_group","Sex","Age","Tube_id","Batch","File_name"]]
meta.columns = meta.columns.str.lower()
meta.rename(columns={'tube_id': 'sample_id'}, inplace=True)
meta['sex'] = meta['sex'].str.lower()
meta['self_reported_ethnicity'] = "Caucasian"
meta['disease'] = "normal"

print(adata.obs.head())
print(meta.head())

print("adata shape:", adata.shape)
print("metadata shape:", meta.shape)

print("adata.obs index example:", adata.obs.index[:5])
print("metadata index example:", meta.index[:5])

meta = meta.loc[adata.obs.index]
print(meta.head())
adata.obs = meta

print("Saving AnnData object...")
adata.write(snakemake.output[0])


