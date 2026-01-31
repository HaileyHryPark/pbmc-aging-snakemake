import scanpy as sc
import pandas as pd
import scipy.sparse
import numpy as np
import mygene

dataset = snakemake.params["dataset"]

SPECIES = "human"
MITO_PREFIX = "MT-"

## Functions
def compute_qc_metrics(adata):
    """Compute nCount_RNA, nFeature_RNA, and percent_mito for an AnnData object."""
    print("Computing QC metrics...")

    X = adata.raw.X if adata.raw is not None else adata.X

    if scipy.sparse.issparse(X):
        adata.obs['nCount_RNA'] = np.ravel(X.sum(axis=1))
        adata.obs['nFeature_RNA'] = np.ravel((X > 0).sum(axis=1))
    else:
        adata.obs['nCount_RNA'] = X.sum(axis=1)
        adata.obs['nFeature_RNA'] = (X > 0).sum(axis=1)

    adata.var_names = adata.var_names.str.upper()
    mito_genes = adata.var_names.str.startswith(MITO_PREFIX)
    if mito_genes.any():
        if scipy.sparse.issparse(X):
            mito_counts = np.ravel(adata[:, mito_genes].X.sum(axis=1))
        else:
            mito_counts = adata[:, mito_genes].X.sum(axis=1)
        adata.obs['percent_mito'] = 100 * mito_counts / adata.obs['nCount_RNA']
    else:
        print("No mitochondrial genes found starting with 'MT-'")
        adata.obs['percent_mito'] = 0.0

    return adata

def map_genes(adata):
    """Map Ensembl IDs to gene symbols using mygene, with caching."""
    print("Mapping Ensembl IDs to gene symbols...")

    ensembl_ids = adata.raw.var_names.tolist() if adata.raw else adata.var_names.tolist()

    print(ensembl_ids)
    print("Querying mygene.info...")
    mg = mygene.MyGeneInfo()
    query = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species=SPECIES)
    map_df = pd.DataFrame(query)[['query', 'symbol']].dropna().drop_duplicates(subset='query')
    map_df.to_csv(snakemake.output[2], index=False)

    ens_to_symbol = dict(zip(map_df['query'], map_df['symbol']))

    mito_rename_dict = {
        'ND1': 'MT-ND1', 'ND2': 'MT-ND2', 'COX1': 'MT-CO1', 'COX2': 'MT-CO2',
        'ATP8': 'MT-ATP8', 'ATP6': 'MT-ATP6', 'COX3': 'MT-CO3', 'ND3': 'MT-ND3',
        'ND4L': 'MT-ND4L', 'ND4': 'MT-ND4', 'ND5': 'MT-ND5', 'ND6': 'MT-ND6',
        'CYTB': 'MT-CYB', 'RNR1': 'MT-RNR1', 'RNR2': 'MT-RNR2'
    }

    new_var_names = []
    for eid in ensembl_ids:
        symbol = ens_to_symbol.get(eid, eid)
        symbol = mito_rename_dict.get(symbol, symbol)
        new_var_names.append(symbol)

    adata.var_names = new_var_names
    adata.var_names_make_unique()

    return adata

def extract_age(adata):
    """Extract numeric age from metadata field."""
    print("Extracting age information from metadata...")
    if 'age' in adata.obs:
        adata.obs['age'] = adata.obs['age']
    elif 'development_stage' in adata.obs:
        adata.obs['age'] = adata.obs['development_stage'].str.extract(r'(\d+)').astype(float)
    else:
        print("No 'development_stage' column found; setting age to NaN.")
        adata.obs['age'] = np.nan
    return adata

## MAIN
if __name__ == "__main__":
    print("Starting preprocessing pipeline...")

    # Load AnnData
    adata = sc.read_h5ad(snakemake.input[0])
    print(f"Loaded AnnData: {adata.shape[0]} cells ~W {adata.shape[1]} genes")
    print(np.min(adata.X), np.max(adata.X))

    if dataset in ["perez"]:
        adata_new = sc.AnnData(
            X=adata.raw.X.copy(),  # use raw count matrix
            var=adata.raw.var.copy(),
            obs=adata.obs.copy()
        )
        adata = adata_new
        print(np.min(adata.X), np.max(adata.X))

    if dataset in ["cima_rna"]:
        adata.obs["donor_id"] = adata.obs["sample"]
        adata.obs["sample_id"] = adata.obs["sample"]
        adata.obs["disease"] = "normal"

    # Gene name mapping
    if dataset not in ["marina", "cima_rna"]:
        adata = map_genes(adata)
    else:
        pd.DataFrame(columns=['query', 'symbol']).to_csv(snakemake.output[2], index=False)

    # Compute QC metrics
    adata = compute_qc_metrics(adata)

    # Extract age metadata
    adata = extract_age(adata)

    # Save processed outputs
    print("Saving processed AnnData object...")
    adata.write(snakemake.output[0])

    print("Saving metadata (obs) table...")
    adata.obs.to_csv(snakemake.output[1])

    print("Preprocessing completed")

