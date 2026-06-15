import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

dfs = [pd.read_csv(f) for f in snakemake.input]

all_columns = set().union(*(df.columns for df in dfs))
common_columns = set(dfs[0].columns)
for df in dfs[1:]:
    common_columns &= set(df.columns)
incomplete_columns = all_columns - common_columns

## Outer join, incomplete columns to NaN
all_dfs = pd.concat(dfs, join="outer")
all_dfs = all_dfs[~all_dfs['rowname'].duplicated(keep=False)]

## Remove onek1k and aida
all_dfs = all_dfs[~all_dfs["dataset"].isin(["onek1k", "aida"])]
all_dfs.to_csv(snakemake.output[0], index=False)

## keep track of all collumns and incomplete columns
with open(snakemake.output[1], "w") as f:
    f.write("=== Column Summary ===\n")
    f.write(f"Total unique columns: {len(all_columns)}\n")
    f.write(f"Columns kept after inner join: {len(common_columns)}\n")
    f.write(f"Columns incomplete: {len(incomplete_columns)}\n\n")

    f.write("Kept columns:\n")
    f.write(", ".join(sorted(common_columns)) + "\n\n")
    
    f.write("Dropped columns:\n")
    f.write(", ".join(sorted(incomplete_columns)) + "\n")

## plot pca
sex_labels = all_dfs['sex']
dataset_labels = all_dfs['dataset']
X_all = all_dfs.drop(columns=['rowname', 'sex', 'age', 'disease', 'dataset', 'ethnicity'], errors='ignore').fillna(0)

scaler = StandardScaler()
scaled_data = scaler.fit_transform(X_all)

pca = PCA(n_components=2)
X_pca = pca.fit_transform(scaled_data)

sex_to_marker = {'male': 'o', 'female': '^'}
datasets = dataset_labels.unique()
colors = plt.cm.tab10.colors

for i, dataset in enumerate(datasets):
    color = colors[i % len(colors)]
    for sex_val in sex_labels.unique():
        idx = (dataset_labels == dataset) & (sex_labels == sex_val)
        plt.scatter(
            X_pca[idx, 0], X_pca[idx, 1],
            label=f'{dataset}, {sex_val}',
            alpha=0.7,
            marker=sex_to_marker.get(sex_val, 'o'),
            color=color,
            edgecolors='k',
            s=50
        )

pc1_var = pca.explained_variance_ratio_[0] * 100
pc2_var = pca.explained_variance_ratio_[1] * 100

plt.title('PCA of pseudobulk expression by dataset (color) and sex (shape)')
plt.xlabel(f'PC1 ({pc1_var:.1f}%)')
plt.ylabel(f'PC2 ({pc2_var:.1f}%)')
plt.grid(True)
plt.legend(bbox_to_anchor=(0.5, -0.1), loc='upper center', ncol = 4)
plt.tight_layout()
plt.savefig(snakemake.output[2], dpi=300)

loadings = pd.Series(pca.components_[0], index=X_all.columns)
top_pc1_genes = loadings.abs().sort_values(ascending=False).head(50)
print(top_pc1_genes)
