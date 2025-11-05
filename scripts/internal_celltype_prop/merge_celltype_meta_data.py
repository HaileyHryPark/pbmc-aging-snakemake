import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

dfs = [pd.read_csv(f) for f in snakemake.input]

all_columns = set().union(*(df.columns for df in dfs))
common_columns = set(dfs[0].columns)
for df in dfs[1:]:
    common_columns &= set(df.columns)
dropped_columns = all_columns - common_columns

all_dfs = pd.concat(dfs, join="inner")
all_dfs.to_csv(snakemake.output[0], index=False)

## keep track of all collumns and dropped columns
with open(snakemake.output[1], "w") as f:
    f.write("=== Column Summary ===\n")
    f.write(f"Total unique columns: {len(all_columns)}\n")
    f.write(f"Columns kept after inner join: {len(common_columns)}\n")
    f.write(f"Columns dropped: {len(dropped_columns)}\n\n")

    f.write("Kept columns:\n")
    f.write(", ".join(sorted(common_columns)) + "\n\n")
    
    f.write("Dropped columns:\n")
    f.write(", ".join(sorted(dropped_columns)) + "\n")

