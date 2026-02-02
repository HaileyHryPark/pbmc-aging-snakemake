import pandas as pd

df1 = pd.read_csv(snakemake.input[0])
df2 = pd.read_csv(snakemake.input[1])

exclude_col = "ethnicity"
cols_to_keep = df2.columns.drop(exclude_col)
print(f"Number of columns to keep: {len(cols_to_keep)}")
print("First 10 columns:", list(cols_to_keep)[:10])

# All cols exist in cima data, want to maintain the columns order so commented below
#cols_to_keep = df1.columns.intersection(cols_to_keep)
#print(f"Number of columns to keep: {len(cols_to_keep)}")
#print("First 10 columns:", list(cols_to_keep)[:10])

df_subset = df1[cols_to_keep]

df_subset.to_csv(snakemake.output[0], index=False)

