import pandas as pd

df1 = pd.read_csv(snakemake.input[0])
df2 = pd.read_csv(snakemake.input[1])

if "self_reported_ethnicity" in df1.columns:
    df1 = df1.rename(
        columns={"self_reported_ethnicity": "ethnicity"}
    )

if "donor_id" in df1.columns:
    df1 = df1.rename(
        columns={"donor_id": "rowname"}
    )

cols_to_keep = df2.columns

print(f"Number of columns to keep: {len(cols_to_keep)}")
print("First 10 columns:", list(cols_to_keep)[:10])

# retain only matching columns
existing_cols = [c for c in cols_to_keep if c in df1.columns]

missing_cols = [c for c in cols_to_keep if c not in df1.columns]

print(f"Existing columns: {len(existing_cols)}")
print(f"Missing columns: {len(missing_cols)}")

if len(missing_cols) > 0:
    print("Missing columns:")
    print(missing_cols)

df_subset = df1[existing_cols]

combined_df = pd.concat(
    [df_subset, df2],
    axis=0,
    ignore_index=True
)

df_subset.to_csv(snakemake.output[0], index=False)
combined_df.to_csv(snakemake.output[1], index=False)

