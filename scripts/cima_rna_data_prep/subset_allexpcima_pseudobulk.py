import pandas as pd

df = pd.read_csv(snakemake.input[0])

non_feature_cols = ["rowname", "age", "sex", "dataset"]

feature_cols = df.columns.difference(non_feature_cols)
non_zero_features = df[feature_cols].loc[:, (df[feature_cols] != 0).all()]

filtered_df = pd.concat([df[non_feature_cols], non_zero_features], axis=1)

filtered_df.to_csv(snakemake.output[0], index=False)
