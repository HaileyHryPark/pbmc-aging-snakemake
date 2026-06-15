import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedShuffleSplit

SEED = 123

# -----------------------------
# Load data
# -----------------------------
df = pd.read_csv(snakemake.input[0]).set_index("rowname")

sex_all = df["sex"].astype(str).values

all_splits = []

# -----------------------------
# Outer subsampling
# -----------------------------
subsampler = StratifiedShuffleSplit(
    n_splits=5,
    train_size=0.5,
    random_state=SEED
)

for fold_idx, (sub_idx, _) in enumerate(
    subsampler.split(df, sex_all),
    start=1
):

    df_sub = df.iloc[sub_idx].copy()

    sex_sub = df_sub["sex"].astype(str).values

    # -------------------------
    # 80/20 split
    # -------------------------
    train_idx, test_idx = next(
        StratifiedShuffleSplit(
            n_splits=1,
            train_size=0.8,
            random_state=SEED
        ).split(df_sub, sex_sub)
    )

    df_train = df_sub.iloc[train_idx].copy()
    df_test = df_sub.iloc[test_idx].copy()

    df_train["outer_fold"] = fold_idx
    df_train["split"] = "train"

    df_test["outer_fold"] = fold_idx
    df_test["split"] = "test"

    all_splits.append(df_train)
    all_splits.append(df_test)

# -----------------------------
# Save
# -----------------------------
split_df = pd.concat(all_splits)
split_df = split_df.reset_index()
split_df = split_df[
    [
        "rowname",
        "age",
        "sex",
        "dataset",
        "outer_fold",
        "split"
    ]
]

split_df.to_csv(
    snakemake.output[0],
    index=False
)
