import os
import json
import joblib
import warnings
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error

warnings.filterwarnings("ignore")

# -----------------------------
# Load External Data
# -----------------------------
def load_external_data(csv_file):
    df = pd.read_csv(csv_file)
    features = [col for col in df.columns if col not in ["rowname","age", "sex", "dataset"]]
    samples = df["rowname"].values.astype(str)
    X = df[features].values
    y = df["age"].values
    sex = df["sex"].values
    dataset = df["dataset"].values
    return X, y, sex, samples, dataset 

# -----------------------------
# Main
# -----------------------------
def main(external_csv, model_paths, prediction_csv):

    # Load external data
    X_ext, y_ext, sex_ext, samples_ext, dataset_ext = load_external_data(external_csv)

    model_list = joblib.load(model_paths)

    all_preds = []

    # -----------------------------
    # Evaluate Each Fold Model
    # -----------------------------
    for fold, model in enumerate(model_list, 1):
        print(f"=== Evaluating Fold {fold} ===")

        preds = model.predict(X_ext)

        df_preds = pd.DataFrame({
            "sample_id": samples_ext,
            "dataset": dataset_ext,
            "actual_age": y_ext,
            "predicted_age": preds,
            "sex": sex_ext,
            "fold": fold
        })
        all_preds.append(df_preds)

    # -----------------------------
    # Combine Results
    # -----------------------------
    preds_df = pd.concat(all_preds)
    preds_df.to_csv(prediction_csv, index=False)

    # -----------------------------
    # Write Summary
    # -----------------------------
    print("Predictions saved:", prediction_csv)

# -----------------------------
# Run
# -----------------------------
if __name__ == "__main__":
    external_csv = snakemake.input[0]
    model_paths = snakemake.input[1]  # you can pass a list of joblib models if multiple folds

    prediction_csv = snakemake.output[0]

    main(external_csv, model_paths, prediction_csv)

