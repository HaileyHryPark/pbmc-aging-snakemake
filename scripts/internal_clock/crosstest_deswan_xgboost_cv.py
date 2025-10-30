import os
import json
import joblib
import warnings
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

warnings.filterwarnings("ignore")

# -----------------------------
# Setup
# -----------------------------
SEED = 123
np.random.seed(SEED)
print(f"Using XGBoost version: {xgb.__version__}")

# -----------------------------
# Load External Data
# -----------------------------
def load_external_data(csv_file):
    df = pd.read_csv(csv_file).set_index("rowname")
    df = df[df['sex'] != snakemake.params["gender"]]
    features = [col for col in df.columns if col not in ["rowname","age", "sex", "dataset"]]
    X = df[features].values
    y = df["age"].values
    sex = df["sex"].values
    donors = df.index.astype(str)
    return X, y, sex, donors

# -----------------------------
# Main
# -----------------------------
def main(external_csv, model_path, param_json_path, original_train_csv, prediction_csv):

    # Load external data
    X_ext, y_ext, sex_ext, donors_ext = load_external_data(external_csv)

    # Load models and params
    models = joblib.load(model_path)  # list of sklearn pipelines
    with open(param_json_path) as f:
        all_params = json.load(f)

    all_predictions = []

    for fold_idx, (pipeline, params) in enumerate(zip(models, all_params), start=1):
        print(f"\n=== Evaluating Fold {fold_idx} ===")

        # No need scaling because the pipeline already includes scaling
        # ----- Predictions -----
        preds = pipeline.predict(X_ext)

        df_pred = pd.DataFrame({
            "donor_id": donors_ext,
            "actual_age": y_ext,
            "predicted_age": preds,
            "sex": sex_ext,
            "fold": fold_idx
        })
        all_predictions.append(df_pred)

    # -----------------------------
    # Save Results
    # -----------------------------
    pd.concat(all_predictions).to_csv(prediction_csv, index=False)

# -----------------------------
# Run
# -----------------------------
if __name__ == "__main__":
    external_csv = snakemake.input[0]
    model_path = snakemake.input[1]
    param_json_path = snakemake.input[2]
    original_train_csv = snakemake.input[3]
    prediction_csv = snakemake.output[0]

    main(external_csv, model_path, param_json_path, original_train_csv, prediction_csv)

