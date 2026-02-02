import os
import json
import joblib
import shap
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
    features = [col for col in df.columns if col not in ["rowname","age", "sex", "dataset"]]
    X = df[features].values
    y = df["age"].values
    sex = df["sex"].values
    samples = df.index.astype(str)
    dataset = df["dataset"].values
    return X, y, sex, samples, dataset

# -----------------------------
# Main
# -----------------------------
def main(external_csv, model_path, param_json_path, original_train_csv, prediction_csv, shap_csv):
    # Load training data to fit same scaler for SHAP
    train_df = pd.read_csv(original_train_csv).set_index("rowname")
    feature_cols = [col for col in train_df.columns if col not in ["age", "sex", "dataset", "disease"]]
    X_train = train_df[feature_cols].values
    
    np.random.seed(SEED)
    background_size = min(100, X_train.shape[0])
    background_idx = np.random.choice(X_train.shape[0], size=background_size, replace=False)

    # Load external data
    X_ext, y_ext, sex_ext, samples_ext, dataset_ext = load_external_data(external_csv)

    # Load models and params
    models = joblib.load(model_path)  # list of sklearn pipelines
    with open(param_json_path) as f:
        all_params = json.load(f)

    all_predictions = []
    all_shap = []

    for fold_idx, (pipeline, params) in enumerate(zip(models, all_params), start=1):
        print(f"\n=== Evaluating Fold {fold_idx} ===")

        # No need scaling because the pipeline already includes scaling
        # ----- Predictions -----
        preds = pipeline.predict(X_ext)

        df_pred = pd.DataFrame({
            "sample_id": samples_ext,
            "dataset": dataset_ext,
            "actual_age": y_ext,
            "predicted_age": preds,
            "sex": sex_ext,
            "fold": fold_idx
        })
        all_predictions.append(df_pred)

        # ----- SHAP values -----
        ## Extract fitted model and scaler
        model = pipeline.named_steps["model"]
        scaler_in_pipeline = pipeline.named_steps["scaler"]
    
        # Prepare reproducible background for SHAP
        X_train_scaled = scaler_in_pipeline.transform(X_train)
        background_data = X_train_scaled[background_idx, :]
        
        explainer = shap.TreeExplainer(model, data=background_data)
        print(f"Using {background_size} samples for SHAP background (fixed seed {SEED}).")

        #X_ext = X_ext[:, [X_train.columns.get_loc(c) for c in feature_cols]]
        X_ext_scaled = scaler_in_pipeline.transform(X_ext)
        shap_values = explainer.shap_values(X_ext_scaled)

        shap_df = pd.DataFrame(shap_values, columns=feature_cols)
        shap_df["sample_id"] = samples_ext
        shap_df["dataset"] = dataset_ext
        shap_df["actual_age"] = y_ext
        shap_df["sex"] = sex_ext
        shap_df["fold"] = fold_idx

        all_shap.append(shap_df)

    # -----------------------------
    # Save Results
    # -----------------------------
    pd.concat(all_predictions).to_csv(prediction_csv, index=False)
    pd.concat(all_shap).to_csv(shap_csv, index=False)

# -----------------------------
# Run
# -----------------------------
if __name__ == "__main__":
    external_csv = snakemake.input[0]
    model_path = snakemake.input[1]
    param_json_path = snakemake.input[2]
    original_train_csv = snakemake.input[3]
    prediction_csv = snakemake.output[0]
    shap_csv = snakemake.output[1]

    main(external_csv, model_path, param_json_path, original_train_csv, prediction_csv, shap_csv)

