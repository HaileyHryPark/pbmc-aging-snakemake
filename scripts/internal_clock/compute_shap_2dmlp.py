import os
import torch
import shap
import json
import joblib
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import StandardScaler
import torch.nn as nn

# -----------------------------
# Setup
# -----------------------------
SEED = 123
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)
torch.backends.cudnn.deterministic = True

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if device.type != "cuda":
    raise SystemExit("CUDA not available. Exiting.")
print(f"Using device: {torch.cuda.get_device_name(0)}")

# -----------------------------
# Model Definition
# -----------------------------
class SimpleMLP(nn.Module):
    def __init__(self, input_dim, hidden_units, dropout, num_layers):
        super().__init__()
        layers = []
        current_dim = input_dim
        for _ in range(num_layers):
            layers.append(nn.Linear(current_dim, hidden_units))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))
            current_dim = hidden_units
        layers.append(nn.Linear(hidden_units, 1))
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)

# -----------------------------
# Load Data
# -----------------------------
def load_internal_data(csv_file):
    df = pd.read_csv(csv_file).set_index("rowname")
    features = [col for col in df.columns if col not in ["age", "sex", "dataset"]]
    X = df[features].values.astype(np.float32)
    y = df["age"].values.astype(np.float32)
    sex = df["sex"].astype(str).values
    donors = df.index.astype(str)
    dataset = df["dataset"].astype(str).values
    return X, y, sex, donors, dataset, features

# -----------------------------
# Load Background
# -----------------------------
def load_background(background_path, all_data, features):
    bg_df = pd.read_csv(background_path)
    bg_df = all_data.loc[all_data.index.isin(bg_df["donor_id"])]
    bg_X = bg_df[features].values.astype(np.float32)
    return bg_X

# -----------------------------
# Main Function
# -----------------------------
def main(input_csv, background_csv, model_path, param_path, scaler_path, shap_output_csv, gender, fold):

    # --------------------------
    # Load data
    # --------------------------
    X, y, sex, donors, dataset, features = load_internal_data(input_csv)
    full_df = pd.read_csv(input_csv).set_index("rowname")

    # --------------------------
    # Select the right fold index
    # --------------------------
    fold_idx = int(fold) - 1  # 1-based to 0-based index

    print(f"\n=== Calculating SHAP for fold {fold} ({gender}) ===")

    # --------------------------
    # Load correct model, params, scaler
    # --------------------------
    print(f"Loading model from {model_path}")
    print(f"Loading params from {param_path}")
    print(f"Loading scaler from {scaler_path}")

    model_list = torch.load(model_path, map_location=device)
    with open(param_path) as f:
        param_list = json.load(f)
    scaler_list = joblib.load(scaler_path)

    model_state = model_list[fold_idx]
    params = param_list[fold_idx]
    scaler = scaler_list[fold_idx]

    # --------------------------
    # Scale features
    # --------------------------
    X_scaled = scaler.transform(X)
    X_tensor = torch.tensor(X_scaled, dtype=torch.float32).to(device)

    # --------------------------
    # Load background reference
    # --------------------------
    bg_X = load_background(background_csv, full_df, features)
    bg_scaled = scaler.transform(bg_X)
    background = torch.tensor(bg_scaled, dtype=torch.float32).to(device)
    print(f"Background reference: {background.shape[0]} samples")

    # --------------------------
    # Model setup
    # --------------------------
    model = SimpleMLP(
        input_dim=X.shape[1],
        hidden_units=params["hidden_units"],
        dropout=params["dropout"],
        num_layers=params["num_layers"]
    ).to(device)
    model.load_state_dict(model_state)
    model.eval()

    # --------------------------
    # SHAP calculation
    # --------------------------
    print("Computing SHAP values...")
    explainer = shap.DeepExplainer(model, background)
    shap_values = explainer.shap_values(X_tensor)

    # --------------------------
    # Save SHAP values
    # --------------------------
    shap_df = pd.DataFrame(shap_values.squeeze(-1), columns=features)
    shap_df["donor_id"] = donors
    shap_df["dataset"] = dataset
    shap_df["sex"] = sex
    shap_df["actual_age"] = y
    shap_df["fold"] = fold
    shap_df["background_n"] = background.shape[0]

    shap_df.to_csv(shap_output_csv, index=False)
    print(f"\nSHAP values for fold {fold} ({gender}) saved to: {shap_output_csv}")
    print(f"Rows: {shap_df.shape[0]} | Columns: {shap_df.shape[1]}")

# -----------------------------
# Run (Snakemake integration)
# -----------------------------
if __name__ == "__main__":
    input_csv = snakemake.input[0]            # Full dataset CSV
    background_csv = snakemake.input[1] # Background donor CSV
    model_path = snakemake.input[2]         # List of 5 model paths
    param_path = snakemake.input[3]         # List of 5 param JSONs
    scaler_path = snakemake.input[4]       # List of 5 scalers (joblib)
    shap_output_csv = snakemake.output[0]
    gender = snakemake.params["gender"]            # "male", "female", or "both"
    if gender == "both":
        fold = 1
    elif gender == "female":
        fold = 4
    elif gender == "male":
        fold = 3

    main(input_csv, background_csv, model_path, param_path, scaler_path, shap_output_csv, gender, fold)

