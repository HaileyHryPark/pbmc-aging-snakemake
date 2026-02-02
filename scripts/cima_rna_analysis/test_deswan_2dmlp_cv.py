import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
import numpy as np
import json
import shap
import joblib
from sklearn.preprocessing import StandardScaler

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
def load_external_data(csv_file):
    df = pd.read_csv(csv_file, dtype={"rowname": str}).set_index("rowname")
    features = [col for col in df.columns if col not in ["rowname", "age", "sex", "dataset"]]
    X = df[features].values.astype(np.float32)
    y = df["age"].values.astype(str)
    sex = df["sex"].values.astype(str)
    samples = df.index.astype(str)
    dataset = df["dataset"].values.astype(str)
    return X, y, sex, samples, dataset

# -----------------------------
# Load Background
# -----------------------------
def load_background(background_path, all_data, features):
    bg_df = pd.read_csv(background_path)
    bg_df = all_data.loc[all_data.index.isin(bg_df["donor_id"])]
    bg_X = bg_df[features].values.astype(np.float32)
    return bg_X

# -----------------------------
# Evaluation Function
# -----------------------------
def evaluate(model, loader):
    model.eval()
    preds = []
    with torch.no_grad():
        for (X_batch, ) in loader:
            X_batch = X_batch.to(device)
            outputs = model(X_batch)
            preds.append(outputs.cpu().numpy())
    preds = np.concatenate(preds).squeeze()
    return preds

# -----------------------------
# Main
# -----------------------------
def main(external_csv, model_path, scaler_dir, param_json_path, original_train_csv, background_csv, prediction_csv, shap_csv):

    # Load external data
    X_ext, y_ext, sex_ext, samples_ext, dataset_ext = load_external_data(external_csv)

    # Load all model states (list of dicts)
    model_states = torch.load(model_path, map_location=device)
    scalers = joblib.load(scaler_dir)

    # Load all best parameters (list of dicts)
    with open(param_json_path) as f:
        all_params = json.load(f)

    all_predictions = []
    all_shap = []

    for fold_idx, (state_dict, params, scaler) in enumerate(zip(model_states, all_params, scalers), start=1):
        print(f"\n=== Evaluating Fold {fold_idx} ===")

        # ----- Model setup -----
        model = SimpleMLP(
            input_dim=X_ext.shape[1],
            hidden_units=params["hidden_units"],
            dropout=params["dropout"],
            num_layers=params["num_layers"]
        ).to(device)
        model.load_state_dict(state_dict)
        model.eval()

        # ----- Scale external data -----
        X_ext_scaled = scaler.transform(X_ext)
        X_ext_tensor = torch.tensor(X_ext_scaled, dtype=torch.float32)
        #y_ext_tensor = torch.tensor(y_ext, dtype=torch.float32)
        #ext_loader = DataLoader(TensorDataset(X_ext_tensor, y_ext_tensor), batch_size=32, shuffle=False)
        ext_loader = DataLoader(TensorDataset(X_ext_tensor), batch_size=32, shuffle=False)

        # ----- Predictions -----
        preds = evaluate(model, ext_loader)

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
        # Load training data to fit same scaler
        train_df = pd.read_csv(original_train_csv).set_index("rowname")
        feature_cols = [col for col in train_df.columns if col not in ["age", "sex", "dataset", "disease"]]

        bg_X = load_background(background_csv, train_df, feature_cols)
        bg_scaled = scaler.transform(bg_X)
        background = torch.tensor(bg_scaled, dtype=torch.float32).to(device)
        print(f"Background reference: {background.shape[0]} samples")

        explainer = shap.DeepExplainer(model, background)
        X_ext_tensor = X_ext_tensor.to(device).requires_grad_()
        shap_values = explainer.shap_values(X_ext_tensor.to(device))

        # Save SHAP values to CSV
        shap_df = pd.DataFrame(shap_values.squeeze(-1), columns=feature_cols)
        #shap_df = pd.DataFrame(shap_values, columns=feature_cols)
        shap_df["sample_id"] = samples_ext
        shap_df["dataset"] = dataset_ext
        shap_df["actual_age"] = y_ext 
        shap_df["sex"] = sex_ext
        shap_df["fold"] = fold_idx

        all_shap.append(shap_df)

    # ----- Save predictions -----
    pd.concat(all_predictions).to_csv(prediction_csv, index=False)
    pd.concat(all_shap).to_csv(shap_csv, index=False)

# -----------------------------
# Run
# -----------------------------
if __name__ == "__main__":
    external_csv = snakemake.input[0]
    model_path = snakemake.input[1]
    param_json_path = snakemake.input[2]
    scaler_path = snakemake.input[3]
    original_train_csv = snakemake.input[4]
    background_csv = snakemake.input[5]
    prediction_csv = snakemake.output[0]
    shap_csv = snakemake.output[1]

    main(external_csv, model_path, scaler_path, param_json_path, original_train_csv, background_csv, prediction_csv, shap_csv)
