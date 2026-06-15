import torch
import torch.nn as nn
import pandas as pd
import numpy as np
import json
import joblib
from sklearn.metrics import mean_squared_error, mean_absolute_error

# -----------------------------
# Setup
# -----------------------------
SEED = 123
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.backends.cudnn.deterministic = True

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# -----------------------------
# Model Definition (MUST MATCH TRAINING)
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
# Load external data
# -----------------------------
def load_external_data(csv_file):
    df = pd.read_csv(csv_file, dtype={"sample_id": str}).set_index("sample_id")
    df["age"] = pd.to_numeric(df["age"], errors="coerce")
    df = df.dropna(subset=["age"]).copy()

    features = [
        c for c in df.columns
        if c not in ["donor_id", "age", "sex", "dataset", "disease"]
    ]

    X = df[features].values.astype(np.float32)
    y = df["age"].values.astype(np.float32)

    return X, y, df, features


# -----------------------------
# Prediction function
# -----------------------------
def predict(model, X):
    model.eval()
    with torch.no_grad():
        X = torch.tensor(X, dtype=torch.float32).to(device)
        preds = model(X).view(-1).cpu().numpy()
    return preds


# -----------------------------
# MAIN
# -----------------------------
def main(
    external_csv,
    model_path,
    param_json_path,
    scaler_path,
    output_pred_csv,
    output_res_csv
):

    # -----------------------------
    # Load external data
    # -----------------------------
    X_ext, y_ext, df_ext, feature_cols = load_external_data(external_csv)

    # -----------------------------
    # Load LODO artifacts
    # -----------------------------
    bundle = torch.load(model_path, map_location=device)
    model_states = bundle["model_states"]
    splits = bundle["splits"]

    scalers = joblib.load(scaler_path)

    with open(param_json_path, "r") as f:
        params_list = json.load(f)

    # -----------------------------
    # storage
    # -----------------------------
    all_preds = []
    summary = []

    dataset_group = df_ext["dataset"].values.copy()

    # merge logic must MATCH training
    dataset_group = np.array(dataset_group)
    dataset_group[np.isin(dataset_group, ["Perez", "Marina"])] = "Perez_Marina"

    # -----------------------------
    # evaluate each LODO model
    # -----------------------------
    for i, held_out in enumerate(splits):

        print(f"\n=== Evaluating LODO model: {held_out} ===")

        model = SimpleMLP(
            input_dim=X_ext.shape[1],
            hidden_units=params_list[i]["hidden_units"],
            dropout=params_list[i]["dropout"],
            num_layers=params_list[i]["num_layers"]
        ).to(device)

        model.load_state_dict(model_states[i])

        scaler = scalers[i]

        # scale external data
        X_scaled = scaler.transform(X_ext)

        # predictions
        preds = predict(model, X_scaled)

        # metrics (GLOBAL external eval for this model)
        eval_datasets = ["immage", "soundlife"]
        mask = ~df_ext["dataset"].isin(eval_datasets)

        rmse = np.sqrt(mean_squared_error(y_ext[mask], preds[mask]))
        mae = mean_absolute_error(y_ext[mask], preds[mask])

        rmse_all = np.sqrt(mean_squared_error(y_ext, preds))
        mae_all = mean_absolute_error(y_ext, preds)

        summary.append({
            "held_out_train_domain": held_out,
            "rmse_external": rmse,
            "mae_external": mae,
            "rmse_external_all": rmse_all,
            "mae_external_all": mae_all,
            "n_external": len(y_ext)
        })

        # store predictions
        df_pred = pd.DataFrame({
            "sample_id": df_ext.index,
            "donor_id": df_ext["donor_id"].values,
            "dataset": df_ext["dataset"].values,
            "disease": df_ext["disease"].values,
            "sex": df_ext["sex"].values,
            "actual_age": y_ext,
            "predicted_age": preds,
            "model_trained_on": held_out
        })

        all_preds.append(df_pred)

    # -----------------------------
    # SAVE OUTPUTS
    # -----------------------------
    pd.concat(all_preds).to_csv(output_pred_csv, index=False)

    pd.DataFrame(summary).to_csv(output_res_csv, index=False)

    print("\n=== DONE ===")
    print(pd.DataFrame(summary))


# -----------------------------
# RUN
# -----------------------------
if __name__ == "__main__":
    main(
        snakemake.input[0],  # external CSV
        snakemake.input[1],  # model .pt
        snakemake.input[2],  # params json
        snakemake.input[3],  # scalers
        snakemake.output[0], # predictions
        snakemake.output[1]  # res summary
    )
