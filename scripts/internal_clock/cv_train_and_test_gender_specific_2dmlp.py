import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.model_selection import RepeatedKFold
import optuna
import json
import random
import joblib

# -----------------------------
# Setup & Reproducibility
# -----------------------------
SEED = 123
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed(SEED)
torch.backends.cudnn.deterministic = True

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
if device.type != "cuda":
    raise SystemExit("CUDA not available. Exiting.")
print(f"CUDA is available. Using GPU: {torch.cuda.get_device_name(0)}")

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
# Data Loading
# -----------------------------
def load_data(csv_file):
    df = pd.read_csv(csv_file).set_index("rowname")
    df = df[df['sex'] == snakemake.params["gender"]]
    features = [col for col in df.columns if col not in ["age", "sex","dataset"]]
    X = torch.tensor(df[features].values, dtype=torch.float32)
    y = torch.tensor(df["age"].values, dtype=torch.float32)
    sex = df["sex"].astype(str).values  
    donors = df.index.astype(str)
    dataset = df["dataset"].astype(str).values
    return X, y, sex, donors, dataset

# -----------------------------
# Training & Evaluation
# -----------------------------
def train_one_epoch(model, loader, optimizer, criterion):
    model.train()
    total_loss = 0
    for X_batch, y_batch in loader:
        X_batch, y_batch = X_batch.to(device), y_batch.to(device)
        optimizer.zero_grad()
        preds = model(X_batch)
        loss = criterion(preds.squeeze(1), y_batch)
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * len(X_batch)
    return total_loss / len(loader.dataset)

def evaluate(model, loader, criterion):
    model.eval()
    preds, labels = [], []
    total_loss = 0
    with torch.no_grad():
        for X_batch, y_batch in loader:
            X_batch, y_batch = X_batch.to(device), y_batch.to(device)
            outputs = model(X_batch)
            loss = criterion(outputs.squeeze(1), y_batch)
            total_loss += loss.item() * len(X_batch)
            preds.append(outputs.cpu().numpy())
            labels.append(y_batch.cpu().numpy())
    preds = np.concatenate(preds).squeeze()
    labels = np.concatenate(labels)
    return total_loss / len(loader.dataset), preds, labels

# -----------------------------
# Main with Optuna
# -----------------------------
def main(input_csv, output_model_dir, output_metrics, output_predictions, param_json, scalers_dir):
    X_all, y_all, sex_all, donor_ids_all, dataset_all = load_data(input_csv)

    outer_cv = RepeatedKFold(n_splits=5, n_repeats=1, random_state=SEED)
    rmse_list = []
    mae_list = []
    all_results = []
    best_params_all_folds = []
    model_states = []
    scalers = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X_all)):
        print(f"\n=== Outer Fold {fold_idx + 1} ===")

        X_train_outer, y_train_outer = X_all[train_idx], y_all[train_idx]
        X_test_outer, y_test_outer = X_all[test_idx], y_all[test_idx]
        donors_test = donor_ids_all[test_idx]
        dataset_test = dataset_all[test_idx]
        sex_test = sex_all[test_idx]

        scaler = StandardScaler()
        X_train_outer = torch.tensor(scaler.fit_transform(X_train_outer), dtype=torch.float32)
        X_test_outer = torch.tensor(scaler.transform(X_test_outer), dtype=torch.float32)
        scalers.append(scaler)

        # Inner CV hyperparameter tuning
        def objective(trial):
            hidden_units = trial.suggest_int("hidden_units", 64, 256)
            num_layers = trial.suggest_int("num_layers", 1, 3)
            dropout = trial.suggest_float("dropout", 0.01, 0.1)
            lr = trial.suggest_float("lr", 1e-4, 1e-2, log=True)
            weight_decay = trial.suggest_float("weight_decay", 1e-6, 1e-2, log=True)

            inner_cv = RepeatedKFold(n_splits=5, n_repeats=1, random_state=SEED)
            val_losses = []

            for inner_train_idx, inner_val_idx in inner_cv.split(X_train_outer, sex_all[train_idx]):
                X_inner_train_raw = X_train_outer[inner_train_idx].numpy()
                X_inner_val_raw = X_train_outer[inner_val_idx].numpy()

                # Fit scaler on inner training data
                scaler_inner = StandardScaler()
                X_inner_train = torch.tensor(scaler_inner.fit_transform(X_inner_train_raw), dtype=torch.float32)
                X_inner_val = torch.tensor(scaler_inner.transform(X_inner_val_raw), dtype=torch.float32)

                y_inner_train = y_train_outer[inner_train_idx]
                y_inner_val = y_train_outer[inner_val_idx]

                train_loader = DataLoader(TensorDataset(X_inner_train, y_inner_train), batch_size=32, shuffle=True)
                val_loader = DataLoader(TensorDataset(X_inner_val, y_inner_val), batch_size=32)

                model = SimpleMLP(X_all.shape[1], hidden_units, dropout, num_layers).to(device)
                optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
                scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=3)
                criterion = nn.MSELoss()

                best_loss = float("inf")
                counter = 0
                for epoch in range(50):
                    train_one_epoch(model, train_loader, optimizer, criterion)
                    val_loss, _, _ = evaluate(model, val_loader, criterion)
                    scheduler.step(val_loss)
                    if val_loss < best_loss:
                        best_loss = val_loss
                        counter = 0
                    else:
                        counter += 1
                        if counter >= 5:
                            break
                val_losses.append(best_loss)
            return np.mean(val_losses)

        study = optuna.create_study(direction="minimize")
        study.optimize(objective, n_trials=25)
        best_params = study.best_params
        best_params_all_folds.append(best_params)

        model = SimpleMLP(X_all.shape[1], best_params["hidden_units"], best_params["dropout"], best_params["num_layers"]).to(device)
        optimizer = optim.Adam(model.parameters(), lr=best_params["lr"], weight_decay=best_params["weight_decay"])
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=3)
        criterion = nn.MSELoss()

        train_loader = DataLoader(TensorDataset(X_train_outer, y_train_outer), batch_size=32, shuffle=True)
        test_loader = DataLoader(TensorDataset(X_test_outer, y_test_outer), batch_size=32)

        best_model_state = None
        best_loss = float("inf")
        counter = 0

        for epoch in range(50):
            train_loss = train_one_epoch(model, train_loader, optimizer, criterion)
            val_loss, _, _ = evaluate(model, train_loader, criterion)
            scheduler.step(val_loss)

            if val_loss < best_loss:
                best_loss = val_loss
                best_model_state = model.state_dict()
                counter = 0
            else:
                counter += 1
                if counter >= 5:
                    break

        model.load_state_dict(best_model_state)
        model_states.append(model.state_dict())

        test_loss, y_pred, y_true = evaluate(model, test_loader, criterion)
        rmse = np.sqrt(mean_squared_error(y_true, y_pred))
        mae = mean_absolute_error(y_true, y_pred)

        rmse_list.append(rmse)
        mae_list.append(mae)

        print(f"Fold {fold_idx + 1} RMSE: {rmse:.4f} | MAE: {mae:.4f}")

        df = pd.DataFrame({
            "donor_id": donors_test,
            "dataset": dataset_test,
            "disease": "healthy",
            "actual_age": y_true,
            "predicted_age": y_pred,
            "sex": sex_test,
            "fold": fold_idx + 1
        })
        all_results.append(df)

    # -----------------------------
    # Save Outputs
    # -----------------------------
    torch.save(model_states, output_model_dir)
    joblib.dump(scalers, scalers_dir)

    pd.concat(all_results).to_csv(output_predictions, index=False)

    results_df = pd.DataFrame({
        "fold": range(1, len(rmse_list) + 1),
        "RMSE": rmse_list,
        "MAE": mae_list
    })
    results_df.to_csv(output_metrics, index=False)

    with open(param_json, "w") as f:
        json.dump(best_params_all_folds, f, indent=2)

if __name__ == "__main__":
    main(snakemake.input[0], snakemake.output[0], snakemake.output[1], snakemake.output[2], snakemake.output[3], snakemake.output[4])

