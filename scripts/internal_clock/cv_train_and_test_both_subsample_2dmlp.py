import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import mean_squared_error, mean_absolute_error
import optuna
#import wandb
import random
import json
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from sklearn.preprocessing import StandardScaler
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
    features = [col for col in df.columns if col not in ["age", "sex","dataset"]]

    X = torch.tensor(df[features].values, dtype=torch.float32)
    y = torch.tensor(df["age"].values, dtype=torch.float32)
    sex = df["sex"].astype(str).values 
    donors = df.index.astype(str)
    dataset = df["dataset"].astype(str).values

    return X, y, sex, donors, dataset, df

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
# Main with Stratified Subsampling
# -----------------------------
def main(input_csv, output_model_dir, output_metrics, output_predictions, param_json, scaler_dir):
    X_all, y_all, sex_all, donor_ids_all, dataset_all, df = load_data(input_csv)

    subsampler = StratifiedShuffleSplit(n_splits=5, train_size=0.5, random_state=SEED)
    rmse_list, mae_list, all_results, best_params_all_folds, model_states, scalers = [], [], [], [], [], []

    for fold_idx, (sub_idx, _) in enumerate(subsampler.split(X_all, sex_all)):
        print(f"\n=== Subsample {fold_idx + 1} (Half size, stratified) ===")

        X_sub = X_all[sub_idx]
        y_sub = y_all[sub_idx]
        sex_sub = sex_all[sub_idx]
        donor_sub = donor_ids_all[sub_idx]
        dataset_sub = dataset_all[sub_idx]

        # 80/20 split of current subsample for training and final evaluation
        sub_train_idx, sub_test_idx = next(StratifiedShuffleSplit(n_splits=1, train_size=0.8, random_state=SEED).split(X_sub, sex_sub))
        
        # Split data
        X_train, y_train = X_sub[sub_train_idx], y_sub[sub_train_idx]
        X_test, y_test = X_sub[sub_test_idx], y_sub[sub_test_idx]

        # Scale features using training set only
        scaler = StandardScaler()
        X_train = torch.tensor(scaler.fit_transform(X_train), dtype=torch.float32)
        X_test = torch.tensor(scaler.transform(X_test), dtype=torch.float32)
        scalers.append(scaler)

        sex_test = sex_sub[sub_test_idx]
        donor_test = donor_sub[sub_test_idx]
        dataset_test = dataset_sub[sub_test_idx]

        # Inner CV for Hyperparameter tuning
        def objective(trial):
            hidden_units = trial.suggest_int("hidden_units", 64, 256)
            num_layers = trial.suggest_int("num_layers", 1, 3)
            dropout = trial.suggest_float("dropout", 0.01, 0.1)
            lr = trial.suggest_float("lr", 1e-4, 1e-2, log=True)
            weight_decay = trial.suggest_float("weight_decay", 1e-6, 1e-2, log=True)

            inner_cv = StratifiedShuffleSplit(n_splits=5, train_size=0.8, random_state=SEED)
            val_losses = []

            for inner_train_idx, inner_val_idx in inner_cv.split(X_sub, sex_sub):
                # Scale with training set of this inner split
                scaler_inner = StandardScaler()
                X_train_inner = torch.tensor(scaler_inner.fit_transform(X_sub[inner_train_idx]), dtype=torch.float32)
                X_val_inner = torch.tensor(scaler_inner.transform(X_sub[inner_val_idx]), dtype=torch.float32)

                y_train_inner = y_sub[inner_train_idx]
                y_val_inner = y_sub[inner_val_idx]

                model = SimpleMLP(X_all.shape[1], hidden_units, dropout, num_layers).to(device)
                optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
                criterion = nn.MSELoss()
                scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=3)

                train_loader = DataLoader(TensorDataset(X_train_inner, y_train_inner), batch_size=32, shuffle=True)
                val_loader = DataLoader(TensorDataset(X_val_inner, y_val_inner), batch_size=32)

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
        #wandb.config.update({f"fold_{fold_idx+1}_{k}": v for k, v in best_params.items()})

        model = SimpleMLP(X_all.shape[1], best_params["hidden_units"], best_params["dropout"], best_params["num_layers"]).to(device)
        optimizer = optim.Adam(model.parameters(), lr=best_params["lr"], weight_decay=best_params["weight_decay"])
        criterion = nn.MSELoss()
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=3)

        train_loader = DataLoader(TensorDataset(X_train, y_train), batch_size=32, shuffle=True)
        test_loader = DataLoader(TensorDataset(X_test, y_test), batch_size=32)

        best_model_state = None
        best_loss = float("inf")
        counter = 0
        for epoch in range(50):
            train_loss = train_one_epoch(model, train_loader, optimizer, criterion)
            val_loss, _, _ = evaluate(model, train_loader, criterion)
            scheduler.step(val_loss)

         #   wandb.log({
          #      f"fold_{fold_idx+1}/train_loss": train_loss,
           #     f"fold_{fold_idx+1}/val_loss": val_loss,
            #    f"fold_{fold_idx+1}/epoch": epoch
            #})

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

        # Evaluate on same set (since no holdout)
        test_loss, y_pred, y_true = evaluate(model, test_loader, criterion)
        rmse = np.sqrt(mean_squared_error(y_true, y_pred))
        mae = mean_absolute_error(y_true, y_pred)

        #wandb.log({
        #    f"fold_{fold_idx+1}/test_rmse": rmse,
        #    f"fold_{fold_idx+1}/test_mae": mae,
        #    f"fold_{fold_idx+1}/test_loss": test_loss
        #})

        rmse_list.append(rmse)
        mae_list.append(mae)

        print(f"Subsample {fold_idx + 1} RMSE: {rmse:.4f} | MAE: {mae:.4f}")

        df_result = pd.DataFrame({
            "donor_id": donor_test,
            "dataset": dataset_test,
            "disease": "healthy",
            "actual_age": y_true,
            "predicted_age": y_pred,
            "sex": sex_test,
            "fold": fold_idx + 1
        })
        all_results.append(df_result)

    # -----------------------------
    # Save Outputs
    # -----------------------------
    torch.save(model_states, output_model_dir)
    joblib.dump(scalers, scaler_dir)

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
