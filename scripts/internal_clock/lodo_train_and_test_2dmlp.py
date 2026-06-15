import os
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error, mean_absolute_error
import random
import json
import joblib

# -----------------------------
# Reproducibility
# -----------------------------
SEED = 123
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed_all(SEED)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


# -----------------------------
# Model (FIXED ARCHITECTURE)
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
# Data loading
# -----------------------------
def load_data(csv_file,  gender="both"):
    df = pd.read_csv(csv_file).set_index("rowname")

    if gender != "both":
        df = df[df["sex"] == gender].copy()

    features = [c for c in df.columns if c not in ["age", "sex", "dataset"]]

    X = df[features].values.astype(np.float32)
    y = df["age"].values.astype(np.float32)

    sex = df["sex"].values
    dataset = df["dataset"].values
    donors = df.index.values

    return X, y, sex, dataset, donors


# -----------------------------
# Robust hyperparameter aggregation
# -----------------------------
def aggregate_hyperparams(param_list):
    df = pd.DataFrame(param_list)

    def gmean(x):
        return np.exp(np.mean(np.log(x)))

    return {
        "hidden_units": int(round(df["hidden_units"].mean())),
        "num_layers": int(round(df["num_layers"].mean())),
        "dropout": df["dropout"].mean(),
        "lr": gmean(df["lr"]),
        "weight_decay": gmean(df["weight_decay"])
    }


# -----------------------------
# Train / eval
# -----------------------------
def train_epoch(model, loader, optimizer, loss_fn):
    model.train()
    total = 0

    for Xb, yb in loader:
        Xb = Xb.to(device)
        yb = yb.to(device).view(-1)

        optimizer.zero_grad()
        pred = model(Xb).view(-1)
        loss = loss_fn(pred, yb)
        loss.backward()
        optimizer.step()

        total += loss.item() * len(Xb)

    return total / len(loader.dataset)


def evaluate(model, loader, loss_fn):
    model.eval()
    preds, trues = [], []
    total = 0

    with torch.no_grad():
        for Xb, yb in loader:
            Xb = Xb.to(device)
            yb = yb.to(device).view(-1)

            pred = model(Xb).view(-1)
            loss = loss_fn(pred, yb)

            total += loss.item() * len(Xb)

            preds.append(pred.cpu().numpy().reshape(-1))
            trues.append(yb.cpu().numpy().reshape(-1))

    preds = np.concatenate(preds)
    trues = np.concatenate(trues)

    return total / len(loader.dataset), preds, trues


# -----------------------------
# Single LODO training run
# -----------------------------
def fit_model(X_train, y_train, X_test, y_test, params):

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    X_train = torch.tensor(X_train, dtype=torch.float32)
    X_test = torch.tensor(X_test, dtype=torch.float32)
    y_train = torch.tensor(y_train, dtype=torch.float32)
    y_test = torch.tensor(y_test, dtype=torch.float32)

    train_loader = DataLoader(
        TensorDataset(X_train, y_train),
        batch_size=32,
        shuffle=True,
        generator=torch.Generator().manual_seed(SEED)
    )

    test_loader = DataLoader(
        TensorDataset(X_test, y_test),
        batch_size=32,
        shuffle=False
    )

    model = SimpleMLP(
        X_train.shape[1],
        params["hidden_units"],
        params["dropout"],
        params["num_layers"]
    ).to(device)

    optimizer = optim.Adam(
        model.parameters(),
        lr=params["lr"],
        weight_decay=params["weight_decay"]
    )

    loss_fn = nn.MSELoss()

    best_state = None
    best_loss = np.inf
    patience = 10
    counter = 0

    for epoch in range(100):
        train_epoch(model, train_loader, optimizer, loss_fn)
        val_loss, _, _ = evaluate(model, train_loader, loss_fn)

        if val_loss < best_loss:
            best_loss = val_loss
            best_state = model.state_dict()
            counter = 0
        else:
            counter += 1
            if counter >= patience:
                break

    model.load_state_dict(best_state)

    test_loss, preds, trues = evaluate(model, test_loader, loss_fn)

    rmse = np.sqrt(mean_squared_error(trues, preds))
    mae = mean_absolute_error(trues, preds)

    return rmse, mae, preds, trues, scaler, model


# -----------------------------
# LODO evaluation
# -----------------------------
def run_lodo(X, y, sex, dataset, donors, params):

    df_results = []

    dataset_group = np.array(dataset)
    dataset_group[np.isin(dataset_group, ["perez", "marina"])] = "perez_marina"

    splits = ["onek1k", "aida", "perez_marina"]

    model_states = []
    scalers = []
    params_log = []
    summary = []

    for held_out in splits:

        test_idx = dataset_group == held_out
        train_idx = ~test_idx

        X_train, y_train = X[train_idx], y[train_idx]
        X_test, y_test = X[test_idx], y[test_idx]

        rmse, mae, preds, trues, scaler, model = fit_model(
            X_train, y_train,
            X_test, y_test,
            params
        )

        print(f"{held_out} | RMSE={rmse:.4f} | MAE={mae:.4f}")

        model_states.append(model.state_dict())
        scalers.append(scaler)

        params_log.append({
            "held_out": held_out,
            **params
        })

        summary.append({
            "held_out": held_out,
            "rmse": rmse,
            "mae": mae,
            "n_train": len(X_train),
            "n_test": len(X_test),
            "train_age": np.mean(y_train),
            "test_age": np.mean(y_test),
            "train_male_pct": np.mean(sex[train_idx] == "male"),
            "test_male_pct": np.mean(sex[test_idx] == "male")
        })

        df_results.append(pd.DataFrame({
            "donor_id": donors[test_idx],
            "dataset": dataset[test_idx],
            "actual_age": y_test,
            "predicted_age": preds,
            "sex": sex[test_idx],
            "held_out": held_out
        }))

    return (
        pd.concat(df_results),
        pd.DataFrame(summary),
        model_states,
        scalers,
        params_log
    )


# -----------------------------
# MAIN
# -----------------------------
def main(input_csv, param_json, output_pred, output_metrics, output_model, params_out, scaler_out):

    gender = snakemake.params["gender"]

    X, y, sex, dataset, donors = load_data(input_csv, gender)

    # load CV-tuned hyperparameters
    with open(param_json, "r") as f:
        param_list = json.load(f)

    params = aggregate_hyperparams(param_list)

    print("\n=== Aggregated Hyperparameters ===")
    print(params)

    preds, summary, models, scalers, params_log = run_lodo(
        X, y, sex, dataset, donors,
        params
    )

    torch.save({
        "model_states": models,
        "architecture": "SimpleMLP",
        "params": params,
        "splits": ["onek1k", "aida", "perez_marina"]
    }, output_model)

    joblib.dump(scalers, scaler_out)

    preds.to_csv(output_pred, index=False)
    summary.to_csv(output_metrics, index=False)

    with open(params_out, "w") as f:
        json.dump(params_log, f, indent=2)


if __name__ == "__main__":
    main(
        snakemake.input[0],
        snakemake.input[1],
        snakemake.output[0],
        snakemake.output[1],
        snakemake.output[2],
        snakemake.output[3],
        snakemake.output[4]
    )

