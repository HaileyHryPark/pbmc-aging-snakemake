import os
import json
import joblib
import warnings
import pandas as pd
import numpy as np
import xgboost as xgb
import optuna
import wandb

from sklearn.model_selection import KFold, RepeatedKFold, StratifiedShuffleSplit, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")

# -----------------------------
# Check GPU Availability
# -----------------------------
print(xgb.config_context())

# Check XGBoost build info
build_info = xgb.build_info()
print(build_info)

# Check if GPU support is enabled
if 'USE_CUDA' in build_info and build_info['USE_CUDA'] == True:
    print("XGBoost is built with GPU support")
else:
    raise SystemExit("XGBoost is not built with GPU support")

# Subsample settings
n_subsamples = 5
sample_fraction = 0.5
random_state = 123

# Load data
data = pd.read_csv(snakemake.input[0], index_col=0)
X_full = data.drop(columns=['age', 'sex', "dataset"]).apply(pd.to_numeric)
y_full = data['age']
gender_full = data['sex']
dataset_full = data['dataset']
donor_ids_full = data.index

# Stratified subsampling by gender
sss = StratifiedShuffleSplit(n_splits=n_subsamples, test_size=sample_fraction, random_state=random_state)

# Output placeholders
all_metrics = []
all_predictions = []
all_models = []
all_params = []

# Loop over stratified subsamples (used as outer folds)
for fold, (subsample_idx, _) in enumerate(sss.split(X_full, gender_full), 1):
    print(f"=== Subsample {fold} ===")

    X_sub = X_full.iloc[subsample_idx]
    y_sub = y_full.iloc[subsample_idx]
    gender_sub = gender_full.iloc[subsample_idx]
    dataset_sub = dataset_full.iloc[subsample_idx]
    donor_ids_sub = donor_ids_full[subsample_idx]

    # Inner CV
    inner_cv = RepeatedKFold(n_splits=5, n_repeats=1, random_state=123)

    # -----------------------------
    # Optuna Inner CV for Subsample
    # -----------------------------
    def objective(trial):
        params = {
            'objective': 'reg:squarederror',
            'eval_metric': 'mae',
            'lambda': trial.suggest_float('lambda', 1e-8, 10.0, log=True),
            'alpha': trial.suggest_float('alpha', 1e-8, 10.0, log=True),
            'colsample_bytree': trial.suggest_float('colsample_bytree', 0.3, 1.0),
            'subsample': trial.suggest_float('subsample', 0.5, 1.0),
            'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
            'max_depth': trial.suggest_int('max_depth', 3, 9),
            'min_child_weight': trial.suggest_int('min_child_weight', 1, 300),
            'tree_method': 'hist',
            'device': 'cuda',
            'random_state': 123,
            "n_jobs": 1,
            "verbosity": 0,
        }

        val_maes = []
        for tr_idx, va_idx in inner_cv.split(X_sub, y_sub):
            X_tr, X_va = X_sub.iloc[tr_idx], X_sub.iloc[va_idx]
            y_tr, y_va = y_sub.iloc[tr_idx], y_sub.iloc[va_idx]

            scaler = StandardScaler()
            X_tr_s = scaler.fit_transform(X_tr)
            X_va_s = scaler.transform(X_va)

            dtrain = xgb.DMatrix(X_tr_s, label=y_tr)
            dval = xgb.DMatrix(X_va_s, label=y_va)

            bst = xgb.train(
                params, dtrain,
                num_boost_round=500,
                evals=[(dtrain, "train"), (dval, "eval")],
                early_stopping_rounds=25,
                verbose_eval=False
            )

            preds = bst.predict(dval, iteration_range=(0, bst.best_iteration))
            val_maes.append(mean_absolute_error(y_va, preds))

        mean_mae = np.mean(val_maes)
        return mean_mae

    study = optuna.create_study(direction="minimize")
    study.optimize(objective, n_trials=25)

    best_params = study.best_params
    best_params.update({
        "tree_method": "gpu_hist",
        "predictor": "gpu_predictor",
        "objective": "reg:squarederror",
        "eval_metric": "mae",
        "n_jobs": 1,
        "random_state": 123
    })


    # -----------------------------
    # Train-Test Split (80/20) within Subsample
    # -----------------------------
    X_train, X_test, y_train, y_test, gender_train, gender_test, donor_train, donor_test, dataset_train, dataset_test = train_test_split(
        X_sub, y_sub, gender_sub, donor_ids_sub, dataset_sub,
        test_size=0.2,
        random_state=fold * 42,
        stratify=gender_sub
    )

    # -----------------------------
    # Train Final Model on Subsample
    # -----------------------------
    model = xgb.XGBRegressor(**best_params)
    pipeline = Pipeline([
        ("scaler", StandardScaler()),
        ("model", model)
    ])

    pipeline.fit(X_train, y_train)
    preds = pipeline.predict(X_test)

    # Save model and params
    all_models.append(pipeline)
    all_params.append(best_params)

    # -----------------------------
    # Evaluate
    # -----------------------------
    rmse = np.sqrt(mean_squared_error(y_test, preds))
    mae = mean_absolute_error(y_test, preds)
    r2 = r2_score(y_test, preds)
    corr, pval = pearsonr(y_test, preds)

    metrics = {"fold": fold, "RMSE": rmse, "MAE": mae, "R2": r2, "r": corr, "p": pval}
    all_metrics.append(metrics)

    fold_preds = pd.DataFrame({
        "donor_id": donor_test,
        "dataset": dataset_test,
        "disease": "healthy",
        "actual_age": y_test,
        "predicted_age": preds,
        "sex": gender_test,
        "fold": fold
    })
    all_predictions.append(fold_preds)

# -----------------------------
# Save Aggregate Results
# -----------------------------
metrics_df = pd.DataFrame(all_metrics)

# Add mean metrics row
mean_metrics = {
    "fold": "mean",
    "RMSE": metrics_df["RMSE"].mean(),
    "MAE": metrics_df["MAE"].mean(),
    "R2": metrics_df["R2"].mean(),
    "r": metrics_df["r"].mean(),
    "p": metrics_df["p"].mean()
}
metrics_df = pd.concat([metrics_df, pd.DataFrame([mean_metrics])], ignore_index=True)

metrics_df.to_csv(snakemake.output[0], index=False)

all_preds_df = pd.concat(all_predictions)
all_preds_df.to_csv(snakemake.output[1], index=False)

# Save full model list and params list
joblib.dump(all_models, snakemake.output[2])
with open(snakemake.output[3], "w") as f:
    json.dump(all_params, f, indent=2)

