import os
import json
import joblib
import warnings
import pandas as pd
import numpy as np

from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.model_selection import KFold, RepeatedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr

warnings.filterwarnings("ignore")
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

# Load data
data = pd.read_csv(snakemake.input[0]).set_index("rowname")
data = data[data['sex'] == snakemake.params["gender"]]
X = data.drop(columns=["age", "sex", "dataset"]).apply(pd.to_numeric)
y = data['age']
gender = data['sex']
donor_ids = data.index

# Cross-validation
outer_cv = KFold(n_splits=5, shuffle=True, random_state=123)
inner_cv = RepeatedKFold(n_splits=5, n_repeats=1, random_state=123)

# Output placeholders
all_metrics = []
all_predictions = []
all_models = []
all_params = []

# Outer CV loop
for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y), 1):
    print(f"=== Outer Fold {fold} ===")

    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
    gender_test = gender.iloc[test_idx]
    donor_test = donor_ids[test_idx]

    # -----------------------------
    # ElasticNetCV for Inner CV
    # -----------------------------
    print("X_train type:", type(X_train))
    print("X_train shape:", getattr(X_train, "shape", None))
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)

    model_cv = ElasticNetCV(
        alphas=np.logspace(-3, 1, 15),
        l1_ratio=np.linspace(0.0, 1.0, 20),
        cv=inner_cv,
        max_iter=5000,
        random_state=123,
        n_jobs=-1
    )
    model_cv.fit(X_train_scaled, y_train)

    best_params = {
        "alpha": model_cv.alpha_,
        "l1_ratio": model_cv.l1_ratio_,
        "max_iter": 5000,
        "random_state": 123
    }

    # -----------------------------
    # Train on Outer Fold
    # -----------------------------
    model = ElasticNet(**best_params)
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
metrics_df.to_csv(snakemake.output[0], index=False)

all_preds_df = pd.concat(all_predictions)
all_preds_df.to_csv(snakemake.output[1], index=False)

# Save full model list and params list
joblib.dump(all_models, snakemake.output[2])
with open(snakemake.output[3], "w") as f:
    json.dump(all_params, f, indent=2)

