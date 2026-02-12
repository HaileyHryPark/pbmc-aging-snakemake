import os
import json
import joblib
import warnings
import pandas as pd
import numpy as np

from sklearn.linear_model import ElasticNet, ElasticNetCV
from sklearn.model_selection import train_test_split, RepeatedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from scipy.stats import pearsonr
from sklearn.utils import resample
from sklearn.model_selection import StratifiedShuffleSplit

warnings.filterwarnings("ignore")
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

# Load data
data = pd.read_csv(snakemake.input[0]).set_index("rowname")
X_full = data.drop(columns=['age', 'sex', 'dataset'])
y_full = data['age']
gender_full = data['sex']
donor_ids_full = data.index

# Subsampling setup
n_subsamples = 5
sample_fraction = 0.5
random_state = 123

# Output placeholders
all_metrics = []
all_predictions = []
all_models = []
all_params = []

# Generate 5 stratified subsamples
strat_splitter = StratifiedShuffleSplit(n_splits=n_subsamples, test_size=sample_fraction, random_state=random_state)

for fold, (sub_idx, _) in enumerate(strat_splitter.split(X_full, gender_full), 1):
    print(f"=== Subsample {fold} ===")

    X = X_full.iloc[sub_idx]
    y = y_full.iloc[sub_idx]
    gender = gender_full.iloc[sub_idx]
    donor_ids = donor_ids_full[sub_idx]

    # Inner CV train/test split (80/20)
    X_train, X_test, y_train, y_test, gender_train, gender_test, donor_train, donor_test = train_test_split(
        X, y, gender, donor_ids,
        test_size=0.2,
        random_state=fold * 42,
        stratify=gender
    )

    # -----------------------------
    # ElasticNetCV for Inner CV
    # -----------------------------
    print("X_train type:", type(X_train))
    print("X_train shape:", getattr(X_train, "shape", None))

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)

    inner_cv = RepeatedKFold(n_splits=5, n_repeats=1, random_state=123)
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
    # Train Final Model on Full Train Set
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

joblib.dump(all_models, snakemake.output[2])
with open(snakemake.output[3], "w") as f:
    json.dump(all_params, f, indent=2)

