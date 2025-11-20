import pandas as pd
import numpy as np
import argparse

SEED = 123

# -----------------------------
# Utility functions
# -----------------------------
def get_seen_donors(pred_df, fold):
    """
    Identify donors that were SEEN during model training.
    In your setup, 5-fold CV means: if fold == k, that donor was unseen in that fold.
    So all donors except those in the test fold of the final model were SEEN.
    """
    all_donors = set(pred_df["donor_id"].unique())
    unseen_donors = set(pred_df.loc[pred_df["fold"] == fold, "donor_id"])
    seen_donors = all_donors - unseen_donors
    return seen_donors, unseen_donors


def stratified_sample(df, n_per_sex=50, random_state=123):
    """
    Sample 50 male + 50 female donors, age-balanced within each sex.
    Returns dataframe of selected donors.
    """
    np.random.seed(SEED)
    sampled_donors = []

    for sex in ["male", "female"]:
        sex_df = df[df["sex"] == sex]
        if len(sex_df) == 0:
            continue

        # Split by median age for balance
        median_age = sex_df["actual_age"].median()
        younger = sex_df[sex_df["actual_age"] <= median_age]
        older = sex_df[sex_df["actual_age"] > median_age]

        n_half = n_per_sex // 2
        younger_sample = younger.sample(n=min(n_half, len(younger)), random_state=random_state)
        older_sample = older.sample(n=min(n_half, len(older)), random_state=random_state)
        sampled_donors.append(pd.concat([younger_sample, older_sample]))

    sampled_df = pd.concat(sampled_donors).drop_duplicates("donor_id")
    return sampled_df


# -----------------------------
# Main function
# -----------------------------
def main(both_csv, both_fold, female_csv, female_fold, male_csv, male_fold, output_csv, n_per_sex=50):
    # Load predictions
    both_df = pd.read_csv(both_csv)
    female_df = pd.read_csv(female_csv)
    male_df = pd.read_csv(male_csv)

    # Identify seen donors for each model
    both_seen, both_unseen = get_seen_donors(both_df, both_fold)
    female_seen, female_unseen = get_seen_donors(female_df, female_fold)
    male_seen, male_unseen = get_seen_donors(male_df, male_fold)

    # Identify overlapping donors
    female_candidates = both_df[
        both_df["donor_id"].isin(both_seen & female_seen)
    ].drop_duplicates("donor_id")

    male_candidates = both_df[
        both_df["donor_id"].isin(both_seen & male_seen)
    ].drop_duplicates("donor_id")

    print(f"Found {len(female_candidates)} female and {len(male_candidates)} male candidate donors for background.")

    # Combine for balanced sampling
    candidates_df = pd.concat([female_candidates, male_candidates])

    # Stratified sampling for background donors
    background_df = stratified_sample(candidates_df, n_per_sex=n_per_sex, random_state=SEED)

    # Save selected donors
    background_df.to_csv(output_csv, index=False)
    print(f"Background donor set saved to: {output_csv}")
    print(f"Total donors: {len(background_df)} (female={sum(background_df['sex']=='female')}, male={sum(background_df['sex']=='male')})")
    print("Sampled donor IDs:")
    print(background_df["donor_id"].tolist())


# -----------------------------
# CLI
# -----------------------------
if __name__ == "__main__":

    both = snakemake.input[0]
    both_fold = snakemake.params["both_fold"]
    female = snakemake.input[1]
    female_fold = snakemake.params["female_fold"]
    male = snakemake.input[2]
    male_fold = snakemake.params["male_fold"]
    output = snakemake.output[0]

    main(both, both_fold, female, female_fold, male, male_fold, output, 50)

