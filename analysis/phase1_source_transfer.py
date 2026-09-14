"""Transfer the ATP response learned in one JAK3 paper to the other.

The 2017 low-ATP assay uses 4 uM and microfluidics; the 2019 assay uses 6 uM
and HTRF. Both high-ATP assays use 1 mM. This is a deliberately difficult
source/method/chemical-series transfer check, not a clean unseen-condition test.
The one shared molecule (PF-06651600) is excluded from each test set.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from phase1_four_arm import OUT, fingerprint, ridge_predict


FILES = {
    "2019_HTRF_6uM_to_1mM": OUT / "jak3_atp_current_20_pairs.csv",
    "2017_microfluidic_4uM_to_1mM": OUT / "jak3_2017_4uM_vs_1mM_matched_pairs.csv",
}


def main() -> None:
    pairs = {name: pd.read_csv(path) for name, path in FILES.items()}
    results = []
    for train_name, train in pairs.items():
        test_name = next(name for name in pairs if name != train_name)
        test_all = pairs[test_name]
        test = test_all[~test_all.molecule_chembl_id.isin(train.molecule_chembl_id)].copy()
        assert len(test) > 0
        x_train = np.vstack([fingerprint(s) for s in train.canonical_smiles_first])
        delta_train = (train.pIC50_second - train.pIC50_first).to_numpy(float)
        constant_delta = float(np.mean(delta_train))
        for _, row in test.iterrows():
            x_test = fingerprint(row.canonical_smiles_first)[None, :]
            actual_delta = float(row.pIC50_second - row.pIC50_first)
            predicted_structure_delta = ridge_predict(x_train, delta_train, x_test)
            for method, predicted_delta in [
                ("zero_delta", 0.0),
                ("source_mean_delta", constant_delta),
                ("source_structure_delta", predicted_structure_delta),
            ]:
                results.append({
                    "train_source": train_name,
                    "test_source": test_name,
                    "method": method,
                    "molecule_chembl_id": row.molecule_chembl_id,
                    "n_training_molecules": len(train),
                    "n_test_molecules_after_shared_id_exclusion": len(test),
                    "actual_pIC50_delta_second_minus_first": actual_delta,
                    "predicted_pIC50_delta_second_minus_first": predicted_delta,
                    "absolute_delta_error_log10": abs(actual_delta - predicted_delta),
                })
    pred = pd.DataFrame(results)
    pred.to_csv(OUT / "phase1_source_transfer_predictions.csv", index=False)
    metrics = (pred.groupby(["train_source", "test_source", "method"], sort=False)
               .agg(n_test_molecules=("molecule_chembl_id", "size"),
                    mae_delta_log10=("absolute_delta_error_log10", "mean"))
               .reset_index())
    metrics.to_csv(OUT / "phase1_source_transfer_metrics.csv", index=False)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
