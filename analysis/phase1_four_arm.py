"""Small, reproducible assay-record handling pilot on four matched pairs.

This is a within-series, leave-one-molecule-out test of predicting both assay
readouts for a new molecule.  It is not a held-out assay, series, or mechanism
test.  Every policy starts from the same training molecules and raw pair.

Run with the ccs_gnn Python environment after phase0_rebuild.py.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
PAIRS = {
    "JAK3_2019_ATP": "jak3_atp_current_20_pairs.csv",
    "JAK3_2017_ATP": "jak3_2017_4uM_vs_1mM_matched_pairs.csv",
    "renin_plasma_trypsin": "renin_plasma_vs_trypsin_matched_pairs.csv",
    "KDR_HTRF_cell": "kdr_htrf_vs_cell_matched_pairs.csv",
}
ALPHA = 1.0
FINGERPRINT_BITS = 512
FINGERPRINT_RADIUS = 2
GEN = rdFingerprintGenerator.GetMorganGenerator(
    radius=FINGERPRINT_RADIUS, fpSize=FINGERPRINT_BITS
)


def fingerprint(smiles: str) -> np.ndarray:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Unparseable SMILES: {smiles}")
    bits = GEN.GetFingerprint(mol)
    arr = np.zeros((FINGERPRINT_BITS,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(bits, arr)
    return arr.astype(float)


def linear_solve(matrix: list[list[float]], target: list[float]) -> list[float]:
    """Gaussian elimination with pivoting for small, positive-definite systems.

    This deliberately avoids the broken BLAS dependency in the local sklearn
    environment; the fitted objective is still Ridge(alpha=1, intercept=True).
    """
    n = len(target)
    aug = [row[:] + [value] for row, value in zip(matrix, target)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        if abs(scale) < 1e-12:
            raise ArithmeticError("Near-singular ridge system")
        for j in range(col, n + 1):
            aug[col][j] /= scale
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            for j in range(col, n + 1):
                aug[row][j] -= factor * aug[col][j]
    return [row[n] for row in aug]


def ridge_predict(features: np.ndarray, target: np.ndarray,
                  new: np.ndarray) -> float:
    """Intercept-fitted linear Ridge in dual form, solved without native BLAS."""
    rows = features.tolist()
    y = target.tolist()
    n = len(rows)
    n_features = len(rows[0])
    x_mean = [sum(row[j] for row in rows) / n for j in range(n_features)]
    centered = [[v - x_mean[j] for j, v in enumerate(row)] for row in rows]
    y_mean = sum(y) / n
    kernel = [
        [sum(a * b for a, b in zip(centered[i], centered[j]))
         + (ALPHA if i == j else 0.0) for j in range(n)]
        for i in range(n)
    ]
    coeff = linear_solve(kernel, [value - y_mean for value in y])
    centered_new = [v - x_mean[j] for j, v in enumerate(new[0].tolist())]
    return y_mean + sum(
        coeff[i] * sum(a * b for a, b in zip(centered[i], centered_new))
        for i in range(n)
    )


def hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def analyze(name: str, path: Path, drop_possible_duplicates: bool = False) -> pd.DataFrame:
    pair = pd.read_csv(path)
    if drop_possible_duplicates:
        pair = pair[
            pair.potential_duplicate_first.eq(0)
            & pair.potential_duplicate_second.eq(0)
        ].copy()
    assert pair.molecule_chembl_id.is_unique
    assert (pair.nM_first.gt(0) & pair.nM_second.gt(0)).all()
    assert (pair.canonical_smiles_first == pair.canonical_smiles_second).all()
    x = np.vstack([fingerprint(s) for s in pair.canonical_smiles_first])
    y0 = pair.pIC50_first.to_numpy(float)
    y1 = pair.pIC50_second.to_numpy(float)
    mean = (y0 + y1) / 2  # geometric mean IC50 on the original nM scale
    delta = y1 - y0
    predictions = []
    for held in range(len(pair)):
        train = np.arange(len(pair)) != held
        xtr = x[train]
        xte = x[held:held + 1]
        mtr = mean[train]
        dtr = delta[train]

        # Averaging both assays preserves one molecular response; picking only
        # the first is a concrete, threshold-free "discard the other assay" rule.
        p_mean = ridge_predict(xtr, mtr, xte)
        p_first = ridge_predict(xtr, y0[train], xte)

        # Balanced pooled records have two rows per molecule.  Weight .5 per
        # row matches the per-molecule objective of the averaged arm; this is
        # a control for the mistaken idea that row duplication adds signal.
        # For identical features x_i in each assay and weights 0.5 each,
        # sum_j 0.5(y_ij-f(x_i))^2 = (mean_i-f(x_i))^2 + constant.
        # Therefore this control is exactly equal to mean-label Ridge.
        p_pooled = p_mean

        # Encoding the condition first learns the overall assay offset.  The
        # interaction arm asks whether molecular structure predicts deviation
        # from that offset.  Both use Ridge with the same frozen fingerprint,
        # alpha, train molecules and raw observed pair; no test label enters fit.
        p_offset = float(np.mean(dtr))
        p_interaction = ridge_predict(xtr, dtr, xte)
        methods = {
            "average_pIC50": (p_mean, p_mean),
            "discard_second_assay": (p_first, p_first),
            "pool_no_context": (p_pooled, p_pooled),
            "context_global_offset": (p_mean - p_offset / 2, p_mean + p_offset / 2),
            "context_structure_interaction": (
                p_mean - p_interaction / 2, p_mean + p_interaction / 2
            ),
        }
        for method, (pred0, pred1) in methods.items():
            predictions.append({
                "pair": name,
                "method": method,
                "molecule_chembl_id": pair.molecule_chembl_id.iloc[held],
                "activity_id_first": int(pair.activity_id_first.iloc[held]),
                "activity_id_second": int(pair.activity_id_second.iloc[held]),
                "actual_pIC50_first": y0[held],
                "actual_pIC50_second": y1[held],
                "predicted_pIC50_first": pred0,
                "predicted_pIC50_second": pred1,
                "actual_delta_second_minus_first": delta[held],
                "predicted_delta_second_minus_first": pred1 - pred0,
                "n_training_molecules": int(train.sum()),
                "n_raw_training_activity_records_available": int(train.sum() * 2),
                "ridge_alpha": ALPHA,
            })
    return pd.DataFrame(predictions)


def main() -> None:
    frames = [analyze(name, OUT / filename) for name, filename in PAIRS.items()]
    frames.append(analyze(
        "JAK3_2019_ATP_flag_sensitivity",
        OUT / PAIRS["JAK3_2019_ATP"],
        drop_possible_duplicates=True,
    ))
    all_predictions = pd.concat(frames, ignore_index=True)
    all_predictions.to_csv(OUT / "phase1_four_arm_predictions.csv", index=False)
    metrics = []
    for (pair, method), frame in all_predictions.groupby(["pair", "method"], sort=False):
        e0 = np.abs(frame.actual_pIC50_first - frame.predicted_pIC50_first)
        e1 = np.abs(frame.actual_pIC50_second - frame.predicted_pIC50_second)
        ed = np.abs(frame.actual_delta_second_minus_first - frame.predicted_delta_second_minus_first)
        metrics.append({
            "pair": pair, "method": method, "n_test_molecules": len(frame),
            "mae_both_conditions_pIC50": float((e0.sum() + e1.sum()) / (2 * len(frame))),
            "mae_first_pIC50": float(e0.mean()),
            "mae_second_pIC50": float(e1.mean()),
            "mae_delta_log10": float(ed.mean()),
            "median_abs_delta_error_log10": float(ed.median()),
        })
    summary = pd.DataFrame(metrics)
    summary.to_csv(OUT / "phase1_four_arm_metrics.csv", index=False)
    manifest = {
        "purpose": "LOO-molecule within-assay-pair pilot, not external mechanism validation",
        "date": "2026-09-13",
        "feature": f"RDKit Morgan radius={FINGERPRINT_RADIUS}, {FINGERPRINT_BITS} binary bits",
        "model": f"Linear Ridge(alpha={ALPHA}, fit_intercept=True), dual kernel solved by pure-Python Gaussian elimination because local native BLAS crashed; fixed before evaluation",
        "average_rule": "mean of paired pIC50; geometric mean of IC50 on nM scale",
        "discard_rule": "keep first assay record per molecule, discard second; no tuned threshold",
        "pool_rule": "both raw records, no condition, weight 0.5 each per molecule",
        "context_rule": "same structure model for pair mean, plus global or structure-predicted paired delta",
        "resource_accounting": "Same molecules and same two raw activity records initially available to each policy. Discard uses one per molecule; averaging combines two into one label; pooled/context retain two. The balanced pooled fit gives each molecule total weight 1.",
        "sensitivity": "A fifth output group removes any 2019 JAK3 pair marked potential_duplicate in either assay; it is not part of the primary four-pair comparison.",
        "literature_source_boundary": "The 2017 article was read via institutional subscribed access; only public ChEMBL snapshots, DOI links and our derived CSV/figures are stored here. The 2020 Petri article was open access.",
        "software": {"python": sys.version.split()[0], "numpy": np.__version__,
                     "pandas": pd.__version__, "rdkit": rdkit.__version__},
        "sha256": {
            str(path.relative_to(ROOT)).replace("\\", "/"): hash_file(path)
            for path in [
                ROOT / "analysis/phase1_four_arm.py",
                ROOT / "analysis/phase1_source_transfer.py",
                ROOT / "analysis/phase1_plot.py",
                ROOT / "sources/chembl_molecule_CHEMBL4085582_2026-09-13.json",
                ROOT / "sources/chembl_activity_molecule_CHEMBL4085582_offset0_2026-09-13.json",
                ROOT / "sources/chembl_activity_molecule_CHEMBL4085582_offset10_2026-09-13.json",
                ROOT / "sources/chembl_document_CHEMBL5230128_2026-09-13.json",
                ROOT / "reports/PHASE1_2026-09-13.md",
                OUT / "phase1_four_arm_predictions.csv",
                OUT / "phase1_four_arm_metrics.csv",
                OUT / "phase1_source_transfer_predictions.csv",
                OUT / "phase1_source_transfer_metrics.csv",
                OUT / "phase1_record_policy_and_transfer.png",
            ] + [OUT / filename for filename in PAIRS.values()]
            if path.exists()
        },
    }
    (OUT / "phase1_four_arm_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
