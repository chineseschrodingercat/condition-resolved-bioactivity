"""Frozen v0 source-level experiment/pair selector (no orthogonal labels as input).

Input CSV: source_id, assay_axis, molecule_id, canonical_smiles,
first_ic50_nM, second_ic50_nM, first_relation, second_relation.
One source/assay pair per invocation. Optional orthogonal outcomes are prohibited
from the input by the fixed schema. Output is deterministic and descriptive.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from itertools import combinations
from pathlib import Path

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator


AXIS_TEST = {
    "cofactor_concentration": "cofactor_competition_and_timecourse",
    "target_site_variant": "allele_specific_engagement_or_site_control",
    "biochemical_to_cellular": "intracellular_exposure_then_target_engagement",
    "biological_matrix": "parent_and_free_exposure_over_time",
    "genetic_cell_context": "orthogonal_process_readout_by_genotype",
}
REQUIRED = {"source_id", "assay_axis", "molecule_id", "canonical_smiles",
            "first_ic50_nM", "second_ic50_nM", "first_relation", "second_relation"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    with args.input.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if set(reader.fieldnames or []) != REQUIRED:
            raise ValueError(f"Input must contain only frozen pre-reveal columns: {sorted(REQUIRED)}")
        raw_rows = list(reader)
    if not raw_rows:
        raise ValueError("Empty challenge source")
    source_ids = {r["source_id"] for r in raw_rows}
    axes = {r["assay_axis"] for r in raw_rows}
    if len(source_ids) != 1 or len(axes) != 1 or next(iter(axes)) not in AXIS_TEST:
        raise ValueError("One source and one registered assay axis are required")
    if len({r["molecule_id"] for r in raw_rows}) != len(raw_rows):
        raise ValueError("Repeated molecule ID within challenge source")
    generator = GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    eligible = []
    excluded = []
    for row in raw_rows:
        mid = row["molecule_id"]
        if row["first_relation"] != "=" or row["second_relation"] != "=":
            excluded.append({"molecule_id": mid, "reason": "censored_or_non_exact_measurement"})
            continue
        try:
            first = float(row["first_ic50_nM"])
            second = float(row["second_ic50_nM"])
        except ValueError:
            excluded.append({"molecule_id": mid, "reason": "missing_numeric_potency"})
            continue
        if not (math.isfinite(first) and math.isfinite(second) and first > 0 and second > 0):
            excluded.append({"molecule_id": mid, "reason": "nonpositive_or_nonfinite_potency"})
            continue
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None:
            excluded.append({"molecule_id": mid, "reason": "structure_not_parseable"})
            continue
        eligible.append({"id": mid, "first_p": 9 - math.log10(first),
                         "mean_p": 9 - math.log10(math.sqrt(first * second)),
                         "delta": math.log10(second / first),
                         "fp": generator.GetFingerprint(mol)})
    result = {
        "selector_version": "v0_frozen_pre_reveal",
        "source_id": next(iter(source_ids)),
        "assay_axis": next(iter(axes)),
        "action_metadata_only": AXIS_TEST[next(iter(axes))],
        "action_contrast_aware": AXIS_TEST[next(iter(axes))],
        "action_comment": "v0 changes which matched molecules to test; an assay-axis-only baseline chooses the same assay family, intentionally exposing any lack of incremental action value",
        "n_raw_molecules": len(raw_rows), "n_eligible": len(eligible), "excluded": excluded,
        "orthogonal_outcomes_seen_by_selector": False,
        "input_sha256": sha256(args.input),
        "script_sha256": sha256(Path(__file__)),
        "rdkit_version": rdBase.rdkitVersion,
    }
    if len(eligible) >= 2:
        med_delta = statistics.median(r["delta"] for r in eligible)
        ranked = []
        for a, b in combinations(eligible, 2):
            similarity = DataStructs.TanimotoSimilarity(a["fp"], b["fp"])
            potency_gap = abs(a["first_p"] - b["first_p"])
            contrast_gap = abs((a["delta"] - med_delta) - (b["delta"] - med_delta))
            score = similarity * contrast_gap / (1 + potency_gap)
            ranked.append((score, contrast_gap, similarity, min(a["id"], b["id"]), max(a["id"], b["id"])))
        ranked.sort(key=lambda x: (-x[0], -x[1], -x[2], x[3], x[4]))
        result["median_log10_second_over_first"] = round(med_delta, 6)
        result["contrast_matched_pair"] = {
            "molecule_ids": [ranked[0][3], ranked[0][4]],
            "priority_score": round(ranked[0][0], 6),
            "response_difference_log10": round(ranked[0][1], 6),
            "morgan_tanimoto": round(ranked[0][2], 6),
        }
        result["first_potency_pair"] = [r["id"] for r in sorted(eligible, key=lambda r: (-r["first_p"], r["id"]))[:2]]
        result["mean_potency_pair"] = [r["id"] for r in sorted(eligible, key=lambda r: (-r["mean_p"], r["id"]))[:2]]
        result["number_of_candidate_pairs"] = len(ranked)
    else:
        result["status"] = "not_enough_exact_structured_molecules_to_select_pair"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k not in {"excluded", "script_sha256", "input_sha256"}}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
