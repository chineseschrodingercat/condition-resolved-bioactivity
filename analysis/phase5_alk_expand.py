"""Validate and expand the original ALK series without rewriting frozen Phase 4.

Each 72 h Ba/F3 MTS allele is kept as its own condition. Missing/censored
pair gaps are not converted to exact values. The patent SRPK bridge is flagged
as a separate source linked by exact ChEMBL molecule identity.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "sources"
OUT = ROOT / "outputs"
DOC = "CHEMBL3745656"
ASSAYS = {
    "WT": "CHEMBL3748820", "C1156Y": "CHEMBL3748821", "F1174L": "CHEMBL3748822",
    "L1196M": "CHEMBL3748823", "L1152R": "CHEMBL3748824", "1151Tins": "CHEMBL3748825",
    "G1202R": "CHEMBL3748826", "G1269A": "CHEMBL3748827", "S1206Y": "CHEMBL3748828",
    "parental_BaF3": "CHEMBL3748829",
}
DOSE_ASSAYS = {"WT_1uM": "CHEMBL3748830", "G1202R_1uM": "CHEMBL3748831",
               "parental_1uM": "CHEMBL3748832"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def records(assay: str) -> tuple[dict[str, dict], list[Path]]:
    paths = list(SRC.glob(f"chembl_activity_{assay}_offset*_phase?_2026-09-13.json"))
    paths.sort(key=lambda p: int(p.stem.split("_offset")[1].split("_")[0]))
    if not paths:
        raise ValueError(f"Missing source pages {assay}")
    out = {}
    offsets, totals = [], set()
    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        offsets.append(payload["page_meta"]["offset"])
        totals.add(payload["page_meta"]["total_count"])
        for rec in payload["activities"]:
            if rec["assay_chembl_id"] != assay or rec["document_chembl_id"] != DOC:
                raise ValueError(f"Bad assay/source at {path}")
            mid = rec["molecule_chembl_id"]
            if mid in out:
                raise ValueError(f"Multiple rows {assay} {mid}")
            out[mid] = rec
    if len(totals) != 1 or len(out) != next(iter(totals)) or offsets != list(range(0, len(out), 10)):
        raise ValueError(f"Incomplete pagination {assay}: {offsets}, {totals}")
    return out, paths


def exact_gap(a: dict, b: dict) -> str:
    if a["standard_relation"] != "=" or b["standard_relation"] != "=":
        return ""
    return f'{abs(math.log10(float(a["standard_value"]) / float(b["standard_value"]))):.6f}'


def main() -> None:
    by_condition, source_paths = {}, []
    for cond, aid in ASSAYS.items():
        by_condition[cond], paths = records(aid)
        source_paths.extend(paths)
        if len(by_condition[cond]) != 18:
            raise ValueError(f"Expected 18 candidate compounds in {aid}")
    mids = set(by_condition["WT"])
    if any(set(v) != mids for v in by_condition.values()):
        raise ValueError("ALK condition molecule identity sets differ")
    smi = {m: by_condition["WT"][m]["canonical_smiles"] for m in mids}
    long = []
    wide = []
    for mid in sorted(mids):
        row = {"molecule_id": mid, "canonical_smiles": smi[mid], "source_document": DOC}
        for cond, aid in ASSAYS.items():
            rec = by_condition[cond][mid]
            if (rec["canonical_smiles"] != smi[mid] or rec["standard_type"] != "IC50"
                or rec["standard_units"] != "nM"):
                raise ValueError(f"Identity or measurement mismatch: {cond}/{mid}")
            row[f"{cond}_relation"] = rec["standard_relation"]
            row[f"{cond}_nM"] = rec["standard_value"]
            long.append({"molecule_id": mid, "condition": cond, "assay_id": aid,
                         "activity_id": rec["activity_id"], "source_document": DOC,
                         "assay_description": rec["assay_description"],
                         "relation": rec["standard_relation"], "value_nM": rec["standard_value"],
                         "potential_duplicate": rec["potential_duplicate"],
                         "data_validity_comment": rec["data_validity_comment"]})
        wide.append(row)
    matrix_path = OUT / "phase5_alk_18x10_condition_matrix.csv"
    long_path = OUT / "phase5_alk_18x10_source_audit.csv"
    write_csv(matrix_path, wide)
    write_csv(long_path, long)

    sel_path = OUT / "phase4_alk_v0_selection.json"
    sel = json.loads(sel_path.read_text(encoding="utf-8"))
    strategies = {"contrast_matched_v0": sel["contrast_matched_pair"]["molecule_ids"],
                  "first_potency": sel["first_potency_pair"],
                  "geometric_mean_potency": sel["mean_potency_pair"]}
    pair_rows = []
    for strategy, pair in strategies.items():
        for cond in ASSAYS:
            a, b = by_condition[cond][pair[0]], by_condition[cond][pair[1]]
            pair_rows.append({"strategy": strategy, "molecule_a": pair[0], "molecule_b": pair[1],
                              "condition": cond, "assay_id": ASSAYS[cond],
                              "a_ic50_nM": a["standard_value"], "a_relation": a["standard_relation"],
                              "b_ic50_nM": b["standard_value"], "b_relation": b["standard_relation"],
                              "within_pair_abs_log10_gap": exact_gap(a, b),
                              "readout": "Ba/F3 72h MTS growth/viability IC50, not direct ALK target engagement",
                              "selection_status": "used_for_pair_selection" if cond in {"WT", "L1196M"}
                              else "not_used_for_pair_selection"})
    pair_path = OUT / "phase5_alk_equal_budget_10_condition_comparison.csv"
    write_csv(pair_path, pair_rows)

    pair = strategies["contrast_matched_v0"]
    dose_rows = []
    for cond, aid in DOSE_ASSAYS.items():
        data, paths = records(aid)
        source_paths.extend(paths)
        if len(data) != 31:
            raise ValueError(f"Expected 31 fixed-dose compounds in {aid}")
        for mid in pair:
            r = data[mid]
            if r["standard_type"] != "Activity" or r["standard_units"] != "%":
                raise ValueError(f"Wrong fixed-dose field in {aid}")
            dose_rows.append({"molecule_id": mid, "condition": cond, "assay_id": aid,
                              "activity_id": r["activity_id"], "relation": r["standard_relation"],
                              "reported_activity_percent_at_1uM": r["standard_value"],
                              "note": "ChEMBL label Activity (%); sign/meaning not inferred from label alone"})
    dose_path = OUT / "phase5_alk_selected_pair_1uM_activity.csv"
    write_csv(dose_path, dose_rows)

    mol_rows = []
    for mid in pair:
        path = SRC / f"chembl_molecule_{mid}_phase5_2026-09-13.json"
        payload = json.loads(path.read_text(encoding="utf-8"))
        obj = payload.get("molecule", payload)
        # Keep a compact, auditable identity/descriptor summary; no exposure is inferred.
        props = obj.get("molecule_properties") or {}
        structures = obj.get("molecule_structures") or {}
        mol_rows.append({"molecule_id": mid, "canonical_smiles": smi[mid],
                         "standard_inchi_key": structures.get("standard_inchi_key", ""),
                         "alogp_calculated": props.get("alogp", ""),
                         "mw_freebase_calculated": props.get("mw_freebase", ""),
                         "psa_calculated": props.get("psa", ""),
                         "exposure_measurement": "not_found_in_current_record"})
        source_paths.append(path)
    mol_path = OUT / "phase5_alk_selected_pair_structure_identity.csv"
    write_csv(mol_path, mol_rows)
    manifest = {"doi": "10.1021/acs.jmedchem.5b01136", "n_molecules": len(mids),
                "n_conditions": len(ASSAYS), "n_activity_rows": len(long),
                "pair_selection_source": "Phase4 frozen WT/L1196M 18-molecule selection",
                "structure_note": "RCSB 3AOX contains an ALK scaffold comparator, not compound 17 or 19",
                "sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha(p)
                           for p in [Path(__file__), sel_path, *source_paths, matrix_path,
                                     long_path, pair_path, dose_path, mol_path]}}
    (OUT / "phase5_alk_expand_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"n_molecules": len(mids), "n_conditions": len(ASSAYS),
                      "selected_profile": [r for r in pair_rows if r["strategy"] == "contrast_matched_v0"],
                      "molecules": mol_rows}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
