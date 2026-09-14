"""Reveal IDH1 endpoints only after phase5_idh1_v0_selection.json exists.

All values remain in their native assay type and ChEMBL relation. The three
two-compound strategies are evaluated with the same candidate set and budget.
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
DOC = "CHEMBL2176947"
ASSAYS = {
    "R132C_enzyme": "CHEMBL2188256",
    "HT1080_2HG": "CHEMBL2188255",
    "R132H_enzyme": "CHEMBL2188262",
    "U87_2HG": "CHEMBL2188258",
    "WT_enzyme": "CHEMBL2188253",
    "HT1080_growth": "CHEMBL2188254",
    "U87_growth": "CHEMBL2188257",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def page_records(assay: str) -> tuple[dict[str, dict], list[Path]]:
    paths = sorted(SRC.glob(f"chembl_activity_{assay}_offset*_phase5_2026-09-13.json"),
                   key=lambda p: int(p.stem.split("_offset")[1].split("_")[0]))
    if not paths:
        raise ValueError(f"No raw pages for {assay}")
    by_id: dict[str, dict] = {}
    total = set()
    offsets = []
    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        offsets.append(payload["page_meta"]["offset"])
        total.add(payload["page_meta"]["total_count"])
        for rec in payload["activities"]:
            if rec["assay_chembl_id"] != assay or rec["document_chembl_id"] != DOC:
                raise ValueError(f"Assay/document mismatch in {path}")
            mid = rec["molecule_chembl_id"]
            if mid in by_id:
                raise ValueError(f"Duplicate molecule within {assay}: {mid}")
            by_id[mid] = rec
    if len(total) != 1 or len(by_id) != next(iter(total)) or offsets != list(range(0, len(by_id), 10)):
        raise ValueError(f"Incomplete activity pages: {assay}, {offsets}, {total}")
    return by_id, paths


def fmt(rec: dict | None) -> str:
    if rec is None:
        return ""
    return f'{rec["standard_relation"]}{rec["standard_value"]}'


def gap(a: dict | None, b: dict | None) -> tuple[str, str]:
    if a is None or b is None or a["standard_units"] != b["standard_units"] or a["standard_type"] != b["standard_type"]:
        return "", "missing_or_incomparable"
    if a["standard_relation"] != "=" or b["standard_relation"] != "=":
        return "", "censored_do_not_rank_exact_gap"
    va, vb = float(a["standard_value"]), float(b["standard_value"])
    return f"{abs(math.log10(va / vb)):.6f}", "exact"


def main() -> None:
    selector_path = OUT / "phase5_idh1_v0_selection.json"
    prereveal = OUT / "phase5_idh1_prereveal_input.csv"
    selector = json.loads(selector_path.read_text(encoding="utf-8"))
    with prereveal.open(encoding="utf-8", newline="") as handle:
        input_rows = list(csv.DictReader(handle))
    mids = {r["molecule_id"] for r in input_rows}
    if len(mids) != 8 or selector["n_eligible"] != len(mids):
        raise ValueError("The frozen 8-molecule input/selection has changed")
    indexed = {}
    raw_paths = []
    for label, assay in ASSAYS.items():
        indexed[label], paths = page_records(assay)
        raw_paths.extend(paths)
        if not mids.issubset(indexed[label]):
            raise ValueError(f"Not every selected-series molecule is in {assay}")

    long = []
    for mid in sorted(mids):
        smiles = next(r["canonical_smiles"] for r in input_rows if r["molecule_id"] == mid)
        for label, assay in ASSAYS.items():
            rec = indexed[label][mid]
            if rec["canonical_smiles"] != smiles:
                raise ValueError(f"SMILES mismatch {mid} {assay}")
            if rec["standard_units"] != "nM" or rec["standard_type"] not in {"IC50", "GI50"}:
                raise ValueError(f"Units or type unexpected {mid} {assay}")
            long.append({"molecule_id": mid, "assay_axis": label, "assay_id": assay,
                         "activity_id": rec["activity_id"], "source_document": DOC,
                         "assay_description": rec["assay_description"],
                         "standard_type": rec["standard_type"], "relation": rec["standard_relation"],
                         "value": rec["standard_value"], "units": rec["standard_units"],
                         "potential_duplicate": rec["potential_duplicate"],
                         "data_validity_comment": rec["data_validity_comment"]})
    long_path = OUT / "phase5_idh1_all_orthogonal_endpoints.csv"
    write_csv(long_path, long)

    strategies = {
        "contrast_matched_v0": selector["contrast_matched_pair"]["molecule_ids"],
        "first_potency": selector["first_potency_pair"],
        "geometric_mean_potency": selector["mean_potency_pair"],
    }
    rows = []
    for strategy, pair in strategies.items():
        if len(pair) != 2 or not set(pair).issubset(mids):
            raise ValueError(f"Bad selection {strategy}: {pair}")
        row = {"strategy": strategy, "molecule_a": pair[0], "molecule_b": pair[1],
               "selection_budget_compounds": 2, "source_document": DOC}
        for label in ASSAYS:
            a, b = indexed[label][pair[0]], indexed[label][pair[1]]
            val, status = gap(a, b)
            row[f"{label}_a_nM"] = fmt(a)
            row[f"{label}_b_nM"] = fmt(b)
            row[f"{label}_endpoint_type"] = a["standard_type"]
            row[f"{label}_within_pair_abs_log10_gap"] = val
            row[f"{label}_gap_status"] = status
        rows.append(row)
    comp_path = OUT / "phase5_idh1_equal_budget_reveal.csv"
    write_csv(comp_path, rows)
    manifest = {"doi": "10.1021/ml300225h", "document_id": DOC,
                "n_molecules": len(mids), "n_assays": len(ASSAYS), "n_rows": len(long),
                "strategies": strategies,
                "limitation": "one publication/chemical series; literature retrospective, not independent source validation",
                "sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha(p)
                           for p in [prereveal, selector_path, Path(__file__), *raw_paths, long_path, comp_path]}}
    manifest_path = OUT / "phase5_idh1_reveal_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"n_molecules": len(mids), "n_rows": len(long),
                      "selected": strategies["contrast_matched_v0"],
                      "pair_comparison": rows}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
