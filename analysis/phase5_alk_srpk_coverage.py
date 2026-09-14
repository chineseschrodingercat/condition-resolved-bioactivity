"""Build a coverage-conditioned ALK/SRPK1 comparison and source audit.

Run with --prepare, then the unmodified frozen Phase 4 selector on its output,
then --reveal. This is explicitly retrospective: SRPK coverage defines the
candidate universe, but SRPK values do not enter the selector.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "sources"
OUT = ROOT / "outputs"
INPUT = OUT / "phase4_alk_prereveal_input.csv"
ELIGIBLE = OUT / "phase5_alk_srpk_coverage_input.csv"
SELECTION = OUT / "phase5_alk_srpk_coverage_v0_selection.json"
PAGES = sorted(SRC.glob("chembl_activity_CHEMBL5736981_offset*_phase5_2026-09-13.json"),
               key=lambda p: int(p.stem.split("_offset")[1].split("_")[0]))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_srpk() -> dict[str, dict]:
    if len(PAGES) != 4:
        raise ValueError("SRPK1 pagination not complete")
    rows = []
    for path in PAGES:
        rows += json.loads(path.read_text(encoding="utf-8"))["activities"]
    if len(rows) != 39:
        raise ValueError("Expected all 39 source assay records")
    out = {}
    for r in rows:
        if r["standard_type"] == "IC50":
            if r["standard_units"] != "nM" or r["molecule_chembl_id"] in out:
                raise ValueError("SRPK1 IC50 schema/duplicate issue")
            out[r["molecule_chembl_id"]] = r
    if len(out) != 13:
        raise ValueError("Expected 13 patent SRPK1 IC50 records")
    return out


def csv_out(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["prepare", "reveal"])
    stage = parser.parse_args().stage
    with INPUT.open(encoding="utf-8", newline="") as f:
        alk = list(csv.DictReader(f))
    srpk = read_srpk()
    common = sorted({r["molecule_id"] for r in alk} & set(srpk))
    if len(common) != 9:
        raise ValueError(f"Expected 9 structurally matched ALK/SRPK compounds, found {len(common)}")
    if stage == "prepare":
        selected = [r for r in alk if r["molecule_id"] in common]
        csv_out(ELIGIBLE, selected)
        print(json.dumps({"n_eligible": len(selected), "input_sha256": sha(ELIGIBLE),
                          "qualification": "SRPK1 outcome availability defines this retrospective subset"}))
        return
    if not ELIGIBLE.exists() or not SELECTION.exists():
        raise ValueError("Run prepare, then the frozen selector, before reveal")
    sel = json.loads(SELECTION.read_text(encoding="utf-8"))
    pairs = {"contrast_matched_v0": sel["contrast_matched_pair"]["molecule_ids"],
             "first_potency": sel["first_potency_pair"],
             "geometric_mean_potency": sel["mean_potency_pair"]}
    source_pair = {"CHEMBL3746665": (">", 1000.0, "patent >10^3; ChEMBL >103 likely exponent transcription"),
                   "CHEMBL3746960": ("=", 2.66, "patent Table 9 exact")}
    table = []
    for mid in common:
        r = srpk[mid]
        rel, val, note = source_pair.get(mid, (r["standard_relation"], float(r["standard_value"]),
                                               "ChEMBL patent transcription, not independently source-checked here"))
        table.append({"molecule_id": mid, "patent": "US11066363B2 Example 14 Table 9",
                      "chembl_document_id": r["document_chembl_id"], "assay_id": r["assay_chembl_id"],
                      "activity_id": r["activity_id"], "chembl_relation": r["standard_relation"],
                      "chembl_value_nM": r["standard_value"], "source_checked_relation": rel,
                      "source_checked_value_or_bound_nM": val, "source_note": note,
                      "standard_type": r["standard_type"]})
    audit_path = OUT / "phase5_alk_srpk_9_compound_audit.csv"
    csv_out(audit_path, table)
    idx = {r["molecule_id"]: r for r in table}
    comparison = []
    for strategy, (a, b) in pairs.items():
        x, y = idx[a], idx[b]
        va, vb = float(x["source_checked_value_or_bound_nM"]), float(y["source_checked_value_or_bound_nM"])
        exact = x["source_checked_relation"] == y["source_checked_relation"] == "="
        lower = ""
        relations = (x["source_checked_relation"], y["source_checked_relation"])
        if not exact and relations.count(">") == 1 and relations.count("=") == 1:
            bounded = x if relations[0] == ">" else y
            measured = y if relations[0] == ">" else x
            if float(bounded["source_checked_value_or_bound_nM"]) >= float(measured["source_checked_value_or_bound_nM"]):
                lower = f'{float(bounded["source_checked_value_or_bound_nM"]) / float(measured["source_checked_value_or_bound_nM"]):.6f}'
        comparison.append({"strategy": strategy, "molecule_a": a, "molecule_b": b,
                           "a_srpk1_relation": x["source_checked_relation"],
                           "a_srpk1_nM_or_bound": va, "b_srpk1_relation": y["source_checked_relation"],
                           "b_srpk1_nM_or_bound": vb,
                           "srpk1_pair_abs_log10_gap_if_exact": f"{abs(math.log10(va / vb)):.6f}" if exact else "",
                           "srpk1_pair_fold_gap_lower_bound_if_valid": lower,
                           "source_level_note": "same chemical series/inventors; patent cross-target assay, not independent lab replication"})
    comp_path = OUT / "phase5_alk_srpk_9_coverage_equal_budget.csv"
    csv_out(comp_path, comparison)
    manifest = {"n_overlap_alk_patent_srpk1": len(common), "candidate_universe_conditioned_on_srpk1_coverage": True,
                "same_budget": 2, "selection_file": str(SELECTION.relative_to(ROOT)),
                "source_correction": "compound 17 >10^3 nM patent, not ChEMBL >103 nM",
                "sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha(p)
                           for p in [Path(__file__), INPUT, ELIGIBLE, SELECTION, *PAGES, audit_path, comp_path]}}
    (OUT / "phase5_alk_srpk_coverage_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"n_common": len(common), "pair_comparison": comparison}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
