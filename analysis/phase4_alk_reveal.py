"""Reveal ALK G1202R and parental-cell controls after the frozen v0 choice.

These are additional endpoints in the same paper; they test specificity of the
retrospective contrast pattern, not independence across source or wet-lab proof.
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
ASSAYS = {"g1202r": "CHEMBL3748826", "parental": "CHEMBL3748829"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(assay_id: str) -> tuple[dict[str, dict], list[Path]]:
    paths = sorted(SRC.glob(f"chembl_activity_{assay_id}_offset*_phase4_2026-09-13.json"))
    rows, offsets, totals = {}, [], set()
    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        page = payload["page_meta"]
        offsets.append(page["offset"])
        totals.add(page["total_count"])
        for rec in payload["activities"]:
            mid = rec["molecule_chembl_id"]
            if mid in rows or rec["assay_chembl_id"] != assay_id or rec["document_chembl_id"] != "CHEMBL3745656":
                raise ValueError(f"Assay/identity mismatch for {mid}")
            if rec["standard_type"] != "IC50" or rec["standard_units"] != "nM":
                raise ValueError(f"Endpoint/unit mismatch for {mid}")
            rows[mid] = rec
    if len(totals) != 1 or len(rows) != next(iter(totals)) or offsets != list(range(0, len(rows), 10)):
        raise ValueError(f"Incomplete pages for {assay_id}")
    return rows, paths


def main() -> None:
    input_path = OUT / "phase4_alk_prereveal_input.csv"
    selection_path = OUT / "phase4_alk_v0_selection.json"
    with input_path.open(newline="", encoding="utf-8") as handle:
        initial = list(csv.DictReader(handle))
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    if selection["input_sha256"] != sha256(input_path) or selection["orthogonal_outcomes_seen_by_selector"] is not False:
        raise ValueError("Frozen pre-reveal input/selector mismatch")
    extras, paths = {}, []
    for alias, assay_id in ASSAYS.items():
        extras[alias], files = load(assay_id)
        paths.extend(files)
        if set(extras[alias]) != {r["molecule_id"] for r in initial}:
            raise ValueError(f"Molecule IDs differ between input and {alias}")
    wide = []
    for row in initial:
        mid = row["molecule_id"]
        out = {"molecule_id": mid, "wt_ic50_nM": row["first_ic50_nM"],
               "l1196m_ic50_nM": row["second_ic50_nM"]}
        for alias in ASSAYS:
            rec = extras[alias][mid]
            if rec["canonical_smiles"] != row["canonical_smiles"]:
                raise ValueError(f"Structure mismatch in reveal for {mid}")
            out[f"{alias}_ic50_nM"] = rec["standard_value"]
            out[f"{alias}_relation"] = rec["standard_relation"]
            out[f"{alias}_activity_id"] = rec["activity_id"]
        wide.append(out)
    output = OUT / "phase4_alk_revealed_extra_endpoints.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(wide[0]))
        writer.writeheader()
        writer.writerows(wide)
    by_id = {r["molecule_id"]: r for r in wide}
    comparisons = []
    for label, ids in [
        ("contrast_matched_v0", selection["contrast_matched_pair"]["molecule_ids"]),
        ("first_potency_baseline", selection["first_potency_pair"]),
        ("mean_potency_baseline", selection["mean_potency_pair"]),
    ]:
        a, b = [by_id[mid] for mid in ids]
        out = {"selection": label, "molecule_A": ids[0], "molecule_B": ids[1]}
        for axis in ("wt", "l1196m", "g1202r", "parental"):
            va, vb = a[f"{axis}_ic50_nM"], b[f"{axis}_ic50_nM"]
            ra = a.get(f"{axis}_relation", "=")
            rb = b.get(f"{axis}_relation", "=")
            out[f"A_{axis}_nM"] = va
            out[f"B_{axis}_nM"] = vb
            out[f"A_{axis}_relation"] = ra
            out[f"B_{axis}_relation"] = rb
            out[f"{axis}_pair_abs_log10_gap"] = (round(abs(math.log10(float(va) / float(vb))), 6)
                                                  if ra == rb == "=" and va not in (None, "") and vb not in (None, "") else "")
        comparisons.append(out)
    comparison_path = OUT / "phase4_alk_pair_baselines_and_reveal.csv"
    with comparison_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparisons[0]))
        writer.writeheader()
        writer.writerows(comparisons)
    manifest = {
        "claim_scope": "Additional G1202R and parental Ba/F3 readouts in the same 2015 ALK study; descriptive replay, not new measurements or independent-source mechanism test",
        "selected_v0_pair": selection["contrast_matched_pair"]["molecule_ids"],
        "selection_frozen_sha256": sha256(selection_path),
        "revealed_assay_ids": ASSAYS,
        "direct_cellular_palk_coverage_for_selected_pair": "No matched numerical phospho-ALK engagement pair in the queried 2015 ChEMBL source; CHEMBL3748833 has one qualitative compound other than both selected IDs",
        "sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p)
                   for p in [*paths, input_path, selection_path, Path(__file__), output, comparison_path]},
    }
    (OUT / "phase4_alk_reveal_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"selected": selection["contrast_matched_pair"]["molecule_ids"],
                      "comparisons": comparisons}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
