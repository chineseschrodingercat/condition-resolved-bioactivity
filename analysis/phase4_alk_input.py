"""Build the pre-reveal ALK WT/L1196M input from two original ChEMBL assays.

Only the paired MTS viability measurements and structural identities enter v0.
Other mutations, biochemical activity, structures and discussion stay outside it.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "sources"
OUT = ROOT / "outputs"
ASSAYS = {"first": "CHEMBL3748820", "second": "CHEMBL3748823"}
SOURCE = "CHEMBL3745656_PMC4907642_BaF3_WT_vs_L1196M_MTS72h"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(assay_id: str) -> tuple[dict[str, dict], list[Path]]:
    paths = sorted(SRC.glob(f"chembl_activity_{assay_id}_offset*_phase4_2026-09-13.json"))
    rows = {}
    totals, offsets = set(), []
    for path in paths:
        page = json.loads(path.read_text(encoding="utf-8"))
        totals.add(page["page_meta"]["total_count"])
        offsets.append(page["page_meta"]["offset"])
        for rec in page["activities"]:
            mid = rec["molecule_chembl_id"]
            if mid in rows or rec["assay_chembl_id"] != assay_id or rec["document_chembl_id"] != "CHEMBL3745656":
                raise ValueError(f"Duplicate or assay/document mismatch for {mid}")
            if rec["standard_type"] != "IC50" or rec["standard_units"] != "nM":
                raise ValueError(f"Unexpected units or endpoint for {mid}")
            if "viability after 72 hrs by MTS assay" not in rec["assay_description"]:
                raise ValueError(f"Not the intended 72h viability format for {mid}")
            rows[mid] = rec
    if len(totals) != 1 or len(rows) != next(iter(totals)) or offsets != list(range(0, len(rows), 10)):
        raise ValueError(f"Pagination incomplete for {assay_id}")
    return rows, paths


def main() -> None:
    first, first_paths = load(ASSAYS["first"])
    second, second_paths = load(ASSAYS["second"])
    if set(first) != set(second):
        raise ValueError(f"Molecule ID mismatch: {set(first) ^ set(second)}")
    public, audit = [], []
    for mid in sorted(first):
        a, b = first[mid], second[mid]
        if a["canonical_smiles"] != b["canonical_smiles"]:
            raise ValueError(f"Molecule structure mismatch: {mid}")
        public.append({
            "source_id": SOURCE, "assay_axis": "target_site_variant",
            "molecule_id": mid, "canonical_smiles": a["canonical_smiles"],
            "first_ic50_nM": a["standard_value"], "second_ic50_nM": b["standard_value"],
            "first_relation": a["standard_relation"], "second_relation": b["standard_relation"],
        })
        audit.append({
            "molecule_chembl_id": mid, "wt_assay_chembl_id": ASSAYS["first"],
            "l1196m_assay_chembl_id": ASSAYS["second"],
            "wt_activity_id": a["activity_id"], "l1196m_activity_id": b["activity_id"],
            "document_chembl_id": a["document_chembl_id"],
            "wt_potential_duplicate": a["potential_duplicate"],
            "l1196m_potential_duplicate": b["potential_duplicate"],
            "wt_data_validity_comment": a["data_validity_comment"],
            "l1196m_data_validity_comment": b["data_validity_comment"],
            "endpoint_note": "Ba/F3 MTS cell viability at 72 h, not a direct phospho-ALK endpoint; original Table 1 caption conflicts with footnote",
        })
    for name, rows in [("phase4_alk_prereveal_input.csv", public), ("phase4_alk_pair_provenance.csv", audit)]:
        path = OUT / name
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    manifest = {
        "source_id": SOURCE, "n_molecules": len(public),
        "known_assay_caption_conflict": "Original Table 1 caption calls phospho-ALK, but its footnote and ChEMBL assay descriptions specify 72 h MTS viability; use the latter",
        "pre_reveal_columns_only": list(public[0]),
        "raw_paths": [str(p.relative_to(ROOT)).replace("\\", "/") for p in first_paths + second_paths],
        "sha256": {str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p)
                   for p in [*first_paths, *second_paths, Path(__file__), OUT / "phase4_alk_prereveal_input.csv", OUT / "phase4_alk_pair_provenance.csv"]},
    }
    path = OUT / "phase4_alk_prereveal_manifest.json"
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"n_molecules": len(public), "input": str(OUT / "phase4_alk_prereveal_input.csv"), "manifest": str(path)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
