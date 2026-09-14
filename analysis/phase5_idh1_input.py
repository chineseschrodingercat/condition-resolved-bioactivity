"""Build the locked-input IDH1 R132C enzyme -> HT1080 2-HG challenge.

This has the same eight-column schema consumed by the Phase 4 frozen v0
selector. Only the two registered IC50 assays enter that algorithm.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "sources"
OUT = ROOT / "outputs"
DOC = "CHEMBL2176947"
ASSAYS = {"first": "CHEMBL2188256", "second": "CHEMBL2188255"}
SOURCE = "CHEMBL2176947_IDH1_R132C_enzyme_vs_HT1080_2HG_48h"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def get(assay: str) -> tuple[dict[str, dict], list[Path]]:
    paths = sorted(SRC.glob(f"chembl_activity_{assay}_offset*_phase5_2026-09-13.json"))
    if not paths:
        raise ValueError(f"Missing original activity pages: {assay}")
    rows = {}
    offsets, totals = [], set()
    for path in paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        offsets.append(payload["page_meta"]["offset"])
        totals.add(payload["page_meta"]["total_count"])
        for rec in payload["activities"]:
            mid = rec["molecule_chembl_id"]
            if mid in rows:
                raise ValueError(f"Duplicate compound in {assay}: {mid}")
            if (rec["assay_chembl_id"] != assay or rec["document_chembl_id"] != DOC
                or rec["standard_type"] != "IC50" or rec["standard_units"] != "nM"):
                raise ValueError(f"Wrong source, assay, type, or units for {mid}")
            rows[mid] = rec
    if len(totals) != 1 or len(rows) != next(iter(totals)) or offsets != list(range(0, len(rows), 10)):
        raise ValueError(f"Incomplete pagination in {assay}")
    return rows, paths


def main() -> None:
    first, p1 = get(ASSAYS["first"])
    second, p2 = get(ASSAYS["second"])
    if set(first) != set(second):
        raise ValueError(f"Molecule sets differ: {set(first) ^ set(second)}")
    if not all("R132C mutant using alpha-ketoglutarate" in r["assay_description"] for r in first.values()):
        raise ValueError("First condition is not R132C enzyme")
    if not all("R132C mutant overexpressed in human HT1080" in r["assay_description"]
               and "2-hydroxyglutarate" in r["assay_description"] for r in second.values()):
        raise ValueError("Second condition is not HT1080 2-HG")
    inp, audit = [], []
    for mid in sorted(first):
        a, b = first[mid], second[mid]
        if a["canonical_smiles"] != b["canonical_smiles"]:
            raise ValueError(f"Structure differs across conditions: {mid}")
        inp.append({"source_id": SOURCE, "assay_axis": "biochemical_to_cellular",
                    "molecule_id": mid, "canonical_smiles": a["canonical_smiles"],
                    "first_ic50_nM": a["standard_value"],
                    "second_ic50_nM": b["standard_value"],
                    "first_relation": a["standard_relation"],
                    "second_relation": b["standard_relation"]})
        audit.append({"molecule_id": mid, "document_chembl_id": DOC,
                      "enzyme_assay_id": ASSAYS["first"], "enzyme_activity_id": a["activity_id"],
                      "cell_2hg_assay_id": ASSAYS["second"], "cell_2hg_activity_id": b["activity_id"],
                      "enzyme_potential_duplicate": a["potential_duplicate"],
                      "cell_potential_duplicate": b["potential_duplicate"],
                      "enzyme_data_validity_comment": a["data_validity_comment"],
                      "cell_data_validity_comment": b["data_validity_comment"],
                      "independent_unit_note": "Same 2012 chemical series, not independent source-level validation"})
    paths = []
    for name, rows in [("phase5_idh1_prereveal_input.csv", inp),
                       ("phase5_idh1_pair_provenance.csv", audit)]:
        path = OUT / name
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        paths.append(path)
    manifest = {
        "source": SOURCE, "paper_doi": "10.1021/ml300225h",
        "n_shared_compounds": len(inp), "allowed_input_columns": list(inp[0]),
        "enzyme_endpoint": "IDH1 R132C biochemical IC50, alpha-ketoglutarate substrate, 60 min, nM",
        "cell_endpoint": "HT1080 IDH1 R132C cellular 2-HG inhibition IC50, 48 h, nM",
        "sha256": {str(path.relative_to(ROOT)).replace("\\", "/"): sha(path)
                   for path in [*p1, *p2, Path(__file__), *paths]},
    }
    (OUT / "phase5_idh1_prereveal_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"n_shared": len(inp), "input_sha256": sha(paths[0])}))


if __name__ == "__main__":
    main()
