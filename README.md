# Condition-resolved bioactivity: public figure and baseline data

This is the source-linked reproducibility package for the five data figures in the **Communications Chemistry manuscript draft** on condition-resolved bioactivity contrasts (snapshot: 14 September 2026). It contains our analysis code, public ChEMBL activity/molecule JSON snapshots needed for the ALK and IDH1 selection replays, the exact processed tables plotted, and selected audit tables. It does **not** contain new wet-lab measurements, publisher PDFs, full-text article XML, credentials, or unpublished personal data.

## Reproduce the analyses and figures

Use Python 3.10 with the package versions in `requirements.txt`. From this repository root, run the following in order. The selector receives only the eight pre-reveal input columns; its frozen code is reused unchanged for ALK, IDH1 and the SRPK-covered ALK subset.

```bash
python analysis/phase1_four_arm.py
python analysis/phase1_source_transfer.py

python analysis/phase4_alk_input.py
python analysis/phase4_v0_pair_selector.py --input outputs/phase4_alk_prereveal_input.csv --output outputs/phase4_alk_v0_selection.json
python analysis/phase4_alk_reveal.py
python analysis/phase5_alk_expand.py

python analysis/phase5_idh1_input.py
python analysis/phase4_v0_pair_selector.py --input outputs/phase5_idh1_prereveal_input.csv --output outputs/phase5_idh1_v0_selection.json
python analysis/phase5_idh1_reveal.py

python analysis/phase5_alk_srpk_coverage.py prepare
python analysis/phase4_v0_pair_selector.py --input outputs/phase5_alk_srpk_coverage_input.csv --output outputs/phase5_alk_srpk_coverage_v0_selection.json
python analysis/phase5_alk_srpk_coverage.py reveal

python analysis/commchem_data_figures_v1.py
```

The first two commands rebuild the Fig. 1 model comparison and cross-paper transfer tables from the four matched-pair CSVs. The ALK and IDH1 blocks recreate their pre-reveal inputs, two-compound decisions and withheld-readout tables from the included ChEMBL snapshots. The SRPK block repeats selection in the nine-compound subset with documented assay coverage, then computes the patent-linked comparison. The last command regenerates Figs. 1–5 as PNG and PDF files in `outputs/`; binary image hashes can vary with the local font/rendering environment. The JAK3/KDR/BTK/PARP contextual tables are preserved source-linked derived snapshots, not rebuilt by these commands. They can be traced using the ChEMBL document, assay, molecule and activity identifiers in the CSVs and source links below. The model uses leave-one-molecule-out fitting within each assay pair; its cross-paper JAK3 transfer changes multiple experimental features and is not a mechanism test.

## Figure-to-data map

| Figure | Main input tables in `outputs/` | Source |
| --- | --- | --- |
| Fig. 1: matched contrasts and simple baselines | `jak3_atp_current_20_pairs.csv`, `jak3_2017_4uM_vs_1mM_matched_pairs.csv`, `renin_plasma_vs_trypsin_matched_pairs.csv`, `kdr_htrf_vs_cell_matched_pairs.csv`, `phase1_four_arm_metrics.csv`, `phase1_source_transfer_metrics.csv` | ChEMBL assay metadata in `assay_sources.csv`; JAK3 2019 [DOI:10.1021/acs.jmedchem.8b01823](https://doi.org/10.1021/acs.jmedchem.8b01823), JAK3 2017 [DOI:10.1021/acs.jmedchem.6b01694](https://doi.org/10.1021/acs.jmedchem.6b01694), KDR 2012 [DOI:10.1016/j.bmcl.2012.05.067](https://doi.org/10.1016/j.bmcl.2012.05.067) |
| Fig. 2: ALK mutant-cell spectrum and SRPK1 | `phase5_alk_18x10_condition_matrix.csv`, `phase5_alk_18x10_source_audit.csv`, `phase5_alk_equal_budget_10_condition_comparison.csv`, `phase5_alk_srpk_9_compound_audit.csv`, `phase5_alk_srpk_9_coverage_equal_budget.csv` | ALK paper [DOI:10.1021/acs.jmedchem.5b01136](https://doi.org/10.1021/acs.jmedchem.5b01136), ChEMBL document `CHEMBL3745656`; SRPK1 patent [US11066363B2](https://patents.google.com/patent/US11066363B2/en), Example 14 Table 9, ChEMBL assay `CHEMBL5736981` |
| Fig. 3: IDH1 enzyme–cell and withheld assay | `phase5_idh1_prereveal_input.csv`, `phase5_idh1_all_orthogonal_endpoints.csv`, `phase5_idh1_equal_budget_reveal.csv` | [DOI:10.1021/ml300225h](https://doi.org/10.1021/ml300225h), ChEMBL document `CHEMBL2176947` |
| Fig. 4: JAK3 and KDR cellular boundaries | `jak3_2017_ATP_PBMC_blood_complete_cases.csv`, `kdr_2012_orthogonal_endpoints_joined.csv` | [JAK3 2017](https://doi.org/10.1021/acs.jmedchem.6b01694), [KDR 2012](https://doi.org/10.1016/j.bmcl.2012.05.067) |
| Fig. 5: BTK and PARP published controls | `btk_2025_corrected_15_binding_pairs.csv`, `btk_crosspaper_8_kinetic_rank.csv`, `parp_2012_cell_catalytic_vs_trapping.csv` | [BTK binding](https://doi.org/10.1021/acsptsci.4c00540), [BTK correction](https://doi.org/10.1021/acsptsci.6c00190), [kinetic ordering](https://doi.org/10.1021/acsptsci.5c00412), [PARP 2012](https://doi.org/10.1158/0008-5472.CAN-12-2753) |

All potency columns labelled `nM` are nanomolar source-reported values. A `relation` column (such as `>` or `<`) marks a bound and must not be treated as an exact measurement. The SRPK1 compound-17 bound is corrected against the patent as `>10^3 nM`; ChEMBL's `>103 nM` is an apparent transcription of the exponent. In the BTK table, mechanism labels follow the published correction. The additional `parp_2024_context_gate.csv` supports discussion of genetic context from [Petropoulos *et al.*](https://doi.org/10.1038/s41586-024-07217-2); it is a qualitative transcription, not a plotted numeric response matrix.

## Interpretation and provenance

These are retrospective analyses of published records. ALK/SRPK1 measurements are from related inventors and are not independent replication; the IDH1 withheld readout comes from the same paper and chemical series as the selection inputs; the KDR polyploidy endpoint is not KDR-specific. Raw activity-level uncertainty was unavailable for some source values. See the `*_manifest.json` files for parameter choices and further source/record notes. The package's plotted values and source identifiers are provided to enable review of those boundaries rather than to imply prospective validation.

The included `sources/` JSON files are snapshots of public [ChEMBL API](https://www.ebi.ac.uk/chembl/api/data/) activity and molecule records. They can be retrieved by assay and molecule IDs in the scripts/CSVs; the relevant publications and patent are linked above. Any reuse of underlying source content is subject to the original providers' terms. Our scripts are shared for inspection and reproduction without an added license grant in this snapshot.
