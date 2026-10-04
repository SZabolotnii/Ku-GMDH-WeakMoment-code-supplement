# Validation-Gated Weak-Moment GMDH — code supplement

Code, derived results and generated tables/figures for the manuscript

> S. Zabolotnii. *Validation-Gated Weak-Moment GMDH for Heavy-Tailed Regression.*

The repository contains the `gmdhpmm` R package (GMDH tournament with PMM inner
estimators, including the weak-windowed WPMM2/WPMM3 estimators and the
validation-gated dispatch `auto-valgate`), the experiment scripts that produced
every number in the paper, the result CSVs those scripts wrote, and the scripts
that turn the CSVs into the paper's tables and figures.

## Layout

The directory names mirror the research repository the results were produced
in, so every script runs unchanged. Scripts find the package by walking up from
the working directory to `paper-1-gmdh-pmm/code/DESCRIPTION`.

```text
.
|-- paper-1-gmdh-pmm/code/            R package `gmdhpmm` (DESCRIPTION, R/, tests/)
|   `-- R/
|       |-- weak_pmm.R                WPMM2, WPMM3 (Section 3.1)
|       |-- cumulants.R               raw and weak cumulant diagnostics
|       |-- dispatch.R                classical and weak dispatch rules
|       |-- inner_estimate.R          entry point: force_method = "WPMM2", "WPMM3",
|       |                             "auto-weak", "auto-valgate" (Section 3.2), ...
|       |-- kg2.R, gmdh.R             KG-2 partial model, MIA tournament
|       `-- external_criterion.R, control.R
|-- paper-4-weak-moment-gmdh/
|   |-- code/experiments/             experiment, evidence and table/figure scripts
|   |-- results/                      derived result CSVs (the paper's evidence)
|   `-- latex/{tables,figures}/       tables and figures as generated from results/
`-- shared/datasets/
    |-- prepare_*.R, screen_cumulants.R   dataset preparation and cumulant pre-screen
    `-- external/
        |-- external_candidates.csv       dataset manifest (path, formula, source, licence)
        `-- processed/                    prepared CSVs (see "Data")
```

## Requirements

- R >= 4.1 (results produced with R 4.5.3).
- R packages: `EstemPMM` (PMM2/PMM3 solvers; CRAN), `data.table`, `pkgload`,
  `testthat`; for the cross-sectional datasets `MASS`, `ISLR2`, `boot`,
  `insuranceData`.
- Python >= 3.10 with `numpy`, `pandas`, `scikit-learn` for the tree-ensemble
  baselines only.

## Quick check

```bash
cd paper-1-gmdh-pmm/code
Rscript -e 'testthat::test_local(".")'          # 118 tests, 0 failures

cd ../../paper-4-weak-moment-gmdh/code
Rscript experiments/build_final_evidence_and_meta.R 10000
Rscript experiments/make_latex_tables.R
Rscript experiments/make_latex_figures.R
```

The last three commands rebuild the final evidence layer, Tables 2–4 and
Figures 1–5 from the shipped CSVs (Table 1, the evaluation design, is written in
the manuscript and has no producer). On the release commit the regenerated
`final_claim_*.csv` files and the three table `.tex` files are byte-identical to
the ones in the repository (`tables_manifest.csv` differs only in its timestamp).

## Map from the paper to the code

Run every script from `paper-4-weak-moment-gmdh/code/`. Scripts taking `R` use
`R = 30` replicate splits by default; an optional second argument appends a
suffix to the output names, so a smoke run does not overwrite the archived CSVs
(`Rscript experiments/run_p4_crosssec_modskew.R 1 _smoke`).

| Paper item | Result CSV(s) in `results/` | Producer script(s) |
|---|---|---|
| Cross-sectional rows (Table 2, Fig. 1) | `crosssec_modskew_*.csv` | `run_p4_crosssec_modskew.R` |
| Insurance severity, freMTPL2 heavy-tail rescue (Table 2, Fig. 2, Section 5.3) | `insurance_severity_*.csv` | `run_insurance_severity.R` |
| Soft-sensor rows (Table 2) | `p4_softsensor_honest_*.csv` | `run_p4_softsensor_honest.R` |
| Platykurtic family | `platykurtic_family_*.csv` | `run_platykurtic_family.R` |
| Tree-ensemble baselines (Table 2) | `ml_baselines_crossdomain*.csv` | `export_crossdomain_claim_rows_for_ml.R`, then `ml_baselines_crossdomain.py --trees 300` |
| Final claim rows, split deltas, meta checkpoint (Tables 2 and 4, Figs. 1 and 4) | `final_claim_evidence_table.csv`, `final_claim_split_deltas.csv`, `final_claim_meta_analysis.csv` | `build_final_evidence_and_meta.R [Bboot]` |
| Grouped inference (Table 4) | `grouped_inference_*.csv` | `run_grouped_inference_hardening.R [Bboot]` |
| Validation-gated selector (Table 3, Fig. 3) | `dispatch_retune.csv` | `run_dispatch_retune.R` |
| Feature-richness boundary (Fig. 5) | `feature_count_sweep.csv` | `run_feature_count_sweep.R` |
| Drilling significance analysis and aggregation sensitivity | `pooled_significance.csv`, `aggregation_sensitivity.csv`, `honest_blocked_cv_*.csv`, `blockboot_ci.csv` | `run_pooled_significance.R`, `run_aggregation_sensitivity.R`, `run_honest_blocked_cv*.R`, `run_blockboot_ci.R` |
| Leakage audit, stacking, drilling tree baselines | `leakage_audit_vibration.csv`, `stacking_blocked.csv`, `ml_baselines_blocked.csv` | `run_leakage_audit.R`, `run_stacking_blocked_cv.R`, `export_blocked_cells_for_ml.R` + `ml_baselines_blocked.py` |
| Tables 2–4, Figures 1–5 | `latex/tables/`, `latex/figures/` | `make_latex_tables.R`, `make_latex_figures.R` |

The remaining `run_weak_*.R` scripts are the development studies behind the
estimator design (proof of concept, regime map, ablation, WPMM3 bandwidth probe,
cascade checks); their CSVs are archived in `results/` as well.

## Data

| Dataset | Used for | In this repository |
|---|---|---|
| freMTPL2 severity (`CASdatasets`, GPL >= 2) | insurance rows | `processed/fremtpl2_severity*.csv`; rebuild with `prepare_fremtpl2.R` |
| UCI Gas Turbine CO/NOx (doi:10.24432/C5WC95, CC BY 4.0) | soft-sensor rows | `processed/gas_turbine_*.csv` |
| SRU soft sensor (Mendeley Data doi:10.17632/kcpnnrn67p.1, CC BY 4.0) | soft-sensor rows | `processed/sru*.csv` |
| UCI Concrete (doi:10.24432/C5PK67, CC BY 4.0) | cross-sectional rows | `processed/concrete.csv` |
| `MASS`, `ISLR2`, `boot`, `insuranceData` R datasets | cross-sectional and insurance rows | loaded from the installed packages |
| Utah FORGE drilling and downhole-vibration logs; Equinor Volve 15/9-F-15 | drilling significance analysis, feature-count sweep | **not redistributed** — download from the Utah FORGE / OpenEI Geothermal Data Repository and the Equinor Volve data village; the loaders in `paper-1-gmdh-pmm/code/R/data_drilling.R` document the expected columns |

The derived results of the drilling scripts are archived in `results/`, so every
table and figure can be rebuilt without the raw drilling files.

## Licence and citation

GNU General Public License v3.0 (see `LICENSE`), the licence of the `gmdhpmm`
package. Dataset files remain under their own licences listed above.

If you use this code, please cite the paper (reference to be added on
publication) and this repository; citation metadata is in `CITATION.cff`.
