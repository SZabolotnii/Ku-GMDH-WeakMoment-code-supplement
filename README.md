# Validation-Gated Weak-Moment GMDH — code supplement

Code, derived results and generated tables/figures for the manuscript

> S. Zabolotnii. *Validation-Gated Weak-Moment Group Method of Data Handling for
> Heavy-Tailed Regression.*

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
|   |-- blocked-power-prereg-2026-10-04.md   decision rule of the ten-fold blocked
|   |                                        check, written before the run
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
- Python >= 3.10 with `numpy`, `pandas`, `scikit-learn`, `joblib` for the
  tree-ensemble baselines only.

## Quick check

```bash
cd paper-1-gmdh-pmm/code
Rscript -e 'testthat::test_local(".")'          # 122 tests, 0 failures

cd ../../paper-4-weak-moment-gmdh/code
Rscript experiments/build_final_evidence_and_meta.R 10000
Rscript experiments/make_blocked_power_summary.R
Rscript experiments/make_claimrow_robustness_summary.R
Rscript experiments/make_latex_tables.R
Rscript experiments/make_latex_figures.R
```

These commands read only the shipped CSVs. They rebuild the per-row evidence
layer, the summaries of the robustness, large-claim, tuned-tree and ten-fold
blocked checks, Tables 2–4 and Figures 1–2 of the paper and Figure S1 of its
supplemental material (Table 1, the evaluation design, is written in the
manuscript and has no producer). On the release commit the regenerated
`final_claim_*.csv`, `blocked_power_summary.csv` and `claimrow_*.csv` files and
the three table `.tex` files are byte-identical to the ones in the repository
(`tables_manifest.csv` differs only in its timestamp, the figure PDFs only in
their creation date).

## Map from the paper to the code

Section numbers refer to the Journal of Applied Statistics version of the paper;
"S" sections are in its supplemental material. Run every script from
`paper-4-weak-moment-gmdh/code/`. Scripts taking `R` use `R = 30` replicate
splits by default; an optional second argument appends a suffix to the output
names, so a smoke run does not overwrite the archived CSVs
(`Rscript experiments/run_p4_crosssec_modskew.R 1 _smoke`). The same holds for
`run_claimrow_robustness.R` and `run_blocked_power.R`, whose first argument is
the number of cores.

| Paper item | Result CSV(s) in `results/` | Producer script(s) |
|---|---|---|
| freMTPL2 heavy-tail rescue (Section 5.1, Fig. 1) | `insurance_severity_*.csv` | `run_insurance_severity.R` |
| Insurance claim severity: gate vs best robust method and vs best tuned tree, MAE, trimmed RMSE, RMSE and large-claim RMSE (Section 5.2, Table 2) | `claimrow_robustness_raw_tail.csv`, `ml_tuned_trees_long_tail.csv`, `claimrow_tail_summary.csv` | `run_claimrow_robustness.R 6 _tail fremtpl2_severity_raw,fremtpl2_severity_log,insurance_autobi_loss,insurance_autoclaims_paid tail`; `export_allrows_for_ml.R`, then `python experiments/ml_tuned_trees.py --tail --suffix _tail --datasets fremtpl2_severity_raw,fremtpl2_severity_log,insurance_autobi_loss,insurance_autoclaims_paid`; `make_claimrow_robustness_summary.R` checks that every refit reproduces the archived trimmed RMSE (it does, to the last digit) |
| Archived per-split fits of the fifteen datasets (seed 0; inputs to Tables 2–3) | `crosssec_modskew_*.csv`, `insurance_severity_*.csv`, `p4_softsensor_honest_*.csv`, `final_claim_*.csv` | `run_p4_crosssec_modskew.R`, `run_insurance_severity.R`, `run_p4_softsensor_honest.R`, then `build_final_evidence_and_meta.R [Bboot]` |
| Validation gate on all fifteen datasets, two extra seeds per method (Sections 5–6, Tables 2–3, Fig. 2) | `claimrow_robustness_raw.csv` (eight datasets), `claimrow_robustness_raw_pool.csv` (the other seven) | `run_claimrow_robustness.R [cores] [suffix] [datasets] [arms]`; the pool file comes from `run_claimrow_robustness.R 9 _pool mass_boston_medv,islr_wage,airquality_ozone,mass_cars93_mpg,insurance_autoclaims_paid,fremtpl2_severity_log,sru_y1_static valgate,seeds,repro`; the Boston and Wage rows were rerun on 2026-10-05 without the race covariates (`black`, `race`) and appended, so they come last |
| Tuned tree ensembles (Sections 5–6, Tables 2–3, Fig. 2) | `ml_tuned_trees_long.csv` (tuned and default fits on every outer split) | `export_allrows_for_ml.R`, then `python experiments/ml_tuned_trees.py` |
| Table 3 and Fig. 2 inputs: gate vs best robust method and vs best tuned tree, untrimmed errors, row-level signed-rank test, gate vs hindsight-best weak estimator | `claimrow_candidate_pool.csv`, `claimrow_tuned_trees_summary.csv`, `claimrow_gate_splits.csv`, `claimrow_untrimmed_summary.csv`, `claimrow_rowlevel_summary.csv`, `claimrow_valgate_summary.csv` | `make_claimrow_robustness_summary.R` |
| Cumulant-keyed vs validation gate on credit balance, twelve splits (Section 6) | `dispatch_retune.csv`, `dispatch_retune_splits.csv` | `run_dispatch_retune.R` |
| Protocol lesson: five folds, seeds 0–2, purge gaps 50 and 100 (Section 7, Table 4) | `claimrow_seeds_summary.csv`, `claimrow_gap_summary.csv`, `claimrow_protocol_lesson.csv` | `run_claimrow_robustness.R`, then `make_claimrow_robustness_summary.R` |
| Ten-fold, five-seed blocked check with a prespecified rule (Section 7, Table 4 last column) | `blocked_power_raw.csv`, `blocked_power_summary.csv` | `run_blocked_power.R`, then `make_blocked_power_summary.R`; the decision rule, written before the run and not registered externally, is in `paper-4-weak-moment-gmdh/blocked-power-prereg-2026-10-04.md` |
| Drilling case study: leakage, aggregation, block bootstrap (Section S1) | `leakage_audit_vibration.csv`, `pooled_significance.csv`, `aggregation_sensitivity.csv`, `honest_blocked_cv_*.csv`, `blockboot_ci.csv` | `run_leakage_audit.R`, `run_pooled_significance.R`, `run_aggregation_sensitivity.R`, `run_honest_blocked_cv*.R`, `run_blockboot_ci.R` |
| Skew-tilt ablation, proof of concept, regime map (Section S2) | `weak_ablation_raw.csv`, `weak_pmm_poc_raw.csv`, `weak_regime_map_raw.csv` | `run_weak_ablation.R`, `run_weak_pmm_poc.R`, `run_weak_regime_map.R` |
| Gate design checks: interleaved inner folds, library without Huber/LAD (Section S3) | `claimrow_ablation_summary.csv` (arms `rinner`, `norobust`) | `run_claimrow_robustness.R`, then `make_claimrow_robustness_summary.R`; interleaved folds are `valgate_inner = "random"` in `gmdh_pmm_control()` |
| Gate design checks: window width 1.5 and 4 times the robust scale (Section S3) | `claimrow_sigma_summary.csv` | as above |
| Convex stack vs discrete selector (Section S3) | `stacking_blocked.csv` | `run_stacking_blocked_cv.R` |
| Feature-richness boundary (Section S4, Fig. S1) | `feature_count_sweep.csv` | `run_feature_count_sweep.R` |
| Tables 2–4, Figures 1–2 and S1 | `latex/tables/`, `latex/figures/` | `make_latex_tables.R`, `make_latex_figures.R` |

The remaining scripts are archived with their CSVs: the other development
studies behind the estimator design (WPMM3 bandwidth probe, cascade and dispatch
checks), the platykurtic family, the untuned tree baselines
(`ml_baselines_crossdomain.py`, whose `ml_baselines_crossdomain_long.csv`
`make_claimrow_robustness_summary.R` still reads, and `ml_baselines_blocked.py`)
and the grouped and meta-analytic inference of an earlier version of the
manuscript (`run_grouped_inference_hardening.R`, `final_claim_meta_analysis.csv`).

## Data

| Dataset | Used for | In this repository |
|---|---|---|
| freMTPL2 severity (`CASdatasets`, GPL >= 2) | insurance rows | `processed/fremtpl2_severity*.csv`; rebuild with `prepare_fremtpl2.R` |
| UCI Gas Turbine CO/NOx (doi:10.24432/C5WC95, CC BY 4.0) | soft-sensor rows | `processed/gas_turbine_*.csv` |
| SRU soft sensor (Mendeley Data doi:10.17632/kcpnnrn67p.1, CC BY 4.0) | soft-sensor rows | `processed/sru*.csv` |
| UCI Concrete (doi:10.24432/C5PK67, CC BY 4.0) | cross-sectional rows | `processed/concrete.csv` |
| `MASS`, `ISLR2`, `boot`, `insuranceData` R datasets; `airquality` from base R | cross-sectional and insurance rows | loaded from the installed packages |
| Utah FORGE drilling and downhole-vibration logs; Equinor Volve 15/9-F-15 | drilling significance analysis, feature-count sweep | **not redistributed** — download from the Utah FORGE / OpenEI Geothermal Data Repository and the Equinor Volve data village; the loaders in `paper-1-gmdh-pmm/code/R/data_drilling.R` document the expected columns |

The derived results of the drilling scripts are archived in `results/`, so every
table and figure can be rebuilt without the raw drilling files.

## Licence and citation

GNU General Public License v3.0 (see `LICENSE`), the licence of the `gmdhpmm`
package. Dataset files remain under their own licences listed above.

If you use this code, please cite the paper (reference to be added on
publication) and this repository; citation metadata is in `CITATION.cff`.
