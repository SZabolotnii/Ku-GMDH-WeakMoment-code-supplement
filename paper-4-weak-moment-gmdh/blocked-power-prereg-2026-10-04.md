# Blocked-row power check — pre-registration (2026-10-04)

Written before `run_blocked_power.R` was run. The rule below decides how the four blocked
claim rows are reported; it is not to be changed after the results are seen.

## Why

`run_claimrow_robustness.R` (2026-10-04) showed that the Table-2 deltas of the three blocked
"strong positive" rows depend on the tournament seed (SRU y1: −15.6 % at seed 0, +0.8 % at
seeds 1 and 2; CO: −10.1 %, +2.4 %, +1.9 %; NOx: −15.0 %, +5.2 %, −15.8 %) and on the purge
gap. With five folds per row the median of five values cannot separate an effect from noise.

## Design

- Rows: `sru_y1_dynamic`, `gas_turbine_co_2015_raw`, `gas_turbine_nox_2015_raw`,
  `sru_y2_static` (control). Same data, the first 3000 observations, same standardisation.
- Outer splits: **10** contiguous blocked folds (was 5), purge gap **20** (primary) and
  **50** (sensitivity).
- Tournament seeds: **5** replicates per fold; each method's per-fold error is the mean over
  the five seeds.
- Methods: the eight Table-2 methods and `auto-valgate` (default library, blocked inner folds,
  sigma_mult 2.5).

## Comparisons (fixed in advance)

For each row and gap:
1. **Gate:** `auto-valgate` vs the row's Table-2 robust method
   (`final_claim_evidence_table.csv`, column `best_robust_stat_method`).
2. **Fixed weak method:** the row's Table-2 weak method (`best_weak_method`) vs the same
   robust method. The method is the one already named in Table 2, not re-selected.
Statistic: the paired median over folds of (weak or gate TRMSE / robust TRMSE − 1), and the
number of folds won.

## Decision rule

A row **supports a gain** if, at gap 20, the gate's paired median is ≤ −3 % **and** it wins
at least 7 of 10 folds, **and** at gap 50 the paired median is still negative. A row
**shows a loss** if the paired median is ≥ +3 % with at most 3 of 10 folds won. Anything
else is reported as **no detectable difference**. The same rule is applied to the fixed
weak method and reported alongside, but the gate decides how the row is described,
because it is the deployable route.

The `sru_y2_static` control is expected to show no detectable difference; if it shows a
gain, the rule is too loose and that is reported as such.

## Outcome (appended after the run, 2026-10-04; the rule above was not changed)

3600 fits, none failed (`results/blocked_power_raw.csv`, summary in
`results/blocked_power_summary.csv`, built by `make_blocked_power_summary.R`).

| Row | Gate verdict | Fixed Table-2 weak method verdict |
|---|---|---|
| sru_y1_dynamic | no detectable difference (−3.1 %, 5/10; gap 50 +0.1 %) | no detectable difference (−4.7 %, 7/10; gap 50 +6.1 %) |
| gas_turbine_co_2015_raw | no detectable difference (+0.9 %, 4/10) | no detectable difference (−0.7 %, 6/10) |
| gas_turbine_nox_2015_raw | no detectable difference (−2.6 %, 6/10) | supports a gain, at the margin (−3.2 %, 7/10; gap 50 −1.5 %) |
| sru_y2_static (control) | no detectable difference (+2.9 %, 4/10) | **supports a gain** (−7.8 %, 7/10; gap 50 −31.5 %) |

The gate, which decides how a row is described, finds no detectable difference on any of
the four blocked rows, and the control behaves as expected. The fixed-method comparison
fails its control: SRU y2 was a tie at five folds, yet passes the rule at ten. The likely
reason is that the comparator is the robust method that happened to be best at five folds
(LSE there), which need not be the best robust method at ten folds. As the pre-registration
says, a control that passes means the rule is too loose for that comparison, so the
fixed-method verdicts, including the marginal NOx gain, are not used as evidence.
