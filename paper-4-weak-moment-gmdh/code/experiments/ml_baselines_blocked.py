#!/usr/bin/env python
"""
EXPERIMENT 1 -- modern-ML baselines under IDENTICAL blocked folds.

Referee objection mirrored: "no comparison to RF/GBM". On the SAME blocked folds
and the SAME 2 raw features as the R honest-CV harness (run_honest_blocked_cv.R /
run_honest_blocked_cv_rop.R), fit RandomForest, HistGradientBoosting, ExtraTrees,
and -- as a PARITY CHECK -- sklearn LinearRegression (== R LSE). Per cell/band we
report the BLOCKED median test trimmed-RMSE for each model.

Splits are byte-identical to R: the R exporter dumped the time/depth-ordered
(f1,f2,y) per cell/band; here we replicate the deterministic blocked-5fold
arithmetic (same gap, folds=5, same train-purge, same guards). Cap=6000 was
already applied in the export (contiguous head), so no cap logic needed here.

trmse() reproduces the R definition EXACTLY:
    q <- quantile(|e|, 0.90, type=7); sqrt(mean(e[|e|<=q]^2))
np.quantile default == R type 7 (linear interpolation), so they coincide.

PARITY MODEL: the R inner estimator fits the order-2 Kolmogorov-Gabor (KG-2)
partial polynomial in the two features:
    P = th0 + th1 f1 + th2 f2 + th3 f1*f2 + th4 f1^2 + th5 f2^2.
The R "LSE" column (bl_LSE) is ridge-OLS (lambda=1e-8 ~ OLS) on THIS 6-term
design (see kg2.R). So the sklearn LSE-parity model must be LinearRegression on
the SAME quadratic design -- otherwise parity will (correctly) fail. The TREE
ensembles take the 2 RAW features (they synthesize their own nonlinearity); this
is the standard fair comparison the referee asks for.

KEY QUESTION: does the best tree ensemble beat best-weak on the cells where weak
wins (gamma4 > ~5)? If RF/GBM dominate there, the contribution is DEFLATED.
"""
import os
import sys
import csv
import numpy as np
import pandas as pd
from sklearn.ensemble import (RandomForestRegressor,
                              HistGradientBoostingRegressor,
                              ExtraTreesRegressor)
from sklearn.linear_model import LinearRegression

IN = "/tmp/p4_ml"
PKG_RESULTS = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "results"))
VIB_CSV = os.path.join(PKG_RESULTS, "honest_blocked_cv_vibration.csv")
ROP_CSV = os.path.join(PKG_RESULTS, "honest_blocked_cv_rop.csv")
OUT_CSV = os.path.join(PKG_RESULTS, "ml_baselines_blocked.csv")

FOLDS = 5


def trmse(e, p=0.90):
    e = np.asarray(e, dtype=float)
    q = np.quantile(np.abs(e), p)            # numpy default == R type 7
    keep = np.abs(e) <= q
    return float(np.sqrt(np.mean(e[keep] ** 2)))


def blocked_indices(n, gap):
    """Yield (train_idx, test_idx) 0-based, mirroring the R 1-based arithmetic."""
    fb = n // FOLDS                          # floor(n/folds)
    out = []
    for k in range(1, FOLDS + 1):
        a = (k - 1) * fb + 1                 # 1-based start
        z = n if k == FOLDS else k * fb      # 1-based end (inclusive)
        tei = np.arange(a, z + 1)            # 1-based test indices
        purge_lo = max(1, a - gap)
        purge_hi = min(n, z + gap)
        purged = set(range(purge_lo, purge_hi + 1))
        tri = np.array([i for i in range(1, n + 1) if i not in purged])
        out.append((tri - 1, tei - 1))       # to 0-based
    return out


def kg2_design(X):
    """6-term KG-2 design WITHOUT intercept (LinearRegression adds intercept).
    Columns: f1, f2, f1*f2, f1^2, f2^2 -- matches kg2.R canonical order."""
    f1 = X[:, 0]
    f2 = X[:, 1]
    return np.column_stack([f1, f2, f1 * f2, f1 ** 2, f2 ** 2])


def fit_models(Xtr, ytr, Xte):
    """Trees use the 2 RAW features; LSE-parity uses the KG-2 quadratic design."""
    preds = {}
    rf = RandomForestRegressor(n_estimators=300, random_state=0, n_jobs=-1)
    rf.fit(Xtr, ytr)
    preds["rf"] = rf.predict(Xte)
    hgb = HistGradientBoostingRegressor(random_state=0)
    hgb.fit(Xtr, ytr)
    preds["hgb"] = hgb.predict(Xte)
    et = ExtraTreesRegressor(n_estimators=300, random_state=0, n_jobs=-1)
    et.fit(Xtr, ytr)
    preds["extra"] = et.predict(Xte)
    lin = LinearRegression()                 # LSE parity on KG-2 quadratic design
    lin.fit(kg2_design(Xtr), ytr)
    preds["lin"] = lin.predict(kg2_design(Xte))
    return preds


def eval_cell(path, gap):
    df = pd.read_csv(path)
    f1 = df["f1"].to_numpy(float)
    f2 = df["f2"].to_numpy(float)
    y = df["y"].to_numpy(float)
    X = np.column_stack([f1, f2])
    n = len(y)
    per_model = {m: [] for m in ("rf", "hgb", "extra", "lin")}
    for tri, tei in blocked_indices(n, gap):
        # mirror R guards
        if gap == 50:
            if len(tri) < 100 or len(tei) < 40:
                continue
        else:  # gap == 30 (ROP)
            if len(tri) < 80 or len(tei) < 30:
                continue
        preds = fit_models(X[tri], y[tri], X[tei])
        for m, p in preds.items():
            per_model[m].append(trmse(y[tei] - p))
    return {m: (float(np.median(v)) if v else float("nan"))
            for m, v in per_model.items()}, n


def load_ref():
    """Best-weak reference + R LSE column from the committed honest CSVs."""
    ref = {}
    vib = pd.read_csv(VIB_CSV)
    for _, row in vib.iterrows():
        name = "vib__%s__%s__%s" % (row["run"], row["sensor"], row["chan"])
        bw = min(row["bl_WPMM2"], row["bl_WPMM3"])
        ref[name] = dict(gamma4=row["g4"], r_lse=row["bl_LSE"], best_weak=bw,
                         bw_name=("WPMM2" if row["bl_WPMM2"] <= row["bl_WPMM3"]
                                  else "WPMM3"))
    rop = pd.read_csv(ROP_CSV)
    for _, row in rop.iterrows():
        name = row["dataset"]
        bw = min(row["bl_WPMM2"], row["bl_WPMM3"])
        ref[name] = dict(gamma4=row["g4"], r_lse=row["bl_LSE"], best_weak=bw,
                         bw_name=("WPMM2" if row["bl_WPMM2"] <= row["bl_WPMM3"]
                                  else "WPMM3"))
    return ref


def main():
    man = pd.read_csv(os.path.join(IN, "manifest.csv"))
    ref = load_ref()
    rows = []
    parity_fails = []
    for _, mrow in man.iterrows():
        dataset = mrow["dataset"]
        cell = mrow["cell"]
        gap = int(mrow["gap"])
        res, n = eval_cell(mrow["file"], gap)
        r = ref.get(cell, {})
        r_lse = r.get("r_lse", float("nan"))
        gamma4 = r.get("gamma4", mrow["gamma4"])
        best_weak = r.get("best_weak", float("nan"))
        # parity check: sklearn LinearRegression blocked median vs R LSE column
        parity_ok = (not np.isnan(r_lse)
                     and abs(res["lin"] - r_lse) <= 5e-3 * max(1.0, abs(r_lse)))
        if not parity_ok and not np.isnan(r_lse):
            parity_fails.append((cell, res["lin"], r_lse))
        rows.append(dict(dataset=dataset, cell=cell, n=n,
                         gamma4=round(gamma4, 4),
                         rf=round(res["rf"], 6),
                         hgb=round(res["hgb"], 6),
                         extra=round(res["extra"], 6),
                         lin=round(res["lin"], 6),
                         r_lse=round(r_lse, 6) if not np.isnan(r_lse) else "",
                         best_weak=round(best_weak, 6)
                         if not np.isnan(best_weak) else "",
                         bw_name=r.get("bw_name", ""),
                         parity_ok=parity_ok))
    # write CSV
    cols = ["dataset", "cell", "n", "gamma4", "rf", "hgb", "extra", "lin",
            "r_lse", "best_weak", "bw_name", "parity_ok"]
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for rrow in rows:
            w.writerow(rrow)

    # ---------- console summary ----------
    print("\n=== PARITY CHECK (sklearn LinearRegression vs R LSE column) ===")
    print("%-40s %10s %10s %8s" % ("cell", "py_lin", "R_LSE", "absdiff"))
    max_absdiff = 0.0
    for rrow in rows:
        if rrow["r_lse"] == "":
            continue
        d = abs(float(rrow["lin"]) - float(rrow["r_lse"]))
        max_absdiff = max(max_absdiff, d)
        print("%-40s %10.5f %10.5f %8.5f%s"
              % (rrow["cell"], rrow["lin"], rrow["r_lse"], d,
                 "  <-- FAIL" if not rrow["parity_ok"] else ""))
    print("max |py_lin - R_LSE| = %.5f over %d cells; parity %s"
          % (max_absdiff, len(rows),
             "PASS (all within 5e-3 rel)"
             if not parity_fails else "FAIL on %d cells" % len(parity_fails)))

    print("\n=== BEST-TREE vs BEST-WEAK per cell ===")
    print("%-40s %6s %9s %9s %9s | %9s %9s | %s"
          % ("cell", "g4", "rf", "hgb", "extra", "best_tree", "best_weak",
             "winner"))
    weak_cells = []   # gamma4 > ~5 cells: where weak is claimed to win
    tree_beats_weak_on_weakcells = 0
    n_weakcells = 0
    for rrow in rows:
        if rrow["best_weak"] == "":
            continue
        bt = min(rrow["rf"], rrow["hgb"], rrow["extra"])
        bw = float(rrow["best_weak"])
        winner = "tree" if bt < bw else "weak"
        g4 = float(rrow["gamma4"])
        flag = ""
        if g4 > 5.0:
            n_weakcells += 1
            flag = " *weak-regime"
            if bt < bw:
                tree_beats_weak_on_weakcells += 1
        print("%-40s %6.2f %9.5f %9.5f %9.5f | %9.5f %9.5f | %s%s"
              % (rrow["cell"], g4, rrow["rf"], rrow["hgb"], rrow["extra"],
                 bt, bw, winner, flag))
        weak_cells.append((rrow["cell"], g4, bt, bw, winner))

    print("\n=== KEY QUESTION (gamma4 > 5 cells: where weak is claimed best) ===")
    print("best tree beats best-weak on %d / %d high-kurtosis cells"
          % (tree_beats_weak_on_weakcells, n_weakcells))
    n_all = len([r for r in rows if r["best_weak"] != ""])
    tree_wins_all = sum(
        1 for r in rows if r["best_weak"] != ""
        and min(r["rf"], r["hgb"], r["extra"]) < float(r["best_weak"]))
    print("best tree beats best-weak on %d / %d cells overall" % (tree_wins_all, n_all))
    print("\nSaved %s" % OUT_CSV)


if __name__ == "__main__":
    main()
