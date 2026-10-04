#!/usr/bin/env python3
"""
Uniform cross-domain ML baseline for Paper 4 claim rows.

Inputs are exported by export_crossdomain_claim_rows_for_ml.R so the Python side
uses the same model matrices and train/test splits as the rebuilt R evidence.
Trees use the full exported feature matrix; the linear model is only a reference
baseline, not a parity check against GMDH-LSE.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import LinearRegression


def find_root() -> Path:
    cur = Path.cwd().resolve()
    for p in [cur, *cur.parents]:
        if (p / "paper-4-weak-moment-gmdh").exists() and (p / "paper-1-gmdh-pmm").exists():
            return p
    raise RuntimeError("PMM-GMDH root not found")


ROOT = find_root()
RESULTS = ROOT / "paper-4-weak-moment-gmdh" / "results"
DEFAULT_IN = Path("/tmp/p4_crossdomain_ml")


def trmse(e: np.ndarray, p: float = 0.90) -> float:
    e = np.asarray(e, dtype=float)
    e = e[np.isfinite(e)]
    if e.size == 0:
        return float("nan")
    q = np.quantile(np.abs(e), p)
    keep = np.abs(e) <= q
    return float(np.sqrt(np.mean(e[keep] ** 2)))


def rmse(e: np.ndarray) -> float:
    e = np.asarray(e, dtype=float)
    e = e[np.isfinite(e)]
    if e.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(e ** 2)))


def parse_idx(s: str) -> np.ndarray:
    if not isinstance(s, str) or not s.strip():
        return np.array([], dtype=int)
    return np.fromstring(s, sep=" ", dtype=int)


def fit_models(xtr: np.ndarray, ytr: np.ndarray, xte: np.ndarray, seed: int, trees: int) -> dict[str, np.ndarray]:
    models = {
        "rf": RandomForestRegressor(n_estimators=trees, random_state=seed, n_jobs=-1),
        "hgb": HistGradientBoostingRegressor(random_state=seed),
        "extra": ExtraTreesRegressor(n_estimators=trees, random_state=seed, n_jobs=-1),
        "linear": LinearRegression(),
    }
    preds: dict[str, np.ndarray] = {}
    for name, model in models.items():
        model.fit(xtr, ytr)
        preds[name] = model.predict(xte)
    return preds


def load_reference() -> dict[tuple[str, str], dict[str, float | str]]:
    refs: dict[tuple[str, str], dict[str, float | str]] = {}

    soft = pd.read_csv(RESULTS / "p4_softsensor_honest_compare.csv")
    for _, r in soft.iterrows():
        refs[(r["dataset"], r["protocol"])] = {
            "best_weak_trmse": float(r["pmm_trmse"]),
            "best_weak_mae": float(r["pmm_mae"]),
            "best_robust_trmse": float(r["rob_trmse"]),
            "best_robust_mae": float(r["rob_mae"]),
            "lse_trmse": float(r["lse_trmse"]),
            "best_weak_method": str(r["best_pmm_trmse_method"]),
            "best_robust_method": str(r["best_rob_trmse_method"]),
        }

    cross = pd.read_csv(RESULTS / "crosssec_modskew_comparison.csv")
    for _, r in cross.iterrows():
        refs[(r["candidate"], "random")] = {
            "best_weak_trmse": float(r["best_pmm_trmse"]),
            "best_weak_mae": np.nan,
            "best_robust_trmse": float(r["best_rob_trmse"]),
            "best_robust_mae": np.nan,
            "lse_trmse": float(r["lse_trmse"]),
            "best_weak_method": str(r["best_pmm_method_tr"]),
            "best_robust_method": str(r["best_robust_tr"]),
        }

    ins = pd.read_csv(RESULTS / "insurance_severity_family.csv")
    for _, r in ins.iterrows():
        refs[(r["candidate"], "random")] = {
            "best_weak_trmse": float(r["pmm_trmse"]),
            "best_weak_mae": float(r["pmm_mae"]),
            "best_robust_trmse": float(r["robust_trmse"]),
            "best_robust_mae": float(r["robust_mae"]),
            "lse_trmse": float(r["lse_trmse"]),
            "best_weak_method": str(r["best_pmm_tr"]),
            "best_robust_method": str(r["best_robust_tr"]),
        }

    return refs


def summarize(long_df: pd.DataFrame, manifest: pd.DataFrame) -> pd.DataFrame:
    refs = load_reference()
    rows = []
    for _, m in manifest.iterrows():
        ds = str(m["dataset"])
        proto = str(m["protocol"])
        sub = long_df[long_df["dataset"] == ds]
        by_model = {}
        for model, g in sub.groupby("model"):
            by_model[model] = {
                "trmse": float(np.median(g["trmse"])),
                "mae": float(np.median(g["mae"])),
                "rmse": float(np.median(g["rmse"])),
            }
        tree_models = ["rf", "hgb", "extra"]
        best_tree = min(tree_models, key=lambda z: by_model[z]["trmse"])
        ref = refs[(ds, proto)]
        best_tree_trmse = by_model[best_tree]["trmse"]
        best_weak_trmse = float(ref["best_weak_trmse"])
        best_robust_trmse = float(ref["best_robust_trmse"])
        rows.append(
            {
                "dataset": ds,
                "protocol": proto,
                "role": m["role"],
                "n": int(m["n"]),
                "p": int(m["p"]),
                "rf_trmse": by_model["rf"]["trmse"],
                "hgb_trmse": by_model["hgb"]["trmse"],
                "extra_trmse": by_model["extra"]["trmse"],
                "linear_trmse": by_model["linear"]["trmse"],
                "best_tree_method": best_tree,
                "best_tree_trmse": best_tree_trmse,
                "best_tree_mae": by_model[best_tree]["mae"],
                "best_weak_method": ref["best_weak_method"],
                "best_weak_trmse": best_weak_trmse,
                "best_robust_method": ref["best_robust_method"],
                "best_robust_trmse": best_robust_trmse,
                "lse_trmse": float(ref["lse_trmse"]),
                "tree_vs_weak_trmse_pct": 100.0 * (best_tree_trmse / best_weak_trmse - 1.0),
                "tree_vs_robust_trmse_pct": 100.0 * (best_tree_trmse / best_robust_trmse - 1.0),
                "weak_vs_tree_trmse_pct": 100.0 * (best_weak_trmse / best_tree_trmse - 1.0),
                "winner_tree_vs_weak": "tree" if best_tree_trmse < best_weak_trmse else "weak",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", default=str(DEFAULT_IN))
    ap.add_argument("--trees", type=int, default=300)
    ap.add_argument("--max-splits", type=int, default=None, help="Smoke-test limit per dataset.")
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()

    indir = Path(args.input_dir)
    manifest = pd.read_csv(indir / "manifest.csv")
    splits = pd.read_csv(indir / "splits.csv")

    long_rows = []
    for _, m in manifest.iterrows():
        ds = str(m["dataset"])
        data = pd.read_csv(m["file"])
        y = data["y"].to_numpy(float)
        X = data.drop(columns=["y"]).to_numpy(float)
        ds_splits = splits[splits["dataset"] == ds].copy()
        if args.max_splits is not None:
            ds_splits = ds_splits.head(args.max_splits)
        print(f"[{ds} | {m['protocol']} | n={len(y)} p={X.shape[1]} splits={len(ds_splits)} role={m['role']}]")
        for _, sp in ds_splits.iterrows():
            train = parse_idx(sp["train_idx"])
            test = parse_idx(sp["test_idx"])
            seed = 91000 + int(sp["split_id"])
            preds = fit_models(X[train], y[train], X[test], seed=seed, trees=args.trees)
            for model, pred in preds.items():
                e = pred - y[test]
                long_rows.append(
                    {
                        "dataset": ds,
                        "protocol": m["protocol"],
                        "role": m["role"],
                        "split_id": int(sp["split_id"]),
                        "model": model,
                        "trmse": trmse(e),
                        "mae": float(np.mean(np.abs(e))),
                        "rmse": rmse(e),
                    }
                )

    long_df = pd.DataFrame(long_rows)
    summary = summarize(long_df, manifest)

    long_path = RESULTS / f"ml_baselines_crossdomain_long{args.suffix}.csv"
    summary_path = RESULTS / f"ml_baselines_crossdomain{args.suffix}.csv"
    long_df.to_csv(long_path, index=False)
    summary.to_csv(summary_path, index=False, quoting=csv.QUOTE_MINIMAL)

    print("\n=== Uniform cross-domain ML baseline ===")
    cols = ["dataset", "role", "best_tree_method", "best_tree_trmse", "best_weak_method",
            "best_weak_trmse", "tree_vs_weak_trmse_pct", "winner_tree_vs_weak"]
    print(summary[cols].to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print(f"\nSaved {summary_path}")
    print(f"Saved {long_path}")


if __name__ == "__main__":
    main()
