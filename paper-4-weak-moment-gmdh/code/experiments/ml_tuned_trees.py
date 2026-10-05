#!/usr/bin/env python3
"""Tree ensembles with inner-CV tuning on the exact outer splits of the R generators.

Answers the 2026-10-04 review point that the tree baselines were untuned defaults while
the GMDH route selects internally. For every outer split, each ensemble is tuned on the
outer training part only, by 4-fold inner CV on trimmed RMSE (90 %) -- the criterion the
validation gate uses. Inner folds are contiguous for blocked (time-series) rows and
shuffled for random-split rows. The untuned defaults of ml_baselines_crossdomain.py are
refitted alongside, so tuned and default trees are compared on the same splits.

Usage: python ml_tuned_trees.py [--input-dir /tmp/p4_allrows_ml] [--max-splits N] [--suffix S]
"""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.model_selection import KFold

from ml_baselines_crossdomain import RESULTS, parse_idx, rmse, trmse

GRIDS = {
    "rf": {"max_features": [1.0, 0.5, "sqrt"], "min_samples_leaf": [1, 5, 20]},
    "extra": {"max_features": [1.0, 0.5, "sqrt"], "min_samples_leaf": [1, 5, 20]},
    "hgb": {"learning_rate": [0.03, 0.1], "max_leaf_nodes": [15, 31], "min_samples_leaf": [20, 50]},
}


def make(model: str, params: dict, seed: int, trees: int):
    if model == "rf":
        return RandomForestRegressor(n_estimators=trees, random_state=seed, n_jobs=1, **params)
    if model == "extra":
        return ExtraTreesRegressor(n_estimators=trees, random_state=seed, n_jobs=1, **params)
    return HistGradientBoostingRegressor(random_state=seed, max_iter=300, **params)


def tune(model, X, y, blocked, seed, trees):
    keys = list(GRIDS[model])
    best, best_err = None, np.inf
    kf = KFold(4, shuffle=not blocked, random_state=None if blocked else seed)
    for combo in itertools.product(*(GRIDS[model][k] for k in keys)):
        params = dict(zip(keys, combo))
        errs = []
        for tr, va in kf.split(X):
            m = make(model, params, seed, trees).fit(X[tr], y[tr])
            errs.append(trmse(m.predict(X[va]) - y[va]))
        e = float(np.mean(errs))
        if e < best_err:
            best, best_err = params, e
    return best


def one_split(ds, proto, role, X, y, sp, trees):
    train, test = parse_idx(sp["train_idx"]), parse_idx(sp["test_idx"])
    seed = 91000 + int(sp["split_id"])
    out = []
    for model in GRIDS:
        variants = {"default": None, "tuned": tune(model, X[train], y[train], proto == "blocked", seed, trees)}
        for kind, params in variants.items():
            if params is None:
                m = (RandomForestRegressor(n_estimators=300, random_state=seed, n_jobs=1) if model == "rf" else
                     ExtraTreesRegressor(n_estimators=300, random_state=seed, n_jobs=1) if model == "extra" else
                     HistGradientBoostingRegressor(random_state=seed))
            else:
                m = make(model, params, seed, 300)
            e = m.fit(X[train], y[train]).predict(X[test]) - y[test]
            out.append({"dataset": ds, "protocol": proto, "role": role, "split_id": int(sp["split_id"]),
                        "model": model, "variant": kind, "params": "" if params is None else repr(params),
                        "trmse": trmse(e), "mae": float(np.mean(np.abs(e))), "rmse": rmse(e)})
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", default="/tmp/p4_allrows_ml")
    ap.add_argument("--tune-trees", type=int, default=100, help="trees per forest during the inner search")
    ap.add_argument("--max-splits", type=int, default=None)
    ap.add_argument("--jobs", type=int, default=9)
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()

    indir = Path(args.input_dir)
    manifest = pd.read_csv(indir / "manifest.csv")
    splits = pd.read_csv(indir / "splits.csv")
    out = RESULTS / f"ml_tuned_trees_long{args.suffix}.csv"
    if out.exists():
        out.unlink()
    for _, m in manifest.iterrows():
        data = pd.read_csv(m["file"])
        y = data["y"].to_numpy(float)
        X = data.drop(columns=["y"]).to_numpy(float)
        sps = splits[splits["dataset"] == m["dataset"]]
        if args.max_splits is not None:
            sps = sps.head(args.max_splits)
        print(f"[{m['dataset']} | {m['protocol']} | n={len(y)} p={X.shape[1]} splits={len(sps)}]", flush=True)
        res = Parallel(n_jobs=args.jobs)(
            delayed(one_split)(m["dataset"], m["protocol"], m["role"], X, y, sp, args.tune_trees)
            for _, sp in sps.iterrows())
        # one dataset at a time, so an interrupted run keeps every finished row
        pd.DataFrame([r for chunk in res for r in chunk]).to_csv(
            out, index=False, mode="a", header=not out.exists())
        print(f"  appended to {out.name}", flush=True)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
