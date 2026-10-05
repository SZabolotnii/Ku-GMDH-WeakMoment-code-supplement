#!/usr/bin/env Rscript
# Summaries of results/claimrow_robustness_raw.csv against the archived claim-row fits.
# Writes results/claimrow_{valgate,ablation,sigma,seeds,gap,untrimmed}_summary.csv and
# results/claimrow_candidate_pool.csv.
# Usage: Rscript experiments/make_claimrow_robustness_summary.R

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[[1]], mustWork = TRUE)) else getwd()
}
res_dir <- normalizePath(file.path(.script_dir(), "..", "..", "results"), mustWork = TRUE)
rd <- function(f) utils::read.csv(file.path(res_dir, f), stringsAsFactors = FALSE)
wr <- function(x, f) { utils::write.csv(x, file.path(res_dir, f), row.names = FALSE); cat("wrote", f, "\n") }

ev <- rd("final_claim_evidence_table.csv")
new <- rd("claimrow_robustness_raw.csv")
pool_file <- file.path(res_dir, "claimrow_robustness_raw_pool.csv")
if (file.exists(pool_file)) new <- rbind(new, rd("claimrow_robustness_raw_pool.csv"))
new <- new[new$ok, , drop = FALSE]
soft <- rd("p4_softsensor_honest_raw.csv")
cross <- rd("crosssec_modskew_raw.csv")
ins <- rd("insurance_severity_raw.csv")
trees <- rd("ml_baselines_crossdomain_long.csv")
PMM <- c("auto", "auto-weak", "WPMM2", "WPMM3")
ROB <- c("LSE", "ridge-LSE", "Huber", "L1")

# Archived fits (tournament seed replicate 0, gap 20) in one long frame.
arch <- rbind(
  data.frame(dataset = soft$dataset, protocol = soft$protocol, split = soft$fold, method = soft$method,
             trmse = soft$trmse, rmse = soft$rmse, ok = soft$ok),
  data.frame(dataset = cross$candidate, protocol = "random", split = cross$rep, method = cross$method,
             trmse = cross$trmse, rmse = cross$rmse, ok = cross$ok),
  data.frame(dataset = ins$candidate, protocol = "random", split = ins$rep, method = ins$method,
             trmse = ins$trmse, rmse = ins$rmse, ok = ins$ok))
arch <- arch[arch$ok %in% c(TRUE, "TRUE"), ]

# Per-split error of one method, averaged over the requested seed replicates.
per_split <- function(d, proto, m, seeds = 0L, arm_new = "seeds", gap = 20L, sigma = 2.5) {
  parts <- list()
  if (0L %in% seeds) {
    if (m == "auto-valgate") {
      z <- new[new$dataset == d & new$arm == "valgate" & new$seed_rep == 0, ]
    } else z <- arch[arch$dataset == d & arch$protocol == proto & arch$method == m, ]
    parts[[1]] <- z[, c("split", "trmse", "rmse")]
  }
  s_new <- setdiff(seeds, 0L)
  if (length(s_new)) {
    arm <- if (m == "auto-valgate") "valgate" else arm_new
    z <- new[new$dataset == d & new$arm == arm & new$method == m & new$seed_rep %in% s_new, ]
    parts[[length(parts) + 1L]] <- z[, c("split", "trmse", "rmse")]
  }
  z <- do.call(rbind, parts)
  stats::aggregate(cbind(trmse, rmse) ~ split, data = z, FUN = mean)
}
pair <- function(a, b, col = "trmse") {
  m <- merge(a, b, by = "split", suffixes = c(".a", ".b"))
  x <- m[[paste0(col, ".a")]]; y <- m[[paste0(col, ".b")]]
  c(n = nrow(m), ratio_med_pct = 100 * (stats::median(x) / stats::median(y) - 1),
    paired_pct = 100 * (stats::median(x / y) - 1), wins = sum(x < y))
}
best_of <- function(d, proto, fam, seeds) {
  meds <- vapply(fam, function(m) stats::median(per_split(d, proto, m, seeds)$trmse), numeric(1))
  names(which.min(meds))
}

SEEDS <- 0:2
vg_rows <- list(); seed_rows <- list(); un_rows <- list(); ab_rows <- list(); sg_rows <- list(); gap_rows <- list()
for (i in seq_len(nrow(ev))) {
  d <- ev$dataset[i]; proto <- ev$protocol[i]; role <- ev$role[i]
  wm <- ev$best_weak_method[i]; rm <- ev$best_robust_stat_method[i]; tm <- ev$best_tree_method[i]
  vg <- per_split(d, proto, "auto-valgate", SEEDS)
  rb <- per_split(d, proto, rm, SEEDS)
  wk <- per_split(d, proto, wm, SEEDS)
  tr <- trees[trees$dataset == d & trees$model == tm, c("split_id", "trmse", "rmse")]
  names(tr)[1] <- "split"
  a <- pair(vg, rb); b <- pair(vg, wk); t <- pair(vg, tr)
  vg_rows[[i]] <- data.frame(dataset = d, role = role, protocol = proto, n_splits = a[["n"]],
    robust_method = rm, best_weak_method = wm, tree_method = tm,
    vg_vs_robust_ratio_med_pct = a[["ratio_med_pct"]], vg_vs_robust_paired_pct = a[["paired_pct"]],
    vg_wins_robust = a[["wins"]],
    vg_vs_bestweak_paired_pct = b[["paired_pct"]],
    vg_vs_tree_ratio_med_pct = t[["ratio_med_pct"]], vg_vs_tree_paired_pct = t[["paired_pct"]],
    vg_wins_tree = t[["wins"]], seeds = "0-2")

  # Seed dispersion of the Table-2 delta: family-best re-selected at each seed replicate.
  for (s in SEEDS) {
    bw <- best_of(d, proto, PMM, s); br <- best_of(d, proto, ROB, s)
    p <- pair(per_split(d, proto, bw, s), per_split(d, proto, br, s))
    seed_rows[[length(seed_rows) + 1L]] <- data.frame(dataset = d, role = role, seed_rep = s,
      best_weak = bw, best_robust = br, weak_vs_robust_ratio_med_pct = p[["ratio_med_pct"]],
      weak_vs_robust_paired_pct = p[["paired_pct"]])
  }

  # Untrimmed tail metric for the Table-2 pair (archived seed 0) and for the gate.
  u <- pair(per_split(d, proto, wm, 0L), per_split(d, proto, rm, 0L), "rmse")
  uv <- pair(per_split(d, proto, "auto-valgate", SEEDS), rb, "rmse")
  un_rows[[i]] <- data.frame(dataset = d, role = role,
    weak_vs_robust_rmse_ratio_med_pct = u[["ratio_med_pct"]], weak_vs_robust_rmse_paired_pct = u[["paired_pct"]],
    weak_rmse_wins = u[["wins"]], vg_vs_robust_rmse_paired_pct = uv[["paired_pct"]], n_splits = u[["n"]])

  # Selection-layer ablations, all at seed replicate 0.
  v0 <- per_split(d, proto, "auto-valgate", 0L)
  for (arm in c("rinner", "norobust")) {
    z <- new[new$dataset == d & new$arm == arm, c("split", "trmse", "rmse")]
    p <- pair(z, v0); q <- pair(z, per_split(d, proto, rm, 0L))
    ab_rows[[length(ab_rows) + 1L]] <- data.frame(dataset = d, role = role, arm = arm,
      vs_valgate_paired_pct = p[["paired_pct"]], vs_valgate_ratio_med_pct = p[["ratio_med_pct"]],
      vs_robust_paired_pct = q[["paired_pct"]], wins_vs_robust = q[["wins"]], n_splits = p[["n"]])
  }

  # Window-width sensitivity, relative to the default sigma_mult = 2.5 at seed replicate 0.
  for (m in c("WPMM2", "WPMM3", "auto-valgate")) for (sg in c(1.5, 4.0)) {
    z <- new[new$dataset == d & new$arm == "sigma" & new$method == m & new$sigma == sg, c("split", "trmse", "rmse")]
    p <- pair(z, per_split(d, proto, m, 0L)); q <- pair(z, per_split(d, proto, rm, 0L))
    sg_rows[[length(sg_rows) + 1L]] <- data.frame(dataset = d, role = role, method = m, sigma_mult = sg,
      vs_default_ratio_med_pct = p[["ratio_med_pct"]], vs_robust_ratio_med_pct = q[["ratio_med_pct"]])
  }

  # Purge-gap sensitivity on the blocked rows: family-best delta and the gate, re-selected per gap.
  if (proto == "blocked") for (g in c(20L, 50L, 100L)) {
    getm <- function(m) {
      if (g == 20L) return(per_split(d, proto, m, 0L))
      z <- new[new$dataset == d & new$arm == "gap" & new$gap == g & new$method == m, c("split", "trmse", "rmse")]
      z
    }
    med <- function(m) stats::median(getm(m)$trmse)
    bw <- PMM[which.min(vapply(PMM, med, numeric(1)))]; br <- ROB[which.min(vapply(ROB, med, numeric(1)))]
    p <- pair(getm(bw), getm(br)); q <- pair(getm("auto-valgate"), getm(br))
    gap_rows[[length(gap_rows) + 1L]] <- data.frame(dataset = d, role = role, gap = g, best_weak = bw,
      best_robust = br, weak_vs_robust_ratio_med_pct = p[["ratio_med_pct"]],
      vg_vs_robust_ratio_med_pct = q[["ratio_med_pct"]], vg_vs_robust_paired_pct = q[["paired_pct"]],
      vg_wins = q[["wins"]], n_splits = q[["n"]])
  }
}
rnd <- function(x) { num <- vapply(x, is.numeric, logical(1)); x[num] <- lapply(x[num], round, 2); x }
wr(rnd(do.call(rbind, vg_rows)), "claimrow_valgate_summary.csv")
wr(rnd(do.call(rbind, seed_rows)), "claimrow_seeds_summary.csv")
wr(rnd(do.call(rbind, un_rows)), "claimrow_untrimmed_summary.csv")
wr(rnd(do.call(rbind, ab_rows)), "claimrow_ablation_summary.csv")
wr(rnd(do.call(rbind, sg_rows)), "claimrow_sigma_summary.csv")
wr(rnd(do.call(rbind, gap_rows)), "claimrow_gap_summary.csv")

# Candidate pool: every dataset the three claim-row generators ran, with its outcome.
sc <- rd("crosssec_modskew_comparison.csv")
isv <- rd("insurance_severity_family.csv")
so <- rd("p4_softsensor_honest_compare.csv")
so <- so[(so$dataset %in% c("sru_y1_static", "sru_y1_dynamic", "sru_y2_static",
                            "gas_turbine_co_2015_raw", "gas_turbine_nox_2015_raw")) & so$protocol == "blocked", ]
pool <- rbind(
  data.frame(dataset = sc$candidate, generator = "run_p4_crosssec_modskew.R", protocol = "random",
             weak_vs_robust_trmse_pct = sc$pmm_vs_robust_trmse_pct),
  data.frame(dataset = isv$candidate, generator = "run_insurance_severity.R", protocol = "random",
             weak_vs_robust_trmse_pct = isv$pmm_vs_robust_trmse_pct),
  data.frame(dataset = so$dataset, generator = "run_p4_softsensor_honest.R", protocol = "blocked",
             weak_vs_robust_trmse_pct = so$pmm_vs_rob_trmse_pct))
pool$outcome <- ifelse(pool$weak_vs_robust_trmse_pct <= -3, "weak ahead (>=3%)",
                ifelse(pool$weak_vs_robust_trmse_pct >= 3, "robust ahead (>=3%)", "within 3%"))
pool$in_table2 <- pool$dataset %in% ev$dataset
pool$table2_role <- ev$role[match(pool$dataset, ev$dataset)]

# The gate on every pool row (Table-2 rows from the summary above, the seven others from
# claimrow_robustness_raw_pool.csv), against the best robust method of the seed-0 pool run.
best_rob <- c(stats::setNames(sc$best_robust_tr, sc$candidate),
              stats::setNames(isv$best_robust_tr, isv$candidate),
              stats::setNames(so$best_rob_trmse_method, so$dataset))
gate <- lapply(seq_len(nrow(pool)), function(i) {
  d <- pool$dataset[i]; proto <- pool$protocol[i]; rm <- best_rob[[d]]
  if (!any(new$dataset == d & new$arm == "valgate")) return(NULL)
  vg <- per_split(d, proto, "auto-valgate", SEEDS); rb <- per_split(d, proto, rm, SEEDS)
  a <- pair(vg, rb); u <- pair(vg, rb, "rmse")
  bw <- best_of(d, proto, PMM, SEEDS)                 # best weak by median test error: hindsight
  h <- pair(vg, per_split(d, proto, bw, SEEDS))
  seeds_delta <- vapply(SEEDS, function(s) {
    bw <- best_of(d, proto, PMM, s); br <- best_of(d, proto, ROB, s)
    pair(per_split(d, proto, bw, s), per_split(d, proto, br, s))[["paired_pct"]]
  }, numeric(1))
  data.frame(dataset = d, robust_method = rm, n_splits = a[["n"]],
             gate_vs_robust_ratio_med_pct = a[["ratio_med_pct"]],
             gate_vs_robust_paired_pct = a[["paired_pct"]], gate_wins = a[["wins"]],
             gate_vs_robust_rmse_paired_pct = u[["paired_pct"]],
             hindsight_weak = bw, gate_vs_hindsight_weak_paired_pct = h[["paired_pct"]],
             family_delta_paired_min = min(seeds_delta), family_delta_paired_max = max(seeds_delta))
})
pool <- merge(pool, do.call(rbind, gate), by = "dataset", all.x = TRUE, sort = FALSE)
wr(rnd(pool), "claimrow_candidate_pool.csv")

# Tuned versus default tree ensembles (ml_tuned_trees.py) and the gate against the best
# tuned ensemble, on every pool row.
tt_file <- file.path(res_dir, "ml_tuned_trees_long.csv")
if (file.exists(tt_file)) {
  tt <- rd("ml_tuned_trees_long.csv")
  trow <- lapply(unique(tt$dataset), function(d) {
    s <- tt[tt$dataset == d, ]
    med <- stats::aggregate(trmse ~ model + variant, data = s, FUN = stats::median)
    bd <- med[med$variant == "default", ][which.min(med$trmse[med$variant == "default"]), ]
    bt <- med[med$variant == "tuned", ][which.min(med$trmse[med$variant == "tuned"]), ]
    tun <- s[s$model == bt$model & s$variant == "tuned", c("split_id", "trmse", "rmse")]
    names(tun)[1] <- "split"
    proto <- s$protocol[[1]]
    a <- pair(per_split(d, proto, "auto-valgate", SEEDS), tun)
    u <- pair(per_split(d, proto, "auto-valgate", SEEDS), tun, "rmse")
    dft <- s[s$model == bd$model & s$variant == "default", c("split_id", "trmse", "rmse")]
    names(dft)[1] <- "split"
    a0 <- pair(per_split(d, proto, "auto-valgate", SEEDS), dft)
    data.frame(dataset = d, best_default_tree = bd$model, best_default_trmse = bd$trmse,
               best_tuned_tree = bt$model, best_tuned_trmse = bt$trmse,
               tuned_vs_default_pct = 100 * (bt$trmse / bd$trmse - 1), n_splits = a[["n"]],
               gate_vs_tuned_ratio_med_pct = a[["ratio_med_pct"]],
               gate_vs_tuned_paired_pct = a[["paired_pct"]], gate_wins = a[["wins"]],
               gate_vs_tuned_rmse_paired_pct = u[["paired_pct"]], gate_rmse_wins = u[["wins"]],
               gate_vs_default_paired_pct = a0[["paired_pct"]], gate_wins_vs_default = a0[["wins"]])
  })
  trees_sum <- do.call(rbind, trow)
  wr(rnd(trees_sum), "claimrow_tuned_trees_summary.csv")

  # Per-split file behind the pool figure: gate, the row's robust method (both averaged
  # over seed replicates 0-2) and the best tuned tree, on the same split.
  sp_rows <- lapply(seq_len(nrow(pool)), function(i) {
    d <- pool$dataset[i]; proto <- pool$protocol[i]
    vg <- per_split(d, proto, "auto-valgate", SEEDS)
    rb <- per_split(d, proto, best_rob[[d]], SEEDS)
    bt <- trees_sum$best_tuned_tree[trees_sum$dataset == d]
    tr <- tt[tt$dataset == d & tt$model == bt & tt$variant == "tuned", c("split_id", "trmse", "rmse")]
    names(tr) <- c("split", "tree_trmse", "tree_rmse")
    m <- merge(merge(vg, rb, by = "split", suffixes = c("_gate", "_robust")), tr, by = "split")
    data.frame(dataset = d, protocol = proto, m)
  })
  wr(do.call(rbind, sp_rows), "claimrow_gate_splits.csv")

  # Row-level summary: one paired delta per row, signed-rank test across rows.
  pr <- merge(pool[, c("dataset", "protocol", "gate_vs_robust_paired_pct", "gate_vs_robust_rmse_paired_pct")],
              trees_sum[, c("dataset", "gate_vs_tuned_paired_pct")], by = "dataset")
  insurance <- c("fremtpl2_severity_raw", "fremtpl2_severity_log", "insurance_autobi_loss",
                 "insurance_autoclaims_paid")
  rl <- function(label, rows, col) {
    x <- pr[[col]][pr$dataset %in% rows]
    p <- tryCatch(stats::wilcox.test(x, exact = FALSE)$p.value, error = function(e) NA_real_)
    data.frame(group = label, comparison = col, rows = length(x), negative = sum(x < 0),
               positive = sum(x > 0), median_pct = stats::median(x), wilcoxon_p = p)
  }
  rlev <- do.call(rbind, lapply(c("gate_vs_robust_paired_pct", "gate_vs_robust_rmse_paired_pct",
                                  "gate_vs_tuned_paired_pct"), function(col)
    rbind(rl("all 15 rows", pr$dataset, col), rl("insurance (4 rows)", insurance, col))))
  rlev$wilcoxon_p <- signif(rlev$wilcoxon_p, 3)
  wr(rnd(rlev), "claimrow_rowlevel_summary.csv")
}

# How the five-fold, single-seed, best-of-family protocol produced the blocked-row gains:
# the Table-2-style delta at seed 0 and at seeds 1-2, the gate on the same five folds, and
# the gate in the pre-registered ten-fold, five-seed check.
bp_file <- file.path(res_dir, "blocked_power_summary.csv")
if (file.exists(bp_file)) {
  bp <- rd("blocked_power_summary.csv"); bp <- bp[bp$compared == "gate", ]
  sd_tab <- rd("claimrow_seeds_summary.csv")
  vg_tab <- rd("claimrow_valgate_summary.csv")
  les <- lapply(bp$dataset, function(d) {
    s <- sd_tab[sd_tab$dataset == d, ]
    v <- vg_tab[vg_tab$dataset == d, ]
    b <- bp[bp$dataset == d, ]
    data.frame(dataset = d, role = b$role,
               hindsight_seed0_ratio_med_pct = s$weak_vs_robust_ratio_med_pct[s$seed_rep == 0],
               hindsight_seed1_ratio_med_pct = s$weak_vs_robust_ratio_med_pct[s$seed_rep == 1],
               hindsight_seed2_ratio_med_pct = s$weak_vs_robust_ratio_med_pct[s$seed_rep == 2],
               gate5_paired_pct = v$vg_vs_robust_paired_pct, gate5_wins = v$vg_wins_robust, gate5_folds = v$n_splits,
               gate10_paired_pct = b$gap20_paired_pct, gate10_wins = b$gap20_wins, gate10_folds = b$folds,
               gate10_gap50_paired_pct = b$gap50_paired_pct, verdict = b$verdict)
  })
  wr(rnd(do.call(rbind, les)), "claimrow_protocol_lesson.csv")
}

# Large-claim error on the insurance rows: the `tail` arm of run_claimrow_robustness.R
# (robust methods and the gate, seed replicates 0-2) and ml_tuned_trees.py --tail. Every
# refit must reproduce the archived trimmed-RMSE of the same split, method and seed.
tail_file <- file.path(res_dir, "claimrow_robustness_raw_tail.csv")
ttail_file <- file.path(res_dir, "ml_tuned_trees_long_tail.csv")
if (file.exists(tail_file) && file.exists(ttail_file) && exists("trees_sum")) {
  tl <- rd("claimrow_robustness_raw_tail.csv"); tl <- tl[tl$ok, ]
  tt_tail <- rd("ml_tuned_trees_long_tail.csv")
  METRICS <- c("trmse", "rmse", "mae", "tail_rmse", "tail_mae")
  repro <- function(d, proto, m) {
    old <- do.call(rbind, lapply(SEEDS, function(s) {
      z <- if (s == 0L && m != "auto-valgate") arch[arch$dataset == d & arch$method == m, c("split", "trmse")]
           else new[new$dataset == d & new$method == m & new$seed_rep == s &
                    new$arm == (if (m == "auto-valgate") "valgate" else "seeds"), c("split", "trmse")]
      data.frame(z, seed_rep = s)
    }))
    x <- merge(tl[tl$dataset == d & tl$method == m, c("split", "seed_rep", "trmse")], old,
               by = c("split", "seed_rep"))
    max(abs(x$trmse.x - x$trmse.y))
  }
  tail_rows <- lapply(unique(tl$dataset), function(d) {
    proto <- "random"; rm <- best_rob[[d]]; bt <- trees_sum$best_tuned_tree[trees_sum$dataset == d]
    avg <- function(m) stats::aggregate(tl[tl$dataset == d & tl$method == m, METRICS],
                                        by = list(split = tl$split[tl$dataset == d & tl$method == m]), FUN = mean)
    vg <- avg("auto-valgate"); rb <- avg(rm)
    tr <- tt_tail[tt_tail$dataset == d & tt_tail$model == bt & tt_tail$variant == "tuned", c("split_id", METRICS)]
    names(tr)[1] <- "split"
    old_tr <- tt[tt$dataset == d & tt$model == bt & tt$variant == "tuned", c("split_id", "trmse")]
    tree_diff <- max(abs(merge(tr, old_tr, by.x = "split", by.y = "split_id")[, c("trmse.x")] -
                         merge(tr, old_tr, by.x = "split", by.y = "split_id")[, c("trmse.y")]))
    out <- data.frame(dataset = d, robust_method = rm, best_tuned_tree = bt, n_splits = nrow(vg),
                      mean_n_tail = mean(tl$n_tail[tl$dataset == d & tl$method == "auto-valgate"]))
    for (k in METRICS) {
      a <- pair(vg, rb, k); b <- pair(vg, tr, k)
      out[[paste0("gate_vs_robust_", k, "_paired_pct")]] <- a[["paired_pct"]]
      out[[paste0("gate_vs_robust_", k, "_wins")]] <- a[["wins"]]
      out[[paste0("gate_vs_tree_", k, "_paired_pct")]] <- b[["paired_pct"]]
      out[[paste0("gate_vs_tree_", k, "_wins")]] <- b[["wins"]]
    }
    out$repro_max_abs_diff <- max(repro(d, proto, "auto-valgate"), repro(d, proto, rm), tree_diff)
    out
  })
  tail_sum <- do.call(rbind, tail_rows)
  cat("tail refits, max |trimmed-RMSE - archived| per dataset:",
      paste(signif(tail_sum$repro_max_abs_diff, 3), collapse = " "), "\n")
  if (any(tail_sum$repro_max_abs_diff > 1e-8)) warning("tail refits do not reproduce the archived fits")
  wr(rnd(tail_sum), "claimrow_tail_summary.csv")
}
