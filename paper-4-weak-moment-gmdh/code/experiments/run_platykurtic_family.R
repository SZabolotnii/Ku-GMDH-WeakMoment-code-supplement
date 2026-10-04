#!/usr/bin/env Rscript
# Rebuilt small platykurtic / light-tailed family generator for Paper 4 controls.
# Usage: Rscript experiments/run_platykurtic_family.R [R=30] [suffix=""]

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[[1]], mustWork = TRUE)) else getwd()
}
source(file.path(.script_dir(), "_p4_rebuilt_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[[1]]) else 30L
suffix <- if (length(args) >= 2) args[[2]] else ""

load_small <- function(id) {
  if (id == "boot_motor_accel") {
    data(motor, package = "boot")
    d <- stats::na.omit(boot::motor)
    return(p4_frame_to_xy(d, accel ~ times + v + strata))
  }
  if (id == "toothgrowth") {
    d <- datasets::ToothGrowth
    d$supp <- factor(d$supp)
    return(p4_frame_to_xy(d, len ~ supp + dose))
  }
  if (id == "iris_versicolor") {
    d <- datasets::iris[datasets::iris$Species == "versicolor", ]
    return(p4_frame_to_xy(d, Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width))
  }
  if (id == "trees_volume") {
    return(p4_frame_to_xy(datasets::trees, Volume ~ Girth + Height))
  }
  if (id == "orange_circ") {
    d <- datasets::Orange
    return(p4_frame_to_xy(d, circumference ~ age + Tree))
  }
  stop("Unknown small dataset: ", id)
}

DATASETS <- c("boot_motor_accel", "toothgrowth", "iris_versicolor", "trees_volume", "orange_circ")

cat(sprintf("=== Rebuilt platykurtic/light family generator | R=%d | suffix='%s' ===\n", R, suffix))
raw <- list()
for (id in DATASETS) {
  dat <- load_small(id)
  cat(sprintf("[%s] n=%d p=%d\n", id, length(dat$y), ncol(dat$X)))
  raw[[length(raw) + 1L]] <- p4_run_random(id, dat, reps = R, seed0 = 75000L, id_col = "dataset")
}
raw <- do.call(rbind, raw)
p4_write(raw[, c("dataset", "rep", "method", "trmse", "mae", "rmse",
                 "pmm2_share", "pmm3_share", "wpmm2_share", "wpmm3_share", "ok")],
         "platykurtic_family_raw", suffix)

per <- p4_permethod(raw, id_col = "dataset", family_upper = TRUE)
methods <- data.frame(dataset = per$dataset, method = per$method, family = per$family,
                      med_trmse = round(per$trmse_med, 4),
                      med_mae = round(per$mae_med, 4),
                      med_rmse = round(per$rmse_med, 4))
p4_write(methods, "platykurtic_family_methods", suffix)

summary <- list()
for (id in DATASETS) {
  z <- p4_compare_one(raw, extra_filter = raw$dataset == id)
  dat <- load_small(id)
  auto <- raw[raw$dataset == id & raw$method == "auto", , drop = FALSE]
  aw <- raw[raw$dataset == id & raw$method == "auto-weak", , drop = FALSE]
  qs <- stats::quantile(z$unit_diff_pct, c(.25, .75), na.rm = TRUE, names = FALSE)
  summary[[length(summary) + 1L]] <- data.frame(
    dataset = id, n = length(dat$y),
    best_pmm_trmse_method = z$best_pmm_tr_method,
    best_robust_trmse_method = z$best_robust_tr_method,
    best_pmm_mae_method = z$best_pmm_mae_method,
    best_robust_mae_method = z$best_robust_mae_method,
    med_pmm_trmse = round(z$pmm_trmse, 4),
    med_robust_trmse = round(z$robust_trmse, 4),
    med_pmm_mae = round(z$pmm_mae, 4),
    med_robust_mae = round(z$robust_mae, 4),
    med_lse_trmse = round(z$lse_trmse, 4),
    pmm_vs_robust_trmse_pct = round(z$pmm_vs_robust_trmse_pct, 2),
    pmm_vs_robust_mae_pct = round(z$pmm_vs_robust_mae_pct, 2),
    pmm_vs_lse_pct = round(z$pmm_vs_lse_pct, 2),
    win_rate_trmse = round(z$win_rate_trmse, 3),
    iqr_diff_lo = round(qs[[1]], 2),
    iqr_diff_hi = round(qs[[2]], 2),
    wilcox_p_optimistic = signif(z$p_wilcox_optimistic, 3),
    auto_pmm3_share = round(stats::median(auto$pmm3_share, na.rm = TRUE), 3),
    autoweak_wpmm3_share = round(stats::median(aw$wpmm3_share, na.rm = TRUE), 3),
    autoweak_wpmm2_share = round(stats::median(aw$wpmm2_share, na.rm = TRUE), 3),
    verdict = p4_verdict(z$pmm_vs_robust_trmse_pct, z$p_wilcox_optimistic, z$win_rate_trmse),
    check.names = FALSE
  )
}
p4_write(do.call(rbind, summary), "platykurtic_family_summary", suffix)
cat("Saved platykurtic_family_*", suffix, ".csv\n", sep = "")
