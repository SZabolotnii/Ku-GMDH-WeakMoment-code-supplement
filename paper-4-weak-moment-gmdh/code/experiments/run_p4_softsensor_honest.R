#!/usr/bin/env Rscript
# Rebuilt honest soft-sensor / emissions generator for Paper 4.
# Usage: Rscript experiments/run_p4_softsensor_honest.R [R=30] [suffix=""]

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[[1]], mustWork = TRUE)) else getwd()
}
source(file.path(.script_dir(), "_p4_rebuilt_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[[1]]) else 30L
suffix <- if (length(args) >= 2) args[[2]] else ""

DATASETS <- c("sru_y1_static", "sru_y1_dynamic", "sru_y2_static",
              "gas_turbine_co_2015_raw", "gas_turbine_nox_2015_raw")

cat(sprintf("=== Rebuilt honest soft-sensor generator | R=%d | suffix='%s' ===\n", R, suffix))
raw <- list()
for (id in DATASETS) {
  dat <- p4_ext_load(id)
  cat(sprintf("[%s] n=%d p=%d\n", id, length(dat$y), ncol(dat$X)))
  raw[[length(raw) + 1L]] <- p4_run_protocols(id, dat, reps = R, seed0 = 76000L,
                                              cap_blocked = 3000L, cap_random = 4000L,
                                              gap = 20L)
}
raw <- do.call(rbind, raw)
p4_write(raw, "p4_softsensor_honest_raw", suffix)

sum_rows <- list()
for (id in DATASETS) for (proto in c("blocked", "random")) {
  for (m in P4_METHODS) {
    sub <- raw[raw$dataset == id & raw$protocol == proto & raw$method == m & raw$ok, , drop = FALSE]
    sum_rows[[length(sum_rows) + 1L]] <- data.frame(
      dataset = id, protocol = proto, method = m,
      trmse_med = stats::median(sub$trmse, na.rm = TRUE),
      mae_med = stats::median(sub$mae, na.rm = TRUE),
      rmse_med = stats::median(sub$rmse, na.rm = TRUE),
      trmse_iqr = stats::IQR(sub$trmse, na.rm = TRUE),
      check.names = FALSE
    )
  }
}
p4_write(do.call(rbind, sum_rows), "p4_softsensor_honest_summary", suffix)

cmp <- list()
for (id in DATASETS) for (proto in c("blocked", "random")) {
  filt <- raw$dataset == id & raw$protocol == proto
  z <- p4_compare_one(raw, extra_filter = filt)
  cmp[[length(cmp) + 1L]] <- data.frame(
    dataset = id, protocol = proto,
    best_pmm_trmse_method = z$best_pmm_tr_method,
    best_rob_trmse_method = z$best_robust_tr_method,
    best_pmm_mae_method = z$best_pmm_mae_method,
    best_rob_mae_method = z$best_robust_mae_method,
    pmm_trmse = round(z$pmm_trmse, 4),
    rob_trmse = round(z$robust_trmse, 4),
    pmm_vs_rob_trmse_pct = round(z$pmm_vs_robust_trmse_pct, 2),
    pmm_mae = round(z$pmm_mae, 4),
    rob_mae = round(z$robust_mae, 4),
    pmm_vs_rob_mae_pct = round(z$pmm_vs_robust_mae_pct, 2),
    lse_trmse = round(z$lse_trmse, 4),
    pmm_vs_lse_trmse_pct = round(z$pmm_vs_lse_pct, 2),
    win_rate_trmse = round(z$win_rate_trmse, 3),
    check.names = FALSE
  )
}
p4_write(do.call(rbind, cmp), "p4_softsensor_honest_compare", suffix)
cat("Saved p4_softsensor_honest_*", suffix, ".csv\n", sep = "")
