#!/usr/bin/env Rscript
# Rebuilt cross-sectional moderate-skew search for Paper 4.
# Usage: Rscript experiments/run_p4_crosssec_modskew.R [R=30] [suffix=""]

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[[1]], mustWork = TRUE)) else getwd()
}
source(file.path(.script_dir(), "_p4_rebuilt_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[[1]]) else 30L
suffix <- if (length(args) >= 2) args[[2]] else ""

load_builtin <- function(id) {
  if (id == "mass_boston_medv") {
    # The screened formula (shared/datasets/screen_cumulants.R): the race-proxy column
    # `black` is excluded.
    data(Boston, package = "MASS")
    return(p4_frame_to_xy(MASS::Boston, medv ~ crim + zn + indus + chas + nox + rm + age + dis +
                            rad + tax + ptratio + lstat))
  }
  if (id == "islr_wage") {
    data(Wage, package = "ISLR2")
    d <- ISLR2::Wage
    # `race` removed (2026-10-05); the screened formula included it.
    return(p4_frame_to_xy(d, wage ~ year + age + maritl + education + jobclass + health + health_ins))
  }
  if (id == "airquality_ozone") {
    d <- stats::na.omit(datasets::airquality)
    return(p4_frame_to_xy(d, Ozone ~ Solar.R + Wind + Temp + Month + Day))
  }
  if (id == "islr_credit_balance") {
    data(Credit, package = "ISLR2")
    d <- stats::na.omit(ISLR2::Credit)
    return(p4_frame_to_xy(d, Balance ~ Income + Limit + Rating + Cards + Age + Education +
                            Own + Student + Married + Region))
  }
  if (id == "mass_cars93_mpg") {
    data(Cars93, package = "MASS")
    d <- stats::na.omit(MASS::Cars93)
    return(p4_frame_to_xy(d, MPG.city ~ EngineSize + Horsepower + RPM + Rev.per.mile +
                            Fuel.tank.capacity + Length + Wheelbase + Width + Weight))
  }
  if (id == "concrete") return(p4_ext_load("concrete"))
  stop("Unknown dataset: ", id)
}

DATASETS <- c("mass_boston_medv", "islr_wage", "airquality_ozone",
              "islr_credit_balance", "concrete", "mass_cars93_mpg")

cat(sprintf("=== Rebuilt cross-sectional moderate-skew search | R=%d | suffix='%s' ===\n", R, suffix))
raw <- list()
for (id in DATASETS) {
  dat <- load_builtin(id)
  cat(sprintf("[%s] n=%d p=%d\n", id, length(dat$y), ncol(dat$X)))
  raw[[length(raw) + 1L]] <- p4_run_random(id, dat, reps = R, seed0 = 73000L, id_col = "candidate")
}
raw <- do.call(rbind, raw)
p4_write(raw[, c("candidate", "rep", "method", "trmse", "mae", "rmse", "bias", "layers", "ok")],
         "crosssec_modskew_raw", suffix)

summary <- p4_permethod(raw, id_col = "candidate")
summary <- summary[, c("candidate", "method", "n_reps", "trmse_med", "trmse_iqr",
                       "mae_med", "mae_iqr", "rmse_med", "bias_med")]
names(summary)[names(summary) == "n_reps"] <- "n_ok"
p4_write(summary, "crosssec_modskew_summary", suffix)

cmp <- list()
for (id in DATASETS) {
  z <- p4_compare_one(raw, extra_filter = raw$candidate == id)
  cmp[[length(cmp) + 1L]] <- data.frame(
    candidate = id,
    best_pmm_method_tr = z$best_pmm_tr_method,
    best_pmm_method_mae = z$best_pmm_mae_method,
    best_robust_tr = z$best_robust_tr_method,
    best_robust_mae = z$best_robust_mae_method,
    pmm_vs_robust_trmse_pct = round(z$pmm_vs_robust_trmse_pct, 2),
    pmm_vs_robust_mae_pct = round(z$pmm_vs_robust_mae_pct, 2),
    pmm_vs_lse_trmse_pct = round(z$pmm_vs_lse_pct, 2),
    pmm_win_rate_trmse = round(z$win_rate_trmse, 3),
    pmm_win_rate_mae = NA_real_,
    wilcox_p_trmse_optimistic = signif(z$p_wilcox_optimistic, 3),
    best_pmm_trmse = round(z$pmm_trmse, 4),
    best_rob_trmse = round(z$robust_trmse, 4),
    lse_trmse = round(z$lse_trmse, 4),
    check.names = FALSE
  )
}
p4_write(do.call(rbind, cmp), "crosssec_modskew_comparison", suffix)
cat("Saved crosssec_modskew_*", suffix, ".csv\n", sep = "")
