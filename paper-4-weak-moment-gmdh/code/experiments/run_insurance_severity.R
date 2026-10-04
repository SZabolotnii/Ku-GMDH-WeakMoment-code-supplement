#!/usr/bin/env Rscript
# Rebuilt insurance/severity generator for Paper 4.
# Usage: Rscript experiments/run_insurance_severity.R [R=30] [suffix=""]

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[[1]], mustWork = TRUE)) else getwd()
}
source(file.path(.script_dir(), "_p4_rebuilt_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[[1]]) else 30L
suffix <- if (length(args) >= 2) args[[2]] else ""

load_ins <- function(id) {
  if (id %in% c("fremtpl2_severity_raw", "fremtpl2_severity_log")) return(p4_ext_load(id))
  if (id == "insurance_autobi_loss") {
    data(AutoBi, package = "insuranceData")
    d <- stats::na.omit(get("AutoBi"))
    d$ATTORNEY <- factor(d$ATTORNEY)
    d$CLMSEX <- factor(d$CLMSEX)
    d$MARITAL <- factor(d$MARITAL)
    d$CLMINSUR <- factor(d$CLMINSUR)
    d$SEATBELT <- factor(d$SEATBELT)
    return(p4_frame_to_xy(d, LOSS ~ ATTORNEY + CLMSEX + MARITAL + CLMINSUR + SEATBELT + CLMAGE))
  }
  if (id == "insurance_autoclaims_paid") {
    data(AutoClaims, package = "insuranceData")
    d <- stats::na.omit(get("AutoClaims"))
    d$CLASS <- factor(d$CLASS)
    d$GENDER <- factor(d$GENDER)
    # STATE has many levels and makes the GMDH tournament prohibitively wide
    # for a routine evidence rebuild. Keep the actuarial risk-profile variables.
    return(p4_frame_to_xy(d, PAID ~ CLASS + GENDER + AGE))
  }
  stop("Unknown insurance dataset: ", id)
}

DATASETS <- c("fremtpl2_severity_raw", "insurance_autobi_loss",
              "insurance_autoclaims_paid", "fremtpl2_severity_log")

cat(sprintf("=== Rebuilt insurance severity generator | R=%d | suffix='%s' ===\n", R, suffix))
raw <- list()
for (id in DATASETS) {
  dat <- load_ins(id)
  cat(sprintf("[%s] n=%d p=%d\n", id, length(dat$y), ncol(dat$X)))
  raw[[length(raw) + 1L]] <- p4_run_random(id, dat, reps = R, seed0 = 74000L, id_col = "candidate")
}
raw <- do.call(rbind, raw)
p4_write(raw[, c("candidate", "rep", "method", "trmse", "mae", "rmse", "bias", "layers", "ok")],
         "insurance_severity_raw", suffix)

per <- p4_permethod(raw, id_col = "candidate")
p4_write(per[, c("candidate", "method", "family", "trmse_med", "trmse_iqr", "mae_med",
                 "mae_iqr", "rmse_med", "bias_med", "layers_med", "n_reps")],
         "insurance_severity_permethod", suffix)

fam <- list()
for (id in DATASETS) {
  z <- p4_compare_one(raw, extra_filter = raw$candidate == id)
  dat <- load_ins(id)
  fam[[length(fam) + 1L]] <- data.frame(
    candidate = id, n = length(dat$y),
    best_pmm_tr = z$best_pmm_tr_method,
    best_robust_tr = z$best_robust_tr_method,
    best_pmm_mae = z$best_pmm_mae_method,
    best_robust_mae = z$best_robust_mae_method,
    pmm_trmse = round(z$pmm_trmse, 4),
    robust_trmse = round(z$robust_trmse, 4),
    pmm_mae = round(z$pmm_mae, 4),
    robust_mae = round(z$robust_mae, 4),
    lse_trmse = round(z$lse_trmse, 4),
    pmm_vs_robust_trmse_pct = round(z$pmm_vs_robust_trmse_pct, 2),
    pmm_vs_robust_mae_pct = round(z$pmm_vs_robust_mae_pct, 2),
    pmm_vs_lse_pct = round(z$pmm_vs_lse_pct, 2),
    win_rate_trmse = round(z$win_rate_trmse, 3),
    p_wilcox_optimistic = signif(z$p_wilcox_optimistic, 3),
    verdict = p4_verdict(z$pmm_vs_robust_trmse_pct, z$p_wilcox_optimistic, z$win_rate_trmse),
    check.names = FALSE
  )
}
p4_write(do.call(rbind, fam), "insurance_severity_family", suffix)
cat("Saved insurance_severity_*", suffix, ".csv\n", sep = "")
