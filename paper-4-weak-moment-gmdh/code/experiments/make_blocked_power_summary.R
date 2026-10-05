#!/usr/bin/env Rscript
# Applies the pre-registered rule of blocked-power-prereg-2026-10-04.md to
# results/blocked_power_raw.csv; writes results/blocked_power_summary.csv.
# Usage: Rscript experiments/make_blocked_power_summary.R

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[[1]], mustWork = TRUE)) else getwd()
}
res_dir <- normalizePath(file.path(.script_dir(), "..", "..", "results"), mustWork = TRUE)
raw <- utils::read.csv(file.path(res_dir, "blocked_power_raw.csv"), stringsAsFactors = FALSE)
ev <- utils::read.csv(file.path(res_dir, "final_claim_evidence_table.csv"), stringsAsFactors = FALSE)
raw <- raw[raw$ok, , drop = FALSE]

# Per-fold error of one method: mean over the tournament-seed replicates.
fold_err <- function(d, g, m, col = "trmse") {
  z <- raw[raw$dataset == d & raw$gap == g & raw$method == m, ]
  stats::aggregate(z[[col]], by = list(fold = z$fold), FUN = mean)
}
cmp <- function(d, g, a, b, col = "trmse") {
  m <- merge(fold_err(d, g, a, col), fold_err(d, g, b, col), by = "fold", suffixes = c(".a", ".b"))
  r <- m$x.a / m$x.b
  c(n = nrow(m), paired_pct = 100 * (stats::median(r) - 1),
    ratio_med_pct = 100 * (stats::median(m$x.a) / stats::median(m$x.b) - 1), wins = sum(r < 1))
}
verdict <- function(p20, w20, n20, p50) {
  if (p20 <= -3 && w20 >= 7 && p50 < 0) return("supports a gain")
  if (p20 >= 3 && w20 <= 3) return("shows a loss")
  "no detectable difference"
}

rows <- list()
for (d in unique(raw$dataset)) {
  e <- ev[ev$dataset == d, ]
  rm <- e$best_robust_stat_method; wm <- e$best_weak_method
  for (who in c("auto-valgate", wm)) {
    g20 <- cmp(d, 20L, who, rm); g50 <- cmp(d, 50L, who, rm)
    u20 <- cmp(d, 20L, who, rm, "rmse")
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = d, role = e$role, compared = if (who == "auto-valgate") "gate" else "fixed weak",
      method = who, robust_method = rm, folds = g20[["n"]],
      gap20_paired_pct = g20[["paired_pct"]], gap20_ratio_med_pct = g20[["ratio_med_pct"]],
      gap20_wins = g20[["wins"]], gap50_paired_pct = g50[["paired_pct"]], gap50_wins = g50[["wins"]],
      gap20_rmse_paired_pct = u20[["paired_pct"]],
      verdict = verdict(g20[["paired_pct"]], g20[["wins"]], g20[["n"]], g50[["paired_pct"]]))
  }
}
out <- do.call(rbind, rows)
num <- vapply(out, is.numeric, logical(1)); out[num] <- lapply(out[num], round, 2)
utils::write.csv(out, file.path(res_dir, "blocked_power_summary.csv"), row.names = FALSE)
print(out, row.names = FALSE)
