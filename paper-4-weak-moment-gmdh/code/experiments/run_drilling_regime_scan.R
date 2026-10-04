#!/usr/bin/env Rscript
# Real-data regime scan for Paper 4: find a STATIONARY heavy-tailed-SKEWED segment
# of real drilling ROP residuals — the regime where weak-PMM2 / the weak dispatch
# is predicted to help. Slides a window along depth, fits an LSE ROP model per
# window, and reports the residual cumulants (gamma3 skew, gamma4 excess kurtosis,
# tail_ratio = sd/(MAD/0.6745)). Flags windows that are heavy AND skewed.
# Usage (from paper-4-weak-moment-gmdh/code): Rscript experiments/run_drilling_regime_scan.R

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
                 else pkgload::load_all(PKG, quiet = TRUE))

scan_one <- function(label, loader_path, win = 500L, step = 250L) {
  d <- tryCatch(load_volve_drilling(file.path(PKG, loader_path)), error = function(e) NULL)
  if (is.null(d)) { cat("  (skip", label, ")\n"); return(NULL) }
  X <- as.data.frame(d$X); y <- as.numeric(d$y)
  cat(sprintf("  %-16s loaded n=%d, cols=[%s]\n", label, length(y), paste(names(X), collapse=",")))
  ord <- order(X$dept); X <- X[ord, ]; y <- y[ord]
  feat <- setdiff(names(X), "dept")
  n <- length(y)
  if (n <= win + step) { cat("    (too few rows for window", win, ")\n"); return(NULL) }
  starts <- seq(1L, n - win, by = step)
  rows <- lapply(starts, function(s) {
    idx <- s:(s + win - 1L); Xi <- X[idx, feat, drop = FALSE]; yi <- y[idx]
    fit <- tryCatch(stats::lm(yi ~ ., data = Xi), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    r <- stats::residuals(fit); cu <- sample_cumulants(r)
    s_eps <- stats::sd(r); mad_s <- stats::mad(r, constant = 1.4826)
    data.frame(label = label, depth0 = X$dept[idx[1]], depth1 = X$dept[idx[win]],
               gamma3 = cu$gamma3, gamma4 = cu$gamma4,
               tail_ratio = if (mad_s > 1e-9) s_eps / mad_s else NA, n = win)
  })
  do.call(rbind, rows)
}

cat("=== Real drilling residual regime scan (window=500, step=250) ===\n")
all <- rbind(
  scan_one("volve_onbottom", "data/volve_onbottom.csv"),
  scan_one("volve_real",     "data/volve_15_9_F15_real.csv")
)
outdir <- "../results"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(all, file.path(outdir, "drilling_regime_scan.csv"), row.names = FALSE)

cat(sprintf("\nScanned %d windows. Global residual character per source:\n", nrow(all)))
for (lb in unique(all$label)) { s <- all[all$label == lb, ]
  cat(sprintf("  %-16s windows=%d | gamma3 med=%+.2f [%.2f,%.2f] | gamma4 med=%.1f [max %.0f] | tail_ratio med=%.2f\n",
      lb, nrow(s), median(s$gamma3), min(s$gamma3), max(s$gamma3),
      median(s$gamma4), max(s$gamma4), median(s$tail_ratio, na.rm=TRUE))) }

# Paper-4 niche = heavy AND skewed: |gamma3|>0.5 and gamma4>4 (asymmetric heavy)
niche <- all[abs(all$gamma3) > 0.5 & all$gamma4 > 4, ]
niche <- niche[order(-niche$gamma4), ]
cat(sprintf("\nHeavy-tailed-skewed windows (|g3|>0.5 & g4>4): %d of %d\n", nrow(niche), nrow(all)))
if (nrow(niche)) { cat("Top candidates (by kurtosis):\n")
  print(utils::head(niche[, c("label","depth0","depth1","gamma3","gamma4","tail_ratio")], 12), row.names = FALSE, digits = 3) }
# also clean-skew (Paper-1 PMM2 regime) for contrast
clean <- all[abs(all$gamma3) > 0.5 & all$gamma4 <= 2 & all$gamma4 > -1, ]
cat(sprintf("\nClean-skew windows (|g3|>0.5 & -1<g4<=2, classical-PMM2 regime): %d\n", nrow(clean)))
cat(sprintf("\nSaved drilling_regime_scan.csv\n"))
