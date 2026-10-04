#!/usr/bin/env Rscript
# Real-data demonstration for Paper 4 (gap A): does the weak inner estimator / weak
# dispatch help on REAL drilling ROP residuals in a stationary heavy-tailed-skewed
# segment, IN-DISTRIBUTION (random train/test splits, so the heavy-tail-skew effect
# is isolated from the extrapolation that confounded Paper 1's reserve windows)?
#
# Per depth band: characterize the LSE-residual cumulants (+ a spike-vs-stationary
# check: gamma4 after trimming the single most extreme residual), then run R random
# 70/30 splits comparing inner estimators on held-out test error (trimmed-RMSE 90%
# and MAE — robust to the test set's own heavy tail). KG-2 on (swob, tqa).
# Usage (from paper-4-weak-moment-gmdh/code): Rscript experiments/run_drilling_realdata.R [R]

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
                 else pkgload::load_all(PKG, quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 40L
METHODS <- c("LSE","Huber","L1","PMM2","WPMM2","WPMM3","auto-weak")
trmse <- function(e, p = 0.90) { q <- stats::quantile(abs(e), p, names = FALSE, na.rm = TRUE); sqrt(mean(e[abs(e) <= q]^2)) }

d <- load_volve_drilling(file.path(PKG, "data/volve_onbottom.csv"))
X <- as.data.frame(d$X); y <- as.numeric(d$y)
ord <- order(X$dept); X <- X[ord, ]; y <- y[ord]

# candidate depth bands (heavy+skew from the regime scan) + a clean-skew contrast
BANDS <- list(
  `heavy_skew_2000_2080` = c(2000, 2080),
  `heavy_skew_1740_1820` = c(1740, 1820),
  `heavy_skew_2110_2200` = c(2110, 2200),
  `mild_1360_1480`       = c(1360, 1480)
)

band_diag <- function(r) {
  cu <- sample_cumulants(r)
  r_trim <- r[-which.max(abs(r))]                 # drop single most extreme residual
  cu_t <- sample_cumulants(r_trim)
  list(g3 = cu$gamma3, g4 = cu$gamma4, g4_trim = cu_t$gamma4,
       tail_ratio = stats::sd(r) / stats::mad(r, constant = 1.4826))
}

fit_pred <- function(method, tr, te) {
  ctrl <- gmdh_pmm_control(B = if (method == "auto-weak") 60L else 0L,
                           force_method = method, max_iter = 60L, weak_sigma_mult = 2.5)
  est <- tryCatch(inner_estimate(tr$swob, tr$tqa, tr$y, ctrl), error = function(e) NULL)
  if (is.null(est)) return(list(pred = rep(NA, nrow(te)), method = NA))
  list(pred = kg2_predict(est$theta, te$swob, te$tqa), method = est$method)
}

cat(sprintf("=== Real Volve drilling, in-distribution estimator comparison | %d splits | KG-2(swob,tqa) ===\n", R))
summary_rows <- list()
for (bn in names(BANDS)) {
  rng <- BANDS[[bn]]; sel <- X$dept >= rng[1] & X$dept < rng[2]
  Xi <- X[sel, ]; yi <- y[sel]; n <- length(yi)
  if (n < 120) { cat(sprintf("\n[%s] depth %g-%g: only n=%d, skip\n", bn, rng[1], rng[2], n)); next }
  fit0 <- stats::lm(yi ~ swob + tqa, data = Xi); bd <- band_diag(stats::residuals(fit0))
  cat(sprintf("\n[%s] depth %g-%g | n=%d | gamma3=%+.2f gamma4=%.1f (trim->%.1f) tail_ratio=%.2f\n",
      bn, rng[1], rng[2], n, bd$g3, bd$g4, bd$g4_trim, bd$tail_ratio))
  spike <- bd$g4 > 8 && bd$g4_trim < 0.4 * bd$g4
  cat(sprintf("   stationarity: %s\n", if (spike) "SPIKE-dominated (g4 collapses on trim) — not cleanly stationary" else "bulk heavy tail (g4 robust to trimming) — usable"))

  band <- data.frame(swob = Xi$swob, tqa = Xi$tqa, y = yi)
  mat <- matrix(NA, R, length(METHODS), dimnames = list(NULL, METHODS))
  methsel <- character(0)
  for (r in seq_len(R)) {
    set.seed(9000 + r)
    tr_idx <- sample.int(n, floor(0.7 * n)); te_idx <- setdiff(seq_len(n), tr_idx)
    tr <- band[tr_idx, ]; te <- band[te_idx, ]
    for (m in METHODS) {
      fp <- fit_pred(m, tr, te); e <- te$y - fp$pred
      mat[r, m] <- trmse(e)
      if (m == "auto-weak" && !is.na(fp$method)) methsel <- c(methsel, fp$method)
    }
  }
  med <- apply(mat, 2, median, na.rm = TRUE); ord_m <- order(med)
  cat("   test trimmed-RMSE median (lower=better):\n")
  for (m in METHODS[ord_m]) cat(sprintf("     %-10s %.4f\n", m, med[m]))
  # paired tests: auto-weak & WPMM2 vs the robust/classical references
  for (pair in list(c("auto-weak","LSE"), c("auto-weak","Huber"), c("auto-weak","L1"),
                    c("auto-weak","PMM2"), c("WPMM2","Huber"))) {
    a <- mat[, pair[1]]; b <- mat[, pair[2]]; ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 5) next
    w <- suppressWarnings(stats::wilcox.test(a[ok], b[ok], paired = TRUE))
    cat(sprintf("     %-9s vs %-5s: %+.1f%% median, win %2.0f%%, p=%.2e\n",
        pair[1], pair[2], 100*(median(a,na.rm=TRUE)/median(b,na.rm=TRUE)-1), 100*mean(b[ok]>a[ok]), w$p.value))
  }
  if (length(methsel)) cat(sprintf("   auto-weak dispatched: %s\n",
      paste(sprintf("%s:%d", names(table(methsel)), as.integer(table(methsel))), collapse=" ")))
  for (m in METHODS) summary_rows[[length(summary_rows)+1L]] <- data.frame(band=bn, method=m, trmse_med=med[m],
      g3=bd$g3, g4=bd$g4, stringsAsFactors=FALSE)
}
outdir <- "../results"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(do.call(rbind, summary_rows), file.path(outdir, "drilling_realdata.csv"), row.names = FALSE)
cat("\nSaved drilling_realdata.csv\n")
