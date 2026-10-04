#!/usr/bin/env Rscript
# Real-data check #2 for Paper 4: FORGE 56-32 (granite geothermal, n~79k, globally
# heavy-tailed-skewed: gamma3=-1.7, gamma4=13.9). Auto-finds stationary heavy+skew
# depth bands, then runs the same in-distribution estimator comparison as the Volve
# check (random 70/30 splits, test trimmed-RMSE, KG-2 on wob+rpm). Question: do
# Huber/LAD win here too (negative robust across two real wells), or is there a
# genuine WPMM2 niche in impulsive granite drilling?
# Usage (from paper-4-weak-moment-gmdh/code): Rscript experiments/run_forge56_realdata.R [R]

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
F1 <- "wob"; F2 <- "rpm"

d <- load_forge_standard_drilling(well = "56-32", cadence = "10sec")
X <- as.data.frame(d$X); y <- as.numeric(d$y)
ord <- order(X$dept); X <- X[ord, ]; y <- y[ord]
n <- length(y); cat(sprintf("FORGE 56-32 10sec loaded: n=%d\n", n))

# --- window scan for heavy+skew bands ---
win <- 800L; step <- 400L; starts <- seq(1L, n - win, by = step)
scan <- do.call(rbind, lapply(starts, function(s) {
  idx <- s:(s + win - 1L); Xi <- X[idx, c(F1, F2), drop = FALSE]; yi <- y[idx]
  fit <- tryCatch(stats::lm(yi ~ ., data = Xi), error = function(e) NULL); if (is.null(fit)) return(NULL)
  r <- stats::residuals(fit); cu <- sample_cumulants(r)
  r_trim <- r[-which.max(abs(r))]; g4t <- sample_cumulants(r_trim)$gamma4
  data.frame(start = s, depth0 = X$dept[idx[1]], depth1 = X$dept[idx[win]],
             g3 = cu$gamma3, g4 = cu$gamma4, g4_trim = g4t,
             tail_ratio = stats::sd(r)/stats::mad(r, constant = 1.4826))
}))
heavy <- scan[abs(scan$g3) > 1 & scan$g4 > 6 & scan$g4_trim > 0.4 * scan$g4, ]   # heavy+skew, not single-spike
heavy <- heavy[order(-heavy$g4), ]
cat(sprintf("scan: %d windows, %d heavy+skew (|g3|>1 & g4>6 & not single-spike)\n", nrow(scan), nrow(heavy)))
if (!nrow(heavy)) { cat("no usable heavy+skew band\n"); quit(save="no") }
print(utils::head(heavy[, c("depth0","depth1","g3","g4","g4_trim","tail_ratio")], 8), row.names = FALSE, digits = 3)

fit_pred <- function(method, tr, te) {
  ctrl <- gmdh_pmm_control(B = if (method == "auto-weak") 60L else 0L,
                           force_method = method, max_iter = 60L, weak_sigma_mult = 2.5)
  est <- tryCatch(inner_estimate(tr[[F1]], tr[[F2]], tr$y, ctrl), error = function(e) NULL)
  if (is.null(est)) return(list(pred = rep(NA, nrow(te)), method = NA))
  list(pred = kg2_predict(est$theta, te[[F1]], te[[F2]]), method = est$method)
}

# take up to 3 distinct heavy+skew bands (widen each window to a band of its rows)
bands <- utils::head(heavy, 3)
summary_rows <- list()
for (bi in seq_len(nrow(bands))) {
  b <- bands[bi, ]; sel <- X$dept >= b$depth0 & X$dept <= b$depth1
  Xi <- X[sel, ]; yi <- y[sel]; nb <- length(yi); if (nb < 120) next
  band <- data.frame(wob = Xi[[F1]], rpm = Xi[[F2]], y = yi); names(band) <- c(F1, F2, "y")
  cat(sprintf("\n[band %d] depth %.0f-%.0f | n=%d | g3=%+.2f g4=%.1f (trim->%.1f) tail_ratio=%.2f\n",
      bi, b$depth0, b$depth1, nb, b$g3, b$g4, b$g4_trim, b$tail_ratio))
  mat <- matrix(NA, R, length(METHODS), dimnames = list(NULL, METHODS)); methsel <- character(0)
  for (r in seq_len(R)) {
    set.seed(13000 + 100*bi + r)
    tr_idx <- sample.int(nb, floor(0.7*nb)); te_idx <- setdiff(seq_len(nb), tr_idx)
    tr <- band[tr_idx, ]; te <- band[te_idx, ]
    for (m in METHODS) { fp <- fit_pred(m, tr, te); mat[r, m] <- trmse(te$y - fp$pred)
      if (m == "auto-weak" && !is.na(fp$method)) methsel <- c(methsel, fp$method) }
  }
  med <- apply(mat, 2, median, na.rm = TRUE)
  cat("   test trimmed-RMSE median:\n"); for (m in METHODS[order(med)]) cat(sprintf("     %-10s %.4f\n", m, med[m]))
  for (pair in list(c("auto-weak","Huber"), c("auto-weak","L1"), c("auto-weak","LSE"), c("WPMM2","Huber"), c("WPMM3","Huber"))) {
    a <- mat[, pair[1]]; b2 <- mat[, pair[2]]; ok <- is.finite(a) & is.finite(b2); if (sum(ok) < 5) next
    w <- suppressWarnings(stats::wilcox.test(a[ok], b2[ok], paired = TRUE))
    cat(sprintf("     %-9s vs %-5s: %+.1f%% median, win %2.0f%%, p=%.2e\n",
        pair[1], pair[2], 100*(median(a,na.rm=TRUE)/median(b2,na.rm=TRUE)-1), 100*mean(b2[ok]>a[ok]), w$p.value))
  }
  if (length(methsel)) cat(sprintf("   auto-weak dispatched: %s\n",
      paste(sprintf("%s:%d", names(table(methsel)), as.integer(table(methsel))), collapse=" ")))
  for (m in METHODS) summary_rows[[length(summary_rows)+1L]] <- data.frame(band=bi, depth0=b$depth0, method=m,
      trmse_med=med[m], g3=b$g3, g4=b$g4, stringsAsFactors=FALSE)
}
outdir <- "../results"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(do.call(rbind, summary_rows), file.path(outdir, "forge56_realdata.csv"), row.names = FALSE)
cat("\nSaved forge56_realdata.csv\n")
