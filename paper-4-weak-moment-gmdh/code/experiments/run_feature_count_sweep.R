#!/usr/bin/env Rscript
# EXPERIMENT 4 -- FEATURE-COUNT SWEEP (commit the mechanism).
# Referee objection: "the weak win is an under-specification artifact; the k=2/k=3
# sweep is not committed." This script answers it on the SAME honest blocked-CV
# protocol as run_honest_blocked_cv.R (contiguous 5 blocks, purge gap=50, cap=6000
# contiguous head, time order preserved, identical trmse). For a vibration cell
# where weak wins at k=2, we add a THIRD same-sensor channel and watch the residual
# skew/kurtosis AND the best-weak-vs-best-baseline advantage move together.
#
# Feature sets (all KG-2 base = StickSlip + GyroXspread; extra channel appended as a
# linear column to the design before the inner estimator):
#   k2            : KG-2(StickSlip, GyroXspread)                  -- the current model
#   k3_gyromed    : + GyroXmed(RPM)        (rotary-speed LEVEL; WEAK residual predictor)
#   k3_rms        : + log <other-axis>rms  (vibration-energy; GENUINELY INFORMATIVE)
# The gyro term is the "richness without information" control; the rms term is the
# "richness WITH information" case. Honest expected pattern: only the INFORMATIVE
# channel collapses residual gamma3 toward 0 and the weak advantage toward 0.
#
# Cells: 14.75-Run1 CSS-007 ShZpeak (clean +skew, weak wins ~27% at k=2) and
#        14.75-Run1 CSS-008 ShZpeak (heavy -skew g4~8.7, weak wins ~19% at k=2).
# Strictly drilling. No estimator-numerics change. Usage: Rscript experiments/run_feature_count_sweep.R

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))

# --- identical metric & inner-estimator harness as run_honest_blocked_cv.R ---
trmse <- function(e, p = 0.90) { q <- stats::quantile(abs(e), p, names = FALSE, na.rm = TRUE); sqrt(mean(e[abs(e) <= q]^2)) }
fitm  <- function(m, b, r, yy) tryCatch(inner_estimate(b, r, yy,
            gmdh_pmm_control(B = 0, force_method = m, max_iter = 60L, weak_sigma_mult = 2.5))$theta,
            error = function(e) rep(NA, 6))
METH  <- c("LSE", "Huber", "L1", "WPMM2", "WPMM3")

# KG-2 inner estimator works on TWO regressors (b, r). To extend the model to k>=3
# WITHOUT changing the inner estimator, we residualize the response on the extra
# channel(s) out-of-fold (linear projection fit on the TRAIN rows only -- no leakage),
# then run the SAME KG-2(b,r) + same estimators on the residualized response. This is
# the standard "partial out the extra term" route: the k>=3 model is
#   y = beta0 + beta_extra*X_extra + KG2(b,r) ,
# estimated by first removing the (train-fit) linear extra part, then fitting KG-2 on
# the remainder, exactly as a richer-feature GMDH node would absorb the extra channel.
# evpt returns, per method, the trmse on the held-out block of the FINAL prediction
# (extra linear part + KG-2 part), so it is directly comparable to the k=2 numbers.
evpt <- function(b, r, yy, X, tri, tei) {
  if (is.null(X)) {                      # k = 2 : no extra channel
    yy_tr <- yy[tri]
    sapply(METH, function(m)
      trmse(yy[tei] - kg2_predict(fitm(m, b[tri], r[tri], yy_tr), b[tei], r[tei])))
  } else {                               # k >= 3 : partial out extra linear channel(s) on TRAIN only
    Xtr <- X[tri, , drop = FALSE]; Xte <- X[tei, , drop = FALSE]
    lf  <- stats::lm.fit(cbind(1, Xtr), yy[tri])           # train-only linear fit of extra part
    cf  <- lf$coefficients; cf[is.na(cf)] <- 0
    pred_extra_tr <- as.numeric(cbind(1, Xtr) %*% cf)
    pred_extra_te <- as.numeric(cbind(1, Xte) %*% cf)
    yres_tr <- yy[tri] - pred_extra_tr                     # KG-2 fits the remainder
    sapply(METH, function(m) {
      th <- fitm(m, b[tri], r[tri], yres_tr)
      pred_te <- pred_extra_te + kg2_predict(th, b[tei], r[tei])
      trmse(yy[tei] - pred_te)
    })
  }
}

cap <- 6000L; gap <- 50L; folds <- 5L; TMP <- "/tmp/forge78b_sensor"
.find1 <- function(pat) { x <- list.files(TMP, pattern = pat, recursive = TRUE, full.names = TRUE); if (length(x)) x[1] else NA }
RUN1 <- .find1("Run 1 - output data\\.csv$")

# Residual cumulants of the FULL-model (in-sample, capped head) residual: this is the
# diagnostic the dispatch sees. For k>=3 the residual is from lm(y ~ b + r + X).
resid_cumulants <- function(b, r, yy, X) {
  if (is.null(X)) cu <- sample_cumulants(resid(stats::lm(yy ~ b + r)))
  else            cu <- sample_cumulants(resid(stats::lm(yy ~ b + r + X)))
  c(g3 = cu$gamma3, g4 = cu$gamma4)
}

# Build one (sensor, chan) cell on the capped contiguous head, returning aligned vectors.
build_cell <- function(d, sensor, chan) {
  col <- function(s) suppressWarnings(as.numeric(d[[paste0(sensor, "_", s)]]))
  ss <- col("StickSlip(%)"); gs <- col("GyroXspread(RPM)"); tg <- col(paste0(chan, "(g)"))
  gmed <- col("GyroXmed(RPM)")
  # informative same-sensor energy channel = a DIFFERENT-axis RMS shock (smooth co-driver
  # of total vibration energy; not the response axis/statistic). ShZpeak -> ShXrms.
  rms_name <- if (grepl("ShZ", chan)) "ShXrms(g)" else if (grepl("ShY", chan)) "ShXrms(g)" else "ShYrms(g)"
  rmsc <- col(rms_name)
  ok <- is.finite(ss) & is.finite(gs) & is.finite(tg) & tg > 0     # preserves time order; same ok-filter as honest CV
  b <- ss[ok]; r <- gs[ok]; y <- log(tg[ok]); gmed <- gmed[ok]; rmsc <- rmsc[ok]
  n <- length(y); if (n < 600) return(NULL)
  if (n > cap) { idx <- 1:cap; b <- b[idx]; r <- r[idx]; y <- y[idx]; gmed <- gmed[idx]; rmsc <- rmsc[idx]; n <- cap }
  # log-transform the rms energy channel (positive), guard non-finite by column mean
  lrms <- log(pmax(rmsc, 1e-6)); lrms[!is.finite(lrms)] <- mean(lrms[is.finite(lrms)])
  gmed[!is.finite(gmed)] <- mean(gmed[is.finite(gmed)])
  list(b = b, r = r, y = y, n = n, gmed = gmed, lrms = lrms, rms_name = rms_name)
}

# blocked 5-fold best-weak-vs-best-baseline for one feature set (X = extra design or NULL)
sweep_blocked <- function(cl, X) {
  b <- cl$b; r <- cl$r; y <- cl$y; n <- cl$n
  fb <- floor(n / folds)
  mb <- t(sapply(1:folds, function(k) { a <- (k - 1) * fb + 1; z <- if (k == folds) n else k * fb
    tei <- a:z; tri <- setdiff(1:n, max(1, a - gap):min(n, z + gap))
    if (length(tri) < 100 || length(tei) < 40) return(rep(NA, length(METH)))
    evpt(b, r, y, X, tri, tei) }))
  bl <- apply(mb, 2, median, na.rm = TRUE); names(bl) <- METH
  bw  <- min(bl["WPMM2"], bl["WPMM3"])
  bwn <- names(which.min(bl[c("WPMM2", "WPMM3")]))
  base_best <- min(bl["LSE"], bl["Huber"], bl["L1"]); base_name <- names(which.min(bl[c("LSE", "Huber", "L1")]))
  cu <- resid_cumulants(b, r, y, X)
  list(bl = bl, bestweak = bwn, base_name = base_name,
       g3 = cu["g3"], g4 = cu["g4"], bestweak_vs_best_pct = 100 * (bw / base_best - 1))
}

cat("=== EXPERIMENT 4: feature-count sweep (commit the mechanism) | blocked-5fold, gap=50, cap=6000 ===\n")
cat(sprintf("%-9s %-7s %-7s %-11s %2s %+7s %6s | %5s %5s %5s %5s %5s | %-8s %-6s %s\n",
            "run","sensor","chan","feat_set","k","g3","g4","LSE","Hub","L1","W2","W3","bweak","base","bweak_vs_best%"))

rows <- list()
if (is.na(RUN1)) stop("14.75-Run1 output CSV not found under ", TMP)
d <- data.table::fread(RUN1, skip = 35L, header = TRUE, fill = TRUE, showProgress = FALSE)

CELLS <- list(c("CSS-007", "ShZpeak"), c("CSS-008", "ShZpeak"))
for (sc in CELLS) {
  sensor <- sc[1]; chan <- sc[2]
  cl <- build_cell(d, sensor, chan); if (is.null(cl)) { cat(sprintf("(cell %s %s skipped)\n", sensor, chan)); next }
  FEATS <- list(
    list(tag = "k2",         k = 2L, X = NULL),
    list(tag = "k3_gyromed", k = 3L, X = matrix(cl$gmed, ncol = 1)),  # richness w/o information (weak predictor)
    list(tag = paste0("k3_", sub("\\(g\\)$", "", cl$rms_name)), k = 3L, X = matrix(cl$lrms, ncol = 1))  # informative
  )
  for (fs in FEATS) {
    s <- sweep_blocked(cl, fs$X)
    cat(sprintf("%-9s %-7s %-7s %-11s %2d %+7.3f %6.2f | %5.3f %5.3f %5.3f %5.3f %5.3f | %-8s %-6s %+6.1f\n",
        "14.75-Run1", sensor, chan, fs$tag, fs$k, s$g3, s$g4,
        s$bl["LSE"], s$bl["Huber"], s$bl["L1"], s$bl["WPMM2"], s$bl["WPMM3"],
        s$bestweak, s$base_name, s$bestweak_vs_best_pct))
    rows[[length(rows) + 1L]] <- data.frame(
      run = "14.75-Run1", sensor = sensor, chan = chan,
      feat_set = fs$tag, k = fs$k,
      gamma3 = unname(s$g3), gamma4 = unname(s$g4),
      bl_LSE = unname(s$bl["LSE"]), bl_Huber = unname(s$bl["Huber"]), bl_L1 = unname(s$bl["L1"]),
      bl_WPMM2 = unname(s$bl["WPMM2"]), bl_WPMM3 = unname(s$bl["WPMM3"]),
      bestweak = s$bestweak, base_best = s$base_name,
      bestweak_vs_best_pct = unname(s$bestweak_vs_best_pct),
      row.names = NULL, stringsAsFactors = FALSE)
  }
}

outdir <- file.path(dirname(dirname(PKG)), "paper-4-weak-moment-gmdh", "results")  # canonical, cwd-independent
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(do.call(rbind, rows), file.path(outdir, "feature_count_sweep.csv"), row.names = FALSE)
cat(sprintf("\nSaved feature_count_sweep.csv to %s\n", outdir))
