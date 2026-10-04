#!/usr/bin/env Rscript
# HONEST BLOCKED-CV across ALL 12 vibration cells (FORGE 78B-32 Sanvean).
# The 12/12 replication and the headline win were RANDOM 70/30 on autocorrelated
# series (residual lag-1 ~0.99) -> leaked. This re-runs every cell under
# leakage-controlled blocked 5-fold (contiguous blocks, purge gap) with time
# order PRESERVED, and -- critically -- compares best-weak against the BEST
# baseline (min of LSE/Huber/L1), not only Huber, so a "win vs Huber" that is
# really "Huber broke, LSE won" is exposed. Random-split medians are reported
# side-by-side to quantify the inflation. Strictly drilling. No re-fit changes
# (a2 reconciliation was comment-only; estimator numerics unchanged).
# Usage: Rscript experiments/run_honest_blocked_cv.R

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))

trmse <- function(e, p = 0.90) { q <- stats::quantile(abs(e), p, names = FALSE, na.rm = TRUE); sqrt(mean(e[abs(e) <= q]^2)) }
fitm  <- function(m, b, r, yy) tryCatch(inner_estimate(b, r, yy,
            gmdh_pmm_control(B = 0, force_method = m, max_iter = 60L, weak_sigma_mult = 2.5))$theta,
            error = function(e) rep(NA, 6))
METH  <- c("LSE", "Huber", "L1", "WPMM2", "WPMM3")
evpt  <- function(b, r, yy, tri, tei) sapply(METH, function(m)
            trmse(yy[tei] - kg2_predict(fitm(m, b[tri], r[tri], yy[tri]), b[tei], r[tei])))

cap <- 6000L; gap <- 50L; folds <- 5L; TMP <- "/tmp/forge78b_sensor"
.find1 <- function(pat) { x <- list.files(TMP, pattern = pat, recursive = TRUE, full.names = TRUE); if (length(x)) x[1] else NA }
RUNS <- list(`10.625` = .find1("Forge 10\\.625 - output data\\.csv$"),
             `14.75-Run1` = .find1("Run 1 - output data\\.csv$"),
             `14.75-Run2` = .find1("Run 2 - output data\\.csv$"))
SENSORS <- c("CSS-008", "CSS-007"); CHANS <- c("ShZpeak", "ShYpeak", "ShZrms")

cell <- function(d, sensor, chan) {
  col <- function(s) suppressWarnings(as.numeric(d[[paste0(sensor, "_", s)]]))
  ss <- col("StickSlip(%)"); gs <- col("GyroXspread(RPM)"); tg <- col(paste0(chan, "(g)"))
  if (is.null(ss) || is.null(gs) || is.null(tg)) return(NULL)
  ok <- is.finite(ss) & is.finite(gs) & is.finite(tg) & tg > 0       # preserves time order
  b <- ss[ok]; r <- gs[ok]; y <- log(tg[ok]); n <- length(y)
  if (n < 600) return(NULL)
  if (n > cap) { b <- b[1:cap]; r <- r[1:cap]; y <- y[1:cap]; n <- cap }  # contiguous head
  cu <- sample_cumulants(resid(lm(y ~ b + r)))
  ac1 <- stats::acf(resid(lm(y ~ b + r)), lag.max = 1, plot = FALSE)$acf[2]
  # BLOCKED 5-fold
  fb <- floor(n / folds)
  mb <- t(sapply(1:folds, function(k) { a <- (k - 1) * fb + 1; z <- if (k == folds) n else k * fb
    tei <- a:z; tri <- setdiff(1:n, max(1, a - gap):min(n, z + gap))
    if (length(tri) < 100 || length(tei) < 40) return(rep(NA, length(METH))); evpt(b, r, y, tri, tei) }))
  bl <- apply(mb, 2, median, na.rm = TRUE)
  # RANDOM 70/30 (for the inflation comparison)
  mr <- t(sapply(1:25, function(s) { set.seed(90000 + s); tri <- sample.int(n, floor(.7 * n)); evpt(b, r, y, tri, setdiff(1:n, tri)) }))
  rn <- apply(mr, 2, median, na.rm = TRUE)
  list(n = n, g3 = cu$gamma3, g4 = cu$gamma4, ac1 = ac1, bl = bl, rn = rn)
}

cat(sprintf("=== HONEST blocked-5fold vs random-70/30 | %d cells | best-weak vs BEST baseline ===\n", 12L))
cat(sprintf("%-10s %-7s %-7s %5s %5s %5s | %-32s | %-26s\n",
            "run","sensor","chan","g3","g4","ac1","BLOCKED  LSE Hub  L1  W2  W3","best-weak vs: LSE  Hub  BEST"))
rows <- list(); surv_best <- 0L; surv_huber <- 0L; tot <- 0L
for (rn0 in names(RUNS)) { f <- RUNS[[rn0]]; if (is.na(f)) { cat(sprintf("  (run %s missing)\n", rn0)); next }
  d <- data.table::fread(f, skip = 35L, header = TRUE, fill = TRUE, showProgress = FALSE)
  for (sn in SENSORS) for (ch in CHANS) {
    res <- tryCatch(cell(d, sn, ch), error = function(e) NULL); if (is.null(res)) next
    tot <- tot + 1L
    bw  <- min(res$bl["WPMM2"], res$bl["WPMM3"]); bwn <- names(which.min(res$bl[c("WPMM2","WPMM3")]))
    base_best <- min(res$bl["LSE"], res$bl["Huber"], res$bl["L1"]); base_name <- names(which.min(res$bl[c("LSE","Huber","L1")]))
    vs_lse  <- 100 * (bw / res$bl["LSE"]   - 1)
    vs_hub  <- 100 * (bw / res$bl["Huber"] - 1)
    vs_best <- 100 * (bw / base_best       - 1)
    if (vs_best < -2) surv_best  <- surv_best  + 1L
    if (vs_hub  < -2) surv_huber <- surv_huber + 1L
    cat(sprintf("%-10s %-7s %-7s %+5.1f %5.1f %5.2f | %5.3f %5.3f %5.3f %5.3f %5.3f | %+5.1f %+5.1f %+5.1f (best=%s)\n",
        rn0, sn, ch, res$g3, res$g4, res$ac1,
        res$bl["LSE"], res$bl["Huber"], res$bl["L1"], res$bl["WPMM2"], res$bl["WPMM3"],
        vs_lse, vs_hub, vs_best, base_name))
    rows[[length(rows)+1L]] <- data.frame(run=rn0, sensor=sn, chan=ch, n=res$n, g3=res$g3, g4=res$g4, ac1=res$ac1,
      bl_LSE=res$bl["LSE"], bl_Huber=res$bl["Huber"], bl_L1=res$bl["L1"], bl_WPMM2=res$bl["WPMM2"], bl_WPMM3=res$bl["WPMM3"],
      rn_LSE=res$rn["LSE"], rn_Huber=res$rn["Huber"], rn_L1=res$rn["L1"], rn_WPMM2=res$rn["WPMM2"], rn_WPMM3=res$rn["WPMM3"],
      bestweak=bwn, blocked_vs_lse=vs_lse, blocked_vs_huber=vs_hub, blocked_vs_best=vs_best,
      random_vs_huber=100*(min(res$rn["WPMM2"],res$rn["WPMM3"])/res$rn["Huber"]-1),
      base_best=base_name, row.names=NULL, stringsAsFactors=FALSE)
  }
}
cat(sprintf("\nUnder BLOCKED CV: best-weak beats the BEST baseline (>2%%) in %d / %d cells; beats Huber in %d / %d.\n",
            surv_best, tot, surv_huber, tot))
outdir <- file.path(dirname(dirname(PKG)), "paper-4-weak-moment-gmdh", "results")  # canonical, cwd-independent
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(do.call(rbind, rows), file.path(outdir, "honest_blocked_cv_vibration.csv"), row.names = FALSE)
cat("Saved honest_blocked_cv_vibration.csv\n")
