#!/usr/bin/env Rscript
# EXPERIMENT 3 -- MOVING-BLOCK-BOOTSTRAP inference for the relative win statistic.
#
# Referee objection: "Paired Wilcoxon over autocorrelated random splits is invalid;
# effective n is far below nominal; you report no CIs." The committed honest-CV
# numbers give a point estimate of the relative win per cell but NO uncertainty
# that respects the residual autocorrelation (lag-1 ~0.98-0.99 vibration; depth-
# autocorrelated ROP). Here, for the KEY win cells ONLY, we:
#   1. Build the SAME leakage-controlled blocked split (contiguous blocks, purge
#      gap=50 vibration / 30 ROP, cap=6000 vibration head; Volve depth-ordered,
#      same 4 bands) and fit each method ONCE on the blocked TRAIN set.
#   2. Evaluate on a CONTIGUOUS held-out test block, recording per-test-point
#      errors in TIME order (so the error series keeps its autocorrelation).
#   3. MOVING-BLOCK-BOOTSTRAP the test-error series with block length L (a few x
#      the residual autocorrelation length; default L=75). Blocks are resampled
#      ONCE per replicate and applied to ALL methods (paired), so the relative
#      win delta = (best-weak trmse / best-baseline trmse - 1) is computed on the
#      same resampled time indices -- honest, autocorrelation-respecting inference.
#   4. Report 95% percentile CI on each method's trmse and on the relative win,
#      and whether the win CI excludes 0.
#
# best-weak  = min over {WPMM2, WPMM3} (point-estimate selection on the FULL test
#              block, frozen across bootstrap replicates -- the method identity is
#              not re-selected per replicate, matching the parent scripts' which.min).
# best-baseline = min over {LSE, Huber, L1} (same freezing rule).
#
# trmse and all data loading / band definitions mirror run_honest_blocked_cv.R and
# run_honest_blocked_cv_rop.R EXACTLY so numbers are directly comparable.
# Strictly drilling. No estimator numeric changes.
# Usage: Rscript experiments/run_blockboot_ci.R [L Bboot]

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))
stopifnot(requireNamespace("boot", quietly = TRUE))

args <- commandArgs(trailingOnly = TRUE)
L     <- if (length(args) >= 1) as.integer(args[1]) else 75L     # moving block length
Bboot <- if (length(args) >= 2) as.integer(args[2]) else 2000L   # bootstrap replicates

trmse <- function(e, p = 0.90) { q <- stats::quantile(abs(e), p, names = FALSE, na.rm = TRUE); sqrt(mean(e[abs(e) <= q]^2)) }
METH  <- c("LSE", "Huber", "L1", "WPMM2", "WPMM3")
WEAK  <- c("WPMM2", "WPMM3"); BASE <- c("LSE", "Huber", "L1")

# residual autocorrelation length: smallest lag where |acf| first drops below 1/e
acf_len <- function(e) {
  a <- stats::acf(e, lag.max = min(400L, length(e) - 1L), plot = FALSE)$acf[-1]
  w <- which(abs(a) < exp(-1)); if (length(w)) w[1] else length(a)
}

# Fold-stratified MOVING-BLOCK-BOOTSTRAP CI on per-method trmse and on the relative win.
# EE is a LIST of per-fold (n_k x 5) TIME-ORDERED error matrices (columns = METH),
# one entry per contiguous blocked test fold. Point estimate = MEDIAN across folds of
# per-method trmse (mirrors run_honest_blocked_cv*.R apply(mb,2,median)), so the win
# point estimate reproduces the committed honest-CV win. The bootstrap resamples
# moving blocks WITHIN each contiguous fold (blocks never cross a fold boundary),
# applied PAIRED to all methods, then takes the per-method median trmse across folds.
mbb_ci <- function(EE, L, Bboot, seed = 12345L) {
  per_fold_tr <- function(Elist) {
    m <- sapply(Elist, function(E) apply(E, 2, trmse))   # 5(methods) x folds
    apply(m, 1, median, na.rm = TRUE)                    # median across folds
  }
  tr_full <- per_fold_tr(EE); names(tr_full) <- METH
  bw_name <- WEAK[which.min(tr_full[WEAK])]
  bb_name <- BASE[which.min(tr_full[BASE])]
  delta_pt <- tr_full[bw_name] / tr_full[bb_name] - 1
  # precompute per-fold block start pools
  pools <- lapply(EE, function(E) { n <- nrow(E); Lk <- min(L, max(1L, floor(n / 4))); nb <- ceiling(n / Lk)
    list(n = n, Lk = Lk, nb = nb, starts = 1:(n - Lk + 1)) })
  set.seed(seed)
  boot_tr <- matrix(NA_real_, Bboot, length(METH)); colnames(boot_tr) <- METH
  boot_delta <- numeric(Bboot)
  for (bI in seq_len(Bboot)) {
    foldtr <- sapply(seq_along(EE), function(k) {
      E <- EE[[k]]; p <- pools[[k]]
      st <- if (length(p$starts) > 1) sample(p$starts, p$nb, replace = TRUE) else rep(1L, p$nb)
      idx <- unlist(lapply(st, function(s) s:(s + p$Lk - 1)))[1:p$n]
      apply(E[idx, , drop = FALSE], 2, trmse)            # 5 methods, this fold
    })                                                    # 5 x folds
    trb <- apply(foldtr, 1, median, na.rm = TRUE); names(trb) <- METH
    boot_tr[bI, ] <- trb
    boot_delta[bI] <- trb[bw_name] / trb[bb_name] - 1     # frozen method identities (paired)
  }
  ci_tr <- apply(boot_tr, 2, stats::quantile, probs = c(.025, .975), na.rm = TRUE)
  ci_d  <- stats::quantile(boot_delta, c(.025, .975), na.rm = TRUE, names = FALSE)
  list(tr_full = tr_full, bw_name = bw_name, bb_name = bb_name,
       delta_pt = as.numeric(delta_pt), ci_tr = ci_tr,
       ci_lo = ci_d[1], ci_hi = ci_d[2], boot_delta = boot_delta)
}

# Fit all 5 methods on (b_tr,r_tr,y_tr); return matrix of time-ordered test errors.
fit_errors <- function(b, r, y, tri, tei) {
  sapply(METH, function(m) {
    th <- tryCatch(inner_estimate(b[tri], r[tri], y[tri],
            gmdh_pmm_control(B = 0, force_method = m, max_iter = 60L, weak_sigma_mult = 2.5))$theta,
            error = function(e) rep(NA, 6))
    y[tei] - kg2_predict(th, b[tei], r[tei])
  })
}

rows <- list()

report <- function(tag, gamma4, b, r, y, gap, seed) {
  n <- length(y); folds <- 5L; fb <- floor(n / folds)
  # Build ALL 5 contiguous blocked folds EXACTLY as run_honest_blocked_cv*.R:
  #   tei = a:z (contiguous), tri = everything except [a-gap, z+gap] (purge).
  # Each fold's test errors stay in TIME order -> a valid block-bootstrap series.
  EE <- list(); n_te_tot <- 0L; acl_acc <- integer(0)
  for (k in 1:folds) {
    a <- (k - 1L) * fb + 1L; z <- if (k == folds) n else k * fb
    tei <- a:z; tri <- setdiff(1:n, max(1, a - gap):min(n, z + gap))
    if (length(tri) < 100 || length(tei) < 40) next
    E <- fit_errors(b, r, y, tri, tei)
    keep <- stats::complete.cases(E); E <- E[keep, , drop = FALSE]
    if (nrow(E) < 40) next
    EE[[length(EE) + 1L]] <- E; n_te_tot <- n_te_tot + nrow(E)
    acl_acc <- c(acl_acc, acf_len(E[, "LSE"]))
  }
  if (length(EE) < 2) { cat(sprintf("%-26s (insufficient folds)\n", tag)); return(invisible(NULL)) }
  acl <- as.integer(round(median(acl_acc)))
  ci <- mbb_ci(EE, L, Bboot, seed = seed)
  win_pct <- 100 * ci$delta_pt
  ci_lo_pct <- 100 * ci$ci_lo; ci_hi_pct <- 100 * ci$ci_hi
  excl <- (ci_lo_pct < 0 && ci_hi_pct < 0) || (ci_lo_pct > 0 && ci_hi_pct > 0)
  cat(sprintf("%-26s g4=%6.2f n_te=%4d nf=%d L=%3d acl=%3d | bw=%-5s(%.4f) bb=%-5s(%.4f) | win=%+6.1f%% CI[%+6.1f, %+6.1f] excl0=%s\n",
      tag, gamma4, n_te_tot, length(EE), L, acl, ci$bw_name, ci$tr_full[ci$bw_name],
      ci$bb_name, ci$tr_full[ci$bb_name], win_pct, ci_lo_pct, ci_hi_pct, excl))
  rows[[length(rows) + 1L]] <<- data.frame(
    cell = tag, gamma4 = gamma4, n_test = n_te_tot, n_folds = length(EE), block_L = L, acf_len = acl,
    best_weak = ci$bw_name, best_base = ci$bb_name,
    trmse_weak = as.numeric(ci$tr_full[ci$bw_name]), trmse_base = as.numeric(ci$tr_full[ci$bb_name]),
    win_pct = win_pct, ci_lo = ci_lo_pct, ci_hi = ci_hi_pct,
    excludes_zero = excl, row.names = NULL, stringsAsFactors = FALSE)
  invisible(NULL)
}

cat(sprintf("=== EXP3 moving-block-bootstrap CI on relative win | L=%d Bboot=%d ===\n", L, Bboot))
cat("win = (best-weak trmse / best-baseline trmse - 1); negative = weak wins; excl0 => significant\n\n")

## ---------- VIBRATION KEY CELLS (FORGE 78B-32) ----------
cap <- 6000L; gapV <- 50L; TMP <- "/tmp/forge78b_sensor"
f_run1 <- { x <- list.files(TMP, pattern = "Run 1 - output data\\.csv$", recursive = TRUE, full.names = TRUE); if (length(x)) x[1] else NA }
vib_cell <- function(d, sensor, chan) {
  col <- function(s) suppressWarnings(as.numeric(d[[paste0(sensor, "_", s)]]))
  ss <- col("StickSlip(%)"); gs <- col("GyroXspread(RPM)"); tg <- col(paste0(chan, "(g)"))
  if (is.null(ss) || is.null(gs) || is.null(tg)) return(NULL)
  ok <- is.finite(ss) & is.finite(gs) & is.finite(tg) & tg > 0   # preserves time order
  b <- ss[ok]; r <- gs[ok]; y <- log(tg[ok]); n <- length(y)
  if (n < 600) return(NULL)
  if (n > cap) { b <- b[1:cap]; r <- r[1:cap]; y <- y[1:cap] }    # contiguous head
  cu <- sample_cumulants(resid(lm(y ~ b + r)))
  list(b = b, r = r, y = y, g4 = cu$gamma4)
}
if (!is.na(f_run1)) {
  d <- data.table::fread(f_run1, skip = 35L, header = TRUE, fill = TRUE, showProgress = FALSE)
  # KEY heavy-tailed win cells (residual gamma4 > ~5) named in the task + neighbour
  KEY_VIB <- list(c("CSS-008","ShZpeak"), c("CSS-007","ShZpeak"),
                  c("CSS-007","ShYpeak"), c("CSS-008","ShYpeak"))
  for (kv in KEY_VIB) {
    cc <- tryCatch(vib_cell(d, kv[1], kv[2]), error = function(e) NULL); if (is.null(cc)) next
    report(sprintf("14.75-Run1_%s_%s", kv[1], kv[2]), cc$g4, cc$b, cc$r, cc$y, gapV,
           seed = 70000L + sum(utf8ToInt(paste0(kv, collapse = ""))))
  }
} else cat("(FORGE Run 1 missing -- vibration cells skipped)\n")

## ---------- VOLVE HEAVY ROP BANDS ----------
gapR <- 30L
dv <- tryCatch(load_volve_drilling(file.path(PKG, "data/volve_onbottom.csv")), error = function(e) NULL)
if (!is.null(dv)) {
  X <- as.data.frame(dv$X); yv <- as.numeric(dv$y); ord <- order(X$dept); X <- X[ord, ]; yv <- yv[ord]
  KEY_ROP <- list(volve_2000_2080 = c(2000,2080), volve_1740_1820 = c(1740,1820),
                  volve_mild_1360_1480 = c(1360,1480))
  for (bn in names(KEY_ROP)) {
    rg <- KEY_ROP[[bn]]; sel <- X$dept >= rg[1] & X$dept < rg[2]
    Xi <- X[sel, ]; yi <- yv[sel]; if (length(yi) < 150) next
    cu <- sample_cumulants(resid(lm(yi ~ Xi$swob + Xi$tqa)))
    report(bn, cu$gamma4, Xi$swob, Xi$tqa, yi, gapR, seed = 40000L + as.integer(rg[1]))
  }
} else cat("(Volve load failed -- ROP bands skipped)\n")

## ---------- save ----------
out <- do.call(rbind, rows)
# columns required by the task: cell, gamma4, win_pct, ci_lo, ci_hi, excludes_zero (+ context)
outdir <- file.path(dirname(dirname(PKG)), "paper-4-weak-moment-gmdh", "results")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(out, file.path(outdir, "blockboot_ci.csv"), row.names = FALSE)
cat(sprintf("\nCells with win CI excluding 0 (genuinely significant): %d / %d\n",
            sum(out$excludes_zero), nrow(out)))
cat("Saved blockboot_ci.csv to", outdir, "\n")
