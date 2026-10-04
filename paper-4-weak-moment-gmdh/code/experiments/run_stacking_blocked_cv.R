#!/usr/bin/env Rscript
# EXPERIMENT 2 -- SUPER-LEARNER / STACKING comparison, referee-proofing.
#
# Referee objection: "the dispatch IS a discrete super-learner; benchmark it
# against the discrete super-learner AND against convex stacking."
#
# Under the SAME honest leakage-controlled blocked 5-fold CV as
#   run_honest_blocked_cv.R      (vibration, 12 FORGE 78B-32 cells)
#   run_honest_blocked_cv_rop.R  (Volve ROP, 4 depth bands)
# we compare, per cell/band, the TEST trimmed-RMSE of:
#   (A) best individual METHOD chosen with hindsight (oracle floor; reference)
#   (B) best-weak  = min(WPMM2, WPMM3)            -- the program's weak family
#   (C) best single baseline = min(LSE, Huber, L1)
#   (D) auto-valgate = the DISCRETE validation-gated super-learner already in
#       the package (picks one candidate by inner held-out trimmed loss)
#   (E) convex-stack = NNLS-on-simplex blend of {LSE,Huber,L1,WPMM2,WPMM3},
#       weights fit on an inner held-out split of the TRAINING fold by
#       minimizing held-out trimmed loss (projected-gradient simplex solver,
#       implemented here -- nnls/quadprog are not installed).
#
# Honesty rules mirrored from the two source scripts:
#   * data loading, cell/band definitions, ok-filter (time order preserved),
#     cap=6000 contiguous head (vibration), gap=50 (vibration) / 30 (ROP),
#     folds=5 contiguous blocks, trmse(p=0.90) -- ALL identical, so numbers
#     are directly comparable to the committed honest baselines.
#   * The stacking inner split and auto-valgate are fit ONLY on training-fold
#     data; the test block is never touched during weight/method selection.
#
# KEY QUESTION: does convex stacking beat the simple validation-gated dispatch?
# If they tie, the honest claim is "our dispatch ~= a discrete super-learner
# inside GMDH, and convex stacking adds little." Report whichever is true.
# Usage: Rscript experiments/run_stacking_blocked_cv.R

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))

set.seed(20260531L)  # only governs the inner sub-split for stacking weight fit

# ---- IDENTICAL metric + single-method fit/predict helpers ----
trmse <- function(e, p = 0.90) { q <- stats::quantile(abs(e), p, names = FALSE, na.rm = TRUE); sqrt(mean(e[abs(e) <= q]^2)) }
METH  <- c("LSE", "Huber", "L1", "WPMM2", "WPMM3")
fitm  <- function(m, b, r, yy) tryCatch(inner_estimate(b, r, yy,
            gmdh_pmm_control(B = 0, force_method = m, max_iter = 60L, weak_sigma_mult = 2.5))$theta,
            error = function(e) rep(NA, 6))
# trimmed-RMSE of every single method on a train/test split
evpt  <- function(b, r, yy, tri, tei) sapply(METH, function(m)
            trmse(yy[tei] - kg2_predict(fitm(m, b[tri], r[tri], yy[tri]), b[tei], r[tei])))

# auto-valgate: the discrete super-learner already in the package
fit_valgate <- function(b, r, yy, tri)
  tryCatch(inner_estimate(b[tri], r[tri], yy[tri],
    gmdh_pmm_control(B = 0, force_method = "auto-valgate", max_iter = 60L,
                     weak_sigma_mult = 2.5,
                     valgate_candidates = METH))$theta,
    error = function(e) rep(NA, 6))

# ---- Convex-stacking solver: weights on the simplex minimizing held-out
#      trimmed loss. Projected-gradient descent with Euclidean simplex
#      projection (Duchi et al. 2008). Self-contained; no external solver.
.proj_simplex <- function(v) {                     # project v onto {w>=0, sum w = 1}
  u <- sort(v, decreasing = TRUE); cs <- cumsum(u)
  rho <- max(which(u + (1 - cs) / seq_along(u) > 0))
  tau <- (cs[rho] - 1) / rho
  pmax(v - tau, 0)
}
# trimmed-squared-error of a blended prediction (matches trmse^2; smooth enough
# for PGD because the trim mask is recomputed each step from current residuals)
.tloss <- function(w, P, y, p = 0.90) {
  e <- y - as.numeric(P %*% w); ae <- abs(e)
  q <- stats::quantile(ae, p, names = FALSE, na.rm = TRUE)
  m <- ae <= q; mean(e[m]^2)
}
nnls_simplex <- function(P, y, steps = 400L, lr = NULL, p = 0.90) {
  # P: (n_val x K) candidate predictions; returns convex weights w (sum 1, >=0).
  K <- ncol(P); ok <- apply(P, 2, function(c) all(is.finite(c)))
  if (sum(ok) == 0L) return(rep(1 / K, K))
  Pf <- P[, ok, drop = FALSE]; Kf <- ncol(Pf)
  w <- rep(1 / Kf, Kf)
  if (is.null(lr)) lr <- 1 / (max(colSums(Pf^2)) / nrow(Pf) + 1e-8)  # ~1/Lipschitz scale
  best_w <- w; best_l <- .tloss(w, Pf, y, p)
  for (it in seq_len(steps)) {
    e <- y - as.numeric(Pf %*% w); ae <- abs(e)
    q <- stats::quantile(ae, p, names = FALSE, na.rm = TRUE); mask <- ae <= q
    g <- -2 * as.numeric(t(Pf[mask, , drop = FALSE]) %*% e[mask]) / max(1, sum(mask))
    w <- .proj_simplex(w - lr * g)
    l <- .tloss(w, Pf, y, p)
    if (is.finite(l) && l < best_l) { best_l <- l; best_w <- w }
  }
  full <- rep(0, K); full[ok] <- best_w; full
}

# Fit a convex stack on the TRAINING fold only, predict the TEST block.
#  1. inner sub-split of the training indices (70/30, contiguous in their own
#     order to avoid leakage between the two training halves);
#  2. fit all candidates on sub-train, predict sub-val -> P (n_subval x K);
#  3. solve simplex weights minimizing held-out trimmed loss;
#  4. refit all candidates on the FULL training fold, blend test predictions.
fit_stack_predict <- function(b, r, yy, tri, tei) {
  ntr <- length(tri)
  if (ntr < 40L) {                      # too small to sub-split: fall back to plain mean blend
    Pte <- sapply(METH, function(m) kg2_predict(fitm(m, b[tri], r[tri], yy[tri]), b[tei], r[tei]))
    okc <- apply(Pte, 2, function(c) all(is.finite(c)))
    if (!any(okc)) return(rep(NA, length(tei)))
    return(rowMeans(Pte[, okc, drop = FALSE]))
  }
  n_si <- floor(0.7 * ntr)
  si  <- tri[seq_len(n_si)]              # sub-train (head of the contiguous training block)
  sv  <- tri[(n_si + 1L):ntr]            # sub-val   (tail)
  Pv <- sapply(METH, function(m) kg2_predict(fitm(m, b[si], r[si], yy[si]), b[sv], r[sv]))
  w  <- nnls_simplex(Pv, yy[sv])
  # refit on FULL training fold, blend on test
  Pte <- sapply(METH, function(m) kg2_predict(fitm(m, b[tri], r[tri], yy[tri]), b[tei], r[tei]))
  okc <- apply(Pte, 2, function(c) all(is.finite(c)))
  if (!any(okc)) return(rep(NA, length(tei)))
  wf <- w; wf[!okc] <- 0; s <- sum(wf); if (s <= 0) wf[okc] <- 1 / sum(okc) else wf <- wf / s
  as.numeric(Pte %*% wf)
}

# ---- Blocked 5-fold evaluation of all four learner families for one cell ----
gap_default <- 50L; folds <- 5L
eval_cell <- function(b, r, y, gap) {
  n <- length(y); fb <- floor(n / folds)
  per_fold <- t(sapply(1:folds, function(k) {
    a <- (k - 1) * fb + 1; z <- if (k == folds) n else k * fb
    tei <- a:z; tri <- setdiff(1:n, max(1, a - gap):min(n, z + gap))
    if (length(tri) < 100 || length(tei) < 40) return(rep(NA, length(METH) + 2L))
    sm <- evpt(b, r, y, tri, tei)                                   # 5 single methods
    vg <- trmse(y[tei] - kg2_predict(fit_valgate(b, r, y, tri), b[tei], r[tei]))
    st <- trmse(y[tei] - fit_stack_predict(b, r, y, tri, tei))
    c(sm, valgate = vg, stack = st)
  }))
  apply(per_fold, 2, median, na.rm = TRUE)
}

emit_row <- function(tag, g3, g4, n, ac1, bl) {
  bw   <- min(bl["WPMM2"], bl["WPMM3"]); bwn <- names(which.min(bl[c("WPMM2","WPMM3")]))
  base <- min(bl["LSE"], bl["Huber"], bl["L1"]); basen <- names(which.min(bl[c("LSE","Huber","L1")]))
  oracle <- min(bl[METH])
  vg <- bl["valgate"]; st <- bl["stack"]
  data.frame(
    cell = tag, g3 = g3, g4 = g4, n = n, ac1 = ac1,
    bl_LSE = bl["LSE"], bl_Huber = bl["Huber"], bl_L1 = bl["L1"],
    bl_WPMM2 = bl["WPMM2"], bl_WPMM3 = bl["WPMM3"],
    oracle_single = oracle,
    best_weak = bw, best_weak_name = bwn,
    best_baseline = base, best_baseline_name = basen,
    valgate = vg, convex_stack = st,
    # head-to-head % (negative => row method is BETTER / lower trmse)
    stack_vs_valgate = 100 * (st / vg - 1),
    stack_vs_bestweak = 100 * (st / bw - 1),
    stack_vs_basebest = 100 * (st / base - 1),
    valgate_vs_bestweak = 100 * (vg / bw - 1),
    valgate_vs_basebest = 100 * (vg / base - 1),
    bestweak_vs_basebest = 100 * (bw / base - 1),
    row.names = NULL, stringsAsFactors = FALSE)
}

rows <- list()

# =================== VIBRATION: 12 FORGE 78B-32 cells ===================
cap <- 6000L; TMP <- "/tmp/forge78b_sensor"
.find1 <- function(pat) { x <- list.files(TMP, pattern = pat, recursive = TRUE, full.names = TRUE); if (length(x)) x[1] else NA }
RUNS <- list(`10.625` = .find1("Forge 10\\.625 - output data\\.csv$"),
             `14.75-Run1` = .find1("Run 1 - output data\\.csv$"),
             `14.75-Run2` = .find1("Run 2 - output data\\.csv$"))
SENSORS <- c("CSS-008", "CSS-007"); CHANS <- c("ShZpeak", "ShYpeak", "ShZrms")

vib_cell <- function(d, sensor, chan) {
  col <- function(s) suppressWarnings(as.numeric(d[[paste0(sensor, "_", s)]]))
  ss <- col("StickSlip(%)"); gs <- col("GyroXspread(RPM)"); tg <- col(paste0(chan, "(g)"))
  if (is.null(ss) || is.null(gs) || is.null(tg)) return(NULL)
  ok <- is.finite(ss) & is.finite(gs) & is.finite(tg) & tg > 0       # preserves time order
  b <- ss[ok]; r <- gs[ok]; y <- log(tg[ok]); n <- length(y)
  if (n < 600) return(NULL)
  if (n > cap) { b <- b[1:cap]; r <- r[1:cap]; y <- y[1:cap]; n <- cap }  # contiguous head
  ee <- resid(lm(y ~ b + r)); cu <- sample_cumulants(ee)
  ac1 <- stats::acf(ee, lag.max = 1, plot = FALSE)$acf[2]
  bl <- eval_cell(b, r, y, gap = 50L)
  list(n = n, g3 = cu$gamma3, g4 = cu$gamma4, ac1 = ac1, bl = bl)
}

cat("=== EXP2 stacking | blocked-5fold | best-weak vs auto-valgate vs convex-stack vs best baseline ===\n")
cat("--- VIBRATION (FORGE 78B-32) ---\n")
for (rn0 in names(RUNS)) { f <- RUNS[[rn0]]; if (is.na(f)) { cat(sprintf("  (run %s missing)\n", rn0)); next }
  d <- data.table::fread(f, skip = 35L, header = TRUE, fill = TRUE, showProgress = FALSE)
  for (sn in SENSORS) for (ch in CHANS) {
    res <- tryCatch(vib_cell(d, sn, ch), error = function(e) NULL); if (is.null(res)) next
    tag <- sprintf("vib:%s/%s/%s", rn0, sn, ch)
    row <- emit_row(tag, res$g3, res$g4, res$n, res$ac1, res$bl)
    rows[[length(rows) + 1L]] <- row
    cat(sprintf("%-30s g4=%5.1f | weak=%.3f base=%.3f vgate=%.3f stack=%.3f | stk-vs-vgate %+5.1f%% stk-vs-weak %+5.1f%%\n",
        tag, res$g4, row$best_weak, row$best_baseline, row$valgate, row$convex_stack,
        row$stack_vs_valgate, row$stack_vs_bestweak))
  }
}

# =================== Volve ROP: 4 depth bands ===================
cat("--- VOLVE ROP (KG-2 swob,tqa) ---\n")
dv <- tryCatch(load_volve_drilling(file.path(PKG, "data/volve_onbottom.csv")), error = function(e) NULL)
if (!is.null(dv)) {
  X <- as.data.frame(dv$X); y <- as.numeric(dv$y); ord <- order(X$dept); X <- X[ord, ]; y <- y[ord]
  BANDS <- list(volve_2000_2080 = c(2000, 2080), volve_1740_1820 = c(1740, 1820),
                volve_2110_2200 = c(2110, 2200), volve_mild_1360_1480 = c(1360, 1480))
  for (bn in names(BANDS)) { rg <- BANDS[[bn]]; sel <- X$dept >= rg[1] & X$dept < rg[2]
    Xi <- X[sel, ]; yi <- y[sel]; if (length(yi) < 150) next
    b <- Xi$swob; r <- Xi$tqa; yy <- yi
    ee <- resid(lm(yi ~ Xi$swob + Xi$tqa)); cu <- sample_cumulants(ee)
    ac1 <- stats::acf(ee, lag.max = 1, plot = FALSE)$acf[2]
    bl <- eval_cell(b, r, yy, gap = 30L)
    tag <- paste0("rop:", bn)
    row <- emit_row(tag, cu$gamma3, cu$gamma4, length(yy), ac1, bl)
    rows[[length(rows) + 1L]] <- row
    cat(sprintf("%-30s g4=%5.1f | weak=%.4f base=%.4f vgate=%.4f stack=%.4f | stk-vs-vgate %+5.1f%% stk-vs-weak %+5.1f%%\n",
        tag, cu$gamma4, row$best_weak, row$best_baseline, row$valgate, row$convex_stack,
        row$stack_vs_valgate, row$stack_vs_bestweak))
  }
} else cat("(Volve load failed)\n")

# =================== summary ===================
DF <- do.call(rbind, rows)
win_thresh <- 2  # %; |delta| < 2% == tie (matches the survival threshold in the source scripts)
n_cells <- nrow(DF)
cat(sprintf("\n=== SUMMARY over %d cells/bands (negative %% => first method better; |.|<%g%% == tie) ===\n", n_cells, win_thresh))
classify <- function(v) { if (v < -win_thresh) "win" else if (v > win_thresh) "lose" else "tie" }
sv <- table(factor(sapply(DF$stack_vs_valgate, classify), levels = c("win","tie","lose")))
cat(sprintf("convex-stack vs auto-valgate : stack better in %d, tie in %d, valgate better in %d\n", sv["win"], sv["tie"], sv["lose"]))
sw <- table(factor(sapply(DF$stack_vs_bestweak, classify), levels = c("win","tie","lose")))
cat(sprintf("convex-stack vs best-weak    : stack better in %d, tie in %d, best-weak better in %d\n", sw["win"], sw["tie"], sw["lose"]))
vw <- table(factor(sapply(DF$valgate_vs_bestweak, classify), levels = c("win","tie","lose")))
cat(sprintf("auto-valgate vs best-weak    : valgate better in %d, tie in %d, best-weak better in %d\n", vw["win"], vw["tie"], vw["lose"]))
sb <- table(factor(sapply(DF$stack_vs_basebest, classify), levels = c("win","tie","lose")))
cat(sprintf("convex-stack vs best baseline: stack better in %d, tie in %d, baseline better in %d\n", sb["win"], sb["tie"], sb["lose"]))
vb <- table(factor(sapply(DF$valgate_vs_basebest, classify), levels = c("win","tie","lose")))
cat(sprintf("auto-valgate vs best baseline: valgate better in %d, tie in %d, baseline better in %d\n", vb["win"], vb["tie"], vb["lose"]))
cat(sprintf("median |stack - valgate| = %.2f%% ; mean stack_vs_valgate = %+.2f%%\n",
            median(abs(DF$stack_vs_valgate), na.rm = TRUE), mean(DF$stack_vs_valgate, na.rm = TRUE)))

outdir <- file.path(dirname(dirname(PKG)), "paper-4-weak-moment-gmdh", "results")  # canonical, cwd-independent
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(DF, file.path(outdir, "stacking_blocked.csv"), row.names = FALSE)
cat(sprintf("\nSaved stacking_blocked.csv (%d rows) to %s\n", nrow(DF), outdir))
