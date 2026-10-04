#!/usr/bin/env Rscript
# HONEST BLOCKED-CV for ROP soft-sensing (Volve + FORGE), to test under
# leakage-controlled folds whether the program's "LSE is uniquely fragile on ROP,
# weak/robust beat it 22-70%" claim survives -- the ROP comparisons were all
# random 70/30 on depth-ordered (autocorrelated) series. Blocked folds = contiguous
# DEPTH blocks (+ purge gap). Same model as run_drilling_realdata.R: KG-2(swob,tqa).
# Reports best-weak AND best-robust vs LSE under blocking. Strictly drilling.
# Usage: Rscript experiments/run_honest_blocked_cv_rop.R

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))

trmse <- function(e, p = 0.90) { q <- stats::quantile(abs(e), p, names = FALSE, na.rm = TRUE); sqrt(mean(e[abs(e) <= q]^2)) }
METH  <- c("LSE", "Huber", "L1", "WPMM2", "WPMM3")
fitp  <- function(m, tr, te) { est <- tryCatch(inner_estimate(tr$swob, tr$tqa, tr$y,
            gmdh_pmm_control(B = 0, force_method = m, max_iter = 60L, weak_sigma_mult = 2.5)), error = function(e) NULL)
          if (is.null(est)) rep(NA, nrow(te)) else kg2_predict(est$theta, te$swob, te$tqa) }
evpt  <- function(band, tri, tei) sapply(METH, function(m) trmse(band$y[tei] - fitp(m, band[tri, ], band[tei, ])))
gap <- 30L; folds <- 5L

run_band <- function(band, tag, g3, g4) {
  n <- nrow(band)
  fb <- floor(n / folds)
  mb <- t(sapply(1:folds, function(k) { a <- (k - 1) * fb + 1; z <- if (k == folds) n else k * fb
    tei <- a:z; tri <- setdiff(1:n, max(1, a - gap):min(n, z + gap))
    if (length(tri) < 80 || length(tei) < 30) return(rep(NA, length(METH))); evpt(band, tri, tei) }))
  bl <- apply(mb, 2, median, na.rm = TRUE)
  mr <- t(sapply(1:25, function(s) { set.seed(40000 + s); tri <- sample.int(n, floor(.7 * n)); evpt(band, tri, setdiff(1:n, tri)) }))
  rnd <- apply(mr, 2, median, na.rm = TRUE)
  bw <- min(bl["WPMM2"], bl["WPMM3"]); rob <- min(bl["Huber"], bl["L1"])
  cat(sprintf("%-28s g3=%+.2f g4=%5.1f | B: LSE=%.4f Hub=%.4f L1=%.4f W2=%.4f W3=%.4f | weak vs LSE %+5.1f%% | robust vs LSE %+5.1f%% | best=%s\n",
      tag, g3, g4, bl["LSE"], bl["Huber"], bl["L1"], bl["WPMM2"], bl["WPMM3"],
      100*(bw/bl["LSE"]-1), 100*(rob/bl["LSE"]-1), names(which.min(bl))))
  data.frame(dataset=tag, g3=g3, g4=g4, n=n,
    bl_LSE=bl["LSE"], bl_Huber=bl["Huber"], bl_L1=bl["L1"], bl_WPMM2=bl["WPMM2"], bl_WPMM3=bl["WPMM3"],
    rn_LSE=rnd["LSE"], rn_Huber=rnd["Huber"], rn_L1=rnd["L1"], rn_WPMM2=rnd["WPMM2"], rn_WPMM3=rnd["WPMM3"],
    blocked_weak_vs_lse=100*(bw/bl["LSE"]-1), blocked_robust_vs_lse=100*(rob/bl["LSE"]-1),
    blocked_best=names(which.min(bl)), row.names=NULL, stringsAsFactors=FALSE)
}

rows <- list()
cat("=== HONEST blocked-CV ROP (depth folds) | KG-2(swob,tqa) | does LSE-fragility survive? ===\n")

# --- Volve ---
dv <- tryCatch(load_volve_drilling(file.path(PKG, "data/volve_onbottom.csv")), error = function(e) NULL)
if (!is.null(dv)) {
  X <- as.data.frame(dv$X); y <- as.numeric(dv$y); ord <- order(X$dept); X <- X[ord, ]; y <- y[ord]
  BANDS <- list(volve_2000_2080=c(2000,2080), volve_1740_1820=c(1740,1820), volve_2110_2200=c(2110,2200), volve_mild_1360_1480=c(1360,1480))
  for (bn in names(BANDS)) { rg <- BANDS[[bn]]; sel <- X$dept >= rg[1] & X$dept < rg[2]
    Xi <- X[sel, ]; yi <- y[sel]; if (length(yi) < 150) next
    band <- data.frame(swob=Xi$swob, tqa=Xi$tqa, y=yi); cu <- sample_cumulants(resid(lm(yi~Xi$swob+Xi$tqa)))
    rows[[length(rows)+1L]] <- run_band(band, bn, cu$gamma3, cu$gamma4) }
} else cat("(Volve load failed)\n")

# --- FORGE (try 56-32; fully guarded — never halt the Volve result/save) ---
tryCatch({
  df <- load_forge_standard_drilling(well = "56-32", cadence = "10sec")
  X <- as.data.frame(df$X); y <- as.numeric(df$y)
  if (is.null(X$swob) || is.null(X$tqa)) stop("FORGE X lacks swob/tqa columns: ", paste(names(X), collapse=","))
  ord <- order(X$dept); X <- X[ord, ]; y <- y[ord]
  n <- length(y); win <- 600L
  starts <- seq(1, n - win, by = win)
  g4s <- sapply(starts, function(s) { idx <- s:(s+win-1); sample_cumulants(resid(lm(y[idx]~X$swob[idx]+X$tqa[idx])))$gamma4 })
  pick <- starts[order(-g4s)][1:min(3, length(starts))]
  for (s in pick) { idx <- s:(s+win-1); Xi <- X[idx, ]; yi <- y[idx]
    band <- data.frame(swob=Xi$swob, tqa=Xi$tqa, y=yi); cu <- sample_cumulants(resid(lm(yi~Xi$swob+Xi$tqa)))
    rows[[length(rows)+1L]] <- run_band(band, sprintf("forge56_d%.0f", Xi$dept[1]), cu$gamma3, cu$gamma4) }
}, error = function(e) cat(sprintf("(FORGE 56-32 ROP skipped: %s)\n", conditionMessage(e))))

outdir <- file.path(dirname(dirname(PKG)), "paper-4-weak-moment-gmdh", "results")  # canonical, cwd-independent
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(do.call(rbind, rows), file.path(outdir, "honest_blocked_cv_rop.csv"), row.names = FALSE)
cat("\nSaved honest_blocked_cv_rop.csv\n")
