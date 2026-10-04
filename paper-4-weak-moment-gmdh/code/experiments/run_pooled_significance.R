#!/usr/bin/env Rscript
# POOLED cross-cell significance (route B, final proofing piece). Block-bootstrap showed
# 0/7 cells individually significant, but the direction was favourable 6/7. This pools
# across ALL pre-specified cells (12 vibration + 4 Volve ROP; NO outcome-based selection)
# to test whether the program-level advantage holds. Three triangulating tests per
# comparison: (a) cross-cell SIGN test (each cell = 1 independent unit), (b) cross-cell
# Wilcoxon signed-rank on per-cell relative wins, (c) POOLED stratified moving-block
# bootstrap (resamples each cell's paired test-error series respecting autocorrelation,
# averages the relative win across cells -> 95% CI on the pooled mean). Reports BOTH an
# oracle comparison (best-weak vs best-baseline / vs LSE) AND honest FIXED/deployable ones
# (WPMM3 vs Huber; deployable auto-valgate vs Huber). Heavy-tail subgroup pre-specified by
# gamma4>4 (a data property). Same blocked CV as run_honest_blocked_cv*.R. Drilling only.
# Usage: Rscript experiments/run_pooled_significance.R

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))

trmse <- function(e, p = 0.90) { q <- stats::quantile(abs(e), p, names = FALSE, na.rm = TRUE); sqrt(mean(e[abs(e) <= q]^2)) }
FITM <- c("LSE", "Huber", "L1", "WPMM2", "WPMM3", "auto-valgate")
fitm <- function(m, b, r, yy) tryCatch(inner_estimate(b, r, yy,
          gmdh_pmm_control(B = 0, force_method = m, max_iter = 60L, weak_sigma_mult = 2.5))$theta,
          error = function(e) rep(NA, 6))
gap_v <- 50L; gap_r <- 30L; folds <- 5L; cap <- 6000L; L <- 75L; B <- 2000L

# ---- collect per-cell, time-ordered, paired test errors per method (blocked folds tile the series) ----
collect <- function(b, r, y, gap) {
  n <- length(y); fb <- floor(n / folds)
  E <- matrix(NA_real_, n, length(FITM), dimnames = list(NULL, FITM)); fid <- integer(n)
  for (k in 1:folds) { a <- (k - 1) * fb + 1; z <- if (k == folds) n else k * fb
    tei <- a:z; tri <- setdiff(1:n, max(1, a - gap):min(n, z + gap))
    if (length(tri) < 100 || length(tei) < 40) next
    fid[tei] <- k
    for (m in FITM) E[tei, m] <- y[tei] - kg2_predict(fitm(m, b[tri], r[tri], y[tri]), b[tei], r[tei]) }
  keep <- fid > 0 & is.finite(rowSums(E[, c("LSE","Huber","L1","WPMM2","WPMM3")]))
  list(E = E[keep, , drop = FALSE], fid = fid[keep])
}
# moving-block resample of row indices WITHIN each fold (blocks never cross fold boundaries)
mbb_idx <- function(fid) {
  out <- integer(0)
  for (k in unique(fid)) { pos <- which(fid == k); m <- length(pos); if (m < 5) { out <- c(out, pos); next }
    nb <- ceiling(m / L); starts <- sample.int(max(1, m - L + 1), nb, replace = TRUE)
    idx <- unlist(lapply(starts, function(s) pos[s:min(m, s + L - 1)]))[1:m]; out <- c(out, idx) }
  out
}
stat_set <- function(E) {  # named vector of comparison ratios from a (resampled) error matrix
  tv <- function(m) trmse(E[, m]); bw <- min(tv("WPMM2"), tv("WPMM3")); bb <- min(tv("LSE"), tv("Huber"), tv("L1"))
  c(weak_vs_base = bw / bb - 1, weak_vs_LSE = bw / tv("LSE") - 1,
    WPMM3_vs_Huber = tv("WPMM3") / tv("Huber") - 1, valgate_vs_Huber = tv("auto-valgate") / tv("Huber") - 1)
}

# ---- build the cell list (mirror honest scripts) ----
cells <- list()
TMP <- "/tmp/forge78b_sensor"
.find1 <- function(pat) { x <- list.files(TMP, pattern = pat, recursive = TRUE, full.names = TRUE); if (length(x)) x[1] else NA }
RUNS <- list(`14.75-Run1` = .find1("Run 1 - output data\\.csv$"), `14.75-Run2` = .find1("Run 2 - output data\\.csv$"))
for (rn in names(RUNS)) { f <- RUNS[[rn]]; if (is.na(f)) next
  d <- data.table::fread(f, skip = 35L, header = TRUE, fill = TRUE, showProgress = FALSE)
  for (sn in c("CSS-008","CSS-007")) for (ch in c("ShZpeak","ShYpeak","ShZrms")) {
    g <- function(s) suppressWarnings(as.numeric(d[[paste0(sn, "_", s)]]))
    ss <- g("StickSlip(%)"); gs <- g("GyroXspread(RPM)"); tg <- g(paste0(ch, "(g)"))
    if (is.null(ss) || is.null(gs) || is.null(tg)) next
    ok <- is.finite(ss) & is.finite(gs) & is.finite(tg) & tg > 0
    b <- ss[ok]; r <- gs[ok]; yy <- log(tg[ok]); nn <- length(yy); if (nn < 600) next
    if (nn > cap) { b <- b[1:cap]; r <- r[1:cap]; yy <- yy[1:cap] }
    cells[[length(cells)+1L]] <- list(tag = sprintf("vib:%s/%s/%s", rn, sn, ch),
      g4 = sample_cumulants(resid(lm(yy ~ b + r)))$gamma4, b = b, r = r, y = yy, gap = gap_v) } }
dv <- tryCatch(load_volve_drilling(file.path(PKG, "data/volve_onbottom.csv")), error = function(e) NULL)
if (!is.null(dv)) { X <- as.data.frame(dv$X); y <- as.numeric(dv$y); o <- order(X$dept); X <- X[o, ]; y <- y[o]
  BANDS <- list(volve_2000_2080=c(2000,2080), volve_1740_1820=c(1740,1820), volve_2110_2200=c(2110,2200), volve_mild_1360_1480=c(1360,1480))
  for (bn in names(BANDS)) { rg <- BANDS[[bn]]; sel <- X$dept >= rg[1] & X$dept < rg[2]
    if (sum(sel) < 150) next; Xi <- X[sel, ]; yi <- y[sel]
    cells[[length(cells)+1L]] <- list(tag = paste0("rop:", bn),
      g4 = sample_cumulants(resid(lm(yi ~ Xi$swob + Xi$tqa)))$gamma4, b = Xi$swob, r = Xi$tqa, y = yi, gap = gap_r) } }

# ---- per-cell point stats + per-cell bootstrap distributions ----
COMPS <- c("weak_vs_base","weak_vs_LSE","WPMM3_vs_Huber","valgate_vs_Huber")
pt <- matrix(NA, length(cells), length(COMPS), dimnames = list(sapply(cells, `[[`, "tag"), COMPS))
g4v <- sapply(cells, `[[`, "g4")
boot_cell <- vector("list", length(cells))   # each: B x COMPS
for (i in seq_along(cells)) { ci <- cells[[i]]
  co <- collect(ci$b, ci$r, ci$y, ci$gap); E <- co$E; fid <- co$fid
  pt[i, ] <- stat_set(E)
  bm <- matrix(NA, B, length(COMPS), dimnames = list(NULL, COMPS))
  for (bb in 1:B) { id <- mbb_idx(fid); bm[bb, ] <- stat_set(E[id, , drop = FALSE]) }
  boot_cell[[i]] <- bm
  cat(sprintf("  %-30s g4=%6.2f | weak/base %+6.1f%%  weak/LSE %+6.1f%%  W3/Hub %+6.1f%%  valgate/Hub %+6.1f%%\n",
      ci$tag, ci$g4, 100*pt[i,1], 100*pt[i,2], 100*pt[i,3], 100*pt[i,4])) }

# ---- pooled inference ----
groups <- list(all = rep(TRUE, length(cells)), heavy_g4gt4 = g4v > 4)
report <- function(rows, sel, comp) {
  v <- pt[sel, comp]; n <- length(v); k <- sum(v < 0)                        # sign: weak better = negative
  sp <- stats::binom.test(k, n, 0.5)$p.value
  wp <- tryCatch(stats::wilcox.test(v, mu = 0)$p.value, error = function(e) NA)
  # pooled MBB: mean across selected cells of the comparison per bootstrap rep
  M <- sapply(which(sel), function(i) boot_cell[[i]][, comp]); pooled <- rowMeans(M)
  ci <- stats::quantile(pooled, c(.025, .975), names = FALSE)
  data.frame(comparison = comp, group = NA, n_cells = n, sign_k = k,
    mean_win_pct = 100*mean(v), sign_p = sp, wilcox_p = wp,
    pooled_mean_pct = 100*mean(pooled), pooled_ci_lo = 100*ci[1], pooled_ci_hi = 100*ci[2],
    pooled_excl0 = (ci[1] < 0 && ci[2] < 0) || (ci[1] > 0 && ci[2] > 0), stringsAsFactors = FALSE)
}
cat("\n=== POOLED cross-cell significance (sign test, Wilcoxon signed-rank, pooled block-bootstrap) ===\n")
out <- list()
for (gn in names(groups)) { sel <- groups[[gn]]
  cat(sprintf("\n[group %s | %d cells]\n", gn, sum(sel)))
  for (comp in COMPS) { row <- report(rows, sel, comp); row$group <- gn; out[[length(out)+1L]] <- row
    cat(sprintf("  %-16s mean=%+6.1f%% | sign %d/%d p=%.4f | wilcox p=%.4f | pooled %+6.1f%% CI[%+.1f,%+.1f] excl0=%s\n",
        comp, row$mean_win_pct, row$sign_k, row$n_cells, row$sign_p, row$wilcox_p,
        row$pooled_mean_pct, row$pooled_ci_lo, row$pooled_ci_hi, row$pooled_excl0)) } }
outdir <- file.path(dirname(dirname(PKG)), "paper-4-weak-moment-gmdh", "results")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(do.call(rbind, out), file.path(outdir, "pooled_significance.csv"), row.names = FALSE)
cat("\nSaved pooled_significance.csv\n")
