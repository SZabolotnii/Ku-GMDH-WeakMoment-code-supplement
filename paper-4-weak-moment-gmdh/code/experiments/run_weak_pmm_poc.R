#!/usr/bin/env Rscript
# Weak-moment PMM proof-of-concept (bridge to Ku_Weak_Moment).
#
# Controlled KG-2 regression where the error is SKEWED *and* HEAVY-TAILED
# (asymmetric contamination) -- the regime that, on real drilling data, made
# classical PMM2 lose: raw cumulants gamma3/gamma4 are destabilized by the heavy
# component, so PMM2's correction is noisy, while Huber/LAD are robust but
# symmetric and cannot exploit the asymmetry. Hypothesis: weak-PMM2 (Gaussian
# window on a robust residual centre + weak-cumulant skew correction) tames the
# heavy tail AND keeps the skew correction, beating both classical PMM2 and the
# robust estimators.
#
# Metrics follow Ku_Weak_Moment MC discipline: MAE and 90%-trimmed RMSE on a
# clean held-out test set (raw RMSE is unstable under heavy tails), plus the
# coefficient L2 error to the known truth.
#
# Usage (run from paper-4-weak-moment-gmdh/code):
#   Rscript experiments/run_weak_pmm_poc.R [N_SEEDS] [N_TRAIN] [CONTAM]
#
# Reuses the shared gmdhpmm package (paper-1-gmdh-pmm/code); the weak-PMM2
# estimator is R/weak_pmm.R there, exposed via force_method = "WPMM2".

.find_pkg <- function() {
  cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm package not found")
    cur <- parent
  }
}
suppressMessages(
  if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
  else pkgload::load_all(.find_pkg(), quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
N_SEEDS <- if (length(args) >= 1) as.integer(args[1]) else 200L
N_TRAIN <- if (length(args) >= 2) as.integer(args[2]) else 200L
CONTAM  <- if (length(args) >= 3) as.numeric(args[3]) else 0.12
N_TEST  <- 5000L
SEED0 <- 75001L
METHODS <- c("WPMM2","PMM2","LSE","Huber","L1")
THETA <- c(0.5, 1.0, -0.8, 0.4, -0.3, 0.2)   # true KG-2 coefficients

# Skewed heavy-tailed error: (1-c) N(0,1) + c * (-|t_df2|*scale) -> left-skewed,
# heavy left tail. Standardized to unit-ish scale via the bulk component.
rerr <- function(n, contam) {
  heavy <- runif(n) < contam
  e <- stats::rnorm(n)
  k <- sum(heavy)
  if (k > 0) e[heavy] <- -abs(stats::rt(k, df = 2)) * 2.5    # one-sided heavy tail
  e
}

gen <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  v1 <- stats::rnorm(n); v2 <- stats::rnorm(n)
  d <- kg2_design(v1, v2)
  Z <- as.matrix(cbind(1, d$b1, d$b2, d$b12, d$b11, d$b22))
  list(v1 = v1, v2 = v2, Z = Z, mu = as.numeric(Z %*% THETA))
}

trmse <- function(e, p = 0.90) { q <- stats::quantile(abs(e), p, names = FALSE); sqrt(mean(e[abs(e) <= q]^2)) }

fit_one <- function(method, d_train) {
  ctrl <- gmdh_pmm_control(B = 0, force_method = method, max_iter = 60L,
                           weak_sigma_mult = 2.5)
  if (method == "WPMM2") ctrl$force_method <- "WPMM2"
  est <- inner_estimate(d_train$b1, d_train$b2, d_train$y, ctrl)
  est$theta
}

run_one <- function(seed) {
  set.seed(seed)
  g <- gen(N_TRAIN); y_tr <- g$mu + rerr(N_TRAIN, CONTAM)
  d_train <- kg2_design(g$v1, g$v2, y_tr)
  gt <- gen(N_TEST, seed + 10000L)            # clean-bulk test (same DGP, eval on truth)
  rows <- list()
  for (m in METHODS) {
    th <- tryCatch(fit_one(m, d_train), error = function(e) rep(NA, 6))
    pred <- kg2_predict(th, gt$v1, gt$v2)
    e <- pred - gt$mu                          # error vs noise-free truth
    rows[[length(rows)+1L]] <- data.frame(seed = seed, method = m,
      mae = mean(abs(e)), trmse = trmse(e), coef_l2 = sqrt(sum((th - THETA)^2)),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

cat(sprintf("=== Weak-PMM PoC | skew+heavy error, contam=%.0f%%, n_train=%d, %d seeds ===\n",
            100*CONTAM, N_TRAIN, N_SEEDS))
t0 <- Sys.time()
all <- do.call(rbind, lapply(seq_len(N_SEEDS), function(s) run_one(SEED0 + s - 1L)))
outdir <- "../results"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(all, file.path(outdir, "weak_pmm_poc_raw.csv"), row.names = FALSE)
agg <- do.call(rbind, by(all, all$method, function(gp) data.frame(method = gp$method[1],
  n = nrow(gp), mae_med = median(gp$mae), mae_iqr = IQR(gp$mae),
  trmse_med = median(gp$trmse), coef_l2_med = median(gp$coef_l2), stringsAsFactors = FALSE)))
agg <- agg[order(agg$mae_med), ]
print(agg, row.names = FALSE, digits = 4)

# paired Wilcoxon: WPMM2 vs each other on MAE
piv <- reshape(all[,c("seed","method","mae")], idvar="seed", timevar="method", direction="wide")
cat("\nPaired (MAE, WPMM2 lower = win):\n")
for (m in setdiff(METHODS, "WPMM2")) {
  a <- piv[["mae.WPMM2"]]; b <- piv[[paste0("mae.", m)]]
  ok <- is.finite(a) & is.finite(b)
  w <- suppressWarnings(stats::wilcox.test(a[ok], b[ok], paired = TRUE))
  cat(sprintf("  WPMM2 vs %-5s: better in %2.0f%% seeds, median dMAE=%+.4f, p=%.2e\n",
      m, 100*mean(b[ok] > a[ok]), median(b[ok]-a[ok]), w$p.value))
}
cat(sprintf("\nSaved weak_pmm_poc_raw.csv | %.0fs\n", as.numeric(difftime(Sys.time(), t0, "secs"))))
