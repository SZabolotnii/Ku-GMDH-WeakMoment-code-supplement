#!/usr/bin/env Rscript
# Mechanism ablation for the weak-PMM2 PoC (Paper 4, Weak-Moment GMDH).
#
# QUESTION. The PoC (weak-pmm-poc-2026-05-27.md) attributes WPMM2's win over
# Huber/LAD to the weak-cumulant a2 term "exploiting asymmetry". But the PoC
# window is ALREADY re-centred on median(r). The Ku_Weak_Moment (WM-1) program
# found that (i) degree-1 weak-PMM = a redescending M-estimator (no polynomial
# content) and (ii) the degree-2 even/skewness term gives ~zero bias correction,
# the real asymmetry fix being window RE-CENTRING. So the PoC win may come from
# re-centred windowing, NOT from the a2 polynomial term.
#
# DESIGN. Identical PoC DGP. A 2x2 ablation of the weak estimator that differs
# from fit_kg2_weak_pmm2() ONLY in two switches:
#   window centre  mu in {0, median(r)}      (recenter off / on)
#   polynomial     a2 in {0, weak-cumulant}  (use_a2 off / on)
# giving W00 (symmetric Welsch), W10 (re-centred Welsch = degree-1 + recenter),
# W01 (symmetric + polynomial), W11 (= WPMM2, the PoC estimator).
# Reference arms LSE/Huber/LAD/PMM2 via the same inner_estimate() as the PoC,
# plus the package WPMM2 to validate that W11 reproduces it.
#
# DECISIVE CONTRASTS.
#   W11 vs W10  -> does the a2 polynomial add anything once re-centred?
#   W10 vs Huber/LAD -> does re-centred degree-1 ALONE already beat them?
# If W11 ~ W10 and W10 already beats Huber/LAD, the "degree>=2 novelty" is empty.
#
# Usage (from paper-4-weak-moment-gmdh/code): Rscript experiments/run_weak_ablation.R [N_SEEDS]

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
N_TRAIN <- 200L; N_TEST <- 4000L; SEED0 <- 75001L
THETA <- c(0.5, 1.0, -0.8, 0.4, -0.3, 0.2)

rerr <- function(n, contam) {           # PoC error: (1-c)N(0,1) + c*(-|t_2|*2.5)
  heavy <- runif(n) < contam
  e <- stats::rnorm(n); k <- sum(heavy)
  if (k > 0) e[heavy] <- -abs(stats::rt(k, df = 2)) * 2.5
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

wridge <- function(Z, y, w, lambda = 1e-8) {   # weighted ridge (mirror of .solve_weighted_ridge)
  w[!is.finite(w) | w < 0] <- 0
  A <- crossprod(Z, Z * w); diag(A) <- diag(A) + lambda
  b <- crossprod(Z, w * y)
  tryCatch(as.numeric(solve(A, b)), error = function(e) rep(NA_real_, ncol(Z)))
}
# Faithful re-implementation of fit_kg2_weak_pmm2 with two switches.
fit_weak <- function(Z, y, sigma_mult = 2.5, recenter = TRUE, use_a2 = TRUE,
                     max_iter = 60L, tol = 1e-7) {
  theta <- wridge(Z, y, rep(1, length(y)))      # LSE warm start
  for (iter in seq_len(max_iter)) {
    r <- as.numeric(y - Z %*% theta)
    mu <- if (recenter) stats::median(r) else 0
    s <- stats::mad(r, constant = 1.4826)
    if (!is.finite(s) || s <= 1e-12) s <- stats::sd(r)
    if (!is.finite(s) || s <= 1e-12) s <- 1
    w <- exp(-0.5 * ((r - mu) / (sigma_mult * s))^2)
    sws <- sum(w); if (!is.finite(sws) || sws <= 0) break
    rc <- r - sum(w * r) / sws
    if (use_a2) {
      c2 <- sum(w * rc^2)/sws; c3 <- sum(w * rc^3)/sws; c4 <- sum(w * rc^4)/sws - 3*c2^2
      denom <- 2*c2^2 + c4
      a2 <- if (!is.finite(c2) || c2 <= 1e-12 || !is.finite(denom) || abs(denom) < 1e-12) 0 else -c3/denom
      a2 <- max(min(a2, 0.5/sqrt(c2 + 1e-12)), -0.5/sqrt(c2 + 1e-12))
    } else a2 <- 0
    W <- w * (1 + a2 * rc); W[!is.finite(W) | W < 0] <- 0
    tn <- wridge(Z, y, W)
    if (any(!is.finite(tn))) break
    if (sqrt(sum((tn - theta)^2))/max(1, sqrt(sum(theta^2))) < tol) { theta <- tn; break }
    theta <- tn
  }
  theta
}
fit_ref <- function(method, d_train) {          # package reference arms (as in PoC)
  ctrl <- gmdh_pmm_control(B = 0, force_method = method, max_iter = 60L, weak_sigma_mult = 2.5)
  inner_estimate(d_train$b1, d_train$b2, d_train$y, ctrl)$theta
}

ARMS_WEAK <- list(W00 = c(recenter=FALSE, use_a2=FALSE), W10 = c(recenter=TRUE,  use_a2=FALSE),
                  W01 = c(recenter=FALSE, use_a2=TRUE ), W11 = c(recenter=TRUE,  use_a2=TRUE ))
ARMS_REF  <- c("LSE","Huber","L1","PMM2","WPMM2")

run_one <- function(seed, contam) {
  set.seed(seed)
  g <- gen(N_TRAIN); y_tr <- g$mu + rerr(N_TRAIN, contam)
  d_train <- kg2_design(g$v1, g$v2, y_tr)
  Ztr <- as.matrix(cbind(1, d_train$b1, d_train$b2, d_train$b12, d_train$b11, d_train$b22))
  gt <- gen(N_TEST, seed + 10000L)
  rows <- list()
  emit <- function(nm, th) {
    pred <- kg2_predict(th, gt$v1, gt$v2); e <- pred - gt$mu
    data.frame(seed=seed, contam=contam, method=nm, mae=mean(abs(e)), trmse=trmse(e),
               coef_l2=sqrt(sum((th-THETA)^2)), b0_abs=abs(th[1]-THETA[1]), stringsAsFactors=FALSE)
  }
  for (nm in names(ARMS_WEAK)) {
    sw <- ARMS_WEAK[[nm]]
    th <- tryCatch(fit_weak(Ztr, d_train$y, recenter=as.logical(sw["recenter"]), use_a2=as.logical(sw["use_a2"])),
                   error=function(e) rep(NA,6))
    rows[[length(rows)+1L]] <- emit(nm, th)
  }
  for (m in ARMS_REF) {
    th <- tryCatch(fit_ref(m, d_train), error=function(e) rep(NA,6))
    rows[[length(rows)+1L]] <- emit(m, th)
  }
  do.call(rbind, rows)
}

cat(sprintf("=== Weak-PMM2 mechanism ablation | %d seeds | n_train=%d ===\n", N_SEEDS, N_TRAIN))
t0 <- Sys.time()
all <- do.call(rbind, lapply(c(0.12, 0.0), function(cc)
  do.call(rbind, lapply(seq_len(N_SEEDS), function(s) run_one(SEED0 + s - 1L, cc)))))
outdir <- "../results"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(all, file.path(outdir, "weak_ablation_raw.csv"), row.names = FALSE)

for (cc in c(0.12, 0.0)) {
  sub <- all[all$contam == cc, ]
  cat(sprintf("\n--- contam = %.0f%% ---\n", 100*cc))
  agg <- do.call(rbind, by(sub, sub$method, function(gp) data.frame(method=gp$method[1],
    mae_med=median(gp$mae), mae_iqr=IQR(gp$mae), trmse_med=median(gp$trmse),
    coef_l2_med=median(gp$coef_l2), b0_abs_med=median(gp$b0_abs), stringsAsFactors=FALSE)))
  print(agg[order(agg$mae_med), ], row.names = FALSE, digits = 4)
  piv <- reshape(sub[,c("seed","method","mae")], idvar="seed", timevar="method", direction="wide")
  contrasts <- list(c("W11","W10"), c("W11","W01"), c("W10","W00"), c("W01","W00"),
                    c("W11","Huber"), c("W11","L1"), c("W10","Huber"), c("W10","L1"),
                    c("W11","WPMM2"))
  cat("paired Wilcoxon (col1 lower MAE = win):\n")
  for (ct in contrasts) {
    a <- piv[[paste0("mae.",ct[1])]]; b <- piv[[paste0("mae.",ct[2])]]
    ok <- is.finite(a) & is.finite(b)
    if (sum(ok) < 5) next
    w <- suppressWarnings(stats::wilcox.test(a[ok], b[ok], paired = TRUE))
    cat(sprintf("  %-5s vs %-5s: win %3.0f%%, median dMAE=%+.4f, p=%.2e\n",
        ct[1], ct[2], 100*mean(b[ok] > a[ok]), median(b[ok]-a[ok]), w$p.value))
  }
}
cat(sprintf("\nSaved weak_ablation_raw.csv | %.0fs\n", as.numeric(difftime(Sys.time(), t0, "secs"))))
