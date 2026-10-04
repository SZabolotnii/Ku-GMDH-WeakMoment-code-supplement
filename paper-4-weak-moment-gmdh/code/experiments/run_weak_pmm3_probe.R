#!/usr/bin/env Rscript
# Weak-PMM3 (degree-3, kurtosis-coupled) probe for Paper 4.
#
# WHY. The a2 (degree-2) ablation/regime-map showed degree-2 wins only under
# ASYMMETRIC heavy contamination and is INERT under symmetric heavy tails. The
# symmetric-heavy regime is the PMM3 route. WM-1 (Ku_Weak_Moment) result for the
# windowed odd-cubic PMM3 score psi(r)=w(r)(kappa^w r - r^3), with
# kappa^w=(m6^w-3 m4^w m2^w)/(m4^w-3 (m2^w)^2): ~MLE-efficient at the optimal
# bandwidth, but its real value is BANDWIDTH INSURANCE -- it holds efficiency at
# wide windows where degree-1 (Welsch) drops. Needs a fast-decay (Schwartz)
# window so m6^w is finite (Gaussian OK even under t2). NB: weak-PMM3 here
# targets SYMMETRIC LEPTOKURTIC heavy tails (raw m6 unstable); it also covers the
# classical platykurtic PMM3 regime (uniform).
#
# This script implements weak-PMM3 (not yet in the package) as a faithful Newton
# solver on the windowed cubic score, validates it, and probes:
#   (1) head-to-head across symmetric regimes at sigma_mult=2.5;
#   (2) a bandwidth sweep (sym_t3) showing the insurance effect vs degree-1.
# Usage (from paper-4-weak-moment-gmdh/code): Rscript experiments/run_weak_pmm3_probe.R [N_SEEDS]

.find_pkg <- function() {
  cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent }
}
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
                 else pkgload::load_all(.find_pkg(), quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
N_SEEDS <- if (length(args) >= 1) as.integer(args[1]) else 150L
N_TRAIN <- 200L; N_TEST <- 4000L; SEED0 <- 75001L
THETA <- c(0.5, 1.0, -0.8, 0.4, -0.3, 0.2)

wridge <- function(Z, y, w, lambda = 1e-8) {
  w[!is.finite(w) | w < 0] <- 0
  A <- crossprod(Z, Z * w); diag(A) <- diag(A) + lambda; b <- crossprod(Z, w * y)
  tryCatch(as.numeric(solve(A, b)), error = function(e) rep(NA_real_, ncol(Z)))
}
# degree-1 weak (recentred Welsch) -- warm start + reference
fit_weak1 <- function(Z, y, sigma_mult = 2.5, max_iter = 60L, tol = 1e-7) {
  theta <- wridge(Z, y, rep(1, length(y)))
  for (i in seq_len(max_iter)) {
    r <- as.numeric(y - Z %*% theta); mu <- stats::median(r)
    s <- stats::mad(r, constant = 1.4826); if (!is.finite(s)||s<=1e-12) s <- stats::sd(r); if (!is.finite(s)||s<=1e-12) s <- 1
    w <- exp(-0.5*((r-mu)/(sigma_mult*s))^2)
    tn <- wridge(Z, y, w); if (any(!is.finite(tn))) break
    if (sqrt(sum((tn-theta)^2))/max(1,sqrt(sum(theta^2))) < tol) { theta <- tn; break }; theta <- tn
  }
  theta
}
# Windowed PMM3 score machinery (frozen kappa^w & weak-mean per evaluation).
.pmm3_eval <- function(theta, Z, y, sigma_mult, meso_tol = 1e-3) {
  r <- as.numeric(y - Z %*% theta); mu <- stats::median(r)
  s <- stats::mad(r, constant = 1.4826); if (!is.finite(s)||s<=1e-12) s <- stats::sd(r); if (!is.finite(s)||s<=1e-12) s <- 1
  w <- exp(-0.5*((r-mu)/(sigma_mult*s))^2); sws <- sum(w)
  if (!is.finite(sws) || sws <= 0) return(NULL)
  rc <- r - sum(w*r)/sws
  m2 <- sum(w*rc^2)/sws; m4 <- sum(w*rc^4)/sws; m6 <- sum(w*rc^6)/sws
  dk <- m4 - 3*m2^2
  meso <- !is.finite(dk) || abs(dk) < meso_tol*max(m2^2, 1e-12)
  if (meso) {                       # mesokurtic -> reduce to degree-1
    g <- crossprod(Z, w*rc); A <- crossprod(Z, Z*w); kappa <- NA_real_
  } else {
    kappa <- (m6 - 3*m4*m2)/dk
    g <- crossprod(Z, w*(kappa*rc - rc^3))            # score U3 = Z'[w(kappa rc - rc^3)]
    A <- crossprod(Z, Z*(w*(kappa - 3*rc^2)))         # -d score/dtheta  (frozen kappa)
  }
  list(g = as.numeric(g), A = A, norm = sqrt(sum(g^2)), kappa = kappa, meso = meso)
}
# Weak-PMM3 by factored IRLS -- the same machinery as the package's WPMM2
# (fit_kg2_weak_pmm2), with the cubic kurtosis weight W = w*(kappa^w - rc^2)/|kappa^w|
# in place of the skew weight w*(1+a2*rc). IRLS solves a sequence of weighted-LS
# problems; its fixed point is the (centred) windowed-PMM3 estimating equation
# sum_i w_i (kappa^w rc_i - rc_i^3) z_i = 0. Weights are redescending (clipped
# >=0: residuals beyond |rc|=sqrt(kappa) get zero weight). Mesokurtic
# (m4^w ~ 3 m2^w^2) -> kappa->inf -> reduces to degree-1 (W=w). Stable because a
# weighted-LS step always solves (no indefinite-Jacobian Newton failure).
fit_weak_pmm3 <- function(Z, y, sigma_mult = 2.5, max_iter = 60L, tol = 1e-7) {
  theta <- wridge(Z, y, rep(1, length(y)))      # LSE warm start
  meso_hits <- 0L; kappa <- NA_real_
  for (iter in seq_len(max_iter)) {
    r <- as.numeric(y - Z %*% theta); mu <- stats::median(r)
    s <- stats::mad(r, constant = 1.4826); if (!is.finite(s)||s<=1e-12) s <- stats::sd(r); if (!is.finite(s)||s<=1e-12) s <- 1
    w <- exp(-0.5*((r-mu)/(sigma_mult*s))^2); sws <- sum(w); if (!is.finite(sws)||sws<=0) break
    rc <- r - sum(w*r)/sws
    m2 <- sum(w*rc^2)/sws; m4 <- sum(w*rc^4)/sws; m6 <- sum(w*rc^6)/sws
    dk <- m4 - 3*m2^2
    if (!is.finite(dk) || abs(dk) < 1e-3*max(m2^2, 1e-12)) {     # mesokurtic -> degree-1
      meso_hits <- meso_hits + 1L; kappa <- NA_real_; W <- w
    } else {
      kappa <- (m6 - 3*m4*m2)/dk
      W <- w * (kappa - rc^2) / max(abs(kappa), 1e-8)            # normalized cubic weight in (-inf, 1]
    }
    W[!is.finite(W) | W < 0] <- 0                               # redescending
    if (sum(W) <= 0) W <- w                                     # safety
    tn <- wridge(Z, y, W); if (any(!is.finite(tn))) break
    if (sqrt(sum((tn - theta)^2))/max(1, sqrt(sum(theta^2))) < tol) { theta <- tn; break }
    theta <- tn
  }
  attr(theta, "meso_hits") <- meso_hits; attr(theta, "kappa") <- kappa
  theta
}
fit_ref <- function(method, d) inner_estimate(d$b1, d$b2, d$y,
  gmdh_pmm_control(B = 0, force_method = method, max_iter = 60L, weak_sigma_mult = 2.5))$theta

# Symmetric error regimes (+ platykurtic), unit-ish scale.
REG <- list(
  gauss        = function(n) stats::rnorm(n),                              # mesokurtic
  uniform_platy= function(n) stats::runif(n, -sqrt(3), sqrt(3)),           # platykurtic g4<0
  sym_t5       = function(n) stats::rt(n, 5)/sqrt(5/3),                    # mild heavy
  sym_t3       = function(n) stats::rt(n, 3)/sqrt(3),                      # moderate heavy
  sym_t2       = function(n) stats::rt(n, 2)/2,                            # heavy (m4 borderline)
  twosided_contam = function(n){e<-stats::rnorm(n);h<-runif(n)<0.12;k<-sum(h);if(k>0)e[h]<-stats::rt(k,2)*2.5;e}
)
gen <- function(n, seed=NULL){ if(!is.null(seed)) set.seed(seed); v1<-stats::rnorm(n); v2<-stats::rnorm(n)
  d<-kg2_design(v1,v2); Z<-as.matrix(cbind(1,d$b1,d$b2,d$b12,d$b11,d$b22)); list(v1=v1,v2=v2,Z=Z,mu=as.numeric(Z%*%THETA)) }

run_one <- function(seed, errfn, sm = 2.5, arms = c("LSE","Huber","L1","W1","WPMM3","PMM3")) {
  set.seed(seed); g <- gen(N_TRAIN); y_tr <- g$mu + errfn(N_TRAIN)
  d_train <- kg2_design(g$v1, g$v2, y_tr)
  Ztr <- as.matrix(cbind(1,d_train$b1,d_train$b2,d_train$b12,d_train$b11,d_train$b22))
  gt <- gen(N_TEST, seed + 10000L)
  emit <- function(nm, th){ e <- kg2_predict(th, gt$v1, gt$v2) - gt$mu
    data.frame(seed=seed, method=nm, mae=mean(abs(e)), coef_l2=sqrt(sum((th-THETA)^2)), stringsAsFactors=FALSE) }
  out <- list()
  if ("LSE"  %in% arms) out[[length(out)+1L]] <- emit("LSE",  tryCatch(fit_ref("LSE",  d_train), error=function(e) rep(NA,6)))
  if ("Huber"%in% arms) out[[length(out)+1L]] <- emit("Huber",tryCatch(fit_ref("Huber",d_train), error=function(e) rep(NA,6)))
  if ("L1"   %in% arms) out[[length(out)+1L]] <- emit("L1",   tryCatch(fit_ref("L1",   d_train), error=function(e) rep(NA,6)))
  if ("PMM3" %in% arms) out[[length(out)+1L]] <- emit("PMM3", tryCatch(fit_ref("PMM3", d_train), error=function(e) rep(NA,6)))
  if ("W1"   %in% arms) out[[length(out)+1L]] <- emit("W1",   tryCatch(fit_weak1(Ztr, d_train$y, sm), error=function(e) rep(NA,6)))
  if ("WPMM3"%in% arms) out[[length(out)+1L]] <- emit("WPMM3",tryCatch(fit_weak_pmm3(Ztr, d_train$y, sm), error=function(e) rep(NA,6)))
  do.call(rbind, out)
}
pair <- function(sub, m1, m2){ piv<-reshape(sub[,c("seed","method","mae")],idvar="seed",timevar="method",direction="wide")
  a<-piv[[paste0("mae.",m1)]]; b<-piv[[paste0("mae.",m2)]]; ok<-is.finite(a)&is.finite(b)
  if(sum(ok)<5) return(c(win=NA,d=NA,p=NA)); w<-suppressWarnings(stats::wilcox.test(a[ok],b[ok],paired=TRUE))
  c(win=100*mean(b[ok]>a[ok]), d=median(b[ok]-a[ok]), p=w$p.value) }

cat(sprintf("=== Weak-PMM3 probe | %d seeds | n_train=%d | sigma_mult=2.5 ===\n", N_SEEDS, N_TRAIN)); t0 <- Sys.time()

cat("\n[1] Head-to-head across symmetric regimes (MAE median):\n")
cat(sprintf("%-16s %7s %7s %7s %7s %7s %7s | %-20s %-20s\n",
    "regime","LSE","Huber","L1","W1","WPMM3","PMM3","WPMM3 vs W1","WPMM3 vs Huber"))
h2h <- list()
for (rg in names(REG)) {
  sub <- do.call(rbind, lapply(seq_len(N_SEEDS), function(s) run_one(SEED0+s-1L, REG[[rg]])))
  sub$regime <- rg; h2h[[rg]] <- sub
  med <- function(m) median(sub$mae[sub$method==m], na.rm=TRUE)
  vW1 <- pair(sub,"WPMM3","W1"); vHu <- pair(sub,"WPMM3","Huber")
  cat(sprintf("%-16s %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f | win%3.0f%% d%+.4f p%.0e | win%3.0f%% d%+.4f p%.0e\n",
      rg, med("LSE"),med("Huber"),med("L1"),med("W1"),med("WPMM3"),med("PMM3"),
      vW1["win"],vW1["d"],vW1["p"], vHu["win"],vHu["d"],vHu["p"]))
}
allh <- do.call(rbind, h2h)
outdir <- "../results"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(allh, file.path(outdir, "weak_pmm3_h2h_raw.csv"), row.names = FALSE)

cat("\n[2] Bandwidth sweep on sym_t3 (insurance: does WPMM3 hold where W1 drops?):\n")
cat(sprintf("%-8s %8s %8s %8s | %-22s\n","sigma","W1_mae","WPMM3_mae","Huber","WPMM3 vs W1"))
sweep <- list()
for (sm in c(1.5, 2.5, 4.0, 6.0, 8.0)) {
  sub <- do.call(rbind, lapply(seq_len(N_SEEDS), function(s) run_one(SEED0+s-1L, REG[["sym_t3"]], sm=sm, arms=c("Huber","W1","WPMM3"))))
  sub$sigma_mult <- sm; sweep[[as.character(sm)]] <- sub
  med <- function(m) median(sub$mae[sub$method==m], na.rm=TRUE); v <- pair(sub,"WPMM3","W1")
  cat(sprintf("%-8.1f %8.4f %8.4f %8.4f | win%3.0f%% d%+.4f p%.0e\n",
      sm, med("W1"), med("WPMM3"), med("Huber"), v["win"], v["d"], v["p"]))
}
utils::write.csv(do.call(rbind, sweep), file.path(outdir, "weak_pmm3_bandwidth_raw.csv"), row.names = FALSE)
cat(sprintf("\nSaved weak_pmm3_h2h_raw.csv + weak_pmm3_bandwidth_raw.csv | %.0fs\n", as.numeric(difftime(Sys.time(), t0, "secs"))))
