#!/usr/bin/env Rscript
# Regime map for the weak-PMM2 a2 (degree-2) term (Paper 4).
#
# The mechanism ablation showed the a2 windowed-PMM term (multiplicative
# asymmetric tilt W = w*(1+a2*rc), NOT WM-1's additive a2*(r^2-s^2)) is the
# DOMINANT driver of the WPMM2 win under one-sided heavy contamination, while
# re-centring is inert. WM-1 reported degree-2 is ~inert for asymmetric bias on
# skew-t / alpha-stable. To reconcile: map the a2 effect (W10 degree-1 recentred
# -> W11 = +a2) across error regimes. Hypothesis: a2 is ACTIVE where the
# (windowed) asymmetry is strong and STABLE (one-sided contamination, light-tail
# skew), and INERT under symmetric heavy tails and weak/continuous skew.
#
# Same KG-2 DGP as the PoC; metric = MAE & intercept bias vs the noise-free truth.
# Usage (from paper-4-weak-moment-gmdh/code): Rscript experiments/run_weak_regime_map.R [N_SEEDS]

.find_pkg <- function() {
  cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent
  }
}
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
                 else pkgload::load_all(.find_pkg(), quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
N_SEEDS <- if (length(args) >= 1) as.integer(args[1]) else 150L
N_TRAIN <- 200L; N_TEST <- 4000L; SEED0 <- 75001L
THETA <- c(0.5, 1.0, -0.8, 0.4, -0.3, 0.2)

# Error regimes (each ~unit bulk scale). g3/g4 character noted.
REGIMES <- list(
  gauss          = function(n) stats::rnorm(n),                                   # symmetric light (clean)
  onesided_contam= function(n){e<-stats::rnorm(n);h<-runif(n)<0.12;k<-sum(h);if(k>0)e[h]<--abs(stats::rt(k,2))*2.5;e}, # STABLE one-sided skew + heavy
  twosided_contam= function(n){e<-stats::rnorm(n);h<-runif(n)<0.12;k<-sum(h);if(k>0)e[h]<-stats::rt(k,2)*2.5;e},        # symmetric heavy contamination
  sym_t3         = function(n) stats::rt(n, df = 3) / sqrt(3),                     # symmetric heavy, light contamination
  exp_skew       = function(n) (stats::rexp(n) - 1),                              # right-skew, LIGHT tail (Paper-1 clean-skew analog)
  lognorm_skew   = function(n){z<-exp(stats::rnorm(n,0,0.6));(z-mean(z))/sd(z)},   # strong stable right-skew, moderate tail
  skewt_heavy    = function(n){z<-stats::rt(n,3);(-abs(z)*0.8 + 0.4*z)}            # continuous skew + heavy (WM-1-like)
)

wridge <- function(Z, y, w, lambda = 1e-8) {
  w[!is.finite(w) | w < 0] <- 0
  A <- crossprod(Z, Z * w); diag(A) <- diag(A) + lambda; b <- crossprod(Z, w * y)
  tryCatch(as.numeric(solve(A, b)), error = function(e) rep(NA_real_, ncol(Z)))
}
fit_weak <- function(Z, y, sigma_mult = 2.5, recenter = TRUE, use_a2 = TRUE, max_iter = 60L, tol = 1e-7) {
  theta <- wridge(Z, y, rep(1, length(y)))
  for (iter in seq_len(max_iter)) {
    r <- as.numeric(y - Z %*% theta)
    mu <- if (recenter) stats::median(r) else 0
    s <- stats::mad(r, constant = 1.4826); if (!is.finite(s) || s <= 1e-12) s <- stats::sd(r); if (!is.finite(s) || s <= 1e-12) s <- 1
    w <- exp(-0.5 * ((r - mu)/(sigma_mult * s))^2); sws <- sum(w); if (!is.finite(sws) || sws <= 0) break
    rc <- r - sum(w * r)/sws
    if (use_a2) {
      c2 <- sum(w*rc^2)/sws; c3 <- sum(w*rc^3)/sws; c4 <- sum(w*rc^4)/sws - 3*c2^2; denom <- 2*c2^2 + c4
      a2 <- if (!is.finite(c2)||c2<=1e-12||!is.finite(denom)||abs(denom)<1e-12) 0 else -c3/denom
      a2 <- max(min(a2, 0.5/sqrt(c2+1e-12)), -0.5/sqrt(c2+1e-12))
    } else a2 <- 0
    W <- w*(1 + a2*rc); W[!is.finite(W)|W<0] <- 0
    tn <- wridge(Z, y, W); if (any(!is.finite(tn))) break
    if (sqrt(sum((tn-theta)^2))/max(1, sqrt(sum(theta^2))) < tol) { theta <- tn; break }
    theta <- tn
  }
  theta
}
fit_ref <- function(method, d) inner_estimate(d$b1, d$b2, d$y,
  gmdh_pmm_control(B = 0, force_method = method, max_iter = 60L, weak_sigma_mult = 2.5))$theta

gen <- function(n, seed=NULL){ if(!is.null(seed)) set.seed(seed); v1<-stats::rnorm(n); v2<-stats::rnorm(n)
  d<-kg2_design(v1,v2); Z<-as.matrix(cbind(1,d$b1,d$b2,d$b12,d$b11,d$b22)); list(v1=v1,v2=v2,Z=Z,mu=as.numeric(Z%*%THETA)) }

run_one <- function(seed, reg, errfn) {
  set.seed(seed); g <- gen(N_TRAIN); y_tr <- g$mu + errfn(N_TRAIN)
  d_train <- kg2_design(g$v1, g$v2, y_tr)
  Ztr <- as.matrix(cbind(1,d_train$b1,d_train$b2,d_train$b12,d_train$b11,d_train$b22))
  gt <- gen(N_TEST, seed + 10000L)
  emit <- function(nm, th){ e <- kg2_predict(th, gt$v1, gt$v2) - gt$mu
    data.frame(seed=seed, regime=reg, method=nm, mae=mean(abs(e)), b0_abs=abs(th[1]-THETA[1]), stringsAsFactors=FALSE) }
  rbind <- list(
    emit("W10", tryCatch(fit_weak(Ztr, d_train$y, recenter=TRUE,  use_a2=FALSE), error=function(e) rep(NA,6))),
    emit("W11", tryCatch(fit_weak(Ztr, d_train$y, recenter=TRUE,  use_a2=TRUE ), error=function(e) rep(NA,6))),
    emit("Huber", tryCatch(fit_ref("Huber", d_train), error=function(e) rep(NA,6))),
    emit("L1",    tryCatch(fit_ref("L1",    d_train), error=function(e) rep(NA,6))),
    emit("LSE",   tryCatch(fit_ref("LSE",   d_train), error=function(e) rep(NA,6))))
  do.call(base::rbind, rbind)
}

cat(sprintf("=== Weak a2 regime map | %d seeds | n_train=%d ===\n", N_SEEDS, N_TRAIN))
t0 <- Sys.time()
all <- do.call(base::rbind, lapply(names(REGIMES), function(rg)
  do.call(base::rbind, lapply(seq_len(N_SEEDS), function(s) run_one(SEED0 + s - 1L, rg, REGIMES[[rg]])))))
outdir <- "../results"; if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(all, file.path(outdir, "weak_regime_map_raw.csv"), row.names = FALSE)

cat(sprintf("\n%-16s %8s %8s %8s %8s | %-22s %-18s\n", "regime","W11_mae","W10_mae","Hub_mae","L1_mae","a2 effect W10->W11","W11 vs best robust"))
for (rg in names(REGIMES)) {
  sub <- all[all$regime == rg, ]
  piv <- reshape(sub[,c("seed","method","mae")], idvar="seed", timevar="method", direction="wide")
  pivb<- reshape(sub[,c("seed","method","b0_abs")], idvar="seed", timevar="method", direction="wide")
  med <- function(m) median(sub$mae[sub$method==m], na.rm=TRUE)
  a<-piv[["mae.W11"]]; b<-piv[["mae.W10"]]; ok<-is.finite(a)&is.finite(b)
  w<-suppressWarnings(stats::wilcox.test(a[ok],b[ok],paired=TRUE))
  db0 <- median(pivb[["b0_abs.W11"]]-pivb[["b0_abs.W10"]], na.rm=TRUE)
  rob <- if (med("Huber") <= med("L1")) "Huber" else "L1"
  c2<-piv[["mae.W11"]]; d2<-piv[[paste0("mae.",rob)]]; ok2<-is.finite(c2)&is.finite(d2)
  w2<-suppressWarnings(stats::wilcox.test(c2[ok2],d2[ok2],paired=TRUE))
  cat(sprintf("%-16s %8.4f %8.4f %8.4f %8.4f | dMAE=%+.4f p=%.1e db0=%+.3f | vs %-5s %+.4f p=%.1e\n",
      rg, med("W11"), med("W10"), med("Huber"), med("L1"),
      median(b[ok]-a[ok]), w$p.value, db0, rob, median(d2[ok2]-c2[ok2]), w2$p.value))
}
cat(sprintf("\nSaved weak_regime_map_raw.csv | %.0fs\n", as.numeric(difftime(Sys.time(), t0, "secs"))))
