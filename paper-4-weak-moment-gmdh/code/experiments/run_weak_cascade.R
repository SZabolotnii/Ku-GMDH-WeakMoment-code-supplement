#!/usr/bin/env Rscript
# Cascade benchmark (H1): does the weak dispatch (force_method="auto-weak") track
# the best-of-family at the FULL GMDH cascade level, across heterogeneous error
# regimes, where no single FIXED inner estimator wins everywhere?
#
# Target: a genuine 2-layer KG-2 composition of 6 inputs (x5,x6 irrelevant), so
# the MIA tournament must build real structure. Error regimes exercise the
# dispatch: gauss (LSE/PMM regime), exp_skew (clean skew -> classical PMM2),
# sym_t2 (symmetric heavy -> WPMM3), onesided_contam (asymmetric heavy -> WPMM2).
# Arms: LSE-GMDH, classical-PMM-GMDH (auto), Huber-GMDH, LAD-GMDH, WPMM2-forced
# (to show that forcing one weak estimator everywhere is NOT enough), weak-dispatch
# (auto-weak). Evaluation: MAE / 90%-trimmed RMSE on a CLEAN held-out test vs the
# noise-free truth.
# Usage (from paper-4-weak-moment-gmdh/code): Rscript experiments/run_weak_cascade.R [N_SEEDS] [B]

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
                 else pkgload::load_all(.find_pkg(), quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
N_SEEDS <- if (length(args) >= 1) as.integer(args[1]) else 30L
BDIAG   <- if (length(args) >= 2) as.integer(args[2]) else 50L
CRIT    <- if (length(args) >= 3) args[3] else "MSE"   # external selection criterion
N_TRAIN <- 320L; N_TEST <- 3000L; SEED0 <- 80001L; D <- 6L

REG <- list(
  gauss          = function(n) rnorm(n),
  exp_skew       = function(n) rexp(n) - 1,
  sym_t2         = function(n) rt(n,2)/2,
  onesided_contam= function(n){e<-rnorm(n);h<-runif(n)<0.12;k<-sum(h);if(k>0)e[h]<--abs(stats::rt(k,2))*2.5;e}
)
ARMS <- list(
  `LSE-GMDH`        = "LSE",
  `classicalPMM`    = "auto",
  `Huber-GMDH`      = "Huber",
  `LAD-GMDH`        = "L1",
  `WPMM2-forced`    = "WPMM2",
  `weak-dispatch`   = "auto-weak"
)

# Noise-free 2-layer KG-2 composition (x5,x6 irrelevant inputs).
truth <- function(X) {
  x1<-X[,1];x2<-X[,2];x3<-X[,3];x4<-X[,4]
  z1 <- 0.5 + 0.8*x1 - 0.6*x2 + 0.4*x1*x2 - 0.2*x1^2
  z2 <- -0.3 + 0.7*x3 + 0.5*x4 - 0.3*x3*x4 + 0.2*x4^2
  0.2 + 1.0*z1 - 0.8*z2 + 0.3*z1*z2
}
gen_X <- function(n) matrix(rnorm(n*D), n, D, dimnames=list(NULL, paste0("x",1:D)))
trmse <- function(e, p=0.90){ q<-stats::quantile(abs(e),p,names=FALSE); sqrt(mean(e[abs(e)<=q]^2)) }

ctrl_for <- function(fm, seed) gmdh_pmm_control(B = if (fm %in% c("auto","auto-weak")) BDIAG else 0L,
  force_method = fm, F = 6L, L_max = 3L, max_iter = 50L, weak_sigma_mult = 2.5,
  criterion = CRIT, seed = seed)

run_one <- function(seed, errfn) {
  set.seed(seed)
  Xtr <- gen_X(N_TRAIN); ytr_clean <- truth(Xtr); ytr <- ytr_clean + errfn(N_TRAIN)
  Xte <- gen_X(N_TEST);  yte_clean <- truth(Xte)
  sapply(names(ARMS), function(arm) {
    fit <- tryCatch(gmdh_pmm(Xtr, ytr, ctrl_for(ARMS[[arm]], seed)), error = function(e) NULL)
    if (is.null(fit)) return(c(mae=NA, trmse=NA))
    pred <- tryCatch(predict(fit, Xte), error = function(e) rep(NA, N_TEST))
    e <- pred - yte_clean
    c(mae = mean(abs(e)), trmse = trmse(e))
  })
}

cat(sprintf("=== Weak cascade benchmark (H1) | %d seeds | B=%d | crit=%s | d=%d, 2-layer KG-2 truth ===\n", N_SEEDS, BDIAG, CRIT, D))
t0 <- Sys.time(); allrows <- list()
for (rg in names(REG)) {
  res <- lapply(seq_len(N_SEEDS), function(s) run_one(SEED0+s-1L, REG[[rg]]))   # each: 2 x arms
  mae_mat <- t(sapply(res, function(r) r["mae", ]))   # seeds x arms
  med <- apply(mae_mat, 2, median, na.rm=TRUE); iqr <- apply(mae_mat, 2, IQR, na.rm=TRUE)
  ord <- order(med)
  cat(sprintf("\n--- %s (MAE median [IQR]) ---\n", rg))
  for (a in names(ARMS)[ord]) cat(sprintf("  %-14s %.4f [%.4f]\n", a, med[a], iqr[a]))
  wd <- mae_mat[, "weak-dispatch"]; bestfix_name <- names(ARMS)[setdiff(ord,which(names(ARMS)=="weak-dispatch"))[1]]
  bf <- mae_mat[, bestfix_name]; ok <- is.finite(wd)&is.finite(bf)
  w <- suppressWarnings(stats::wilcox.test(wd[ok], bf[ok], paired=TRUE))
  cat(sprintf("  weak-dispatch vs best-fixed (%s): %+.1f%% median, win %.0f%%, p=%.2e\n",
      bestfix_name, 100*(median(wd,na.rm=TRUE)/median(bf,na.rm=TRUE)-1), 100*mean(bf[ok]>wd[ok]), w$p.value))
  for (a in names(ARMS)) allrows[[length(allrows)+1L]] <- data.frame(regime=rg, arm=a,
    mae_med=med[a], mae_iqr=iqr[a], stringsAsFactors=FALSE)
}
outdir <- "../results"; if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
fn <- sprintf("weak_cascade_summary_%s.csv", CRIT)
utils::write.csv(do.call(rbind, allrows), file.path(outdir, fn), row.names=FALSE)
cat(sprintf("\nSaved %s | %.0fs\n", fn, as.numeric(difftime(Sys.time(), t0, "secs"))))
