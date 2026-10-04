#!/usr/bin/env Rscript
# H3 validation: does the calibrated weak dispatch (force_method="auto-weak") pick
# a near-best inner estimator in EACH error regime, and never the catastrophic one?
#
# Per regime, compare the dispatched estimator against the fixed candidates
# {LSE, Huber, L1, WPMM2, WPMM3, PMM2}. H3 holds if auto-weak tracks the per-regime
# best-of-family (the oracle) and, in particular, avoids WPMM2 under light-tailed
# bulk skew (exp/lognorm) where it is harmful.
# Usage (from paper-4-weak-moment-gmdh/code): Rscript experiments/run_weak_dispatch_validation.R [N_SEEDS]

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
                 else pkgload::load_all(.find_pkg(), quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
N_SEEDS <- if (length(args) >= 1) as.integer(args[1]) else 150L
N_TRAIN <- 200L; N_TEST <- 4000L; SEED0 <- 75001L; BDIAG <- 80L
THETA <- c(0.5, 1.0, -0.8, 0.4, -0.3, 0.2)

REG <- list(
  gauss          = function(n) rnorm(n),
  uniform_platy  = function(n) runif(n, -sqrt(3), sqrt(3)),
  exp_skew       = function(n) rexp(n) - 1,
  lognorm_skew   = function(n){z<-exp(rnorm(n,0,0.6));(z-mean(z))/sd(z)},
  sym_t3         = function(n) rt(n,3)/sqrt(3),
  sym_t2         = function(n) rt(n,2)/2,
  twosided_contam= function(n){e<-rnorm(n);h<-runif(n)<0.12;k<-sum(h);if(k>0)e[h]<-rt(k,2)*2.5;e},
  onesided_contam= function(n){e<-rnorm(n);h<-runif(n)<0.12;k<-sum(h);if(k>0)e[h]<--abs(rt(k,2))*2.5;e},
  skewt_heavy    = function(n){z<-rt(n,3);(-abs(z)*0.8+0.4*z)}
)
FIXED <- c("LSE","Huber","L1","WPMM2","WPMM3","PMM2")
gen <- function(n, seed=NULL){ if(!is.null(seed)) set.seed(seed); v1<-rnorm(n); v2<-rnorm(n)
  d<-kg2_design(v1,v2); list(v1=v1,v2=v2,mu=as.numeric(cbind(1,d$b1,d$b2,d$b12,d$b11,d$b22)%*%THETA)) }

run_one <- function(seed, errfn) {
  set.seed(seed); g <- gen(N_TRAIN); y_tr <- g$mu + errfn(N_TRAIN); d_tr <- kg2_design(g$v1,g$v2,y_tr)
  gt <- gen(N_TEST, seed+10000L)
  mae <- function(th) mean(abs(kg2_predict(th, gt$v1, gt$v2) - gt$mu))
  out <- list()
  aw <- inner_estimate(g$v1, g$v2, y_tr, gmdh_pmm_control(B=BDIAG, force_method="auto-weak"))
  out[["auto-weak"]] <- list(mae=mae(aw$theta), disp=aw$method)
  for (m in FIXED) {
    th <- tryCatch(inner_estimate(g$v1, g$v2, y_tr, gmdh_pmm_control(B=0, force_method=m))$theta,
                   error=function(e) rep(NA,6))
    out[[m]] <- list(mae=mae(th), disp=m)
  }
  out
}

cat(sprintf("=== Weak dispatch validation (H3) | %d seeds ===\n", N_SEEDS)); t0 <- Sys.time()
rows <- list()
for (rg in names(REG)) {
  res <- lapply(seq_len(N_SEEDS), function(s) run_one(SEED0+s-1L, REG[[rg]]))
  med <- function(m) median(sapply(res, function(r) r[[m]]$mae), na.rm=TRUE)
  maes <- sapply(c("auto-weak", FIXED), med)
  best_fixed <- FIXED[which.min(maes[FIXED])]
  disp_tab <- table(sapply(res, function(r) r[["auto-weak"]]$disp))
  disp_str <- paste(sprintf("%s:%d", names(disp_tab), as.integer(disp_tab)), collapse=" ")
  cat(sprintf("\n%-16s dispatch={ %s }\n", rg, disp_str))
  cat(sprintf("  MAE: auto-weak=%.4f | best-fixed=%s(%.4f) | LSE=%.4f Huber=%.4f L1=%.4f WPMM2=%.4f WPMM3=%.4f PMM2=%.4f\n",
      maes["auto-weak"], best_fixed, maes[best_fixed], maes["LSE"], maes["Huber"], maes["L1"], maes["WPMM2"], maes["WPMM3"], maes["PMM2"]))
  cat(sprintf("  auto-weak vs best-fixed: %+.1f%% | vs LSE: %+.1f%%\n",
      100*(maes["auto-weak"]/maes[best_fixed]-1), 100*(maes["auto-weak"]/maes["LSE"]-1)))
  rows[[rg]] <- data.frame(regime=rg, auto_weak=maes["auto-weak"], best_fixed=best_fixed,
                           best_fixed_mae=maes[best_fixed], lse=maes["LSE"], dispatch=disp_str, stringsAsFactors=FALSE)
}
outdir <- "../results"; if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
utils::write.csv(do.call(rbind, rows), file.path(outdir, "weak_dispatch_validation.csv"), row.names=FALSE)
cat(sprintf("\nSaved weak_dispatch_validation.csv | %.0fs\n", as.numeric(difftime(Sys.time(), t0, "secs"))))
