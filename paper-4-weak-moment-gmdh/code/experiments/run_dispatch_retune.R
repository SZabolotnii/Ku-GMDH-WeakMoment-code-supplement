#!/usr/bin/env Rscript
# DISPATCH RE-TUNING test: does the DEPLOYABLE dispatch capture the oracle forced-method
# wins on the cross-domain positives (credit_balance, sru_y2_static, gas_turbine_nox)?
# Compares per dataset, under the honest protocol: forced {LSE,Huber,L1,WPMM2,WPMM3},
# the two dispatches {auto-weak (cumulant gate), auto-valgate (validation gate)}, and the
# oracle reference best-PMM=min(WPMM2,WPMM3) / best-robust=min(LSE,Huber,L1). Run this
# BEFORE and AFTER editing the auto-valgate inner CV to measure the improvement.
# Usage: Rscript experiments/run_dispatch_retune.R

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg(); ROOT <- dirname(dirname(PKG))
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))
trmse <- function(e, p = 0.90) { e <- e[is.finite(e)]; q <- stats::quantile(abs(e), p, names = FALSE); sqrt(mean(e[abs(e) <= q]^2)) }
mae <- function(e) mean(abs(e[is.finite(e)]))

# --- loaders (mirror run_realworld.R) ---
model_x <- function(formula, data) { X <- stats::model.matrix(formula, data = data); keep <- colnames(X) != "(Intercept)"
  X <- X[, keep, drop = FALSE]; X <- X[, apply(X,2,function(z) all(is.finite(z)) && stats::sd(z) > 1e-12), drop = FALSE]
  storage.mode(X) <- "double"; X }
ext_load <- function(id) { man <- utils::read.csv(file.path(ROOT,"shared/datasets/external/external_candidates.csv"), stringsAsFactors = FALSE)
  row <- man[match(id, man$id), ]; p <- row$path[[1]]; if (!grepl("^/", p)) p <- file.path(ROOT, p)
  d <- stats::na.omit(utils::read.csv(p, stringsAsFactors = TRUE)); f <- stats::as.formula(row$formula[[1]])
  mf <- stats::model.frame(f, data = d, na.action = stats::na.omit)
  list(y = as.numeric(stats::model.response(mf)), X = model_x(f, mf)) }
credit_load <- function() { d <- stats::na.omit(ISLR2::Credit)
  list(y = d$Balance, X = model_x(~ Income+Limit+Rating+Cards+Age+Education+Own+Student+Married+Region, d)) }

DATASETS <- list(
  list(id="islr_credit_balance", proto="random", get=credit_load),
  list(id="sru_y2_static",       proto="blocked", get=function() ext_load("sru_y2_static")),
  list(id="gas_turbine_nox_2015_raw", proto="blocked", get=function() ext_load("gas_turbine_nox_2015_raw"))
)
METHODS <- c("LSE","Huber","L1","WPMM2","WPMM3","auto-weak","auto-valgate")
ctrl_for <- function(m, seed) gmdh_pmm_control(L_max=3L, F=6L, epsilon=-Inf, seed=seed,
  B = if (m %in% c("auto-weak")) 60L else 0L, force_method = m)
std <- function(X, tr, te) { mu <- colMeans(X[tr,,drop=FALSE]); sg <- apply(X[tr,,drop=FALSE],2,stats::sd); sg[!is.finite(sg)|sg<=1e-12] <- 1
  list(tr = sweep(sweep(X[tr,,drop=FALSE],2,mu,"-"),2,sg,"/"), te = sweep(sweep(X[te,,drop=FALSE],2,mu,"-"),2,sg,"/")) }
# average each (method,split) over SEEDS to damp tournament + bootstrap-gate stochasticity
SEEDS <- 1:4
fit_eval <- function(X, y, tr, te) { xs <- std(X, tr, te)
  sapply(METHODS, function(m) {
    vals <- sapply(SEEDS, function(s) { f <- try(gmdh_pmm(xs$tr, y[tr], ctrl_for(m, 1000L+s)), silent = TRUE)
      if (inherits(f,"try-error")) return(c(NA,NA)); e <- predict(f, xs$te) - y[te]; c(trmse(e), mae(e)) })
    c(mean(vals[1,], na.rm=TRUE), mean(vals[2,], na.rm=TRUE)) }) }

cat("=== Dispatch re-tuning: does the deployable dispatch capture the oracle win? ===\n")
rows <- list()
for (ds in DATASETS) {
  dat <- ds$get(); X <- dat$X; y <- dat$y; n <- length(y); cap <- 3000L
  if (ds$proto == "blocked" && n > cap) { X <- X[1:cap,,drop=FALSE]; y <- y[1:cap]; n <- cap }   # contiguous head, time order
  TR <- matrix(NA, 0, length(METHODS)); MA <- matrix(NA, 0, length(METHODS))
  if (ds$proto == "random") {
    for (s in 1:12) { set.seed(5000+s); te <- sort(sample.int(n, ceiling(0.25*n))); tr <- setdiff(seq_len(n), te)
      r <- fit_eval(X,y,tr,te); TR <- rbind(TR, r[1,]); MA <- rbind(MA, r[2,]) }
  } else {
    folds <- 5L; gap <- 20L; fb <- floor(n/folds)
    for (k in 1:folds) { a<-(k-1)*fb+1; z<-if(k==folds) n else k*fb; te<-a:z; tr<-setdiff(1:n, max(1,a-gap):min(n,z+gap))
      if (length(tr)<80||length(te)<30) next; r <- fit_eval(X,y,tr,te); TR <- rbind(TR, r[1,]); MA <- rbind(MA, r[2,]) }
  }
  colnames(TR) <- METHODS; colnames(MA) <- METHODS
  mt <- apply(TR,2,median,na.rm=TRUE); mm <- apply(MA,2,median,na.rm=TRUE)
  bestpmm <- min(mt["WPMM2"],mt["WPMM3"]); bestrob <- min(mt["LSE"],mt["Huber"],mt["L1"])
  cat(sprintf("\n[%s | %s | n=%d]\n", ds$id, ds$proto, n))
  for (m in METHODS) cat(sprintf("  %-13s trmse=%.4f  mae=%.4f\n", m, mt[m], mm[m]))
  cat(sprintf("  ORACLE best-PMM=%.4f  best-robust=%.4f\n", bestpmm, bestrob))
  cat(sprintf("  >> auto-weak vs best-robust   %+6.1f%% trmse\n", 100*(mt["auto-weak"]/bestrob-1)))
  cat(sprintf("  >> auto-valgate vs best-robust %+6.1f%% trmse | vs oracle-best-PMM %+6.1f%%\n",
              100*(mt["auto-valgate"]/bestrob-1), 100*(mt["auto-valgate"]/bestpmm-1)))
  rows[[length(rows)+1L]] <- data.frame(dataset=ds$id, proto=ds$proto, n=n,
    t(round(mt,4)), oracle_bestpmm=round(bestpmm,4), best_robust=round(bestrob,4),
    autoweak_vs_rob=round(100*(mt["auto-weak"]/bestrob-1),2),
    autovalgate_vs_rob=round(100*(mt["auto-valgate"]/bestrob-1),2),
    autovalgate_vs_oracle=round(100*(mt["auto-valgate"]/bestpmm-1),2),
    check.names=FALSE, row.names=NULL) }
outdir <- file.path(ROOT,"paper-4-weak-moment-gmdh","results"); if(!dir.exists(outdir)) dir.create(outdir,recursive=TRUE)
utils::write.csv(do.call(rbind,rows), file.path(outdir,"dispatch_retune.csv"), row.names=FALSE)
cat("\nSaved dispatch_retune.csv\n")
