#!/usr/bin/env Rscript
# Full-cascade test of the validation-gated dispatch on impulsive drilling vibration.
# Earlier the vibration win + valgate were shown at the single-KG-2 (node) level; here
# we run the FULL gmdh_pmm() tournament on multi-feature vibration data and confirm the
# win holds end-to-end. Target log(ShZpeak); predictors = 5 downhole dynamics channels
# (StickSlip, GyroXspread, GyroXmed, GyroXmax, Temperature) so the MIA tournament builds
# real structure. Arms: forced LSE/Huber/L1/WPMM2/WPMM3-GMDH, cumulant auto-weak, and
# validation-gated auto-valgate. Outer 70/30 split (fit on train, score held-out test by
# trimmed-RMSE); MAE external criterion (heavy-tail-safe selection). Drilling domain only.
# Usage: Rscript experiments/run_vibration_cascade.R [R]

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args) >= 1) as.integer(args[1]) else 20L
trmse <- function(e,p=0.90){ q<-stats::quantile(abs(e),p,names=FALSE,na.rm=TRUE); sqrt(mean(e[abs(e)<=q]^2)) }
cap <- 3000L; SENSOR <- "CSS-008"
PRED <- paste0(SENSOR, "_", c("StickSlip(%)","GyroXspread(RPM)","GyroXmed(RPM)","GyroXmax(RPM)","Temperature(C)"))
TARGET <- paste0(SENSOR, "_ShZpeak(g)")
ARMS <- list(`LSE-GMDH`="LSE",`Huber-GMDH`="Huber",`LAD-GMDH`="L1",`WPMM2-GMDH`="WPMM2",
             `WPMM3-GMDH`="WPMM3",`auto-weak`="auto-weak",`auto-valgate`="auto-valgate")

f <- list.files("/tmp/forge78b_sensor", pattern="Run 2 - output data.csv$", recursive=TRUE, full.names=TRUE)[1]
d <- data.table::fread(f, skip=35L, header=TRUE, fill=TRUE, showProgress=FALSE)
X <- sapply(PRED, function(c) suppressWarnings(as.numeric(d[[c]])))
y <- log(suppressWarnings(as.numeric(d[[TARGET]])))
ok <- is.finite(y) & apply(is.finite(X),1,all)
X <- X[ok,,drop=FALSE]; y <- y[ok]; colnames(X) <- c("stick","gspread","gmed","gmax","temp")
n <- nrow(X); cat(sprintf("=== Vibration full-cascade (%s, log ShZpeak) | n=%d | %d preds | %d seeds ===\n", SENSOR, n, ncol(X), R))
if (n > cap) { set.seed(1); ii<-sort(sample.int(n,cap)); X<-X[ii,];y<-y[ii]; n<-cap }
cu <- sample_cumulants(resid(lm(y ~ X))); cat(sprintf("global resid g3=%+.2f g4=%.1f\n", cu$gamma3, cu$gamma4))

ctrl_for <- function(fm) gmdh_pmm_control(B = if (fm %in% c("auto-weak")) 50L else 0L,
  force_method = fm, F = 6L, L_max = 3L, max_iter = 60L, weak_sigma_mult = 2.5, criterion = "MAE")
mat <- matrix(NA, R, length(ARMS), dimnames=list(NULL, names(ARMS)))
for (rr in seq_len(R)) {
  set.seed(95000 + rr); tri <- sample.int(n, floor(0.7*n)); tei <- setdiff(seq_len(n), tri)
  for (a in names(ARMS)) {
    fit <- tryCatch(gmdh_pmm(X[tri,], y[tri], ctrl_for(ARMS[[a]])), error=function(e) NULL)
    mat[rr,a] <- if (is.null(fit)) NA else trmse(y[tei] - tryCatch(predict(fit, X[tei,]), error=function(e) rep(NA,length(tei))))
  }
}
med <- apply(mat,2,median,na.rm=TRUE); best <- names(which.min(med))
cat("\ntest trimmed-RMSE median (full cascade):\n")
for (a in names(med)[order(med)]) cat(sprintf("   %-13s %.4f%s\n", a, med[a], if(a==best)"  <--" else ""))
for (cmp in c("Huber-GMDH","LAD-GMDH","auto-weak")) {
  av<-mat[,"auto-valgate"]; b<-mat[,cmp]; o<-is.finite(av)&is.finite(b); if(sum(o)<5) next
  w<-suppressWarnings(stats::wilcox.test(av[o],b[o],paired=TRUE))
  cat(sprintf("   auto-valgate vs %-11s: %+.1f%% median, win %2.0f%%, p=%.2e\n", cmp, 100*(median(av,na.rm=TRUE)/median(b,na.rm=TRUE)-1), 100*mean(b[o]>av[o]), w$p.value))
}
outdir<-"../results"; if(!dir.exists(outdir)) dir.create(outdir,recursive=TRUE)
utils::write.csv(data.frame(arm=names(med), trmse_med=as.numeric(med)), file.path(outdir,"vibration_cascade.csv"), row.names=FALSE)
cat("\nSaved vibration_cascade.csv\n")
