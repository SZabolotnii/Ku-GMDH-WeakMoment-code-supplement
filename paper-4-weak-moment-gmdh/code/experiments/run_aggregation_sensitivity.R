#!/usr/bin/env Rscript
# AGGREGATION-SENSITIVITY audit (route B, the second material factor after leakage).
# The honest_blocked_cv.R "7/12 vibration wins" used MEDIAN-of-per-fold trimmed-RMSE;
# the pooled-significance test used a single POOLED-error trimmed-RMSE over all out-of-fold
# predictions. On the non-stationary drilling series these DISAGREE — and the per-fold
# matrix shows why: the weak estimators blow up on some folds and win on others, so the
# median masks high cross-fold variance. This tabulates, per cell, the weak-vs-best-baseline
# advantage under BOTH aggregations (+ the standard deviation of weak's per-fold trmse as an
# instability flag). Same blocked CV as run_honest_blocked_cv*.R. Drilling only. No bootstrap.
# Usage: Rscript experiments/run_aggregation_sensitivity.R

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))

trmse <- function(e, p = 0.90) { q <- stats::quantile(abs(e), p, names = FALSE, na.rm = TRUE); sqrt(mean(e[abs(e) <= q]^2)) }
M <- c("LSE","Huber","L1","WPMM2","WPMM3"); folds <- 5L; cap <- 6000L
fitm <- function(m,b,r,yy) tryCatch(inner_estimate(b,r,yy, gmdh_pmm_control(B=0,force_method=m,max_iter=60L,weak_sigma_mult=2.5))$theta, error=function(e) rep(NA,6))

one_cell <- function(b, r, y, gap) {
  n <- length(y); fb <- floor(n/folds)
  pf <- matrix(NA, folds, length(M), dimnames=list(NULL,M))      # per-fold trmse
  E  <- matrix(NA, n, length(M), dimnames=list(NULL,M))          # pooled errors
  for (k in 1:folds) { a<-(k-1)*fb+1; z<-if(k==folds) n else k*fb; tei<-a:z; tri<-setdiff(1:n, max(1,a-gap):min(n,z+gap))
    if (length(tri)<100||length(tei)<40) next
    for (m in M) { e <- y[tei]-kg2_predict(fitm(m,b[tri],r[tri],y[tri]),b[tei],r[tei]); pf[k,m]<-trmse(e); E[tei,m]<-e } }
  med <- apply(pf,2,median,na.rm=TRUE); pl <- apply(E,2,function(c) trmse(c[is.finite(c)]))
  wb_med <- 100*(min(med["WPMM2"],med["WPMM3"])/min(med["LSE"],med["Huber"],med["L1"])-1)
  wb_pl  <- 100*(min(pl["WPMM2"],pl["WPMM3"]) /min(pl["LSE"],pl["Huber"],pl["L1"]) -1)
  wsd <- max(sd(pf[,"WPMM2"],na.rm=TRUE)/mean(pf[,"WPMM2"],na.rm=TRUE),
             sd(pf[,"WPMM3"],na.rm=TRUE)/mean(pf[,"WPMM3"],na.rm=TRUE))           # weak cross-fold CV
  bsd <- min(sd(pf[,"Huber"],na.rm=TRUE)/mean(pf[,"Huber"],na.rm=TRUE),
             sd(pf[,"L1"],na.rm=TRUE)/mean(pf[,"L1"],na.rm=TRUE))                  # robust cross-fold CV
  list(wb_med=wb_med, wb_pl=wb_pl, weak_foldCV=wsd, robust_foldCV=bsd)
}

cells <- list(); TMP <- "/tmp/forge78b_sensor"
.f1 <- function(p){ x<-list.files(TMP,pattern=p,recursive=TRUE,full.names=TRUE); if(length(x)) x[1] else NA }
for (rn in c("Run 1 - output data\\.csv$","Run 2 - output data\\.csv$")) { f<-.f1(rn); if(is.na(f)) next
  tag0 <- if (grepl("1 -",rn)) "R1" else "R2"; d<-data.table::fread(f,skip=35L,header=TRUE,fill=TRUE,showProgress=FALSE)
  for (sn in c("CSS-008","CSS-007")) for (ch in c("ShZpeak","ShYpeak","ShZrms")) {
    g<-function(s) suppressWarnings(as.numeric(d[[paste0(sn,"_",s)]])); ss<-g("StickSlip(%)"); gs<-g("GyroXspread(RPM)"); tg<-g(paste0(ch,"(g)"))
    if(is.null(ss)||is.null(gs)||is.null(tg)) next; ok<-is.finite(ss)&is.finite(gs)&is.finite(tg)&tg>0
    b<-ss[ok]; r<-gs[ok]; y<-log(tg[ok]); if(length(y)<600) next; if(length(y)>cap){b<-b[1:cap];r<-r[1:cap];y<-y[1:cap]}
    cells[[length(cells)+1L]]<-list(tag=sprintf("vib:%s/%s/%s",tag0,sn,ch), g4=sample_cumulants(resid(lm(y~b+r)))$gamma4, b=b,r=r,y=y,gap=50L) } }
dv<-tryCatch(load_volve_drilling(file.path(PKG,"data/volve_onbottom.csv")),error=function(e) NULL)
if(!is.null(dv)){ X<-as.data.frame(dv$X); y<-as.numeric(dv$y); o<-order(X$dept); X<-X[o,]; y<-y[o]
  for(bn in list(c(2000,2080),c(1740,1820),c(2110,2200),c(1360,1480))){ sel<-X$dept>=bn[1]&X$dept<bn[2]; if(sum(sel)<150) next
    Xi<-X[sel,]; yi<-y[sel]; cells[[length(cells)+1L]]<-list(tag=sprintf("rop:volve_%g_%g",bn[1],bn[2]),
      g4=sample_cumulants(resid(lm(yi~Xi$swob+Xi$tqa)))$gamma4, b=Xi$swob,r=Xi$tqa,y=yi,gap=30L) } }

cat("=== Aggregation sensitivity: weak vs best-baseline under MEDIAN-of-fold vs POOLED-error ===\n")
cat(sprintf("%-26s %6s | %10s %10s | %9s | flip? | weak-foldCV robust-foldCV\n","cell","g4","median%","pooled%","|diff|"))
rows<-list()
for(ci in cells){ res<-one_cell(ci$b,ci$r,ci$y,ci$gap)
  flip <- sign(res$wb_med) != sign(res$wb_pl)
  cat(sprintf("%-26s %6.2f | %+9.1f %+9.1f | %8.1f | %-5s | %6.2f %6.2f\n",
      ci$tag, ci$g4, res$wb_med, res$wb_pl, abs(res$wb_med-res$wb_pl), ifelse(flip,"FLIP",""), res$weak_foldCV, res$robust_foldCV))
  rows[[length(rows)+1L]]<-data.frame(cell=ci$tag, g4=ci$g4, weak_vs_base_median=res$wb_med, weak_vs_base_pooled=res$wb_pl,
    sign_flip=flip, weak_foldCV=res$weak_foldCV, robust_foldCV=res$robust_foldCV, stringsAsFactors=FALSE) }
R<-do.call(rbind,rows)
cat(sprintf("\nMEDIAN agg: weak beats baseline (>2%%) in %d/%d ; POOLED agg: %d/%d ; sign flips: %d/%d\n",
  sum(R$weak_vs_base_median< -2), nrow(R), sum(R$weak_vs_base_pooled< -2), nrow(R), sum(R$sign_flip), nrow(R)))
cat(sprintf("weak cross-fold CV (median over cells) = %.2f vs robust = %.2f  (higher => weak less stable across folds)\n",
  median(R$weak_foldCV,na.rm=TRUE), median(R$robust_foldCV,na.rm=TRUE)))
outdir<-file.path(dirname(dirname(PKG)),"paper-4-weak-moment-gmdh","results"); if(!dir.exists(outdir)) dir.create(outdir,recursive=TRUE)
utils::write.csv(R, file.path(outdir,"aggregation_sensitivity.csv"), row.names=FALSE)
cat("Saved aggregation_sensitivity.csv\n")
