#!/usr/bin/env Rscript
# LEAKAGE AUDIT. All earlier real-data comparisons used RANDOM 70/30 splits, but the
# data are autocorrelated time series (1 Hz vibration; depth-ordered ROP). Random
# splitting puts temporally-adjacent points in train and test -> optimistic, possibly
# inflating the win. This checks whether the headline vibration win survives
# leakage-controlled (blocked / contiguous-holdout) splits. Time order is PRESERVED
# (no random subsampling). Drilling domain. Usage: Rscript experiments/run_leakage_audit.R

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))
trmse <- function(e,p=0.90){ q<-stats::quantile(abs(e),p,names=FALSE,na.rm=TRUE); sqrt(mean(e[abs(e)<=q]^2)) }
fitm <- function(m,b,r,yy) tryCatch(inner_estimate(b,r,yy, gmdh_pmm_control(B=0,force_method=m,max_iter=60L,weak_sigma_mult=2.5))$theta, error=function(e) rep(NA,6))
METH <- c("LSE","Huber","L1","WPMM2","WPMM3")
ev <- function(b,r,yy,tri,tei){ sapply(METH, function(m) trmse(yy[tei]-kg2_predict(fitm(m,b[tri],r[tri],yy[tri]),b[tei],r[tei]))) }

# load vibration, PRESERVE time order, take a contiguous segment (no random subsample)
f <- list.files("/tmp/forge78b_sensor", pattern="Run 2 - output data.csv$", recursive=TRUE, full.names=TRUE)[1]
d <- data.table::fread(f, skip=35L, header=TRUE, fill=TRUE, showProgress=FALSE)
S<-"CSS-008"; g<-function(c) suppressWarnings(as.numeric(d[[paste0(S,"_",c)]]))
ss<-g("StickSlip(%)"); gs<-g("GyroXspread(RPM)"); tg<-g("ShZpeak(g)")
ok<-is.finite(ss)&is.finite(gs)&is.finite(tg)&tg>0            # ok preserves row (time) order
b<-ss[ok]; r<-gs[ok]; y<-log(tg[ok]); n<-length(y)
if (n>8000) { b<-b[1:8000]; r<-r[1:8000]; y<-y[1:8000]; n<-8000 }   # contiguous head segment
cat(sprintf("=== Leakage audit: vibration log(ShZpeak) | n=%d (contiguous, time-ordered) | g3=%+.2f g4=%.1f ===\n",
            n, sample_cumulants(resid(lm(y~b+r)))$gamma3, sample_cumulants(resid(lm(y~b+r)))$gamma4))
report <- function(tag, mat){ md<-apply(mat,2,median,na.rm=TRUE); wb<-min(md["WPMM2"],md["WPMM3"]); wn<-names(which.min(md[c("WPMM2","WPMM3")]))
  cat(sprintf("  %-18s LSE=%.4f Huber=%.4f L1=%.4f WPMM2=%.4f WPMM3=%.4f | best-weak(%s) vs Huber %+.1f%%\n",
      tag, md["LSE"],md["Huber"],md["L1"],md["WPMM2"],md["WPMM3"], wn, 100*(wb/md["Huber"]-1))) }

# (1) RANDOM 70/30 (the original, leakage-prone protocol)
m1 <- t(sapply(1:25, function(s){ set.seed(80000+s); tri<-sample.int(n,floor(.7*n)); ev(b,r,y,tri,setdiff(1:n,tri)) }))
report("random-70/30", m1)
# (2) BLOCKED 5-fold: each interior contiguous 20% block is test, the rest is train (+gap)
gap <- 50L; folds <- 5L; fb <- floor(n/folds)
m2 <- t(sapply(1:folds, function(k){ a<-(k-1)*fb+1; z<-if(k==folds) n else k*fb; tei<-a:z
  tri<-setdiff(1:n, max(1,a-gap):min(n,z+gap)); if(length(tri)<50||length(tei)<20) return(rep(NA,length(METH))); ev(b,r,y,tri,tei) }))
report("blocked-5fold", m2)
# (3) CONTIGUOUS holdout: train first 70% (time), test last 30%
report("contig-holdout", matrix(ev(b,r,y, 1:floor(.7*n), (floor(.7*n)+1):n), 1, dimnames=list(NULL,METH)))
# autocorrelation of residuals (why random leaks)
r0 <- resid(lm(y~b+r)); ac1 <- acf(r0, lag.max=1, plot=FALSE)$acf[2]
cat(sprintf("\n  residual lag-1 autocorrelation = %.3f (high => random split leaks)\n", ac1))
outdir<-"../results"; if(!dir.exists(outdir)) dir.create(outdir,recursive=TRUE)
utils::write.csv(rbind(cbind(scheme="random",as.data.frame(t(apply(m1,2,median,na.rm=TRUE)))),
                       cbind(scheme="blocked5",as.data.frame(t(apply(m2,2,median,na.rm=TRUE))))),
                 file.path(outdir,"leakage_audit_vibration.csv"), row.names=FALSE)
cat("Saved leakage_audit_vibration.csv\n")
