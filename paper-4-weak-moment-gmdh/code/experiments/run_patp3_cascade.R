#!/usr/bin/env Rscript
# Tests the curator's idea: at the well-specified cascade the vibration residual is
# symmetric-heavy (gamma3~0, gamma4~13) — the regime for the FRACTIONAL-power PMM3
# (PATP i=3), which needs only nu_{2p}<nu_4 (stable) where the integer windowed WPMM3
# needs m6^w. Both a synthetic symmetric-heavy control (estimator level) and the full
# vibration cascade (FORGE 78B-32, k=5 dynamics features) compare PATP3 vs Huber/WPMM3,
# and auto-valgate with PATP3 added to the candidate set. Drilling domain only.
# Requires the Run 2 sensor CSV extracted under /tmp/forge78b_sensor.
# Usage: Rscript experiments/run_patp3_cascade.R

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))
trmse <- function(e,p=0.90){ q<-stats::quantile(abs(e),p,names=FALSE,na.rm=TRUE); sqrt(mean(e[abs(e)<=q]^2)) }
rows <- list()

## (1) synthetic symmetric-heavy control (estimator level)
TH <- c(0.5,1,-0.8,0.4,-0.3,0.2)
sym_heavy <- function(n){ e<-rnorm(n); h<-runif(n)<0.08; k<-sum(h); if(k>0) e[h]<-rt(k,2)*3; e }
M <- c("LSE","Huber","L1","WPMM3","PATP3"); mat <- matrix(NA,40,length(M),dimnames=list(NULL,M)); g4 <- numeric(40)
for(s in 1:40){ set.seed(600+s); v1<-rnorm(300); v2<-rnorm(300); mu<-as.numeric(cbind(1,v1,v2,v1*v2,v1^2,v2^2)%*%TH); y<-mu+sym_heavy(300)
  gt1<-rnorm(4000); gt2<-rnorm(4000); gmu<-as.numeric(cbind(1,gt1,gt2,gt1*gt2,gt1^2,gt2^2)%*%TH)
  g4[s]<-sample_cumulants(resid(lm(y~v1+v2+I(v1*v2)+I(v1^2)+I(v2^2))))$gamma4
  for(m in M) mat[s,m]<-trmse(kg2_predict(tryCatch(inner_estimate(v1,v2,y,gmdh_pmm_control(B=0,force_method=m,max_iter=60L))$theta,error=function(e)rep(NA,6)),gt1,gt2)-gmu) }
md1<-apply(mat,2,median,na.rm=TRUE)
cat(sprintf("[1] synthetic symmetric-heavy (g4~%.0f): %s | PATP3 vs Huber %+.1f%%, vs WPMM3 %+.1f%%\n",
  median(g4), paste(sprintf("%s=%.4f",names(md1),md1),collapse=" "), 100*(md1["PATP3"]/md1["Huber"]-1), 100*(md1["PATP3"]/md1["WPMM3"]-1)))
for(m in M) rows[[length(rows)+1L]]<-data.frame(test="syn_sym_heavy",arm=m,trmse_med=md1[m],stringsAsFactors=FALSE)

## (2) full vibration cascade (k=5 dynamics features), out-of-sample
f <- list.files("/tmp/forge78b_sensor", pattern="Run 2 - output data.csv$", recursive=TRUE, full.names=TRUE)
if(length(f)){ d<-data.table::fread(f[1],skip=35L,header=TRUE,fill=TRUE,showProgress=FALSE)
  S<-"CSS-008"; g<-function(c) suppressWarnings(as.numeric(d[[paste0(S,"_",c)]]))
  PR<-c("StickSlip(%)","GyroXspread(RPM)","GyroXmed(RPM)","GyroXmax(RPM)","Temperature(C)")
  X<-sapply(PR,g); y<-log(g("ShZpeak(g)")); ok<-is.finite(y)&apply(is.finite(X),1,all); X<-X[ok,,drop=FALSE];y<-y[ok]
  n<-nrow(X); set.seed(1); if(n>3000){ii<-sort(sample.int(n,3000));X<-X[ii,];y<-y[ii];n<-3000}
  cu<-sample_cumulants(resid(lm(y~X)))
  ARMS<-c("Huber","WPMM3","PATP3","valgate+PATP3")
  ctrl<-function(a) if(a=="valgate+PATP3") gmdh_pmm_control(B=0,force_method="auto-valgate",F=6L,L_max=3L,criterion="MAE",valgate_candidates=c("LSE","Huber","L1","WPMM2","WPMM3","PATP3")) else gmdh_pmm_control(B=0,force_method=a,F=6L,L_max=3L,criterion="MAE")
  m2<-matrix(NA,18,length(ARMS),dimnames=list(NULL,ARMS))
  for(rr in 1:18){set.seed(97000+rr);tri<-sample.int(n,floor(.7*n));tei<-setdiff(1:n,tri)
    for(a in ARMS){fit<-tryCatch(gmdh_pmm(X[tri,,drop=FALSE],y[tri],ctrl(a)),error=function(e)NULL)
      m2[rr,a]<-if(is.null(fit))NA else trmse(y[tei]-tryCatch(predict(fit,X[tei,,drop=FALSE]),error=function(e)rep(NA,length(tei))))}}
  md2<-apply(m2,2,median,na.rm=TRUE)
  cat(sprintf("[2] vibration full cascade (k=5, resid g3=%+.2f g4=%.1f): %s\n", cu$gamma3, cu$gamma4, paste(sprintf("%s=%.4f",names(md2),md2),collapse=" ")))
  for(cmp in c("Huber","WPMM3")){p3<-m2[,"PATP3"];b<-m2[,cmp];o<-is.finite(p3)&is.finite(b)
    w<-suppressWarnings(stats::wilcox.test(p3[o],b[o],paired=TRUE)); cat(sprintf("    PATP3 vs %-6s: %+.1f%% win %2.0f%% p=%.2e\n",cmp,100*(median(p3)/median(b)-1),100*mean(b[o]>p3[o]),w$p.value))}
  for(a in ARMS) rows[[length(rows)+1L]]<-data.frame(test="vib_cascade_k5",arm=a,trmse_med=md2[a],stringsAsFactors=FALSE)
} else cat("[2] skipped (sensor CSV not extracted)\n")
outdir<-"../results"; if(!dir.exists(outdir)) dir.create(outdir,recursive=TRUE)
utils::write.csv(do.call(rbind,rows), file.path(outdir,"patp3_cascade.csv"), row.names=FALSE)
cat("Saved patp3_cascade.csv\n")
