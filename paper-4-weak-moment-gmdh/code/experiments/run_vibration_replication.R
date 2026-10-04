#!/usr/bin/env Rscript
# Drilling-domain replication of the impulsive-vibration win (FORGE 78B-32 Sanvean).
# Independent replications WITHIN the same well: 3 bit runs (10.625", 14.75 Run 1,
# 14.75 Run 2) x 2 sensors (CSS-008, CSS-007) x 3 shock channels (ShZpeak, ShYpeak,
# ShZrms). Each = log(shock) ~ KG-2(StickSlip, GyroXspread) of the same sensor,
# in-distribution (random 70/30 splits, test trimmed-RMSE). Reports whether
# WPMM2/WPMM3 beat Huber, and tallies replication across cells. Strictly drilling.
# Usage: Rscript experiments/run_vibration_replication.R [R]

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE); R <- if (length(args) >= 1) as.integer(args[1]) else 25L
trmse <- function(e,p=0.90){ q<-stats::quantile(abs(e),p,names=FALSE,na.rm=TRUE); sqrt(mean(e[abs(e)<=q]^2)) }
cap <- 4000L; TMP <- "/tmp/forge78b_sensor"

.find1 <- function(pat){ x <- list.files(TMP, pattern=pat, recursive=TRUE, full.names=TRUE); if(length(x)) x[1] else NA }
RUNS <- list(
  `10.625`     = .find1("Forge 10\\.625 - output data\\.csv$"),
  `14.75-Run1` = .find1("Run 1 - output data\\.csv$"),
  `14.75-Run2` = .find1("Run 2 - output data\\.csv$")
)
SENSORS <- c("CSS-008","CSS-007"); CHANS <- c("ShZpeak","ShYpeak","ShZrms")
fit_base <- function(m,b,r,yy) tryCatch(inner_estimate(b,r,yy, gmdh_pmm_control(B=0,force_method=m,max_iter=60L,weak_sigma_mult=2.5))$theta, error=function(e) rep(NA,6))
METH <- c("LSE","Huber","L1","WPMM2","WPMM3")

cmp_cell <- function(d, sensor, chan){
  col <- function(suffix) d[[paste0(sensor, "_", suffix)]]
  ss <- suppressWarnings(as.numeric(col("StickSlip(%)")))
  gs <- suppressWarnings(as.numeric(col("GyroXspread(RPM)")))
  tg <- suppressWarnings(as.numeric(col(paste0(chan,"(g)"))))
  if (is.null(ss) || is.null(gs) || is.null(tg)) return(NULL)
  ok <- is.finite(ss) & is.finite(gs) & is.finite(tg) & tg > 0
  b <- ss[ok]; r <- gs[ok]; y <- log(tg[ok]); n <- length(y)
  if (n < 300) return(NULL)
  if (n > cap){ ii<-sort(sample.int(n,cap)); b<-b[ii];r<-r[ii];y<-y[ii]; n<-cap }
  cu <- sample_cumulants(resid(lm(y~b+r)))
  mat <- matrix(NA,R,length(METH),dimnames=list(NULL,METH))
  for(rr in seq_len(R)){ set.seed(70000+rr); tri<-sample.int(n,floor(0.7*n)); tei<-setdiff(seq_len(n),tri)
    for(m in METH) mat[rr,m]<-trmse(y[tei]-kg2_predict(fit_base(m,b[tri],r[tri],y[tri]),b[tei],r[tei])) }
  md <- apply(mat,2,median,na.rm=TRUE)
  w3 <- mat[,"WPMM3"]; hu <- mat[,"Huber"]; okk <- is.finite(w3)&is.finite(hu)
  p <- if(sum(okk)>=5) suppressWarnings(stats::wilcox.test(w3[okk],hu[okk],paired=TRUE))$p.value else NA
  list(n=n, g3=cu$gamma3, g4=cu$gamma4, med=md,
       wpmm3_vs_huber=100*(md["WPMM3"]/md["Huber"]-1), wpmm2_vs_huber=100*(md["WPMM2"]/md["Huber"]-1),
       win=100*mean(hu[okk]>w3[okk]), p=p, best=names(which.min(md)))
}

cat(sprintf("=== FORGE 78B-32 vibration replication (drilling) | %d splits | 3 runs x 2 sensors x 3 channels ===\n", R))
cat(sprintf("%-10s %-7s %-8s %6s %6s | %8s %8s %8s %6s %5s %s\n","run","sensor","chan","g3","g4","WPMM2","WPMM3","Huber","vs Hu","win%","best"))
rows <- list(); win_cells <- 0L; tot_cells <- 0L
for(rn in names(RUNS)){ f <- RUNS[[rn]]; if(is.na(f)){ cat(sprintf("  (run %s missing)\n", rn)); next }
  d <- data.table::fread(f, skip=35L, header=TRUE, fill=TRUE, showProgress=FALSE)
  for(sn in SENSORS) for(ch in CHANS){
    res <- tryCatch(cmp_cell(d, sn, ch), error=function(e) NULL); if(is.null(res)) next
    tot_cells <- tot_cells + 1L; beat <- res$wpmm3_vs_huber < -2 || res$wpmm2_vs_huber < -2; if(beat) win_cells <- win_cells + 1L
    cat(sprintf("%-10s %-7s %-8s %+6.1f %6.1f | %8.4f %8.4f %8.4f %+5.1f%% %5.0f %s\n",
        rn, sn, ch, res$g3, res$g4, res$med["WPMM2"], res$med["WPMM3"], res$med["Huber"], res$wpmm3_vs_huber, res$win, res$best))
    rows[[length(rows)+1L]] <- data.frame(run=rn, sensor=sn, chan=ch, n=res$n, g3=res$g3, g4=res$g4,
      wpmm2=res$med["WPMM2"], wpmm3=res$med["WPMM3"], huber=res$med["Huber"], lse=res$med["LSE"],
      wpmm3_vs_huber=res$wpmm3_vs_huber, win=res$win, p=res$p, best=res$best, stringsAsFactors=FALSE)
  }
}
cat(sprintf("\nReplication: weak beats Huber (>2%%) in %d / %d cells\n", win_cells, tot_cells))
outdir <- "../results"; if(!dir.exists(outdir)) dir.create(outdir,recursive=TRUE)
utils::write.csv(do.call(rbind,rows), file.path(outdir,"vibration_replication.csv"), row.names=FALSE)
cat("Saved vibration_replication.csv\n")
