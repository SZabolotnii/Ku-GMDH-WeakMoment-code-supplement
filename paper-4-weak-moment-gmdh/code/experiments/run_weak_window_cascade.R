#!/usr/bin/env Rscript
# H4: window-scale plateau at the CASCADE level. Sweep weak_sigma_mult for the
# weak-dispatch GMDH and show (a) a stability plateau where cascade MAE varies
# little, and (b) the MAD-tied default (sigma_mult = 2.5) lies in it. Run on the
# heavy regimes (onesided_contam, sym_t2) where the window scale matters.
# Usage (from paper-4-weak-moment-gmdh/code): Rscript experiments/run_weak_window_cascade.R [N_SEEDS] [B]

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm)
                 else pkgload::load_all(.find_pkg(), quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
N_SEEDS <- if (length(args) >= 1) as.integer(args[1]) else 30L
BDIAG   <- if (length(args) >= 2) as.integer(args[2]) else 50L
N_TRAIN <- 320L; N_TEST <- 3000L; SEED0 <- 80001L; D <- 6L
SIGMAS  <- c(1.5, 2.0, 2.5, 3.5, 5.0)

REG <- list(
  onesided_contam= function(n){e<-rnorm(n);h<-runif(n)<0.12;k<-sum(h);if(k>0)e[h]<--abs(stats::rt(k,2))*2.5;e},
  sym_t2         = function(n) stats::rt(n,2)/2
)
truth <- function(X){ x1<-X[,1];x2<-X[,2];x3<-X[,3];x4<-X[,4]
  z1 <- 0.5 + 0.8*x1 - 0.6*x2 + 0.4*x1*x2 - 0.2*x1^2
  z2 <- -0.3 + 0.7*x3 + 0.5*x4 - 0.3*x3*x4 + 0.2*x4^2
  0.2 + 1.0*z1 - 0.8*z2 + 0.3*z1*z2 }
gen_X <- function(n) matrix(rnorm(n*D), n, D, dimnames=list(NULL, paste0("x",1:D)))

run_one <- function(seed, errfn, sm) {
  set.seed(seed); Xtr<-gen_X(N_TRAIN); ytr<-truth(Xtr)+errfn(N_TRAIN); Xte<-gen_X(N_TEST); yte<-truth(Xte)
  ctrl <- gmdh_pmm_control(B=BDIAG, force_method="auto-weak", F=6L, L_max=3L,
                           weak_sigma_mult=sm, criterion="MAE", seed=seed)
  fit <- tryCatch(gmdh_pmm(Xtr, ytr, ctrl), error=function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  mean(abs(predict(fit, Xte) - yte))
}

cat(sprintf("=== Weak-dispatch cascade window-scale plateau (H4) | %d seeds | B=%d ===\n", N_SEEDS, BDIAG))
t0 <- Sys.time(); rows <- list()
for (rg in names(REG)) {
  cat(sprintf("\n--- %s : cascade MAE median by sigma_mult ---\n", rg))
  med <- numeric(length(SIGMAS))
  for (i in seq_along(SIGMAS)) {
    v <- sapply(seq_len(N_SEEDS), function(s) run_one(SEED0+s-1L, REG[[rg]], SIGMAS[i]))
    med[i] <- median(v, na.rm=TRUE)
    rows[[length(rows)+1L]] <- data.frame(regime=rg, sigma_mult=SIGMAS[i], mae_med=med[i], stringsAsFactors=FALSE)
  }
  best <- min(med); within10 <- SIGMAS[med <= 1.10*best]
  for (i in seq_along(SIGMAS))
    cat(sprintf("  sigma=%.1f : MAE=%.4f%s\n", SIGMAS[i], med[i], if (SIGMAS[i]==2.5) "  <- MAD-tied default" else ""))
  cat(sprintf("  plateau (within +10%% of best): sigma in [%.1f, %.1f]; default 2.5 in plateau: %s\n",
      min(within10), max(within10), 2.5 %in% within10))
}
outdir <- "../results"; if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
utils::write.csv(do.call(rbind, rows), file.path(outdir, "weak_window_cascade.csv"), row.names=FALSE)
cat(sprintf("\nSaved weak_window_cascade.csv | %.0fs\n", as.numeric(difftime(Sys.time(), t0, "secs"))))
