#!/usr/bin/env Rscript
# EXPORTER for Experiment 1 (modern-ML baselines under identical blocked folds).
# Dumps, per vibration cell (12) and per Volve ROP band (4), the EXACT
# TIME/DEPTH-ORDERED (f1, f2, y) triples used by the honest blocked-CV harnesses
# (run_honest_blocked_cv.R / run_honest_blocked_cv_rop.R). The Python side reads
# these and replicates the deterministic blocked-5fold arithmetic (same gap/cap).
# No estimator is run here -- this only mirrors data loading + cell/band defs so
# the sklearn comparison is on byte-identical splits. Strictly drilling/vibration.
# Usage: Rscript experiments/export_blocked_cells_for_ml.R

.find_pkg <- function() { cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat { cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur); if (identical(parent, cur)) stop("gmdhpmm not found"); cur <- parent } }
PKG <- .find_pkg()
suppressMessages(if (requireNamespace("gmdhpmm", quietly = TRUE)) library(gmdhpmm) else pkgload::load_all(PKG, quiet = TRUE))

OUT <- "/tmp/p4_ml"
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
manifest <- list()

## ------------------------------------------------------------------ VIBRATION
# Mirror run_honest_blocked_cv.R: cap=6000 contiguous head, ok-filter preserves
# time order, f1=StickSlip(%), f2=GyroXspread(RPM), y=log(<chan>(g)).
cap <- 6000L; TMP <- "/tmp/forge78b_sensor"
.find1 <- function(pat) { x <- list.files(TMP, pattern = pat, recursive = TRUE, full.names = TRUE); if (length(x)) x[1] else NA }
RUNS <- list(`10.625` = .find1("Forge 10\\.625 - output data\\.csv$"),
             `14.75-Run1` = .find1("Run 1 - output data\\.csv$"),
             `14.75-Run2` = .find1("Run 2 - output data\\.csv$"))
SENSORS <- c("CSS-008", "CSS-007"); CHANS <- c("ShZpeak", "ShYpeak", "ShZrms")

cell <- function(d, sensor, chan) {
  col <- function(s) suppressWarnings(as.numeric(d[[paste0(sensor, "_", s)]]))
  ss <- col("StickSlip(%)"); gs <- col("GyroXspread(RPM)"); tg <- col(paste0(chan, "(g)"))
  if (is.null(ss) || is.null(gs) || is.null(tg)) return(NULL)
  ok <- is.finite(ss) & is.finite(gs) & is.finite(tg) & tg > 0       # preserves time order
  b <- ss[ok]; r <- gs[ok]; y <- log(tg[ok]); n <- length(y)
  if (n < 600) return(NULL)
  if (n > cap) { b <- b[1:cap]; r <- r[1:cap]; y <- y[1:cap]; n <- cap }  # contiguous head
  cu <- sample_cumulants(resid(lm(y ~ b + r)))
  list(f1 = b, f2 = r, y = y, n = n, g3 = cu$gamma3, g4 = cu$gamma4)
}

for (rn0 in names(RUNS)) { f <- RUNS[[rn0]]; if (is.na(f)) { cat(sprintf("  (run %s missing)\n", rn0)); next }
  d <- data.table::fread(f, skip = 35L, header = TRUE, fill = TRUE, showProgress = FALSE)
  for (sn in SENSORS) for (ch in CHANS) {
    res <- tryCatch(cell(d, sn, ch), error = function(e) NULL); if (is.null(res)) next
    name <- sprintf("vib__%s__%s__%s", rn0, sn, ch)
    fn <- file.path(OUT, paste0(name, ".csv"))
    utils::write.csv(data.frame(f1 = res$f1, f2 = res$f2, y = res$y), fn, row.names = FALSE)
    manifest[[length(manifest)+1L]] <- data.frame(dataset = "vibration", cell = name,
      file = fn, n = res$n, gamma4 = res$g4, gap = 50L, cap = 6000L, stringsAsFactors = FALSE)
    cat(sprintf("exported %-40s n=%5d g4=%6.2f\n", name, res$n, res$g4))
  }
}

## --------------------------------------------------------------------- VOLVE ROP
# Mirror run_honest_blocked_cv_rop.R: depth-ordered, same 4 bands, gap=30,
# f1=swob, f2=tqa, y=ROP.
dv <- tryCatch(load_volve_drilling(file.path(PKG, "data/volve_onbottom.csv")), error = function(e) NULL)
if (!is.null(dv)) {
  X <- as.data.frame(dv$X); y <- as.numeric(dv$y); ord <- order(X$dept); X <- X[ord, ]; y <- y[ord]
  BANDS <- list(volve_2000_2080=c(2000,2080), volve_1740_1820=c(1740,1820),
                volve_2110_2200=c(2110,2200), volve_mild_1360_1480=c(1360,1480))
  for (bn in names(BANDS)) { rg <- BANDS[[bn]]; sel <- X$dept >= rg[1] & X$dept < rg[2]
    Xi <- X[sel, ]; yi <- y[sel]; if (length(yi) < 150) next
    cu <- sample_cumulants(resid(lm(yi~Xi$swob+Xi$tqa)))
    name <- bn
    fn <- file.path(OUT, paste0("rop__", name, ".csv"))
    utils::write.csv(data.frame(f1 = Xi$swob, f2 = Xi$tqa, y = yi), fn, row.names = FALSE)
    manifest[[length(manifest)+1L]] <- data.frame(dataset = "volve_rop", cell = name,
      file = fn, n = length(yi), gamma4 = cu$gamma4, gap = 30L, cap = NA_integer_, stringsAsFactors = FALSE)
    cat(sprintf("exported %-40s n=%5d g4=%6.2f\n", paste0("rop__", name), length(yi), cu$gamma4))
  }
} else cat("(Volve load failed)\n")

man <- do.call(rbind, manifest)
utils::write.csv(man, file.path(OUT, "manifest.csv"), row.names = FALSE)
cat(sprintf("\nWrote %d cell/band files + manifest.csv to %s\n", nrow(man), OUT))
