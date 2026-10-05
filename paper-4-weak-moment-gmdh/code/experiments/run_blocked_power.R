#!/usr/bin/env Rscript
# Blocked-row power check (pre-registered in blocked-power-prereg-2026-10-04.md):
# the four blocked claim rows with 10 contiguous folds, purge gaps 20 and 50, and five
# tournament seeds per fold, for the eight Table-2 methods and the validation gate.
# Usage: Rscript experiments/run_blocked_power.R [cores=9] [suffix=""]

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[[1]], mustWork = TRUE)) else getwd()
}
source(file.path(.script_dir(), "_p4_rebuilt_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
CORES <- if (length(args) >= 1) as.integer(args[[1]]) else 9L
suffix <- if (length(args) >= 2) args[[2]] else ""

ROWS <- c("sru_y1_dynamic", "gas_turbine_co_2015_raw", "gas_turbine_nox_2015_raw", "sru_y2_static")
FOLDS <- 10L
GAPS <- c(20L, 50L)
SEEDS <- 0:4
SEED0 <- 77000L                         # fresh seed range, distinct from the 5-fold archive
SEED_STEP <- 1000000L
METHODS <- c(P4_METHODS, "auto-valgate")

out_path <- file.path(P4$results, paste0("blocked_power_raw", suffix, ".csv"))
if (file.exists(out_path)) file.remove(out_path)
t0 <- Sys.time()
for (id in ROWS) {
  dat <- p4_ext_load(id)
  n <- min(length(dat$y), 3000L)
  X <- dat$X[seq_len(n), , drop = FALSE]; y <- dat$y[seq_len(n)]
  jobs <- list()
  for (g in GAPS) for (sp in p4_blocked_splits(n, folds = FOLDS, gap = g))
    for (j in seq_along(METHODS)) for (s in SEEDS)
      jobs[[length(jobs) + 1L]] <- list(gap = g, sp = sp, j = j, s = s)
  cat(sprintf("[%s] n=%d p=%d: %d fits\n", id, n, ncol(X), length(jobs)))
  res <- parallel::mclapply(jobs, function(jb) {
    seed <- SEED0 + jb$sp$fold * 100L + jb$j + jb$s * SEED_STEP
    ev <- p4_fit_eval(X, y, jb$sp$train, jb$sp$test, METHODS[[jb$j]], seed = seed)
    data.frame(dataset = id, gap = jb$gap, fold = jb$sp$fold, method = METHODS[[jb$j]],
               seed_rep = jb$s, trmse = ev$trmse, mae = ev$mae, rmse = ev$rmse, ok = ev$ok)
  }, mc.cores = CORES)
  bad <- vapply(res, function(r) !is.data.frame(r), logical(1))
  if (any(bad)) stop(sum(bad), " worker failures on ", id)
  res <- do.call(rbind, res)
  utils::write.table(res, out_path, sep = ",", row.names = FALSE, qmethod = "double",
                     col.names = !file.exists(out_path), append = file.exists(out_path))
  cat(sprintf("  done, %.1f min elapsed\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("Saved", out_path, "\n")
