#!/usr/bin/env Rscript
# Robustness layer for the eight final claim rows (Table 2), answering the
# 2026-10-04 PaperMentor review:
#   valgate    the deployable validation gate on every claim row, 4 tournament seeds
#   rinner     ablation: the same gate with random (interleaved) inner folds
#   norobust   ablation: the gate without Huber/LAD in its candidate library
#   sigma      window-width sensitivity, sigma_mult 1.5 and 4.0 (WPMM2, WPMM3, gate)
#   seeds      two extra tournament seeds for all eight Table-2 methods
#   gap        purge-gap sensitivity on the blocked rows, gap 50 and 100
#   repro      re-runs two archived cells per row to confirm the splits reproduce
# The outer splits and tournament seeds follow the producer scripts exactly
# (run_p4_softsensor_honest.R, run_p4_crosssec_modskew.R, run_insurance_severity.R),
# so every new fit is paired with the archived fits on the same split.
# Usage: Rscript experiments/run_claimrow_robustness.R [cores=9] [suffix=""]

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[[1]], mustWork = TRUE)) else getwd()
}
source(file.path(.script_dir(), "_p4_rebuilt_helpers.R"))

args <- commandArgs(trailingOnly = TRUE)
CORES <- if (length(args) >= 1) as.integer(args[[1]]) else 9L
suffix <- if (length(args) >= 2) args[[2]] else ""
ONLY_ROWS <- if (length(args) >= 3 && nzchar(args[[3]])) strsplit(args[[3]], ",")[[1]] else NULL
ONLY_ARMS <- if (length(args) >= 4 && nzchar(args[[4]])) strsplit(args[[4]], ",")[[1]] else NULL

source(file.path(.script_dir(), "_p4_rows.R"))   # ROWS, POOL, load_row()
VG_J <- length(P4_METHODS) + 1L          # seed slot of the gate, after the eight methods
SEED_STEP <- 1000000L                    # offset between tournament-seed replicates

ctrl <- function(method, seed, sigma = 2.5, cands = NULL, inner = "blocked") {
  gmdh_pmm_control(
    L_max = 3L, F = 6L, epsilon = -Inf, seed = seed,
    B = if (method == "auto-weak") 60L else 0L,
    force_method = method, criterion = "MSE", max_iter = 60L,
    weak_sigma_mult = sigma, valgate_folds = 4L, valgate_inner = inner,
    valgate_candidates = cands %||% c("LSE", "Huber", "L1", "WPMM2", "WPMM3")
  )
}
`%||%` <- function(a, b) if (is.null(a)) b else a

fit_one <- function(X, y, train, test, task) {
  xs <- p4_std_split(X, train, test)
  fit <- try(gmdh_pmm(xs$train, y[train], ctrl(task$method, task$seed, task$sigma,
                                               task$cands, task$inner)), silent = TRUE)
  pred <- if (inherits(fit, "try-error")) NULL else try(stats::predict(fit, xs$test), silent = TRUE)
  if (is.null(pred) || inherits(pred, "try-error")) {
    return(data.frame(trmse = NA_real_, mae = NA_real_, rmse = NA_real_, q90 = NA_real_, ok = FALSE))
  }
  e <- as.numeric(pred - y[test])
  data.frame(trmse = p4_trmse(e), mae = mean(abs(e)), rmse = p4_rmse(e),
             q90 = stats::quantile(abs(e), 0.9, names = FALSE), ok = all(is.finite(e)))
}

splits_for <- function(row, n, gap) {
  if (row$protocol == "random") {
    lapply(seq_len(30L), function(r) c(p4_random_split(n, row$seed0 + r), split = r,
                                       base = row$seed0 + r * 100L))
  } else {
    lapply(p4_blocked_splits(n, gap = gap), function(sp) c(sp, split = sp$fold,
                                                            base = row$seed0 + 5000L + sp$fold * 100L))
  }
}

# One task per (arm, method, seed replicate); a vector `s` expands into several tasks.
task <- function(arm, method, j, s, sigma = 2.5, cands = NULL, inner = "blocked", gap = 20L)
  lapply(s, function(si) list(arm = arm, method = method, j = j, s = si, sigma = sigma,
                              cands = cands, inner = inner, gap = gap))

out_path <- file.path(P4$results, paste0("claimrow_robustness_raw", suffix, ".csv"))
if (file.exists(out_path)) file.remove(out_path)
t0 <- Sys.time()
if (!is.null(ONLY_ROWS)) ROWS <- rbind(ROWS, POOL)[rbind(ROWS, POOL)$dataset %in% ONLY_ROWS, , drop = FALSE]
for (i in seq_len(nrow(ROWS))) {
  row <- ROWS[i, ]
  dat <- load_row(row$dataset)
  X <- dat$X; y <- dat$y
  if (row$protocol == "blocked" && length(y) > 3000L) { X <- X[1:3000, , drop = FALSE]; y <- y[1:3000] }
  arms <- c(task("valgate", "auto-valgate", VG_J, 0:3),
            task("rinner", "auto-valgate", VG_J, 0L, inner = "random"),
            task("norobust", "auto-valgate", VG_J, 0L, cands = c("LSE", "WPMM2", "WPMM3")))
  for (sg in c(1.5, 4.0)) for (m in c("WPMM2", "WPMM3", "auto-valgate"))
    arms <- c(arms, task("sigma", m, if (m == "auto-valgate") VG_J else match(m, P4_METHODS), 0L, sigma = sg))
  for (j in seq_along(P4_METHODS)) arms <- c(arms, task("seeds", P4_METHODS[[j]], j, 1:2))
  if (row$protocol == "blocked") for (g in c(50L, 100L)) {
    for (j in seq_along(P4_METHODS)) arms <- c(arms, task("gap", P4_METHODS[[j]], j, 0L, gap = g))
    arms <- c(arms, task("gap", "auto-valgate", VG_J, 0L, gap = g))
  }
  for (m in c("Huber", "WPMM2")) arms <- c(arms, task("repro", m, match(m, P4_METHODS), 0L))
  if (!is.null(ONLY_ARMS)) arms <- Filter(function(a) a$arm %in% ONLY_ARMS, arms)
  jobs <- list()
  for (a in arms) {
    sps <- splits_for(row, length(y), a$gap)
    if (a$arm == "repro") sps <- sps[1:2]
    for (sp in sps) jobs[[length(jobs) + 1L]] <- list(a = a, sp = sp)
  }
  cat(sprintf("[%s] %s n=%d p=%d: %d fits\n", row$dataset, row$protocol, length(y), ncol(X), length(jobs)))
  res <- parallel::mclapply(jobs, function(jb) {
    a <- jb$a; sp <- jb$sp
    a$seed <- sp$base + a$j + a$s * SEED_STEP
    ev <- fit_one(X, y, sp$train, sp$test, a)
    data.frame(dataset = row$dataset, protocol = row$protocol, split = sp$split,
               arm = a$arm, method = a$method, seed_rep = a$s, sigma = a$sigma,
               library = if (is.null(a$cands)) "full" else paste(a$cands, collapse = "+"),
               inner = a$inner, gap = if (row$protocol == "blocked") a$gap else NA_integer_,
               ev, stringsAsFactors = FALSE)
  }, mc.cores = CORES)
  bad <- vapply(res, function(r) inherits(r, "try-error") || !is.data.frame(r), logical(1))
  if (any(bad)) stop(sum(bad), " worker failures on ", row$dataset)
  res <- do.call(rbind, res)
  utils::write.table(res, out_path, sep = ",", row.names = FALSE, qmethod = "double",
                     col.names = !file.exists(out_path), append = file.exists(out_path))
  cat(sprintf("  done, %.1f min elapsed\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("Saved", out_path, "\n")
