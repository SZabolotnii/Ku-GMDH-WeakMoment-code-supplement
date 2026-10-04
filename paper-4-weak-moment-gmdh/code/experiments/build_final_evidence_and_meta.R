#!/usr/bin/env Rscript
# Build the final Paper 4 claim evidence table and a dataset-level meta-summary.
#
# This script does not fit new models. It reads already rebuilt/frozen CSV files,
# maps each manuscript claim row to its exact source, and summarizes paired
# split-level deltas at the dataset-row level to avoid split pseudo-replication.
#
# Usage: Rscript experiments/build_final_evidence_and_meta.R [Bboot]

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[[1]], mustWork = TRUE)) else getwd()
}

.find_root <- function() {
  cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(cur, "paper-4-weak-moment-gmdh")) &&
        dir.exists(file.path(cur, "paper-1-gmdh-pmm"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("PMM-GMDH root not found")
    cur <- parent
  }
}

args <- commandArgs(trailingOnly = TRUE)
Bboot <- if (length(args) >= 1) as.integer(args[1]) else 10000L
set.seed(20260603L)

ROOT <- .find_root()
RES <- file.path(ROOT, "paper-4-weak-moment-gmdh", "results")

read_result <- function(name) {
  utils::read.csv(file.path(RES, name), stringsAsFactors = FALSE, check.names = FALSE)
}

soft_cmp <- read_result("p4_softsensor_honest_compare.csv")
soft_raw <- read_result("p4_softsensor_honest_raw.csv")
cross_cmp <- read_result("crosssec_modskew_comparison.csv")
cross_raw <- read_result("crosssec_modskew_raw.csv")
ins_cmp <- read_result("insurance_severity_family.csv")
ins_raw <- read_result("insurance_severity_raw.csv")
ml_sum <- read_result("ml_baselines_crossdomain.csv")
ml_long <- read_result("ml_baselines_crossdomain_long.csv")

claim_order <- c(
  "sru_y1_dynamic",
  "gas_turbine_co_2015_raw",
  "gas_turbine_nox_2015_raw",
  "islr_credit_balance",
  "insurance_autobi_loss",
  "fremtpl2_severity_raw",
  "sru_y2_static",
  "concrete"
)

claim_domain <- c(
  sru_y1_dynamic = "soft_sensor",
  gas_turbine_co_2015_raw = "emissions",
  gas_turbine_nox_2015_raw = "emissions",
  islr_credit_balance = "finance",
  insurance_autobi_loss = "insurance",
  fremtpl2_severity_raw = "insurance",
  sru_y2_static = "soft_sensor",
  concrete = "cross_section_control"
)

claim_source <- c(
  sru_y1_dynamic = "p4_softsensor_honest_compare.csv",
  gas_turbine_co_2015_raw = "p4_softsensor_honest_compare.csv",
  gas_turbine_nox_2015_raw = "p4_softsensor_honest_compare.csv",
  islr_credit_balance = "crosssec_modskew_comparison.csv",
  insurance_autobi_loss = "insurance_severity_family.csv",
  fremtpl2_severity_raw = "insurance_severity_family.csv",
  sru_y2_static = "p4_softsensor_honest_compare.csv",
  concrete = "crosssec_modskew_comparison.csv"
)

raw_source <- c(
  sru_y1_dynamic = "p4_softsensor_honest_raw.csv",
  gas_turbine_co_2015_raw = "p4_softsensor_honest_raw.csv",
  gas_turbine_nox_2015_raw = "p4_softsensor_honest_raw.csv",
  islr_credit_balance = "crosssec_modskew_raw.csv",
  insurance_autobi_loss = "insurance_severity_raw.csv",
  fremtpl2_severity_raw = "insurance_severity_raw.csv",
  sru_y2_static = "p4_softsensor_honest_raw.csv",
  concrete = "crosssec_modskew_raw.csv"
)

stat_lookup <- function(dataset, protocol) {
  if (dataset %in% soft_cmp$dataset) {
    r <- soft_cmp[soft_cmp$dataset == dataset & soft_cmp$protocol == protocol, , drop = FALSE]
    return(list(
      best_weak_method = r$best_pmm_trmse_method,
      best_robust_method = r$best_rob_trmse_method,
      best_weak_trmse = r$pmm_trmse,
      best_robust_trmse = r$rob_trmse,
      weak_vs_robust_pct = r$pmm_vs_rob_trmse_pct,
      lse_trmse = r$lse_trmse,
      weak_vs_lse_pct = r$pmm_vs_lse_trmse_pct
    ))
  }
  if (dataset %in% cross_cmp$candidate) {
    r <- cross_cmp[cross_cmp$candidate == dataset, , drop = FALSE]
    return(list(
      best_weak_method = r$best_pmm_method_tr,
      best_robust_method = r$best_robust_tr,
      best_weak_trmse = r$best_pmm_trmse,
      best_robust_trmse = r$best_rob_trmse,
      weak_vs_robust_pct = r$pmm_vs_robust_trmse_pct,
      lse_trmse = r$lse_trmse,
      weak_vs_lse_pct = r$pmm_vs_lse_trmse_pct
    ))
  }
  if (dataset %in% ins_cmp$candidate) {
    r <- ins_cmp[ins_cmp$candidate == dataset, , drop = FALSE]
    return(list(
      best_weak_method = r$best_pmm_tr,
      best_robust_method = r$best_robust_tr,
      best_weak_trmse = r$pmm_trmse,
      best_robust_trmse = r$robust_trmse,
      weak_vs_robust_pct = r$pmm_vs_robust_trmse_pct,
      lse_trmse = r$lse_trmse,
      weak_vs_lse_pct = r$pmm_vs_lse_pct
    ))
  }
  stop("No stat lookup for dataset: ", dataset)
}

raw_values <- function(dataset, protocol, method) {
  if (dataset %in% soft_raw$dataset) {
    r <- soft_raw[soft_raw$dataset == dataset & soft_raw$protocol == protocol &
                    soft_raw$method == method, c("fold", "trmse", "mae"), drop = FALSE]
    names(r)[1] <- "split_id"
    return(r)
  }
  if (dataset %in% cross_raw$candidate) {
    r <- cross_raw[cross_raw$candidate == dataset & cross_raw$method == method,
                   c("rep", "trmse", "mae"), drop = FALSE]
    names(r)[1] <- "split_id"
    return(r)
  }
  if (dataset %in% ins_raw$candidate) {
    r <- ins_raw[ins_raw$candidate == dataset & ins_raw$method == method,
                 c("rep", "trmse", "mae"), drop = FALSE]
    names(r)[1] <- "split_id"
    return(r)
  }
  stop("No raw values for dataset: ", dataset)
}

paired_deltas <- function(dataset, protocol, weak_method, robust_method, tree_method) {
  w <- raw_values(dataset, protocol, weak_method)
  rb <- raw_values(dataset, protocol, robust_method)
  names(w)[2:3] <- c("weak_trmse", "weak_mae")
  names(rb)[2:3] <- c("robust_trmse", "robust_mae")
  m <- merge(w, rb, by = "split_id", all = FALSE)
  tree <- ml_long[ml_long$dataset == dataset & ml_long$model == tree_method,
                  c("split_id", "trmse", "mae"), drop = FALSE]
  names(tree)[2:3] <- c("tree_trmse", "tree_mae")
  m <- merge(m, tree, by = "split_id", all = FALSE)
  m$dataset <- dataset
  m$protocol <- protocol
  m$weak_method <- weak_method
  m$robust_method <- robust_method
  m$tree_method <- tree_method
  m$weak_vs_robust_stat_pct <- 100 * (m$weak_trmse / m$robust_trmse - 1)
  m$weak_vs_tree_pct <- 100 * (m$weak_trmse / m$tree_trmse - 1)
  m
}

rows <- list()
delta_rows <- list()

for (dataset in claim_order) {
  ms <- ml_sum[ml_sum$dataset == dataset, , drop = FALSE]
  if (nrow(ms) != 1L) stop("Expected one ML summary row for ", dataset)
  st <- stat_lookup(dataset, ms$protocol)
  if (!identical(st$best_weak_method, ms$best_weak_method)) {
    stop("Best weak method mismatch for ", dataset, ": stat=", st$best_weak_method,
         " ml=", ms$best_weak_method)
  }
  if (!identical(st$best_robust_method, ms$best_robust_method)) {
    stop("Best robust method mismatch for ", dataset, ": stat=", st$best_robust_method,
         " ml=", ms$best_robust_method)
  }
  d <- paired_deltas(dataset, ms$protocol, st$best_weak_method,
                     st$best_robust_method, ms$best_tree_method)
  delta_rows[[length(delta_rows) + 1L]] <- d
  weak_vs_tree_pct <- ms$weak_vs_tree_trmse_pct
  claim_reading <- if (st$weak_vs_robust_pct < 0 && weak_vs_tree_pct < 0) {
    "positive_vs_robust_stat_and_tree"
  } else if (st$weak_vs_robust_pct < 0 && weak_vs_tree_pct >= 0) {
    "positive_vs_robust_stat_only"
  } else if (st$weak_vs_robust_pct >= 0 && weak_vs_tree_pct < 0) {
    "tree_positive_control_or_rescue"
  } else {
    "negative_or_loss_control"
  }
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = dataset,
    domain = unname(claim_domain[[dataset]]),
    role = ms$role,
    protocol = ms$protocol,
    n = ms$n,
    p = ms$p,
    splits_used = nrow(d),
    best_weak_method = st$best_weak_method,
    best_weak_trmse = st$best_weak_trmse,
    best_robust_stat_method = st$best_robust_method,
    best_robust_stat_trmse = st$best_robust_trmse,
    weak_vs_robust_stat_trmse_pct = st$weak_vs_robust_pct,
    lse_trmse = st$lse_trmse,
    weak_vs_lse_trmse_pct = st$weak_vs_lse_pct,
    best_tree_method = ms$best_tree_method,
    best_tree_trmse = ms$best_tree_trmse,
    weak_vs_tree_trmse_pct = weak_vs_tree_pct,
    canonical_stat_source = unname(claim_source[[dataset]]),
    canonical_raw_source = unname(raw_source[[dataset]]),
    canonical_ml_source = "ml_baselines_crossdomain.csv",
    canonical_ml_long_source = "ml_baselines_crossdomain_long.csv",
    claim_reading = claim_reading,
    stringsAsFactors = FALSE
  )
}

evidence <- do.call(rbind, rows)
deltas <- do.call(rbind, delta_rows)
deltas <- deltas[order(match(deltas$dataset, claim_order), deltas$split_id), ]

meta_one <- function(evidence, deltas, group_name, sel, comparison, B = Bboot) {
  ds <- evidence$dataset[sel]
  v <- sapply(ds, function(z) {
    col <- if (comparison == "weak_vs_robust_stat") {
      "weak_vs_robust_stat_pct"
    } else if (comparison == "weak_vs_tree") {
      "weak_vs_tree_pct"
    } else {
      stop("Unknown comparison: ", comparison)
    }
    stats::median(deltas[deltas$dataset == z, col], na.rm = TRUE)
  })
  names(v) <- ds
  n <- length(v)
  weak_wins <- sum(v < 0)
  sign_p <- stats::binom.test(weak_wins, n, 0.5)$p.value
  wilcox_p <- tryCatch(stats::wilcox.test(v, mu = 0, exact = FALSE)$p.value,
                       error = function(e) NA_real_)
  boot_mean <- replicate(B, mean(sample(v, n, replace = TRUE)))
  ci <- stats::quantile(boot_mean, c(0.025, 0.975), names = FALSE)
  data.frame(
    group = group_name,
    comparison = comparison,
    n_dataset_rows = n,
    weak_win_count = weak_wins,
    sign_test_p = sign_p,
    wilcox_p = wilcox_p,
    median_delta_pct = stats::median(v),
    mean_delta_pct = mean(v),
    bootstrap_mean_pct = mean(boot_mean),
    bootstrap_ci_lo_pct = ci[1],
    bootstrap_ci_hi_pct = ci[2],
    bootstrap_excludes_zero = (ci[1] < 0 && ci[2] < 0) || (ci[1] > 0 && ci[2] > 0),
    dataset_rows = paste(ds, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

groups <- list(
  all_claim_rows = rep(TRUE, nrow(evidence)),
  non_loss_rows = evidence$role != "loss_control",
  positive_rows = evidence$role %in% c("strong_positive", "modest_positive")
)

meta_rows <- list()
for (gn in names(groups)) {
  for (cmp in c("weak_vs_robust_stat", "weak_vs_tree")) {
    meta_rows[[length(meta_rows) + 1L]] <- meta_one(evidence, deltas, gn, groups[[gn]], cmp)
  }
}
meta <- do.call(rbind, meta_rows)

utils::write.csv(evidence, file.path(RES, "final_claim_evidence_table.csv"), row.names = FALSE)
utils::write.csv(deltas, file.path(RES, "final_claim_split_deltas.csv"), row.names = FALSE)
utils::write.csv(meta, file.path(RES, "final_claim_meta_analysis.csv"), row.names = FALSE)

cat("=== Paper 4 final claim evidence table ===\n")
print(evidence[, c("dataset", "role", "protocol", "best_weak_method",
                   "weak_vs_robust_stat_trmse_pct", "best_tree_method",
                   "weak_vs_tree_trmse_pct", "claim_reading")],
      row.names = FALSE)
cat("\n=== Dataset-level meta-summary (negative delta = weak/PMM better) ===\n")
print(meta[, c("group", "comparison", "n_dataset_rows", "weak_win_count",
               "median_delta_pct", "mean_delta_pct", "sign_test_p",
               "bootstrap_ci_lo_pct", "bootstrap_ci_hi_pct",
               "bootstrap_excludes_zero")],
      row.names = FALSE)
cat("\nSaved final_claim_evidence_table.csv, final_claim_split_deltas.csv, final_claim_meta_analysis.csv\n")
