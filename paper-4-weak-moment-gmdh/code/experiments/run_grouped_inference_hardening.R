#!/usr/bin/env Rscript
# Paper 4 pre-submission grouped-inference hardening.
#
# This script does not refit models and does not change the frozen claim table.
# It reads the final claim-row evidence and split-level paired deltas, collapses
# splits to dataset rows, and then checks how fragile the programme-level claim
# is under protocol/domain/role grouping. The purpose is editorial hygiene:
# avoid using split rows as independent evidence and avoid a broad claim that
# disappears after domain clustering.
#
# Usage:
#   Rscript paper-4-weak-moment-gmdh/code/experiments/run_grouped_inference_hardening.R [Bboot]

.find_root <- function() {
  cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(cur, "paper-4-weak-moment-gmdh")) &&
        dir.exists(file.path(cur, "paper-1-gmdh-pmm"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("PMM-GMDH root not found")
    cur <- parent
  }
}

args <- commandArgs(trailingOnly = TRUE)
Bboot <- if (length(args) >= 1) as.integer(args[[1]]) else 10000L
set.seed(20260603L)

ROOT <- .find_root()
P4 <- file.path(ROOT, "paper-4-weak-moment-gmdh")
RES <- file.path(P4, "results")

read_result <- function(name) {
  utils::read.csv(file.path(RES, name), stringsAsFactors = FALSE, check.names = FALSE)
}

evidence <- read_result("final_claim_evidence_table.csv")
deltas <- read_result("final_claim_split_deltas.csv")

comparisons <- c("weak_vs_robust_stat", "weak_vs_tree")
delta_col <- c(
  weak_vs_robust_stat = "weak_vs_robust_stat_pct",
  weak_vs_tree = "weak_vs_tree_pct"
)

dataset_delta <- function(cmp) {
  col <- unname(delta_col[[cmp]])
  out <- lapply(split(deltas, deltas$dataset), function(d) {
    data.frame(
      dataset = d$dataset[[1]],
      dataset_delta_pct = stats::median(d[[col]], na.rm = TRUE),
      split_count = sum(is.finite(d[[col]])),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  merge(
    evidence[, c("dataset", "domain", "role", "protocol", "claim_reading")],
    out,
    by = "dataset",
    all.x = TRUE,
    sort = FALSE
  )
}

scope_rows <- function(df, scope) {
  if (scope == "all_claim_rows") return(rep(TRUE, nrow(df)))
  if (scope == "non_loss_rows") return(df$role != "loss_control")
  if (scope == "positive_rows") return(df$role %in% c("strong_positive", "modest_positive"))
  if (scope == "blocked_rows") return(df$protocol == "blocked")
  if (scope == "random_rows") return(df$protocol == "random")
  stop("Unknown scope: ", scope)
}

unit_values <- function(df, unit_col) {
  sp <- split(df$dataset_delta_pct, df[[unit_col]])
  rows <- lapply(names(sp), function(unit) {
    v <- sp[[unit]]
    data.frame(
      unit = unit,
      n_dataset_rows = length(v),
      unit_delta_pct = stats::median(v, na.rm = TRUE),
      unit_mean_delta_pct = mean(v, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

boot_ci_mean <- function(v, B) {
  n <- length(v)
  if (n == 0L || any(!is.finite(v))) return(c(NA_real_, NA_real_, NA_real_))
  if (n == 1L) return(c(mean(v), NA_real_, NA_real_))
  boot <- replicate(B, mean(sample(v, n, replace = TRUE)))
  qs <- stats::quantile(boot, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  c(mean(boot), qs)
}

summarize_units <- function(uv, cmp, scope, unit_col) {
  v <- uv$unit_delta_pct
  n <- length(v)
  wins <- sum(v < 0, na.rm = TRUE)
  sign_p <- if (n > 0) stats::binom.test(wins, n, 0.5)$p.value else NA_real_
  wilcox_p <- if (n >= 2) {
    tryCatch(stats::wilcox.test(v, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
  } else {
    NA_real_
  }
  ci <- boot_ci_mean(v, Bboot)
  data.frame(
    scope = scope,
    comparison = cmp,
    unit = unit_col,
    n_units = n,
    weak_win_units = wins,
    sign_test_p = sign_p,
    wilcox_p = wilcox_p,
    median_unit_delta_pct = stats::median(v, na.rm = TRUE),
    mean_unit_delta_pct = mean(v, na.rm = TRUE),
    bootstrap_mean_pct = ci[[1]],
    bootstrap_ci_lo_pct = ci[[2]],
    bootstrap_ci_hi_pct = ci[[3]],
    bootstrap_excludes_zero = is.finite(ci[[2]]) && is.finite(ci[[3]]) &&
      ((ci[[2]] < 0 && ci[[3]] < 0) || (ci[[2]] > 0 && ci[[3]] > 0)),
    unit_values = paste(sprintf("%s:%+.3f", uv$unit, uv$unit_delta_pct), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

scopes <- c("all_claim_rows", "non_loss_rows", "positive_rows", "blocked_rows", "random_rows")
unit_cols <- c("dataset", "domain", "protocol", "role")

summary_rows <- list()
unit_detail_rows <- list()

for (cmp in comparisons) {
  dd <- dataset_delta(cmp)
  for (scope in scopes) {
    ds <- dd[scope_rows(dd, scope), , drop = FALSE]
    if (nrow(ds) == 0L) next
    for (unit_col in unit_cols) {
      uv <- unit_values(ds, unit_col)
      uv$scope <- scope
      uv$comparison <- cmp
      uv$unit_type <- unit_col
      unit_detail_rows[[length(unit_detail_rows) + 1L]] <- uv[
        , c("scope", "comparison", "unit_type", "unit", "n_dataset_rows",
            "unit_delta_pct", "unit_mean_delta_pct")
      ]
      summary_rows[[length(summary_rows) + 1L]] <- summarize_units(uv, cmp, scope, unit_col)
    }
  }
}

summary <- do.call(rbind, summary_rows)
unit_details <- do.call(rbind, unit_detail_rows)

protocol_rows <- list()
for (cmp in comparisons) {
  dd <- dataset_delta(cmp)
  for (protocol in sort(unique(dd$protocol))) {
    sub <- dd[dd$protocol == protocol, , drop = FALSE]
    v <- sub$dataset_delta_pct
    ci <- boot_ci_mean(v, Bboot)
    protocol_rows[[length(protocol_rows) + 1L]] <- data.frame(
      comparison = cmp,
      protocol = protocol,
      n_dataset_rows = nrow(sub),
      weak_win_rows = sum(v < 0),
      median_delta_pct = stats::median(v),
      mean_delta_pct = mean(v),
      sign_test_p = stats::binom.test(sum(v < 0), length(v), 0.5)$p.value,
      bootstrap_mean_pct = ci[[1]],
      bootstrap_ci_lo_pct = ci[[2]],
      bootstrap_ci_hi_pct = ci[[3]],
      bootstrap_excludes_zero = is.finite(ci[[2]]) && is.finite(ci[[3]]) &&
        ((ci[[2]] < 0 && ci[[3]] < 0) || (ci[[2]] > 0 && ci[[3]] > 0)),
      dataset_rows = paste(sub$dataset, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}
protocol_summary <- do.call(rbind, protocol_rows)

loo_rows <- list()
for (cmp in comparisons) {
  dd <- dataset_delta(cmp)
  domains <- sort(unique(dd$domain))
  for (omit in domains) {
    sub <- dd[dd$domain != omit, , drop = FALSE]
    v <- sub$dataset_delta_pct
    loo_rows[[length(loo_rows) + 1L]] <- data.frame(
      comparison = cmp,
      omitted_domain = omit,
      n_dataset_rows = nrow(sub),
      weak_win_rows = sum(v < 0),
      median_delta_pct = stats::median(v),
      mean_delta_pct = mean(v),
      dataset_rows = paste(sub$dataset, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}
loo <- do.call(rbind, loo_rows)

utils::write.csv(summary, file.path(RES, "grouped_inference_summary.csv"), row.names = FALSE)
utils::write.csv(unit_details, file.path(RES, "grouped_inference_unit_details.csv"), row.names = FALSE)
utils::write.csv(protocol_summary, file.path(RES, "grouped_inference_protocol_summary.csv"), row.names = FALSE)
utils::write.csv(loo, file.path(RES, "grouped_inference_leave_one_domain.csv"), row.names = FALSE)

cat("=== Paper 4 grouped-inference hardening ===\n")
cat("Negative delta means weak/PMM has lower trimmed-RMSE.\n\n")
show_cols <- c("scope", "comparison", "unit", "n_units", "weak_win_units",
               "median_unit_delta_pct", "mean_unit_delta_pct", "sign_test_p",
               "bootstrap_ci_lo_pct", "bootstrap_ci_hi_pct", "bootstrap_excludes_zero")
print(summary[summary$scope %in% c("all_claim_rows", "positive_rows") &
                summary$unit %in% c("dataset", "domain", "protocol"),
              show_cols],
      row.names = FALSE)
cat("\nSaved grouped_inference_*.csv in ", RES, "\n", sep = "")
