#!/usr/bin/env Rscript
# Export final Paper 4 cross-domain claim rows for a uniform sklearn ML baseline.
#
# The exporter writes model-matrix X, response y, and the exact train/test splits
# used by the rebuilt R generators. Python then fits RF / HistGB / ExtraTrees on
# the same rows and validation discipline.
#
# Usage: Rscript experiments/export_crossdomain_claim_rows_for_ml.R

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[[1]], mustWork = TRUE)) else getwd()
}
source(file.path(.script_dir(), "_p4_rebuilt_helpers.R"))

OUT <- "/tmp/p4_crossdomain_ml"
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

load_claim_dataset <- function(id) {
  if (id %in% c("sru_y1_dynamic", "sru_y2_static",
                "gas_turbine_co_2015_raw", "gas_turbine_nox_2015_raw",
                "concrete", "fremtpl2_severity_raw")) {
    return(p4_ext_load(id))
  }
  if (id == "islr_credit_balance") {
    data(Credit, package = "ISLR2")
    d <- stats::na.omit(ISLR2::Credit)
    return(p4_frame_to_xy(d, Balance ~ Income + Limit + Rating + Cards + Age + Education +
                            Own + Student + Married + Region))
  }
  if (id == "insurance_autobi_loss") {
    data(AutoBi, package = "insuranceData")
    d <- stats::na.omit(get("AutoBi"))
    d$ATTORNEY <- factor(d$ATTORNEY)
    d$CLMSEX <- factor(d$CLMSEX)
    d$MARITAL <- factor(d$MARITAL)
    d$CLMINSUR <- factor(d$CLMINSUR)
    d$SEATBELT <- factor(d$SEATBELT)
    return(p4_frame_to_xy(d, LOSS ~ ATTORNEY + CLMSEX + MARITAL + CLMINSUR + SEATBELT + CLMAGE))
  }
  stop("Unknown claim dataset: ", id)
}

claim_rows <- data.frame(
  dataset = c("sru_y1_dynamic", "gas_turbine_co_2015_raw", "gas_turbine_nox_2015_raw",
              "islr_credit_balance", "insurance_autobi_loss", "fremtpl2_severity_raw",
              "sru_y2_static", "concrete"),
  protocol = c("blocked", "blocked", "blocked", "random", "random", "random", "blocked", "random"),
  role = c("strong_positive", "strong_positive", "strong_positive",
           "modest_positive", "modest_positive", "rescue_tie",
           "tie_control", "loss_control"),
  seed0 = c(76000L, 76000L, 76000L, 73000L, 74000L, 74000L, 76000L, 73000L),
  reps = c(NA_integer_, NA_integer_, NA_integer_, 30L, 30L, 30L, NA_integer_, 30L),
  folds = c(5L, 5L, 5L, NA_integer_, NA_integer_, NA_integer_, 5L, NA_integer_),
  gap = c(20L, 20L, 20L, NA_integer_, NA_integer_, NA_integer_, 20L, NA_integer_),
  source_csv = c("p4_softsensor_honest_compare.csv", "p4_softsensor_honest_compare.csv",
                 "p4_softsensor_honest_compare.csv", "crosssec_modskew_comparison.csv",
                 "insurance_severity_family.csv", "insurance_severity_family.csv",
                 "p4_softsensor_honest_compare.csv", "crosssec_modskew_comparison.csv"),
  stringsAsFactors = FALSE
)

manifest <- list()
splits <- list()

for (i in seq_len(nrow(claim_rows))) {
  row <- claim_rows[i, ]
  dat <- load_claim_dataset(row$dataset)
  X <- as.matrix(dat$X)
  y <- as.numeric(dat$y)
  if (row$protocol == "blocked" && length(y) > 3000L) {
    X <- X[seq_len(3000L), , drop = FALSE]
    y <- y[seq_len(3000L)]
  }
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  data_file <- file.path(OUT, paste0(row$dataset, ".csv"))
  utils::write.csv(data.frame(y = y, X, check.names = FALSE), data_file, row.names = FALSE)

  manifest[[length(manifest) + 1L]] <- data.frame(
    dataset = row$dataset, protocol = row$protocol, role = row$role,
    file = data_file, n = length(y), p = ncol(X), seed0 = row$seed0,
    reps = row$reps, folds = row$folds, gap = row$gap, source_csv = row$source_csv,
    stringsAsFactors = FALSE
  )

  if (row$protocol == "blocked") {
    sp <- p4_blocked_splits(length(y), folds = row$folds, gap = row$gap)
    for (j in seq_along(sp)) {
      splits[[length(splits) + 1L]] <- data.frame(
        dataset = row$dataset, split_id = sp[[j]]$fold, protocol = row$protocol,
        train_idx = paste(sp[[j]]$train - 1L, collapse = " "),
        test_idx = paste(sp[[j]]$test - 1L, collapse = " "),
        stringsAsFactors = FALSE
      )
    }
  } else {
    for (rep in seq_len(row$reps)) {
      sp <- p4_random_split(length(y), seed = row$seed0 + rep)
      splits[[length(splits) + 1L]] <- data.frame(
        dataset = row$dataset, split_id = rep, protocol = row$protocol,
        train_idx = paste(sp$train - 1L, collapse = " "),
        test_idx = paste(sp$test - 1L, collapse = " "),
        stringsAsFactors = FALSE
      )
    }
  }
  cat(sprintf("exported %-28s protocol=%-7s n=%4d p=%2d role=%s\n",
              row$dataset, row$protocol, length(y), ncol(X), row$role))
}

utils::write.csv(do.call(rbind, manifest), file.path(OUT, "manifest.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, splits), file.path(OUT, "splits.csv"), row.names = FALSE)
cat(sprintf("\nWrote %d datasets and %d splits to %s\n", length(manifest), length(splits), OUT))
