#!/usr/bin/env Rscript
# Export all fifteen candidate-pool rows (the eight Table-2 rows and the seven left out)
# with the exact outer splits of the R generators, for ml_tuned_trees.py.
# Same file format as export_crossdomain_claim_rows_for_ml.R.
# Usage: Rscript experiments/export_allrows_for_ml.R [out_dir=/tmp/p4_allrows_ml]

.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) dirname(normalizePath(f[[1]], mustWork = TRUE)) else getwd()
}
source(file.path(.script_dir(), "_p4_rebuilt_helpers.R"))
source(file.path(.script_dir(), "_p4_rows.R"))

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) >= 1) args[[1]] else "/tmp/p4_allrows_ml"
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)

ev <- utils::read.csv(file.path(P4$results, "final_claim_evidence_table.csv"), stringsAsFactors = FALSE)
all_rows <- rbind(ROWS, POOL)
manifest <- list(); splits <- list()
for (i in seq_len(nrow(all_rows))) {
  row <- all_rows[i, ]
  dat <- load_row(row$dataset)
  X <- as.matrix(dat$X); y <- as.numeric(dat$y)
  if (row$protocol == "blocked" && length(y) > 3000L) { X <- X[1:3000, , drop = FALSE]; y <- y[1:3000] }
  colnames(X) <- paste0("x", seq_len(ncol(X)))
  data_file <- file.path(OUT, paste0(row$dataset, ".csv"))
  utils::write.csv(data.frame(y = y, X, check.names = FALSE), data_file, row.names = FALSE)
  role <- ev$role[match(row$dataset, ev$dataset)]
  manifest[[i]] <- data.frame(dataset = row$dataset, protocol = row$protocol,
                              role = if (is.na(role)) "pool" else role, file = data_file,
                              n = length(y), p = ncol(X), stringsAsFactors = FALSE)
  sps <- if (row$protocol == "blocked") {
    lapply(p4_blocked_splits(length(y), gap = 20L), function(s) list(id = s$fold, train = s$train, test = s$test))
  } else {
    lapply(seq_len(30L), function(r) { s <- p4_random_split(length(y), row$seed0 + r)
      list(id = r, train = s$train, test = s$test) })
  }
  for (s in sps) splits[[length(splits) + 1L]] <- data.frame(
    dataset = row$dataset, split_id = s$id, protocol = row$protocol,
    train_idx = paste(s$train - 1L, collapse = " "), test_idx = paste(s$test - 1L, collapse = " "),
    stringsAsFactors = FALSE)
  cat(sprintf("exported %-28s %-7s n=%5d p=%2d\n", row$dataset, row$protocol, length(y), ncol(X)))
}
utils::write.csv(do.call(rbind, manifest), file.path(OUT, "manifest.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, splits), file.path(OUT, "splits.csv"), row.names = FALSE)
cat("Wrote", nrow(all_rows), "datasets and", length(splits), "splits to", OUT, "\n")
