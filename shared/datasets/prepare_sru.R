#!/usr/bin/env Rscript
# Prepare normalized MIMO sulfur recovery unit (SRU) soft-sensor data for C5 external screening.
#
# Source:
#   Mendeley Data: Data for: Echo-state networks for Soft Sensor design in an SRU process
#   DOI: 10.17632/kcpnnrn67p.1, CC BY 4.0
#   CSV mirror:
#   https://raw.githubusercontent.com/softsensors/soft-sensor-data/main/SRU.csv
#
# Run from repository root:
#   Rscript shared/datasets/prepare_sru.R

find_root <- function() {
  here <- normalizePath(getwd(), mustWork = TRUE)
  cur <- here
  repeat {
    if (file.exists(file.path(cur, "ROADMAP.md")) &&
        file.exists(file.path(cur, "paper-1-gmdh-pmm", "code", "DESCRIPTION"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("Cannot locate repository root from ", here)
    cur <- parent
  }
}

root <- find_root()
ext_dir <- file.path(root, "shared", "datasets", "external")
raw_dir <- file.path(ext_dir, "raw")
processed_dir <- file.path(ext_dir, "processed")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

url <- "https://raw.githubusercontent.com/softsensors/soft-sensor-data/main/SRU.csv"
raw_path <- file.path(raw_dir, "SRU.csv")
needs_download <- !file.exists(raw_path) || file.info(raw_path)$size <= 0
if (needs_download) {
  message("Downloading SRU.csv mirror...")
  utils::download.file(url, raw_path, mode = "wb", quiet = FALSE)
}

d <- utils::read.csv(raw_path, stringsAsFactors = FALSE, check.names = FALSE,
                     fileEncoding = "UTF-8-BOM")
names(d) <- sub("^\ufeff", "", enc2utf8(names(d)))
expected <- c("time", paste0("u", 1:5), "y1", "y2")
missing <- setdiff(expected, names(d))
if (length(missing)) stop("Missing expected columns: ", paste(missing, collapse = ", "))

d <- d[, expected]
for (nm in expected) d[[nm]] <- as.numeric(d[[nm]])
d <- d[stats::complete.cases(d), ]
d <- d[order(d$time), ]

lag_vec <- function(x, k) c(rep(NA_real_, k), head(x, -k))
lagged <- d
for (target in c("y1", "y2")) {
  lagged[[paste0(target, "_lag1")]] <- lag_vec(lagged[[target]], 1L)
  lagged[[paste0(target, "_lag2")]] <- lag_vec(lagged[[target]], 2L)
}
for (input in paste0("u", 1:5)) {
  lagged[[paste0(input, "_lag1")]] <- lag_vec(lagged[[input]], 1L)
}
lagged <- lagged[stats::complete.cases(lagged), ]

static_path <- file.path(processed_dir, "sru.csv")
lagged_path <- file.path(processed_dir, "sru_lagged.csv")
utils::write.csv(d, static_path, row.names = FALSE)
utils::write.csv(lagged, lagged_path, row.names = FALSE)

manifest_path <- file.path(ext_dir, "external_candidates.csv")
manifest <- if (file.exists(manifest_path)) {
  utils::read.csv(manifest_path, stringsAsFactors = FALSE)
} else {
  data.frame(id = character(), label = character(), path = character(),
             formula = character(), source = character(), license = character(),
             note = character(), stringsAsFactors = FALSE)
}

required_cols <- c("id", "label", "path", "formula", "source", "license", "note")
for (nm in setdiff(required_cols, names(manifest))) manifest[[nm]] <- character(nrow(manifest))
manifest <- manifest[, required_cols, drop = FALSE]

base_formula <- "TARGET ~ u1 + u2 + u3 + u4 + u5"
dynamic_formula <- paste(
  "TARGET ~ u1 + u2 + u3 + u4 + u5 +",
  "y1_lag1 + y1_lag2 + y2_lag1 + y2_lag2 +",
  "u1_lag1 + u2_lag1 + u3_lag1 + u4_lag1 + u5_lag1"
)
source_text <- paste(
  "Mendeley Data DOI 10.17632/kcpnnrn67p.1;",
  "GitHub mirror softsensors/soft-sensor-data/SRU.csv"
)
license_text <- "CC BY 4.0 via Mendeley Data record; GitHub mirror has no explicit repository license"

new_rows <- data.frame(
  id = c("sru_y1_static", "sru_y2_static", "sru_y1_dynamic", "sru_y2_dynamic"),
  label = c(
    "SRU soft sensor: y1 ~ current inputs",
    "SRU soft sensor: y2 ~ current inputs",
    "SRU soft sensor: y1 ~ current inputs + output/input lags",
    "SRU soft sensor: y2 ~ current inputs + output/input lags"
  ),
  path = c(
    "shared/datasets/external/processed/sru.csv",
    "shared/datasets/external/processed/sru.csv",
    "shared/datasets/external/processed/sru_lagged.csv",
    "shared/datasets/external/processed/sru_lagged.csv"
  ),
  formula = c(
    sub("TARGET", "y1", base_formula),
    sub("TARGET", "y2", base_formula),
    sub("TARGET", "y1", dynamic_formula),
    sub("TARGET", "y2", dynamic_formula)
  ),
  source = source_text,
  license = license_text,
  note = c(
    "Normalized MIMO SRU soft-sensor data; static intake for DeepResearch backup candidate.",
    "Normalized MIMO SRU soft-sensor data; static intake for DeepResearch backup candidate.",
    "Lagged dynamic soft-sensor intake; lags respect row order within the single SRU sequence.",
    "Lagged dynamic soft-sensor intake; lags respect row order within the single SRU sequence."
  ),
  stringsAsFactors = FALSE
)

manifest <- manifest[!(manifest$id %in% new_rows$id), , drop = FALSE]
manifest <- rbind(manifest, new_rows)
utils::write.csv(manifest, manifest_path, row.names = FALSE)

message("Wrote:")
message("  ", static_path, " (n=", nrow(d), ")")
message("  ", lagged_path, " (n=", nrow(lagged), ")")
message("Updated manifest: ", manifest_path)
