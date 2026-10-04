#!/usr/bin/env Rscript
# Prepare UCI Gas Turbine CO/NOx emission data for C5 external screening.
#
# Source:
#   Gas Turbine CO and NOx Emission Data Set, UCI ML Repository
#   DOI: 10.24432/C5WC95, CC BY 4.0
#
# Run from repository root:
#   Rscript shared/datasets/prepare_gasturbine.R

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

url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/00551/pp_gas_emission.zip"
zip_path <- file.path(raw_dir, "pp_gas_emission.zip")
if (!file.exists(zip_path)) {
  message("Downloading UCI Gas Turbine archive...")
  utils::download.file(url, zip_path, mode = "wb", quiet = FALSE)
}

files <- utils::unzip(zip_path, list = TRUE)
csv_files <- files$Name[grepl("^gt_20[0-9]{2}\\.csv$", files$Name)]
if (!length(csv_files)) stop("No gt_YYYY.csv files found in archive: ", zip_path)

read_year <- function(name) {
  d <- utils::read.csv(unz(zip_path, name), stringsAsFactors = FALSE)
  year <- as.integer(sub("^gt_([0-9]{4})\\.csv$", "\\1", basename(name)))
  d$year <- year
  d
}
all_data <- do.call(rbind, lapply(csv_files, read_year))
all_data <- all_data[stats::complete.cases(all_data), ]

expected <- c("AT", "AP", "AH", "AFDP", "GTEP", "TIT", "TAT", "TEY", "CDP", "CO", "NOX", "year")
missing <- setdiff(expected, names(all_data))
if (length(missing)) stop("Missing expected columns: ", paste(missing, collapse = ", "))

all_path <- file.path(processed_dir, "gas_turbine_2011_2015.csv")
csv_2015_path <- file.path(processed_dir, "gas_turbine_2015.csv")
utils::write.csv(all_data, all_path, row.names = FALSE)
utils::write.csv(all_data[all_data$year == 2015, ], csv_2015_path, row.names = FALSE)

manifest_path <- file.path(ext_dir, "external_candidates.csv")
manifest <- if (file.exists(manifest_path)) {
  utils::read.csv(manifest_path, stringsAsFactors = FALSE)
} else {
  data.frame(id = character(), label = character(), path = character(),
             formula = character(), source = character(), license = character(), note = character(),
             stringsAsFactors = FALSE)
}

new_rows <- data.frame(
  id = c("gas_turbine_co_2015_raw", "gas_turbine_nox_2015_raw",
         "gas_turbine_co_all_raw", "gas_turbine_nox_all_raw"),
  label = c("UCI Gas Turbine 2015: CO emissions",
            "UCI Gas Turbine 2015: NOx emissions",
            "UCI Gas Turbine 2011-2015: CO emissions",
            "UCI Gas Turbine 2011-2015: NOx emissions"),
  path = c("shared/datasets/external/processed/gas_turbine_2015.csv",
           "shared/datasets/external/processed/gas_turbine_2015.csv",
           "shared/datasets/external/processed/gas_turbine_2011_2015.csv",
           "shared/datasets/external/processed/gas_turbine_2011_2015.csv"),
  formula = c(
    "CO ~ AT + AP + AH + AFDP + GTEP + TIT + TAT + CDP + TEY",
    "NOX ~ AT + AP + AH + AFDP + GTEP + TIT + TAT + CDP + TEY",
    "CO ~ year + AT + AP + AH + AFDP + GTEP + TIT + TAT + CDP + TEY",
    "NOX ~ year + AT + AP + AH + AFDP + GTEP + TIT + TAT + CDP + TEY"
  ),
  source = "UCI Gas Turbine CO and NOx Emission Data Set, DOI 10.24432/C5WC95, CC BY 4.0",
  license = "CC BY 4.0",
  note = c(
    "Industrial emissions candidate from DeepResearch PDF; single-year subset for fast C5 intake.",
    "Industrial emissions candidate from DeepResearch PDF; single-year subset for fast C5 intake.",
    "Full 2011-2015 version; use for robust follow-up after fast intake.",
    "Full 2011-2015 version; use for robust follow-up after fast intake."
  ),
  stringsAsFactors = FALSE
)

manifest <- manifest[!(manifest$id %in% new_rows$id), , drop = FALSE]
manifest <- rbind(manifest, new_rows)
utils::write.csv(manifest, manifest_path, row.names = FALSE)

message("Wrote:")
message("  ", all_path, " (n=", nrow(all_data), ")")
message("  ", csv_2015_path, " (n=", sum(all_data$year == 2015), ")")
message("Updated manifest: ", manifest_path)
