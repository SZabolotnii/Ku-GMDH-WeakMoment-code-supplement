#!/usr/bin/env Rscript
# Prepare the first external C5 insurance-severity candidate.
#
# Source: CASdatasets / freMTPL2 French Motor Third-Party Liability datasets.
# The script downloads the two official .rda files, builds a deterministic
# positive-claim severity table, and writes an external manifest consumed by
# screen_cumulants.R and experiments/run_realworld.R.
#
# Run from the repository root:
#   Rscript shared/datasets/prepare_fremtpl2.R [n_sample]

args <- commandArgs(trailingOnly = TRUE)
N_SAMPLE <- if (length(args) >= 1) as.integer(args[1]) else 4000L

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
base_dir <- file.path(root, "shared", "datasets", "external")
raw_dir <- file.path(base_dir, "raw")
processed_dir <- file.path(base_dir, "processed")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

sources <- data.frame(
  name = c("freMTPL2freq", "freMTPL2sev"),
  url = c(
    "https://raw.githubusercontent.com/dutangc/CASdatasets/master/data/freMTPL2freq.rda",
    "https://raw.githubusercontent.com/dutangc/CASdatasets/master/data/freMTPL2sev.rda"
  ),
  stringsAsFactors = FALSE
)

download_one <- function(name, url) {
  dest <- file.path(raw_dir, paste0(name, ".rda"))
  if (!file.exists(dest) || file.info(dest)$size <= 0) {
    message("Downloading ", url)
    utils::download.file(url, dest, mode = "wb", quiet = FALSE)
  }
  dest
}

paths <- Map(download_one, sources$name, sources$url)
env <- new.env(parent = emptyenv())
invisible(lapply(paths, load, envir = env))
freMTPL2freq <- get("freMTPL2freq", envir = env)
freMTPL2sev <- get("freMTPL2sev", envir = env)

sev <- stats::aggregate(ClaimAmount ~ IDpol, data = freMTPL2sev, FUN = sum)
names(sev)[names(sev) == "ClaimAmount"] <- "ClaimAmount"
dat <- merge(freMTPL2freq, sev, by = "IDpol", all = FALSE)
dat <- dat[dat$ClaimAmount > 0 & dat$ClaimNb > 0 & dat$Exposure > 0, ]
dat$log_claim_amount <- log(dat$ClaimAmount)
dat$log_density <- log(dat$Density)
dat$Area <- factor(dat$Area)
dat$VehGas <- factor(dat$VehGas)

keep <- c("IDpol", "ClaimAmount", "log_claim_amount", "Exposure", "VehPower",
          "VehAge", "DrivAge", "BonusMalus", "Density", "log_density",
          "Area", "VehGas", "ClaimNb")
dat <- stats::na.omit(dat[, keep])

set.seed(20260525L)
if (nrow(dat) > N_SAMPLE) {
  bins <- cut(dat$ClaimAmount,
              breaks = unique(stats::quantile(dat$ClaimAmount, probs = seq(0, 1, 0.1), na.rm = TRUE)),
              include.lowest = TRUE)
  per_bin <- max(1L, ceiling(N_SAMPLE / length(levels(bins))))
  idx <- unlist(tapply(seq_len(nrow(dat)), bins, function(ii) sample(ii, min(length(ii), per_bin))),
                use.names = FALSE)
  if (length(idx) > N_SAMPLE) idx <- sample(idx, N_SAMPLE)
  dat_sample <- dat[sort(idx), ]
} else {
  dat_sample <- dat
}

full_path <- file.path(processed_dir, "fremtpl2_severity.csv")
sample_path <- file.path(processed_dir, "fremtpl2_severity_sample.csv")
utils::write.csv(dat, full_path, row.names = FALSE)
utils::write.csv(dat_sample, sample_path, row.names = FALSE)

rel_sample <- "shared/datasets/external/processed/fremtpl2_severity_sample.csv"
manifest <- data.frame(
  id = c("fremtpl2_severity_raw", "fremtpl2_severity_log"),
  label = c(
    "freMTPL2 severity: claim amount ~ policy risk factors",
    "freMTPL2 severity: log claim amount ~ policy risk factors"
  ),
  path = c(rel_sample, rel_sample),
  formula = c(
    "ClaimAmount ~ Exposure + VehPower + VehAge + DrivAge + BonusMalus + log_density + Area + VehGas",
    "log_claim_amount ~ Exposure + VehPower + VehAge + DrivAge + BonusMalus + log_density + Area + VehGas"
  ),
  source = c("CASdatasets::freMTPL2freq/freMTPL2sev", "CASdatasets::freMTPL2freq/freMTPL2sev"),
  license = c("GPL-2 via CASdatasets package metadata", "GPL-2 via CASdatasets package metadata"),
  note = c(
    "Positive-claim policy-level severity sample; claim amounts aggregated by policy ID.",
    "Log-severity contrast for the same positive-claim policy-level sample."
  ),
  stringsAsFactors = FALSE
)
manifest_path <- file.path(base_dir, "external_candidates.csv")
utils::write.csv(manifest, manifest_path, row.names = FALSE)

cat(sprintf("Prepared freMTPL2 severity: full n=%d, sample n=%d\n", nrow(dat), nrow(dat_sample)))
cat(sprintf("Saved %s\n", full_path))
cat(sprintf("Saved %s\n", sample_path))
cat(sprintf("Saved %s\n", manifest_path))
