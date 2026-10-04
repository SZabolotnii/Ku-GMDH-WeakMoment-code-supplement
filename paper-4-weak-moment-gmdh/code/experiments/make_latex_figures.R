#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
script_arg <- args[grepl("^--file=", args)]
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "paper-4-weak-moment-gmdh/code/experiments/make_latex_figures.R"
script_path <- normalizePath(script_path, mustWork = TRUE)
p4_dir <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = TRUE)
results_dir <- file.path(p4_dir, "results")
fig_dir <- file.path(p4_dir, "latex", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_result <- function(name) {
  path <- file.path(results_dir, name)
  if (!file.exists(path)) stop("Missing result file: ", path, call. = FALSE)
  read.csv(path, check.names = FALSE)
}

short_label <- function(x) {
  map <- c(
    sru_y1_dynamic = "SRU y1",
    gas_turbine_co_2015_raw = "GT CO",
    gas_turbine_nox_2015_raw = "GT NOx",
    islr_credit_balance = "Credit",
    insurance_autobi_loss = "AutoBi",
    fremtpl2_severity_raw = "freMTPL2",
    sru_y2_static = "SRU y2",
    concrete = "Concrete"
  )
  out <- unname(map[as.character(x)])
  out[is.na(out)] <- as.character(x)[is.na(out)]
  out
}

open_pdf <- function(name, width = 7.2, height = 4.8) {
  pdf(file.path(fig_dir, name), width = width, height = height, useDingbats = FALSE)
  par(family = "Helvetica", mar = c(4.4, 7.2, 2.2, 1.2), las = 1, cex = 0.9)
}

close_pdf <- function() invisible(dev.off())

claim <- read_result("final_claim_evidence_table.csv")
meta <- read_result("final_claim_meta_analysis.csv")

claim$dataset_label <- short_label(claim$dataset)
claim$dataset_label <- factor(claim$dataset_label, levels = rev(claim$dataset_label))

open_pdf("fig_final_claim_delta_map.pdf", width = 7.4, height = 5.4)
vals <- rbind(
  `vs robust stat` = claim$weak_vs_robust_stat_trmse_pct,
  `vs tree` = claim$weak_vs_tree_trmse_pct
)
colnames(vals) <- as.character(claim$dataset_label)
ylim <- c(0.5, ncol(vals) + 0.5)
bp <- barplot(vals[, ncol(vals):1], beside = TRUE, horiz = TRUE,
              col = c("#2b7a78", "#d95f02"), border = NA,
              xlab = "Weak/PMM trimmed-RMSE delta, % (negative is better)",
              names.arg = rev(as.character(claim$dataset_label)),
              cex.names = 0.82, xlim = range(c(vals, 0), na.rm = TRUE) * c(1.08, 1.08))
abline(v = 0, lty = 2, col = "grey35")
legend("bottomright", legend = rownames(vals), fill = c("#2b7a78", "#d95f02"),
       bty = "n", cex = 0.82)
title("Final claim rows: statistical and tree boundaries", cex.main = 1)
close_pdf()

open_pdf("fig_meta_bootstrap_ci.pdf", width = 7.2, height = 4.6)
meta$label <- paste(meta$group, meta$comparison, sep = "\n")
ord <- seq_len(nrow(meta))
xlim <- range(c(meta$bootstrap_ci_lo_pct, meta$bootstrap_ci_hi_pct, 0), na.rm = TRUE)
plot(meta$bootstrap_mean_pct, rev(ord), pch = 19, col = "#2b7a78",
     xlim = xlim * c(1.08, 1.08), yaxt = "n",
     xlab = "Dataset-row bootstrap mean delta, % (negative is better)",
     ylab = "")
segments(meta$bootstrap_ci_lo_pct, rev(ord), meta$bootstrap_ci_hi_pct, rev(ord),
         col = "#2b7a78", lwd = 2)
abline(v = 0, lty = 2, col = "grey35")
axis(2, at = rev(ord), labels = meta$label, cex.axis = 0.68)
title("Dataset-row meta checkpoint", cex.main = 1)
close_pdf()

dispatch_path <- file.path(results_dir, "dispatch_retune.csv")
if (file.exists(dispatch_path)) {
  disp <- read.csv(dispatch_path, check.names = FALSE)
  disp$label <- c("Credit", "SRU y2", "GT NOx")[seq_len(nrow(disp))]
  open_pdf("fig_validation_gate_selector.pdf", width = 7.1, height = 3.8)
  dvals <- rbind(
    `cumulant gate` = disp$autoweak_vs_rob,
    `validation gate` = disp$autovalgate_vs_rob
  )
  colnames(dvals) <- disp$label
  barplot(dvals, beside = TRUE, horiz = TRUE, names.arg = disp$label,
          col = c("#7570b3", "#1b9e77"), border = NA,
          xlab = "Delta vs best robust baseline, % (negative is better)",
          xlim = range(c(dvals, 0), na.rm = TRUE) * c(1.08, 1.08),
          cex.names = 0.82)
  abline(v = 0, lty = 2, col = "grey35")
  legend("bottomright", legend = rownames(dvals), fill = c("#7570b3", "#1b9e77"),
         bty = "n", cex = 0.82)
  title("Validation-gated dispatch selector check", cex.main = 1)
  close_pdf()
}

ins_path <- file.path(results_dir, "insurance_severity_permethod.csv")
if (file.exists(ins_path)) {
  ins <- read.csv(ins_path, check.names = FALSE)
  rescue <- subset(ins, candidate == "fremtpl2_severity_raw")
  keep <- c("auto", "auto-weak", "WPMM2", "WPMM3", "LSE", "Huber", "L1")
  rescue <- rescue[match(keep, rescue$method), ]
  rescue <- rescue[!is.na(rescue$method), ]
  open_pdf("fig_fremtpl2_heavytail_rescue.pdf", width = 7.1, height = 4.2)
  barplot(rescue$trmse_med, names.arg = rescue$method, col = ifelse(rescue$family == "pmm", "#2b7a78", "#d95f02"),
          border = NA, ylab = "Median test trimmed-RMSE", las = 2, cex.names = 0.72)
  legend("topright", legend = c("PMM/weak family", "robust stat"), fill = c("#2b7a78", "#d95f02"),
         bty = "n", cex = 0.82)
  title("freMTPL2 raw severity: heavy-tail rescue", cex.main = 1)
  close_pdf()
}

fs_path <- file.path(results_dir, "feature_count_sweep.csv")
if (file.exists(fs_path)) {
  fs <- read.csv(fs_path, check.names = FALSE)
  open_pdf("fig_feature_richness_boundary.pdf", width = 7.1, height = 4.4)
  cols <- ifelse(grepl("^k2", fs$feat_set), "#1b9e77", "#d95f02")
  plot(fs$gamma3, fs$bestweak_vs_best_pct, pch = 19, col = cols,
       xlab = "Residual skewness gamma3", ylab = "Weak vs best robust, %",
       ylim = range(c(fs$bestweak_vs_best_pct, 0), na.rm = TRUE) * c(1.12, 1.12))
  abline(h = 0, lty = 2, col = "grey35")
  text(fs$gamma3, fs$bestweak_vs_best_pct, labels = fs$feat_set, pos = 3, cex = 0.62)
  legend("topright", legend = c("low feature count", "added feature"), fill = c("#1b9e77", "#d95f02"),
         bty = "n", cex = 0.82)
  title("Feature richness boundary", cex.main = 1)
  close_pdf()
}

cat("Generated figures in ", fig_dir, "\n", sep = "")
