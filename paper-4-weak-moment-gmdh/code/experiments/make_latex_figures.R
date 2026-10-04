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

open_pdf("fig_final_claim_delta_map.pdf", width = 7.4, height = 4.5)
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
legend("topright", legend = rownames(vals), fill = c("#2b7a78", "#d95f02"),
       bty = "n", cex = 0.82)
title("Final claim rows: statistical and tree boundaries", cex.main = 1)
close_pdf()

# Dataset-row meta checkpoint. One panel per comparison, each on its own scale (the tree
# deltas span -62 to +141 %). Grey points: the per-row deltas (median over splits of the
# paired relative difference); dark point and bar: bootstrap mean and 95 % interval.
grouped <- read_result("grouped_inference_summary.csv")
scope_lab <- c(all_claim_rows = "all 8 rows", non_loss_rows = "without the\nloss control",
               positive_rows = "positive rows")
comp_lab <- c(weak_vs_robust_stat = "vs best robust statistical baseline",
              weak_vs_tree = "vs best tree ensemble")
open_pdf("fig_meta_bootstrap_ci.pdf", width = 7.4, height = 3.9)
layout(matrix(1:2, 1), widths = c(1.15, 1))
for (cmp in names(comp_lab)) {
  m <- meta[meta$comparison == cmp, ]
  m <- m[match(names(scope_lab), m$group), ]
  rows <- lapply(m$group, function(s) {
    v <- grouped$unit_values[grouped$scope == s & grouped$comparison == cmp & grouped$unit == "dataset"]
    as.numeric(sub("^[^:]+:", "", strsplit(v, ";")[[1]]))
  })
  xr <- range(c(unlist(rows), m$bootstrap_ci_lo_pct, m$bootstrap_ci_hi_pct, 0))
  xr <- xr + c(-0.06, 0.06) * diff(xr)
  par(mar = c(4.2, if (cmp == "weak_vs_robust_stat") 7.4 else 1.2, 2.4, 0.8))
  k <- nrow(m)
  plot(NA, xlim = xr, ylim = c(0.4, k + 0.6), yaxt = "n", xlab = "", ylab = "",
       main = comp_lab[[cmp]], cex.main = 0.85, font.main = 1, cex.axis = 0.8)
  abline(v = 0, lty = 2, col = "grey35")
  for (i in seq_len(k)) {
    y <- k + 1 - i
    points(rows[[i]], rep(y + 0.22, length(rows[[i]])), pch = 16, cex = 0.7,
           col = adjustcolor("grey40", 0.6))
    segments(m$bootstrap_ci_lo_pct[i], y - 0.08, m$bootstrap_ci_hi_pct[i], y - 0.08,
             col = "#2b7a78", lwd = 2.2)
    points(m$bootstrap_mean_pct[i], y - 0.08, pch = 19, col = "#2b7a78", cex = 1.1)
    text(xr[2], y - 0.36, sprintf("better on %d of %d", m$weak_win_count[i], m$n_dataset_rows[i]),
         pos = 2, cex = 0.68, col = "grey25")
  }
  if (cmp == "weak_vs_robust_stat") axis(2, at = rev(seq_len(k)), labels = scope_lab[m$group],
                                         tick = FALSE, cex.axis = 0.8)
  mtext("dataset-row delta, % (negative is better)", side = 1, line = 2.6, cex = 0.75)
}
close_pdf()

# Per-split selector check: for every split, the trimmed-RMSE of each route divided by the
# best robust baseline of that row (the robust method with the lowest median, as in the
# table), on the same split. Points are splits (seed-averaged), bars are medians; the right
# margin prints the paired median and the number of splits won.
dispatch_path <- file.path(results_dir, "dispatch_retune_splits.csv")
if (file.exists(dispatch_path)) {
  sp <- read.csv(dispatch_path, check.names = FALSE)
  rows_ds <- c("islr_credit_balance", "sru_y2_static", "gas_turbine_nox_2015_raw")
  routes <- c(`cumulant gate (auto-weak)` = "auto-weak",
              `validation gate (auto-valgate)` = "auto-valgate",
              `best forced weak (hindsight)` = "oracle")
  rcol <- c("#7570b3", "#1b9e77", "grey45")
  open_pdf("fig_validation_gate_selector.pdf", width = 7.4, height = 4.3)
  par(mar = c(4.2, 6.2, 2.6, 7.6), xpd = FALSE)
  pct_at <- c(-30, -15, 0, 25, 50, 100, 200)
  plot(NA, xlim = c(0.68, 3.0), ylim = c(0.4, length(rows_ds) * 4 - 0.6), log = "x",
       yaxt = "n", xaxt = "n", xlab = "", ylab = "")
  axis(1, at = 1 + pct_at / 100, labels = ifelse(pct_at == 0, "0", sprintf("%+d %%", pct_at)),
       cex.axis = 0.78)
  mtext("trimmed-RMSE vs best robust baseline, same split (negative is better)",
        side = 1, line = 2.6, cex = 0.8)
  abline(v = 1, lty = 2, col = "grey35")
  mtext("paired median, splits won", side = 4, line = 0.4, at = length(rows_ds) * 4 - 0.3,
        las = 1, adj = 0, cex = 0.62, col = "grey25", padj = 0)
  for (j in seq_along(rows_ds)) {
    d <- sp[sp$dataset == rows_ds[j], ]
    rob <- c("LSE", "Huber", "L1")[which.min(sapply(c("LSE", "Huber", "L1"), function(m) median(d[[m]])))]
    wk <- c("WPMM2", "WPMM3")[which.min(sapply(c("WPMM2", "WPMM3"), function(m) median(d[[m]])))]
    d$oracle <- d[[wk]]
    base <- (length(rows_ds) - j) * 4
    for (i in seq_along(routes)) {
      y <- base + 4 - i
      r <- d[[routes[[i]]]] / d[[rob]]
      set.seed(10 * j + i)
      points(r, y + runif(length(r), -0.15, 0.15), pch = 16, cex = 0.7, col = adjustcolor(rcol[i], 0.7))
      segments(median(r), y - 0.32, median(r), y + 0.32, lwd = 2.5, col = rcol[i])
      mtext(sprintf("%+.1f %%, %d of %d", 100 * (median(r) - 1), sum(r < 1), length(r)),
            side = 4, line = 0.4, at = y, las = 1, adj = 0, cex = 0.66, col = rcol[i])
    }
    axis(2, at = base + 2, labels = sprintf("%s\n(vs %s)", short_label(rows_ds[j]),
                                            if (rob == "L1") "LAD" else rob),
         tick = FALSE, cex.axis = 0.8)
    if (j > 1) abline(h = base + 3.5, col = "grey85")
  }
  par(xpd = NA)
  legend(x = 0.68, y = length(rows_ds) * 4 + 0.9, legend = names(routes), col = rcol, pch = 16,
         bty = "n", cex = 0.72, horiz = TRUE, xjust = 0, yjust = 0.5)
  close_pdf()
}

# Per-split strip chart: each point is one of the 30 random splits, trimmed-RMSE divided by
# LAD (L1) on the same split; the bar is the median. Left panel, log scale: the failures.
# Right panel, linear zoom: the windowed estimators against Huber, i.e. the tie.
strip_panel <- function(vals, cols, xlim, log_x, at, lab, main) {
  k <- length(vals)
  plot(NA, xlim = xlim, ylim = c(0.5, k + 0.5), log = if (log_x) "x" else "",
       yaxt = "n", xaxt = "n", xlab = "", ylab = "", main = main, cex.main = 0.9, font.main = 1)
  abline(v = 1, lty = 2, col = "grey35")
  axis(1, at = at, labels = lab, cex.axis = 0.8)
  for (i in seq_len(k)) {
    y <- k + 1 - i
    if (median(vals[[i]]) > xlim[2]) {
      text(xlim[2], y, sprintf("off scale (median %.1f)", median(vals[[i]])), cex = 0.72,
           col = cols[i], pos = 2)
      next
    }
    boxplot(vals[[i]], at = y, horizontal = TRUE, add = TRUE, axes = FALSE, outline = FALSE,
            boxwex = 0.6, col = adjustcolor(cols[i], 0.25), border = cols[i], lwd = 1.2,
            medlwd = 2.5, staplewex = 0.4)
    set.seed(i)
    points(vals[[i]], y + runif(length(vals[[i]]), -0.12, 0.12), pch = 16, cex = 0.45,
           col = adjustcolor(cols[i], 0.75))
  }
}

ins_path <- file.path(results_dir, "insurance_severity_raw.csv")
if (file.exists(ins_path)) {
  ins <- read.csv(ins_path, check.names = FALSE)
  ins <- subset(ins, candidate == "fremtpl2_severity_raw")
  w <- reshape(ins[, c("rep", "method", "trmse")], idvar = "rep", timevar = "method", direction = "wide")
  ratio <- function(m) w[[paste0("trmse.", m)]] / w[["trmse.L1"]]
  labs <- c(LSE = "LSE", auto = "classical PMM (auto)", WPMM2 = "WPMM2", WPMM3 = "WPMM3",
            `auto-weak` = "windowed dispatch (auto-weak)", Huber = "Huber")
  cols <- c(LSE = "#d95f02", auto = "#7570b3", WPMM2 = "#2b7a78", WPMM3 = "#2b7a78",
            `auto-weak` = "#2b7a78", Huber = "#d95f02")
  vals <- lapply(names(labs), ratio)
  open_pdf("fig_fremtpl2_heavytail_rescue.pdf", width = 6.6, height = 3.3)
  layout(matrix(1:2, 1), widths = c(1.45, 1))
  par(mar = c(4.2, 11.8, 2.2, 0.6))
  strip_panel(vals, cols, xlim = c(0.9, 6.5), log_x = TRUE, at = c(1, 1.5, 2, 3, 4, 6),
              lab = c("1", "1.5", "2", "3", "4", "6"), main = "all estimators (log scale)")
  axis(2, at = rev(seq_along(labs)), labels = labs, cex.axis = 0.8, tick = FALSE)
  mtext("trimmed-RMSE / LAD, same split", side = 1, line = 2.6, cex = 0.8)
  par(mar = c(4.2, 0.8, 2.2, 0.8))
  strip_panel(vals, cols, xlim = c(0.965, 1.035), log_x = FALSE,
              at = c(0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03),
              lab = c("0.97", "", "0.99", "1", "1.01", "", "1.03"), main = "zoom near 1 (linear)")
  mtext("trimmed-RMSE / LAD, same split", side = 1, line = 2.6, cex = 0.8)
  close_pdf()
}

# Feature-count sweep on two drilling sensors: each arrow starts at the two-feature node
# (k2) and ends at one of the two three-feature nodes. Direct labels, no legend, so nothing
# covers a point.
fs_path <- file.path(results_dir, "feature_count_sweep.csv")
if (file.exists(fs_path)) {
  fs <- read.csv(fs_path, check.names = FALSE)
  open_pdf("fig_feature_richness_boundary.pdf", width = 7.1, height = 4.4)
  par(mar = c(4.4, 4.6, 1.2, 1.2))
  scol <- c(`CSS-007` = "#1b9e77", `CSS-008` = "#d95f02")
  plot(NA, xlim = range(fs$gamma3) + c(-0.15, 0.25),
       ylim = range(c(fs$bestweak_vs_best_pct, 0)) + c(-5, 6),
       xlab = expression("residual skewness " * gamma[3]),
       ylab = "best weak vs best robust, % (negative is better)", cex.axis = 0.85)
  rect(par("usr")[1], par("usr")[3], par("usr")[2], 0, col = adjustcolor("#2b7a78", 0.06), border = NA)
  abline(h = 0, lty = 2, col = "grey35")
  abline(v = 0, lty = 3, col = "grey70")
  text(par("usr")[1], -1.5, "weak/PMM better", adj = c(-0.05, 1), cex = 0.7, col = "grey30")
  text(par("usr")[1], 1.5, "robust baseline better", adj = c(-0.05, 0), cex = 0.7, col = "grey30")
  for (s in names(scol)) {
    d <- fs[fs$sensor == s, ]
    b <- d[d$feat_set == "k2", ]
    for (f in setdiff(d$feat_set, "k2")) {
      e <- d[d$feat_set == f, ]
      arrows(b$gamma3, b$bestweak_vs_best_pct, e$gamma3, e$bestweak_vs_best_pct,
             length = 0.08, col = adjustcolor(scol[[s]], 0.7), lwd = 1.4)
    }
    points(d$gamma3, d$bestweak_vs_best_pct, pch = ifelse(d$feat_set == "k2", 17, 19),
           col = scol[[s]], cex = 1.2)
    lab <- ifelse(d$feat_set == "k2", paste0(s, ", 2 features"),
                  c(k3_gyromed = "+ GyroXmed (control)", k3_ShXrms = "+ ShXrms")[d$feat_set])
    text(d$gamma3, d$bestweak_vs_best_pct, lab, pos = ifelse(d$feat_set == "k2", 1, 4),
         cex = 0.72, col = scol[[s]], offset = 0.6)
  }
  close_pdf()
}

cat("Generated figures in ", fig_dir, "\n", sep = "")
