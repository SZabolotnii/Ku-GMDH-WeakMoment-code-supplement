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

open_pdf <- function(name, width = 7.2, height = 4.8) {
  .cur_pdf <<- file.path(fig_dir, name)
  pdf(.cur_pdf, width = width, height = height, useDingbats = FALSE)
  par(family = "Helvetica", mar = c(4.4, 7.2, 2.2, 1.2), las = 1, cex = 0.9)
}

# Journals reject PDFs with unembedded fonts, and pdf() does not embed the base-14
# fonts; Ghostscript (embedFonts) embeds them after the device is closed.
close_pdf <- function() {
  invisible(dev.off())
  tmp <- tempfile(fileext = ".pdf")
  grDevices::embedFonts(.cur_pdf, outfile = tmp,
                        options = "-dPDFSETTINGS=/prepress -dEmbedAllFonts=true -dSubsetFonts=true")
  invisible(file.copy(tmp, .cur_pdf, overwrite = TRUE))
}

# Validation gate on all fifteen candidate rows. Each point is one split: the gate's
# trimmed-RMSE divided by the comparator's on the same split (gate and robust baseline
# averaged over three tournament seeds); the bar is the median. Left: against the row's
# best robust baseline. Right: against the best inner-CV-tuned tree ensemble.
gs <- read_result("claimrow_gate_splits.csv")
groups <- list(
  `insurance severity` = c("fremtpl2_severity_raw", "fremtpl2_severity_log", "insurance_autobi_loss",
                           "insurance_autoclaims_paid"),
  `other cross-sectional` = c("islr_credit_balance", "islr_wage", "mass_boston_medv", "concrete",
                              "airquality_ozone", "mass_cars93_mpg"),
  `soft sensors, blocked folds` = c("sru_y1_dynamic", "sru_y1_static", "sru_y2_static",
                                    "gas_turbine_co_2015_raw", "gas_turbine_nox_2015_raw"))
row_lab <- c(fremtpl2_severity_raw = "freMTPL2, raw", fremtpl2_severity_log = "freMTPL2, log",
             insurance_autobi_loss = "AutoBi", insurance_autoclaims_paid = "AutoClaims",
             islr_credit_balance = "credit balance", islr_wage = "Wage", mass_boston_medv = "Boston",
             concrete = "concrete", airquality_ozone = "air quality", mass_cars93_mpg = "Cars93",
             sru_y1_dynamic = "SRU y1, dynamic", sru_y1_static = "SRU y1, static",
             sru_y2_static = "SRU y2, static", gas_turbine_co_2015_raw = "gas-turbine CO",
             gas_turbine_nox_2015_raw = "gas-turbine NOx")
gcol <- c(`insurance severity` = "#2b7a78", `other cross-sectional` = "#7570b3",
          `soft sensors, blocked folds` = "#d95f02")
ys <- list(); y <- 0
for (g in rev(names(groups))) { for (d in rev(groups[[g]])) { y <- y + 1; ys[[d]] <- y }; y <- y + 0.8 }
pool_panel <- function(num, den, xlim, at, main, labels_left) {
  plot(NA, xlim = xlim, ylim = c(0.4, y - 0.4), log = "x", yaxt = "n", xaxt = "n",
       xlab = "", ylab = "", main = main, cex.main = 0.9, font.main = 1)
  abline(v = 1, lty = 2, col = "grey35")
  axis(1, at = 1 + at / 100, labels = ifelse(at == 0, "0", sprintf("%+d %%", at)), cex.axis = 0.75)
  for (g in names(groups)) for (d in groups[[g]]) {
    z <- gs[gs$dataset == d, ]
    raw <- z[[num]] / z[[den]]
    r <- pmin(pmax(raw, xlim[1]), xlim[2])
    out <- raw < xlim[1] | raw > xlim[2]   # off-axis splits: open triangle at the edge
    set.seed(match(d, names(row_lab)))
    points(r, ys[[d]] + runif(length(r), -0.18, 0.18), pch = ifelse(out, 2, 16),
           cex = ifelse(out, 0.6, 0.5), col = adjustcolor(gcol[[g]], 0.6))
    segments(median(z[[num]] / z[[den]]), ys[[d]] - 0.32, median(z[[num]] / z[[den]]), ys[[d]] + 0.32,
             lwd = 2.4, col = gcol[[g]])
  }
  if (labels_left) {
    axis(2, at = unlist(ys), labels = row_lab[names(ys)], tick = FALSE, cex.axis = 0.72)
    for (g in names(groups)) mtext(g, side = 2, line = 0.3, at = max(unlist(ys[groups[[g]]])) + 0.62,
                                   adj = 1, cex = 0.62, font = 3, col = gcol[[g]])
  }
  mtext("gate vs comparator, same split (negative is better)", side = 1, line = 2.5, cex = 0.7)
}
open_pdf("fig_pool_rows.pdf", width = 7.4, height = 5.6)
layout(matrix(1:2, 1), widths = c(1.35, 1))
par(mar = c(4, 7.6, 2.2, 0.6), xpd = FALSE)
pool_panel("trmse_gate", "trmse_robust", c(0.6, 1.6), c(-40, -20, 0, 20, 50),
           "vs best robust baseline", TRUE)
par(mar = c(4, 0.8, 2.2, 0.8))
pool_panel("trmse_gate", "tree_trmse", c(0.3, 3.3), c(-60, -30, 0, 50, 100, 200),
           "vs tuned tree ensemble", FALSE)
close_pdf()

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
    # A right-hand label that would run past the plot frame goes to the left of its point
    # ("+ GyroXmed (control)" on CSS-008).
    near_right <- d$gamma3 + strwidth(lab, cex = 0.72) + 2 * strwidth("m", cex = 0.72) >
      par("usr")[2]
    text(d$gamma3, d$bestweak_vs_best_pct, lab,
         pos = ifelse(d$feat_set == "k2", 1, ifelse(near_right, 2, 4)),
         cex = 0.72, col = scol[[s]], offset = 0.6)
  }
  close_pdf()
}

cat("Generated figures in ", fig_dir, "\n", sep = "")
