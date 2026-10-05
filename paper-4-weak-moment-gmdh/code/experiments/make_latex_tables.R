#!/usr/bin/env Rscript
# Generates the result tables of the manuscript from results/*.csv:
#   table_insurance.tex       Table 2: the gate on the four insurance rows, trimmed, untrimmed
#                             and large-claim error, against robust methods and tuned trees
#   table_pool_rows.tex       Table 3: the gate on the other eleven candidate rows
#   table_protocol_lesson.tex Table 4: how a five-fold, single-seed, best-of-family protocol
#                             produced the blocked-row gains
# Usage: Rscript experiments/make_latex_tables.R

script_path <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_path, value = TRUE)
if (length(file_arg) == 0) {
  stop("Cannot determine script path from commandArgs(). Run this script via Rscript.")
}
script_file <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_file), "../../.."), mustWork = TRUE)
p4_dir <- file.path(root, "paper-4-weak-moment-gmdh")
results_dir <- file.path(p4_dir, "results")
tables_dir <- file.path(p4_dir, "latex", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(name) read.csv(file.path(results_dir, name), stringsAsFactors = FALSE, check.names = FALSE)
fmt_pct <- function(x) {
  v <- round(as.numeric(x), 1)
  if (v == 0) return("0.0")
  gsub("-", "$-$", sprintf("%+.1f", v), fixed = TRUE)
}
num_word <- function(k) {
  w <- c("none", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten",
         "eleven", "twelve", "thirteen", "fourteen", "fifteen")
  if (k >= 0 && k <= 15) w[[k + 1]] else as.character(k)
}

dataset_label <- c(
  fremtpl2_severity_raw = "freMTPL2 severity, raw",
  fremtpl2_severity_log = "freMTPL2 severity, log",
  insurance_autobi_loss = "AutoBi loss",
  insurance_autoclaims_paid = "AutoClaims paid",
  islr_credit_balance = "credit balance",
  islr_wage = "Wage",
  mass_boston_medv = "Boston housing",
  concrete = "concrete strength",
  airquality_ozone = "air quality (ozone)",
  mass_cars93_mpg = "Cars93 fuel economy",
  sru_y1_dynamic = "SRU \\(y_1\\), dynamic",
  sru_y1_static = "SRU \\(y_1\\), static",
  sru_y2_static = "SRU \\(y_2\\), static",
  gas_turbine_co_2015_raw = "gas-turbine CO",
  gas_turbine_nox_2015_raw = "gas-turbine NOx"
)
robust_label <- c(L1 = "LAD", Huber = "Huber", LSE = "LSE", `ridge-LSE` = "ridge-LSE")
groups <- list(
  "Insurance severity (random splits)" = c("fremtpl2_severity_raw", "fremtpl2_severity_log",
                                           "insurance_autobi_loss", "insurance_autoclaims_paid"),
  "Other cross-sectional (random splits)" = c("islr_credit_balance", "islr_wage", "mass_boston_medv",
                                              "concrete", "airquality_ozone", "mass_cars93_mpg"),
  "Soft sensors and emissions (blocked folds)" = c("sru_y1_dynamic", "sru_y1_static", "sru_y2_static",
                                                   "gas_turbine_co_2015_raw", "gas_turbine_nox_2015_raw")
)

# A cell is bold when the gate is at least 2 % better in the paired median, wins at least
# 70 % of the splits, and the row has at least ten splits (the five-fold blocked rows did not
# hold up at ten folds, Table 3); nothing else is bold.
WIN_PCT <- -2
WIN_SHARE <- 0.7
MIN_SPLITS <- 10
cell <- function(pct, wins, n) {
  s <- sprintf("%s (%d/%d)", fmt_pct(pct), as.integer(wins), as.integer(n))
  if (as.numeric(pct) <= WIN_PCT && as.numeric(wins) >= WIN_SHARE * as.numeric(n) &&
      as.numeric(n) >= MIN_SPLITS) s <- paste0("\\textbf{", s, "}")
  s
}

tree_label <- c(rf = "RF", extra = "ET", hgb = "HGB")
# Interact tables: \tbl{caption}{body} inside a table float.
float <- function(label, caption, spec, head, body) c(
  "\\begin{table}[t]",
  sprintf("\\tbl{%s}", caption),
  "{\\small\\setlength{\\tabcolsep}{4pt}",
  sprintf("\\begin{tabular}{%s}", spec), "\\toprule", head, "\\midrule", body, "\\bottomrule",
  "\\end{tabular}}", sprintf("\\label{%s}", label), "\\end{table}")
PAIRED <- paste0("Each cell is the median over splits of the paired relative difference ",
                 "\\(100(e_{\\text{gate}}/e_{\\text{comparator}}-1)\\), in \\%%, with the number of splits on which the gate is better; ",
                 "negative means the gate is better. ")
BOLD <- "Bold: at least 2\\,\\%% better and better on at least 70\\,\\%% of at least ten splits."

# Table 2: the four insurance rows, three errors, against the robust baseline and the tuned tree.
ins <- groups[[1]]
tl <- read_csv("claimrow_tail_summary.csv")
ins_body <- character(0)
for (cmp in c("robust", "tree")) {
  ins_body <- c(ins_body, sprintf("\\multicolumn{6}{@{}l}{\\emph{%s}} \\\\",
                                  if (cmp == "robust") "Gate against the robust baseline" else "Gate against the tuned tree ensemble"))
  for (d in ins) {
    r <- tl[tl$dataset == d, ]
    who <- if (cmp == "robust") robust_label[[r$robust_method]] else tree_label[[r$best_tuned_tree]]
    cells <- vapply(c("mae", "trmse", "rmse", "tail_rmse"), function(k)
      cell(r[[sprintf("gate_vs_%s_%s_paired_pct", cmp, k)]], r[[sprintf("gate_vs_%s_%s_wins", cmp, k)]], r$n_splits),
      character(1))
    ins_body <- c(ins_body, paste0(paste(sprintf("\\quad %s", dataset_label[[d]]), who, paste(cells, collapse = " & "),
                                         sep = " & "), " \\\\"))
  }
}
big_win <- tl$gate_vs_tree_trmse_wins == tl$n_splits
ins_tex <- c(
  "% Generated by code/experiments/make_latex_tables.R from results/claimrow_tail_summary.csv.",
  float("tab:insurance",
        # no sprintf here, so the doubled %% of the shared caption pieces is undone by hand
        gsub("%%", "%", fixed = TRUE, x = paste0("On the three raw-scale insurance datasets the gate has a lower mean absolute and trimmed error than tuned tree ensembles on every split, ",
               "while on the largest claims the two are within about 6\\,\\%%; it is never more than 1.5\\,\\%% worse than the best robust method on any error. ",
               "The gate is GMDH with the inner estimator of each partial model chosen by inner cross-validation (Section~\\ref{the-dispatch-cumulant-keyed-under-fires-validation-gated-blocked-inner-folds-delivers}); freMTPL2 log is a log-scale control. ", PAIRED,
               "Errors: mean absolute error (MAE), trimmed RMSE (TRMSE, the 10\\,\\%% largest absolute errors discarded), RMSE, and the RMSE on the test claims above the ",
               "90\\,\\%% quantile of the training response (large claims). Gate and robust methods are averaged over three tournament seeds; 30 random splits. ",
               "Robust baseline: the robust method (LSE, ridge-LSE, Huber, LAD) with the lowest median test error. Tree: random forest (RF), ExtraTrees (ET) or ",
               "histogram gradient boosting (HGB), whichever has the lowest median test error after inner-CV tuning on each outer training set. ",
               "Both comparators are thus chosen by test error, an oracle choice that favours them. ", BOLD)),
        "@{}llrrrr@{}",
        "dataset & comparator & MAE (wins) & TRMSE (wins) & RMSE (wins) & large claims (wins) \\\\",
        ins_body))
writeLines(ins_tex, file.path(tables_dir, "table_insurance.tex"))

# Table 3: the other eleven rows (scope).
pool <- read_csv("claimrow_candidate_pool.csv")
trees <- read_csv("claimrow_tuned_trees_summary.csv")
pool <- merge(pool, trees[, c("dataset", "gate_vs_tuned_paired_pct", "gate_wins", "best_tuned_tree",
                             "gate_vs_tuned_rmse_paired_pct")],
              by = "dataset", suffixes = c("", "_tree"))
oth <- pool[!pool$dataset %in% ins, ]
body <- character(0)
for (g in names(groups)[-1]) {
  body <- c(body, sprintf("\\multicolumn{6}{@{}l}{\\emph{%s}} \\\\", g))
  for (d in groups[[g]]) {
    r <- pool[pool$dataset == d, ]
    body <- c(body, paste0(paste(
      sprintf("\\quad %s", dataset_label[[d]]),
      robust_label[[r$robust_method]],
      cell(r$gate_vs_robust_paired_pct, r$gate_wins, r$n_splits),
      fmt_pct(r$gate_vs_robust_rmse_paired_pct),
      sprintf("%s, %s", cell(r$gate_vs_tuned_paired_pct, r$gate_wins_tree, r$n_splits), tree_label[[r$best_tuned_tree]]),
      fmt_pct(r$gate_vs_tuned_rmse_paired_pct),
      sep = " & "), " \\\\"))
  }
}
ties <- sum(abs(oth$gate_vs_robust_paired_pct) <= 2)
oth_tree_losses <- sum(oth$gate_vs_tuned_paired_pct > 0)
pool_tex <- c(
  "% Generated by code/experiments/make_latex_tables.R from results/claimrow_candidate_pool.csv",
  "% and results/claimrow_tuned_trees_summary.csv.",
  float("tab:pool-rows",
        sprintf(paste0("On the %s other datasets the gate is within 2\\,\\%% of the best robust baseline on %s and behind the tuned tree ensemble on %s. ",
                       PAIRED, "TRMSE: trimmed RMSE; RMSE columns: the same statistic without trimming (paired median only). ",
                       "Gate and robust baseline are averaged over three tournament seeds; 30 random splits, or five blocked folds for the soft-sensor and emission data. ",
                       "Robust baseline: the robust method (LSE, ridge-LSE, Huber, LAD) with the lowest median test error. Tree: random forest (RF), ExtraTrees (ET) or histogram gradient boosting (HGB), ",
                       "whichever has the lowest median test error after inner-CV tuning on each outer training set; both choices use the test errors and favour the comparator. ", BOLD),
                num_word(nrow(oth)), num_word(ties), num_word(oth_tree_losses)),
        "@{}lllrlr@{}",
        c(" & & \\multicolumn{2}{c}{gate vs robust baseline} & \\multicolumn{2}{c}{gate vs tuned tree} \\\\",
          "\\cmidrule(lr){3-4}\\cmidrule(l){5-6}",
          "dataset & robust & TRMSE (wins) & RMSE & TRMSE (wins), tree & RMSE \\\\"),
        body))
writeLines(pool_tex, file.path(tables_dir, "table_pool_rows.tex"))

# Table 4: the protocol lesson on the four blocked rows.
les <- read_csv("claimrow_protocol_lesson.csv")
les_body <- vapply(seq_len(nrow(les)), function(i) {
  r <- les[i, ]
  paste0(paste(
    dataset_label[[r$dataset]],
    fmt_pct(r$hindsight_seed0_ratio_med_pct),
    sprintf("%s, %s", fmt_pct(r$hindsight_seed1_ratio_med_pct), fmt_pct(r$hindsight_seed2_ratio_med_pct)),
    sprintf("%s (%d/%d)", fmt_pct(r$gate5_paired_pct), r$gate5_wins, r$gate5_folds),
    sprintf("%s (%d/%d)", fmt_pct(r$gate10_paired_pct), r$gate10_wins, r$gate10_folds),
    sep = " & "), " \\\\")
}, character(1))
lesson_tex <- c(
  "% Generated by code/experiments/make_latex_tables.R from results/claimrow_protocol_lesson.csv.",
  float("tab:protocol-lesson",
        paste0(
          sprintf("Gains of %.0f--%.0f\\,\\%% on the blocked datasets came from five folds, one tournament seed and a best-of-family choice; ",
                  floor(min(-les$hindsight_seed0_ratio_med_pct[les$role != "tie_control"])),
                  ceiling(max(-les$hindsight_seed0_ratio_med_pct[les$role != "tie_control"]))),
          "with ten folds and five seeds none of them meets the prespecified decision rule of Section~\\ref{protocol-lesson}. ",
          "Columns 2--3: best weak method against best robust method, each chosen by its median error on the test folds, ",
          "as a ratio of medians minus one (\\%), at tournament seed 0 and at two further seeds. ",
          "Columns 4--5: the validation gate against the dataset's robust baseline, paired median (\\%) and folds won; ",
          "five folds averaged over three seeds, and ten folds averaged over five seeds (decision rule fixed before the run). ",
          "SRU: sulfur recovery unit, outputs \\(y_1\\) and \\(y_2\\); \\(y_2\\) is the control."),
        "@{}lrrrr@{}",
        c(" & \\multicolumn{2}{c}{best of family, 5 folds} & \\multicolumn{2}{c}{validation gate} \\\\",
          "\\cmidrule(lr){2-3}\\cmidrule(l){4-5}",
          "dataset & seed 0 & seeds 1, 2 & 5 folds & 10 folds \\\\"),
        les_body))
writeLines(lesson_tex, file.path(tables_dir, "table_protocol_lesson.tex"))

manifest <- data.frame(
  table = c("table_insurance.tex", "table_pool_rows.tex", "table_protocol_lesson.tex"),
  source_csv = c("results/claimrow_tail_summary.csv",
                 "results/claimrow_candidate_pool.csv;results/claimrow_tuned_trees_summary.csv",
                 "results/claimrow_protocol_lesson.csv"),
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(tables_dir, "tables_manifest.csv"), row.names = FALSE)
cat("Generated LaTeX tables in ", tables_dir, "\n", sep = "")
