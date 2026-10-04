#!/usr/bin/env Rscript
# C5 dataset screening for Paper 1 (GMDH-PMM).
#
# This script is intentionally conservative. It screens candidate regression
# datasets by fitting a simple OLS model, bootstrapping residual cumulants, and
# applying the same PMM dispatch logic used inside gmdhpmm. Passing this screen
# does not prove a flagship dataset; it only says the dataset deserves a full
# GMDH-PMM run.
#
# Run from the repository root:
#   Rscript shared/datasets/screen_cumulants.R [B]

args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args) >= 1) as.integer(args[1]) else 500L

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
pkg_dir <- file.path(root, "paper-1-gmdh-pmm", "code")

load_gmdhpmm <- function() {
  if (requireNamespace("gmdhpmm", quietly = TRUE)) {
    library(gmdhpmm)
    return(invisible(TRUE))
  }
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(pkg_dir, quiet = TRUE)
    return(invisible(TRUE))
  }
  rdir <- file.path(pkg_dir, "R")
  for (f in c("gmdhpmm-package.R", "kg2.R", "cumulants.R", "dispatch.R")) {
    sys.source(file.path(rdir, f), envir = .GlobalEnv)
  }
  invisible(TRUE)
}
load_gmdhpmm()

candidate <- function(id, label, data, formula, source, note = "") {
  list(id = id, label = label, data = data, formula = formula, source = source, note = note)
}

has_package <- function(pkg) requireNamespace(pkg, quietly = TRUE)

load_data <- function(pkg, name) {
  env <- new.env(parent = globalenv())
  suppressWarnings(utils::data(list = name, package = pkg, envir = env))
  if (!exists(name, envir = env, inherits = FALSE)) {
    stop("Dataset ", pkg, "::", name, " is unavailable in this R installation.")
  }
  get(name, envir = env, inherits = FALSE)
}

make_external_candidates <- function() {
  manifest_path <- file.path(root, "shared", "datasets", "external", "external_candidates.csv")
  if (!file.exists(manifest_path)) return(list())
  manifest <- utils::read.csv(manifest_path, stringsAsFactors = FALSE)
  required <- c("id", "label", "path", "formula", "source", "note")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) stop("External manifest is missing columns: ", paste(missing, collapse = ", "))

  rows <- split(manifest, seq_len(nrow(manifest)))
  out <- lapply(rows, function(row) {
    data_path <- row$path[[1]]
    if (!grepl("^(/|[A-Za-z]:)", data_path)) data_path <- file.path(root, data_path)
    if (!file.exists(data_path)) {
      warning("Skipping external candidate ", row$id[[1]], ": file not found: ", data_path)
      return(NULL)
    }
    dat <- utils::read.csv(data_path, stringsAsFactors = TRUE)
    candidate(row$id[[1]], row$label[[1]], dat, stats::as.formula(row$formula[[1]]),
              row$source[[1]], row$note[[1]])
  })
  Filter(Negate(is.null), out)
}

make_candidates <- function() {
  data(iris, package = "datasets")
  data(trees, package = "datasets")
  data(cars, package = "datasets")
  data(faithful, package = "datasets")
  data(mtcars, package = "datasets")
  data(airquality, package = "datasets")
  data(swiss, package = "datasets")
  data(women, package = "datasets")
  data(pressure, package = "datasets")
  data(ToothGrowth, package = "datasets")

  cands <- list(
    candidate("iris_versicolor_sl_sw", "Iris versicolor: Sepal.Length ~ Sepal.Width",
              subset(iris, Species == "versicolor"), Sepal.Length ~ Sepal.Width,
              "datasets::iris", "PMM3 sanity candidate from prior PMM screening."),
    candidate("iris_setosa_sl_sw", "Iris setosa: Sepal.Length ~ Sepal.Width",
              subset(iris, Species == "setosa"), Sepal.Length ~ Sepal.Width,
              "datasets::iris", "Small n; useful negative/contrast cell."),
    candidate("iris_virginica_sl_sw", "Iris virginica: Sepal.Length ~ Sepal.Width",
              subset(iris, Species == "virginica"), Sepal.Length ~ Sepal.Width,
              "datasets::iris", "Small n; useful negative/contrast cell."),
    candidate("trees_volume", "Trees: Volume ~ Girth + Height",
              trees, Volume ~ Girth + Height, "datasets::trees"),
    candidate("cars_dist", "Cars: dist ~ speed",
              cars, dist ~ speed, "datasets::cars"),
    candidate("faithful_eruptions", "Old Faithful: eruptions ~ waiting",
              faithful, eruptions ~ waiting, "datasets::faithful"),
    candidate("mtcars_mpg", "mtcars: mpg ~ wt + hp + disp",
              mtcars, mpg ~ wt + hp + disp, "datasets::mtcars"),
    candidate("airquality_ozone", "Airquality: Ozone ~ Solar.R + Wind + Temp",
              stats::na.omit(airquality), Ozone ~ Solar.R + Wind + Temp, "datasets::airquality"),
    candidate("swiss_fertility", "Swiss: Fertility ~ socio-economic covariates",
              swiss, Fertility ~ Agriculture + Examination + Education + Catholic + Infant.Mortality,
              "datasets::swiss"),
    candidate("women_weight", "Women: weight ~ height",
              women, weight ~ height, "datasets::women", "Very small n; included as nonlinear negative control."),
    candidate("pressure", "Pressure: pressure ~ temperature",
              pressure, pressure ~ temperature, "datasets::pressure", "Nonlinear physical curve; negative control."),
    candidate("toothgrowth", "ToothGrowth: len ~ dose + supp",
              ToothGrowth, len ~ dose + supp, "datasets::ToothGrowth")
  )

  if (has_package("MASS")) {
    Insurance <- load_data("MASS", "Insurance")
    Insurance$claim_rate <- Insurance$Claims / pmax(Insurance$Holders, 1)
    birthwt <- load_data("MASS", "birthwt")
    birthwt$race <- factor(birthwt$race)
    mcycle <- load_data("MASS", "mcycle")
    Cars93 <- stats::na.omit(load_data("MASS", "Cars93"))
    Boston <- load_data("MASS", "Boston")
    Rubber <- load_data("MASS", "Rubber")
    cpus <- stats::na.omit(load_data("MASS", "cpus"))
    petrol <- load_data("MASS", "petrol")

    cands <- c(cands, list(
      candidate("mass_insurance_claims", "MASS Insurance: Claims ~ Holders + District + Group + Age",
                Insurance, Claims ~ Holders + District + Group + Age,
                "MASS::Insurance", "Car-insurance count data; screening only, not yet a severity model."),
      candidate("mass_insurance_rate", "MASS Insurance: Claims/Holders ~ District + Group + Age",
                Insurance, claim_rate ~ District + Group + Age,
                "MASS::Insurance", "Rate version of the insurance screen."),
      candidate("mass_birthwt_bwt", "Birth weight: bwt ~ clinical risk factors",
                birthwt, bwt ~ age + lwt + race + smoke + ptl + ht + ui + ftv,
                "MASS::birthwt", "Medical/clinical continuous response."),
      candidate("mass_mcycle_accel_poly3", "Motorcycle crash: accel ~ cubic(time)",
                mcycle, accel ~ times + I(times^2) + I(times^3),
                "MASS::mcycle", "Engineering dynamics; nonlinear baseline before full GMDH."),
      candidate("mass_cars93_price", "Cars93: Price ~ vehicle covariates",
                Cars93, Price ~ MPG.city + Horsepower + EngineSize + RPM + Weight + Type + Origin,
                "MASS::Cars93", "Automotive price data; moderate n with mixed covariates."),
      candidate("mass_cars93_mpg", "Cars93: MPG.city ~ vehicle covariates",
                Cars93, MPG.city ~ Horsepower + EngineSize + RPM + Weight + Type + Origin,
                "MASS::Cars93", "Automotive fuel-economy contrast."),
      candidate("mass_boston_medv", "Boston housing: medv ~ structural covariates",
                Boston, medv ~ crim + zn + indus + chas + nox + rm + age + dis + rad + tax + ptratio + lstat,
                "MASS::Boston", "Screening only: overused dataset; exclude sensitive race proxy from formula."),
      candidate("mass_rubber_loss", "Rubber: loss ~ hard + tens",
                Rubber, loss ~ hard + tens,
                "MASS::Rubber", "Materials/tyre testing, small n."),
      candidate("mass_cpus_perf", "Computer CPUs: perf ~ hardware characteristics",
                cpus, perf ~ syct + mmin + mmax + cach + chmin + chmax,
                "MASS::cpus", "Engineering/computer performance data."),
      candidate("mass_petrol_y", "Petrol refinery: Y ~ SG + VP + V10 + EP",
                petrol, Y ~ SG + VP + V10 + EP,
                "MASS::petrol", "Small industrial process dataset.")
    ))
  }

  if (has_package("insuranceData")) {
    AutoBi <- stats::na.omit(load_data("insuranceData", "AutoBi"))
    AutoBi$ATTORNEY <- factor(AutoBi$ATTORNEY)
    AutoBi$CLMSEX <- factor(AutoBi$CLMSEX)
    AutoBi$MARITAL <- factor(AutoBi$MARITAL)
    AutoBi$CLMINSUR <- factor(AutoBi$CLMINSUR)
    AutoBi$SEATBELT <- factor(AutoBi$SEATBELT)
    AutoBi$log_loss <- log1p(AutoBi$LOSS)

    AutoClaims <- stats::na.omit(load_data("insuranceData", "AutoClaims"))
    AutoClaims$log_paid <- log1p(AutoClaims$PAID)

    cands <- c(cands, list(
      candidate("insurance_autobi_loss", "AutoBi bodily injury: LOSS ~ claimant/legal covariates",
                AutoBi, LOSS ~ ATTORNEY + CLMSEX + MARITAL + CLMINSUR + SEATBELT + CLMAGE,
                "insuranceData::AutoBi", "DeepResearch candidate: bodily-injury claim severity, raw loss scale."),
      candidate("insurance_autobi_log_loss", "AutoBi bodily injury: log1p(LOSS) ~ claimant/legal covariates",
                AutoBi, log_loss ~ ATTORNEY + CLMSEX + MARITAL + CLMINSUR + SEATBELT + CLMAGE,
                "insuranceData::AutoBi", "Log-scale control for the bodily-injury severity candidate."),
      candidate("insurance_autoclaims_paid", "AutoClaims property claims: PAID ~ rating covariates",
                AutoClaims, PAID ~ STATE + CLASS + GENDER + AGE,
                "insuranceData::AutoClaims", "DeepResearch candidate: closed automobile claim severity, raw paid scale."),
      candidate("insurance_autoclaims_log_paid", "AutoClaims property claims: log1p(PAID) ~ rating covariates",
                AutoClaims, log_paid ~ STATE + CLASS + GENDER + AGE,
                "insuranceData::AutoClaims", "Log-scale control for the automobile claim severity candidate.")
    ))
  }

  if (has_package("boot")) {
    acme <- load_data("boot", "acme")
    motor <- load_data("boot", "motor")
    nuclear <- stats::na.omit(load_data("boot", "nuclear"))
    salinity <- load_data("boot", "salinity")
    calcium <- load_data("boot", "calcium")

    cands <- c(cands, list(
      candidate("boot_acme_returns", "Acme monthly excess returns: acme ~ market",
                acme, acme ~ market,
                "boot::acme", "Financial returns; small n but strong asymmetric residual screen."),
      candidate("boot_motor_accel_poly3", "Motorcycle crash (boot): accel ~ cubic(time)",
                motor, accel ~ times + I(times^2) + I(times^3),
                "boot::motor", "Engineering dynamics; independent source of mcycle-like data."),
      candidate("boot_nuclear_cost", "Nuclear plant cost ~ construction covariates",
                nuclear, cost ~ date + t1 + t2 + cap + pr + ne + ct + bw + cum.n + pt,
                "boot::nuclear", "Industrial construction cost, small n."),
      candidate("boot_salinity", "Water salinity ~ lag + trend + discharge",
                salinity, sal ~ lag + trend + dis,
                "boot::salinity", "Environmental hydrology, small n."),
      candidate("boot_calcium", "Calcium uptake ~ time",
                calcium, cal ~ time,
                "boot::calcium", "Biochemical uptake, small n.")
    ))
  }

  if (has_package("survival")) {
    solder <- load_data("survival", "solder")
    cands <- c(cands, list(
      candidate("survival_solder_skips", "Soldering experiment: skips ~ process factors",
                solder, skips ~ Opening + Solder + Mask + PadType + Panel,
                "survival::solder", "Manufacturing process count response; strong PMM2 screening candidate.")
    ))
  }

  if (has_package("ISLR2")) {
    Credit <- load_data("ISLR2", "Credit")
    Wage <- load_data("ISLR2", "Wage")
    Auto <- load_data("ISLR2", "Auto")
    Smarket <- load_data("ISLR2", "Smarket")
    Weekly <- load_data("ISLR2", "Weekly")
    NYSE <- load_data("ISLR2", "NYSE")
    Carseats <- load_data("ISLR2", "Carseats")

    cands <- c(cands, list(
      candidate("islr_credit_balance", "Credit balance ~ account covariates",
                Credit, Balance ~ Income + Limit + Rating + Cards + Age + Education + Own + Student + Married + Region,
                "ISLR2::Credit", "Credit/finance continuous response with asymmetric residuals."),
      candidate("islr_wage", "Wage ~ demographic and job covariates",
                Wage, wage ~ year + age + maritl + race + education + jobclass + health + health_ins,
                "ISLR2::Wage", "Large social/economic dataset; use with care in manuscript framing."),
      candidate("islr_auto_mpg", "Auto MPG ~ vehicle covariates",
                Auto, mpg ~ cylinders + displacement + horsepower + weight + acceleration + year + origin,
                "ISLR2::Auto", "Automotive engineering/efficiency benchmark."),
      candidate("islr_smarket_today", "Smarket daily return ~ lagged returns + volume",
                Smarket, Today ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume,
                "ISLR2::Smarket", "Financial-return decoy / stress test."),
      candidate("islr_weekly_today", "Weekly S&P return ~ lagged returns + volume",
                Weekly, Today ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume,
                "ISLR2::Weekly", "Financial-return decoy / stress test."),
      candidate("islr_nyse_return", "NYSE DJ return ~ volume/volatility/day",
                NYSE, DJ_return ~ log_volume + log_volatility + day_of_week,
                "ISLR2::NYSE", "Long financial-return decoy; expected heavy tails but not necessarily PMM-favorable."),
      candidate("islr_nyse_volatility", "NYSE log-volatility ~ return/volume/day",
                NYSE, log_volatility ~ DJ_return + log_volume + day_of_week,
                "ISLR2::NYSE", "Financial volatility contrast."),
      candidate("islr_carseats_sales", "Carseats sales ~ store/product covariates",
                Carseats, Sales ~ CompPrice + Income + Advertising + Population + Price + ShelveLoc + Age + Education + Urban + US,
                "ISLR2::Carseats", "Retail sales regression; likely negative control.")
    ))
  }

  c(cands, make_external_candidates())
}

decision_from <- function(method, diag, n) {
  are <- switch(method,
                PMM2 = if (is.finite(diag$g2) && diag$g2 > 0) 1 / diag$g2 else 1,
                PMM3 = if (is.finite(diag$g3) && diag$g3 > 0) 1 / diag$g3 else 1,
                1)
  if (method %in% c("PMM2", "PMM3") && is.finite(are) && are >= 1.10 && n >= 40) {
    return(c(decision = "candidate", expected_are = sprintf("%.3f", are)))
  }
  if (method %in% c("PMM2", "PMM3") && is.finite(are) && are >= 1.03) {
    return(c(decision = "weak-candidate", expected_are = sprintf("%.3f", are)))
  }
  c(decision = "reject", expected_are = sprintf("%.3f", are))
}

screen_one <- function(x, idx) {
  mf <- stats::model.frame(x$formula, data = x$data, na.action = stats::na.omit)
  fit <- stats::lm(x$formula, data = mf)
  eps <- stats::residuals(fit)
  diag <- bootstrap_cumulant_diag(eps, B = B, robust = TRUE, seed = 20260525L + idx)
  method <- dispatch_method(diag)
  dec <- decision_from(method, diag, length(eps))
  data.frame(
    id = x$id,
    label = x$label,
    source = x$source,
    formula = paste(deparse(x$formula), collapse = " "),
    n = length(eps),
    p = length(stats::coef(fit)),
    gamma3 = round(diag$gamma3, 3),
    se_gamma3 = round(diag$se_gamma3, 3),
    gamma4 = round(diag$gamma4, 3),
    se_gamma4 = round(diag$se_gamma4, 3),
    gamma6 = round(diag$gamma6, 3),
    g2 = round(diag$g2, 3),
    g3 = round(diag$g3, 3),
    dispatch = method,
    expected_are = as.numeric(dec[["expected_are"]]),
    decision = dec[["decision"]],
    note = x$note,
    stringsAsFactors = FALSE
  )
}

cat(sprintf("=== C5 residual-cumulant screening (B=%d) ===\n\n", B))
cands <- make_candidates()
res <- do.call(rbind, Map(screen_one, cands, seq_along(cands)))
res <- res[order(match(res$decision, c("candidate", "weak-candidate", "reject")),
                 -res$expected_are), ]
print(res[, c("id", "n", "gamma3", "gamma4", "dispatch", "expected_are", "decision")],
      row.names = FALSE)

outdir <- file.path(root, "shared", "datasets")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
utils::write.csv(res, file.path(outdir, "screening_results.csv"), row.names = FALSE)

cat(sprintf("\nSaved %s\n", file.path(outdir, "screening_results.csv")))
