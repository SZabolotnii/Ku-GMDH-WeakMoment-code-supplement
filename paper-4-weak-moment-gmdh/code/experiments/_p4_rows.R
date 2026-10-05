# Row tables and data loaders shared by run_claimrow_robustness.R and
# export_allrows_for_ml.R. Source after _p4_rebuilt_helpers.R.

# Loaders copied from the producer scripts (they cannot be sourced without running them).
load_row <- function(id) {
  if (id == "islr_credit_balance") {
    data(Credit, package = "ISLR2")
    d <- stats::na.omit(ISLR2::Credit)
    return(p4_frame_to_xy(d, Balance ~ Income + Limit + Rating + Cards + Age + Education +
                            Own + Student + Married + Region))
  }
  if (id == "mass_boston_medv") {
    # The screened formula: the race-proxy column `black` is excluded.
    data(Boston, package = "MASS")
    return(p4_frame_to_xy(MASS::Boston, medv ~ crim + zn + indus + chas + nox + rm + age + dis +
                            rad + tax + ptratio + lstat))
  }
  if (id == "islr_wage") {
    data(Wage, package = "ISLR2")
    return(p4_frame_to_xy(ISLR2::Wage, wage ~ year + age + maritl + race + education + jobclass +
                            health + health_ins))
  }
  if (id == "airquality_ozone") {
    d <- stats::na.omit(datasets::airquality)
    return(p4_frame_to_xy(d, Ozone ~ Solar.R + Wind + Temp + Month + Day))
  }
  if (id == "mass_cars93_mpg") {
    data(Cars93, package = "MASS")
    d <- stats::na.omit(MASS::Cars93)
    return(p4_frame_to_xy(d, MPG.city ~ EngineSize + Horsepower + RPM + Rev.per.mile +
                            Fuel.tank.capacity + Length + Wheelbase + Width + Weight))
  }
  if (id == "insurance_autoclaims_paid") {
    data(AutoClaims, package = "insuranceData")
    d <- stats::na.omit(get("AutoClaims"))
    d$CLASS <- factor(d$CLASS)
    d$GENDER <- factor(d$GENDER)
    return(p4_frame_to_xy(d, PAID ~ CLASS + GENDER + AGE))
  }
  if (id == "insurance_autobi_loss") {
    data(AutoBi, package = "insuranceData")
    d <- stats::na.omit(get("AutoBi"))
    for (v in c("ATTORNEY", "CLMSEX", "MARITAL", "CLMINSUR", "SEATBELT")) d[[v]] <- factor(d[[v]])
    return(p4_frame_to_xy(d, LOSS ~ ATTORNEY + CLMSEX + MARITAL + CLMINSUR + SEATBELT + CLMAGE))
  }
  p4_ext_load(id)
}

ROWS <- data.frame(
  dataset = c("sru_y1_dynamic", "gas_turbine_co_2015_raw", "gas_turbine_nox_2015_raw",
              "islr_credit_balance", "insurance_autobi_loss", "fremtpl2_severity_raw",
              "sru_y2_static", "concrete"),
  protocol = c("blocked", "blocked", "blocked", "random", "random", "random", "blocked", "random"),
  seed0 = c(76000L, 76000L, 76000L, 73000L, 74000L, 74000L, 76000L, 73000L),
  stringsAsFactors = FALSE
)
# The seven candidate-pool rows the generators ran but Table 2 left out. In
# run_claimrow_robustness.R they run only when named explicitly, so its default run
# reproduces the claim-row file.
POOL <- data.frame(
  dataset = c("mass_boston_medv", "islr_wage", "airquality_ozone", "mass_cars93_mpg",
              "insurance_autoclaims_paid", "fremtpl2_severity_log", "sru_y1_static"),
  protocol = c("random", "random", "random", "random", "random", "random", "blocked"),
  seed0 = c(73000L, 73000L, 73000L, 73000L, 74000L, 74000L, 76000L),
  stringsAsFactors = FALSE
)
