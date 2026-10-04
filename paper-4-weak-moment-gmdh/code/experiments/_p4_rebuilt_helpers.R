# Shared helpers for rebuilt Paper 4 evidence generators.

p4_find_pkg <- function() {
  cur <- normalizePath(getwd(), mustWork = TRUE)
  repeat {
    cand <- file.path(cur, "paper-1-gmdh-pmm", "code")
    if (file.exists(file.path(cand, "DESCRIPTION"))) return(cand)
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("gmdhpmm package not found")
    cur <- parent
  }
}

p4_init <- function() {
  pkg <- p4_find_pkg()
  root <- dirname(dirname(pkg))
  suppressMessages(
    if (requireNamespace("gmdhpmm", quietly = TRUE)) {
      library(gmdhpmm)
    } else {
      pkgload::load_all(pkg, quiet = TRUE)
    }
  )
  out <- file.path(root, "paper-4-weak-moment-gmdh", "results")
  if (!dir.exists(out)) dir.create(out, recursive = TRUE)
  list(pkg = pkg, root = root, results = out)
}

P4 <- p4_init()
P4_METHODS <- c("auto", "auto-weak", "WPMM2", "WPMM3", "LSE", "ridge-LSE", "Huber", "L1")
P4_PMM <- c("auto", "auto-weak", "WPMM2", "WPMM3")
P4_ROBUST <- c("LSE", "ridge-LSE", "Huber", "L1")

p4_trmse <- function(e, p = 0.90) {
  e <- e[is.finite(e)]
  if (!length(e)) return(NA_real_)
  q <- stats::quantile(abs(e), p, names = FALSE, na.rm = TRUE)
  sqrt(mean(e[abs(e) <= q]^2, na.rm = TRUE))
}

p4_rmse <- function(e) {
  e <- e[is.finite(e)]
  if (!length(e)) return(NA_real_)
  sqrt(mean(e^2, na.rm = TRUE))
}

p4_model_x <- function(formula, data) {
  X <- stats::model.matrix(formula, data = data)
  keep <- colnames(X) != "(Intercept)"
  X <- X[, keep, drop = FALSE]
  good <- apply(X, 2, function(z) all(is.finite(z)) && stats::sd(z) > 1e-12)
  X <- X[, good, drop = FALSE]
  storage.mode(X) <- "double"
  X
}

p4_frame_to_xy <- function(data, formula) {
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  list(y = as.numeric(stats::model.response(mf)), X = p4_model_x(formula, mf))
}

p4_ext_load <- function(id) {
  man <- utils::read.csv(file.path(P4$root, "shared/datasets/external/external_candidates.csv"),
                         stringsAsFactors = FALSE)
  row <- man[match(id, man$id), , drop = FALSE]
  if (!nrow(row) || is.na(row$id)) stop("Unknown external candidate: ", id)
  p <- row$path[[1]]
  if (!grepl("^/", p)) p <- file.path(P4$root, p)
  d <- stats::na.omit(utils::read.csv(p, stringsAsFactors = TRUE))
  p4_frame_to_xy(d, stats::as.formula(row$formula[[1]]))
}

p4_std_split <- function(X, train, test) {
  X <- as.matrix(X)
  mu <- colMeans(X[train, , drop = FALSE])
  sg <- apply(X[train, , drop = FALSE], 2, stats::sd)
  sg[!is.finite(sg) | sg <= 1e-12] <- 1
  list(
    train = sweep(sweep(X[train, , drop = FALSE], 2, mu, "-"), 2, sg, "/"),
    test = sweep(sweep(X[test, , drop = FALSE], 2, mu, "-"), 2, sg, "/")
  )
}

p4_method_shares <- function(fit) {
  z <- c(pmm2_share = 0, pmm3_share = 0, wpmm2_share = 0, wpmm3_share = 0)
  if (is.null(fit) || is.null(fit$nodes)) return(z)
  methods <- vapply(Filter(function(nd) !is.null(nd) && identical(nd$type, "model"), fit$nodes),
                    function(nd) nd$method, character(1))
  if (!length(methods)) return(z)
  z["pmm2_share"] <- mean(methods == "PMM2")
  z["pmm3_share"] <- mean(methods == "PMM3")
  z["wpmm2_share"] <- mean(methods == "WPMM2")
  z["wpmm3_share"] <- mean(methods == "WPMM3")
  z
}

p4_ctrl <- function(method, seed, criterion = "MSE", L_max = 3L, F = 6L) {
  gmdh_pmm_control(
    L_max = L_max, F = F, epsilon = -Inf, seed = seed,
    B = if (method == "auto-weak") 60L else 0L,
    force_method = method, criterion = criterion, max_iter = 60L,
    weak_sigma_mult = 2.5, valgate_folds = 4L
  )
}

p4_fit_eval <- function(X, y, train, test, method, seed, criterion = "MSE") {
  xs <- p4_std_split(X, train, test)
  fit <- try(gmdh_pmm(xs$train, y[train], p4_ctrl(method, seed, criterion = criterion)), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(data.frame(method = method, trmse = NA_real_, mae = NA_real_, rmse = NA_real_,
                      bias = NA_real_, layers = NA_integer_, ok = FALSE,
                      as.list(p4_method_shares(NULL)), check.names = FALSE))
  }
  pred <- try(stats::predict(fit, xs$test), silent = TRUE)
  if (inherits(pred, "try-error")) {
    return(data.frame(method = method, trmse = NA_real_, mae = NA_real_, rmse = NA_real_,
                      bias = NA_real_, layers = length(fit$layers), ok = FALSE,
                      as.list(p4_method_shares(fit)), check.names = FALSE))
  }
  e <- as.numeric(pred - y[test])
  data.frame(method = method, trmse = p4_trmse(e), mae = mean(abs(e), na.rm = TRUE),
             rmse = p4_rmse(e), bias = mean(e, na.rm = TRUE),
             layers = length(fit$layers), ok = all(is.finite(e)),
             as.list(p4_method_shares(fit)), check.names = FALSE)
}

p4_random_split <- function(n, seed, test_frac = 0.25) {
  set.seed(seed)
  test <- sort(sample.int(n, max(1L, ceiling(test_frac * n))))
  list(train = setdiff(seq_len(n), test), test = test)
}

p4_blocked_splits <- function(n, folds = 5L, gap = 20L, min_train = 80L, min_test = 30L) {
  fb <- floor(n / folds)
  out <- list()
  for (k in seq_len(folds)) {
    a <- (k - 1L) * fb + 1L
    z <- if (k == folds) n else k * fb
    test <- a:z
    train <- setdiff(seq_len(n), max(1L, a - gap):min(n, z + gap))
    if (length(train) >= min_train && length(test) >= min_test) {
      out[[length(out) + 1L]] <- list(fold = k, train = train, test = test)
    }
  }
  out
}

p4_run_random <- function(dataset_id, dat, reps = 30L, seed0 = 70000L, criterion = "MSE",
                          id_col = "candidate") {
  rows <- list()
  n <- length(dat$y)
  for (rep in seq_len(reps)) {
    sp <- p4_random_split(n, seed0 + rep)
    for (j in seq_along(P4_METHODS)) {
      method <- P4_METHODS[[j]]
      ev <- p4_fit_eval(dat$X, dat$y, sp$train, sp$test, method,
                        seed = seed0 + rep * 100L + j, criterion = criterion)
      rows[[length(rows) + 1L]] <- cbind(setNames(data.frame(dataset_id, rep), c(id_col, "rep")),
                                         ev, row.names = NULL)
    }
  }
  do.call(rbind, rows)
}

p4_run_protocols <- function(dataset_id, dat, reps = 30L, seed0 = 81000L, cap_blocked = 3000L,
                             cap_random = 4000L, gap = 20L, criterion = "MSE") {
  rows <- list()
  Xr <- dat$X
  yr <- dat$y
  if (!is.null(cap_random) && length(yr) > cap_random) {
    Xr <- Xr[seq_len(cap_random), , drop = FALSE]
    yr <- yr[seq_len(cap_random)]
  }
  for (rep in seq_len(reps)) {
    sp <- p4_random_split(length(yr), seed0 + rep)
    for (j in seq_along(P4_METHODS)) {
      method <- P4_METHODS[[j]]
      ev <- p4_fit_eval(Xr, yr, sp$train, sp$test, method,
                        seed = seed0 + rep * 100L + j, criterion = criterion)
      rows[[length(rows) + 1L]] <- cbind(data.frame(dataset = dataset_id, protocol = "random",
                                                    fold = rep),
                                         ev[, c("method", "trmse", "mae", "rmse", "layers", "ok")],
                                         row.names = NULL)
    }
  }
  Xb <- dat$X
  yb <- dat$y
  if (!is.null(cap_blocked) && length(yb) > cap_blocked) {
    Xb <- Xb[seq_len(cap_blocked), , drop = FALSE]
    yb <- yb[seq_len(cap_blocked)]
  }
  for (sp in p4_blocked_splits(length(yb), gap = gap)) {
    for (j in seq_along(P4_METHODS)) {
      method <- P4_METHODS[[j]]
      ev <- p4_fit_eval(Xb, yb, sp$train, sp$test, method,
                        seed = seed0 + 5000L + sp$fold * 100L + j, criterion = criterion)
      rows[[length(rows) + 1L]] <- cbind(data.frame(dataset = dataset_id, protocol = "blocked",
                                                    fold = sp$fold),
                                         ev[, c("method", "trmse", "mae", "rmse", "layers", "ok")],
                                         row.names = NULL)
    }
  }
  do.call(rbind, rows)
}

p4_permethod <- function(raw, id_col = "candidate", family_upper = FALSE) {
  ids <- unique(raw[[id_col]])
  rows <- list()
  for (id in ids) {
    for (m in P4_METHODS) {
      sub <- raw[raw[[id_col]] == id & raw$method == m & raw$ok, , drop = FALSE]
      fam <- if (m %in% P4_PMM) "pmm" else "robust"
      if (family_upper && fam == "pmm") fam <- "PMM"
      rows[[length(rows) + 1L]] <- data.frame(
        x_id = id, method = m, family = fam,
        trmse_med = stats::median(sub$trmse, na.rm = TRUE),
        trmse_iqr = stats::IQR(sub$trmse, na.rm = TRUE),
        mae_med = stats::median(sub$mae, na.rm = TRUE),
        mae_iqr = stats::IQR(sub$mae, na.rm = TRUE),
        rmse_med = stats::median(sub$rmse, na.rm = TRUE),
        bias_med = if ("bias" %in% names(sub)) stats::median(sub$bias, na.rm = TRUE) else NA_real_,
        layers_med = if ("layers" %in% names(sub)) stats::median(sub$layers, na.rm = TRUE) else NA_real_,
        n_reps = nrow(sub),
        check.names = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  names(out)[names(out) == "x_id"] <- id_col
  out
}

p4_compare_one <- function(raw, id_col = "candidate", extra_filter = rep(TRUE, nrow(raw))) {
  sub <- raw[extra_filter & raw$ok, , drop = FALSE]
  med <- stats::aggregate(cbind(trmse, mae) ~ method, data = sub, FUN = stats::median)
  getv <- function(methods, metric) {
    mm <- med[med$method %in% methods, c("method", metric), drop = FALSE]
    mm[which.min(mm[[metric]]), , drop = FALSE]
  }
  bp_tr <- getv(P4_PMM, "trmse")
  br_tr <- getv(P4_ROBUST, "trmse")
  bp_ma <- getv(P4_PMM, "mae")
  br_ma <- getv(P4_ROBUST, "mae")
  lse <- med$trmse[match("LSE", med$method)]
  units <- sort(unique(if ("rep" %in% names(sub)) sub$rep else sub$fold))
  pmm_u <- robust_u <- numeric(0)
  for (u in units) {
    su <- if ("rep" %in% names(sub)) sub[sub$rep == u, , drop = FALSE] else sub[sub$fold == u, , drop = FALSE]
    if (!all(P4_PMM %in% su$method) || !all(P4_ROBUST %in% su$method)) next
    pmm_u <- c(pmm_u, min(su$trmse[su$method %in% P4_PMM], na.rm = TRUE))
    robust_u <- c(robust_u, min(su$trmse[su$method %in% P4_ROBUST], na.rm = TRUE))
  }
  p <- tryCatch(stats::wilcox.test(pmm_u, robust_u, paired = TRUE, alternative = "less")$p.value,
                error = function(e) NA_real_)
  list(
    best_pmm_tr_method = bp_tr$method[[1]], best_robust_tr_method = br_tr$method[[1]],
    best_pmm_mae_method = bp_ma$method[[1]], best_robust_mae_method = br_ma$method[[1]],
    pmm_trmse = bp_tr$trmse[[1]], robust_trmse = br_tr$trmse[[1]],
    pmm_mae = bp_ma$mae[[1]], robust_mae = br_ma$mae[[1]], lse_trmse = lse,
    pmm_vs_robust_trmse_pct = 100 * (bp_tr$trmse[[1]] / br_tr$trmse[[1]] - 1),
    pmm_vs_robust_mae_pct = 100 * (bp_ma$mae[[1]] / br_ma$mae[[1]] - 1),
    pmm_vs_lse_pct = 100 * (bp_tr$trmse[[1]] / lse - 1),
    win_rate_trmse = mean(pmm_u < robust_u, na.rm = TRUE),
    p_wilcox_optimistic = p,
    unit_diff_pct = 100 * (pmm_u / robust_u - 1)
  )
}

p4_verdict <- function(pct, p, win_rate) {
  if (is.finite(pct) && pct < -2 && is.finite(p) && p < 0.05 && is.finite(win_rate) && win_rate >= 0.6) {
    "genuine-positive"
  } else if (is.finite(pct) && pct < -1.5 && is.finite(win_rate) && win_rate >= 0.55) {
    "marginal"
  } else if (is.finite(pct) && pct > 2) {
    "loses"
  } else {
    "tie"
  }
}

p4_write <- function(x, stem, suffix = "") {
  utils::write.csv(x, file.path(P4$results, paste0(stem, suffix, ".csv")), row.names = FALSE)
}
