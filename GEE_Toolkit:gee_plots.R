############################################################
# GEE Plotting Toolkit  Dr. Fariborz Aref
# Works with geepack::geeglm models  binary or Gaussian
# Produces  forest plot  diagnostics  marginal effects
# ROC AUC for binomial  cluster influence visualization
############################################################

# 0) Packages
required <- c("geepack","broom","dplyr","ggplot2","stringr","rlang","purrr","readr")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages(invisible(lapply(required, library, character.only = TRUE)))

# Optional for ROC  binomial only
if (!("pROC" %in% rownames(installed.packages()))) {
  install.packages("pROC", repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages(library(pROC))

# 1) Helpers

.ensure_dir <- function(path = "figs") {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  invisible(path)
}

# Tidy robust estimates from geeglm summary
.tidy_gee <- function(model) {
  s <- summary(model)
  out <- tibble::tibble(
    term     = rownames(s$coefficients),
    estimate = s$coefficients[, "Estimate"],
    std_err  = s$coefficients[, "Std.err"],
    p_value  = s$coefficients[, "Pr(>|W|)"]
  )
  out |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::mutate(
      conf_low  = estimate - 1.96 * std_err,
      conf_high = estimate + 1.96 * std_err
    )
}

# Predict with simple CI on the linear predictor
.predict_ci <- function(model, newdata, type = c("link","response"), level = 0.95) {
  type <- match.arg(type)
  X    <- model.matrix(formula(model), data = newdata)
  beta <- coef(model)

  # Robust covariance
  smry <- summary(model)
  V    <- smry$cov.unscaled * smry$scale

  eta     <- as.vector(X %*% beta)
  se_eta  <- sqrt(rowSums((X %*% V) * X))
  z       <- qnorm(1 - (1 - level) / 2)
  eta_lo  <- eta - z * se_eta
  eta_hi  <- eta + z * se_eta

  if (type == "link") {
    return(list(fit = eta, lo = eta_lo, hi = eta_hi))
  } else {
    inv <- model$family$linkinv
    return(list(fit = inv(eta), lo = inv(eta_lo), hi = inv(eta_hi)))
  }
}

# 2) Plots

# 2.1 Coefficient forest plot
gee_forest_plot <- function(model,
                            save_path = "figs/gee_forest.png",
                            title = "GEE coefficients with robust 95 percent CI") {
  .ensure_dir(dirname(save_path))
  df <- .tidy_gee(model) |>
    dplyr::mutate(term = stringr::str_replace_all(term, ":", " × "))

  g <- ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = reorder(term, estimate))) +
    ggplot2::geom_vline(xintercept = 0, linetype = 3) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = conf_low, xmax = conf_high), width = 0.15) +
    ggplot2::labs(x = "Estimate  link scale", y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 13)

  ggplot2::ggsave(save_path, g, width = 7.2, height = 4.8, dpi = 300)
  invisible(g)
}

# 2.2 Residuals vs fitted
gee_resid_plot <- function(model,
                           save_path = "figs/gee_residuals.png",
                           title = "GEE  Pearson residuals vs fitted") {
  .ensure_dir(dirname(save_path))
  df <- data.frame(
    fitted = fitted(model),
    resid  = residuals(model, type = "pearson")
  )
  g <- ggplot2::ggplot(df, ggplot2::aes(x = fitted, y = resid)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 3) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::geom_smooth(method = "loess", se = FALSE) +
    ggplot2::labs(x = "Fitted values", y = "Pearson residuals", title = title) +
    ggplot2::theme_minimal(base_size = 13)

  ggplot2::ggsave(save_path, g, width = 7.2, height = 4.8, dpi = 300)
  invisible(g)
}

# 2.3 Marginal curve for a focal variable
gee_marginal_curve <- function(model, data, var, by = NULL, grid_n = 100,
                               save_path = "figs/gee_marginal.png",
                               title = "Marginal effect curve with 95 percent CI") {
  .ensure_dir(dirname(save_path))
  stopifnot(var %in% names(data))

  fam  <- model$family$family
  link <- model$family$link

  # Reference row at medians or first level
  ref <- data |>
    dplyr::summarise(dplyr::across(dplyr::everything(), \(x) {
      if (is.numeric(x)) stats::median(x, na.rm = TRUE)
      else if (is.factor(x)) levels(x)[1]
      else if (is.character(x)) x[1]
      else x[1]
    }))

  # Grid over var
  if (is.numeric(data[[var]])) {
    grid <- tibble::tibble(!!var := seq(min(data[[var]], na.rm = TRUE),
                                        max(data[[var]], na.rm = TRUE),
                                        length.out = grid_n))
  } else {
    grid <- tibble::tibble(!!var := sort(unique(data[[var]])))
  }

  newdata <- tidyr::crossing(ref, grid)

  # Optional facet by
  if (!is.null(by)) {
    stopifnot(by %in% names(data))
    levs <- if (is.factor(data[[by]])) levels(data[[by]]) else sort(unique(data[[by]]))
    newdata <- tidyr::crossing(newdata, tibble::tibble(!!by := levs))
  }

  preds   <- .predict_ci(model, newdata, type = "response")
  plot_df <- newdata |>
    dplyr::mutate(fit = preds$fit, lo = preds$lo, hi = preds$hi)

  g <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[var]], y = fit)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      x = var,
      y = ifelse(fam == "binomial", "Predicted probability", "Predicted mean"),
      title = title,
      subtitle = paste0("Family  ", fam, "   Link  ", link)
    ) +
    ggplot2::theme_minimal(base_size = 13)

  if (!is.null(by)) g <- g + ggplot2::facet_wrap(stats::as.formula(paste("~", by)))

  ggplot2::ggsave(save_path, g, width = 7.8, height = 5.0, dpi = 300)
  invisible(g)
}

# 2.4 ROC AUC for binomial models
gee_roc_plot <- function(model, data, yvar,
                         save_path = "figs/gee_roc.png",
                         title = "ROC curve  population average GEE") {
  .ensure_dir(dirname(save_path))
  stopifnot(model$family$family == "binomial")
  stopifnot(yvar %in% names(data))

  pr      <- predict(model, type = "response")
  roc_obj <- pROC::roc(response = data[[yvar]], predictor = pr, quiet = TRUE)
  auc_val <- pROC::auc(roc_obj)

  g <- ggplot2::ggplot(data.frame(
    tpr = roc_obj$sensitivities,
    fpr = 1 - roc_obj$specificities
  ), ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_abline(linetype = 3) +
    ggplot2::geom_path() +
    ggplot2::labs(
      x = "False positive rate",
      y = "True positive rate",
      title = title,
      subtitle = paste0("AUC  ", round(as.numeric(auc_val), 3))
    ) +
    ggplot2::theme_minimal(base_size = 13)

  ggplot2::ggsave(save_path, g, width = 6.0, height = 5.6, dpi = 300)
  invisible(g)
}

# 2.5 Cluster influence  leave one cluster out on a target term
gee_cluster_influence <- function(model, data, cluster_id, target_term,
                                  save_path = "figs/gee_cluster_influence.png",
                                  title = "Cluster influence on coefficient  delta when dropped") {
  .ensure_dir(dirname(save_path))
  stopifnot(cluster_id %in% names(data))
  base <- coef(model)[target_term]
  if (is.na(base)) stop("target_term not found in coefficients")

  levs <- unique(data[[cluster_id]])
  deltas <- purrr::map_dbl(levs, function(cl) {
    subdat <- data[data[[cluster_id]] != cl, , drop = FALSE]
    m_sub <- geepack::geeglm(
      formula(model),
      id     = subdat[[cluster_id]],
      data   = subdat,
      family = model$family,
      corstr = model$corstr
    )
    coef(m_sub)[target_term] - base
  })

  df <- tibble::tibble(cluster = levs, delta = deltas)
  g <- ggplot2::ggplot(df, ggplot2::aes(x = reorder(cluster, delta), y = delta)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 3) +
    ggplot2::geom_point(size = 2) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = cluster_id,
      y = expression(Delta ~ "coefficient"),
      title = title,
      subtitle = paste0("Target term  ", target_term)
    ) +
    ggplot2::theme_minimal(base_size = 13)

  ggplot2::ggsave(save_path, g, width = 7.2, height = 6.2, dpi = 300)
  invisible(g)
}

# 3) Small demo  safe to delete in your repo
if (identical(Sys.getenv("RUN_GEE_PLOT_DEMO", "1"), "1")) {
  set.seed(2025)
  n_clusters <- 40
  n_i <- 30
  N <- n_clusters * n_i
  df <- data.frame(
    id  = rep(seq_len(n_clusters), each = n_i),
    x   = rnorm(N, 0, 1),
    z   = sample(c(0,1), N, replace = TRUE),
    age = rnorm(N, 40, 10)
  )
  lp <- -0.6 + 0.9*df$x + 0.5*df$z + 0.02*(df$age - 40)
  p  <- 1 / (1 + exp(-lp))
  df$y <- rbinom(N, 1, p)

  m <- geepack::geeglm(y ~ x + z + age, id = df$id, data = df,
                       family = binomial("logit"), corstr = "exchangeable")

  gee_forest_plot(m, "figs/demo_forest.png")
  gee_resid_plot(m,  "figs/demo_resid.png")
  gee_marginal_curve(m, df, var = "x", by = "z",
                     save_path = "figs/demo_marginal.png",
                     title = "Predicted probability by x  faceted by z")
  gee_roc_plot(m, df, yvar = "y", save_path = "figs/demo_roc.png")
  gee_cluster_influence(m, df, cluster_id = "id",
                        target_term = "x", save_path = "figs/demo_influence.png")
  message("Demo plots saved to figs")
}
