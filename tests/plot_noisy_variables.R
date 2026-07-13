library(dplyr)
library(tidyr)
library(ggplot2)
source(here::here("tests", "sim_helpers.R"))

# Pivot results_noisy_variables to long format with model and metric columns.
# Each (N, V, noisy, model) combination now has N_REP rows (one per replicate).
pivot_noisy_results <- function(results) {
  results |>
    tidyr::pivot_longer(
      cols      = c(plt_mse, plt_crps, plt_noise,
                    uglt_mse, uglt_crps, uglt_noise,
                    splt_mse, splt_crps, splt_noise),
      names_to  = c("model", "metric"),
      names_sep = "_",
      values_to = "value"
    ) |>
    dplyr::mutate(
      model = dplyr::case_when(
        model == "plt"  ~ "Non-sparse PLT",
        model == "splt" ~ "Sparse PLT",
        model == "uglt" ~ "UGLT"
      ),
      metric = dplyr::case_when(
        metric == "mse"   ~ "Goodness of Fit",
        metric == "crps"  ~ "CRPS",
        metric == "noise" ~ "Noise Ratio"
      ),
      model  = factor(model,  levels = names(MODEL_COLORS)),
      metric = factor(metric, levels = c("CRPS", "Goodness of Fit", "Noise Ratio")),
      N      = factor(N),
      noisy  = ifelse(noisy, "Noisy Variables", "Non-noisy Variables"),
      V      = paste0("V: ", V)
    )
}

# Boxplot across replicates, faceted by noisy condition (columns) and metric (rows), split by V
plot_noisy_variables <- function(results) {
  long <- pivot_noisy_results(results)

  ggplot(long, aes(x = N, y = value, fill = model)) +
    geom_boxplot(
      position     = position_dodge2(width = 0.8, preserve = "single"),
      width        = 0.7,
      outlier.size = 0.8
    ) +
    facet_grid(metric ~ noisy + V, scales = "free_y") +
    scale_fill_manual(values = MODEL_COLORS, name = "Model") +
    labs(x = "Number of Observations (N)", y = NULL) +
    bsfa_boxplot_theme()
}

# Separate plots per noisy condition (to match the two-panel layout in the thesis)
plot_noisy_variables_split <- function(results, noisy_val = FALSE) {
  label <- if (noisy_val) "Noisy Variables" else "Non-noisy Variables"
  long  <- pivot_noisy_results(results) |>
    dplyr::filter(noisy == label)

  ggplot(long, aes(x = N, y = value, fill = model)) +
    geom_boxplot(
      position     = position_dodge2(width = 0.8, preserve = "single"),
      width        = 0.7,
      outlier.size = 0.8
    ) +
    facet_grid(metric ~ V, scales = "free_y") +
    scale_fill_manual(values = MODEL_COLORS, name = "Model") +
    labs(x = "Number of Observations (N)", y = NULL, title = label) +
    bsfa_boxplot_theme()
}

# Combine both panels side by side using patchwork
plot_noisy_variables_patchwork <- function(results) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Install patchwork: install.packages('patchwork')")

  p_non_noisy <- plot_noisy_variables_split(results, noisy_val = FALSE) +
    theme(legend.position = "none")
  p_noisy     <- plot_noisy_variables_split(results, noisy_val = TRUE) +
    theme(legend.position = "none")

  legend <- cowplot::get_legend(
    plot_noisy_variables_split(results, noisy_val = FALSE)
  )

  (p_non_noisy | p_noisy) /
    patchwork::wrap_elements(legend) +
    patchwork::plot_layout(heights = c(10, 1))
}

# Same as pivot_noisy_results(), but the "Noise Ratio" metric is populated
# from noise_error (|true - estimated| noise ratio) instead of the raw
# estimated noise ratio. Requires results_ne (see add_noise_error_noisy_variables()).
pivot_noisy_results_ne <- function(results_ne) {
  results_ne |>
    dplyr::mutate(
      plt_noise  = plt_noise_error,
      uglt_noise = uglt_noise_error,
      splt_noise = splt_noise_error
    ) |>
    pivot_noisy_results()
}

plot_noisy_variables_split_ne <- function(results_ne, noisy_val = FALSE) {
  label <- if (noisy_val) "Noisy Variables" else "Non-noisy Variables"
  long  <- pivot_noisy_results_ne(results_ne) |>
    dplyr::filter(noisy == label)

  ggplot(long, aes(x = N, y = value, fill = model)) +
    geom_boxplot(
      position     = position_dodge2(width = 0.8, preserve = "single"),
      width        = 0.7,
      outlier.size = 0.8
    ) +
    facet_grid(metric ~ V, scales = "free_y") +
    scale_fill_manual(values = MODEL_COLORS, name = "Model") +
    labs(x = "Number of Observations (N)", y = NULL, title = label) +
    bsfa_boxplot_theme()
}

plot_noisy_variables_patchwork_ne <- function(results_ne) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Install patchwork: install.packages('patchwork')")

  p_non_noisy <- plot_noisy_variables_split_ne(results_ne, noisy_val = FALSE) +
    theme(legend.position = "none")
  p_noisy     <- plot_noisy_variables_split_ne(results_ne, noisy_val = TRUE) +
    theme(legend.position = "none")

  legend <- cowplot::get_legend(
    plot_noisy_variables_split_ne(results_ne, noisy_val = FALSE)
  )

  (p_non_noisy | p_noisy) /
    patchwork::wrap_elements(legend) +
    patchwork::plot_layout(heights = c(10, 1))
}

# ---- Noise Error (true vs. estimated signal ratio) --------------------------
# Noise Error = |tr(Lambda Lambda^T)/tr(Omega) - tr(Lambda_est Lambda_est^T)/tr(cov(y))|
# See true_noise_ratio() in sim_helpers.R for the caveat on the estimated term.

add_noise_error_noisy_variables <- function(results, samples) {
  true_ratio <- purrr::map_dbl(samples, ~ true_noise_ratio(.x$Lambda, .x$Sigma))
  results |>
    dplyr::mutate(
      true_ratio       = true_ratio,
      plt_noise_error  = abs(true_ratio - plt_noise),
      uglt_noise_error = abs(true_ratio - uglt_noise),
      splt_noise_error = abs(true_ratio - splt_noise)
    )
}

plot_noisy_variables_noise_error <- function(results_ne) {
  long <- results_ne |>
    tidyr::pivot_longer(
      cols      = c(plt_noise_error, uglt_noise_error, splt_noise_error),
      names_to  = "model",
      values_to = "value"
    ) |>
    dplyr::mutate(
      model = dplyr::case_when(
        model == "plt_noise_error"  ~ "Non-sparse PLT",
        model == "splt_noise_error" ~ "Sparse PLT",
        model == "uglt_noise_error" ~ "UGLT"
      ),
      model = factor(model, levels = names(MODEL_COLORS)),
      N     = factor(N),
      noisy = ifelse(noisy, "Noisy Variables", "Non-noisy Variables"),
      V     = paste0("V: ", V)
    )

  ggplot(long, aes(x = N, y = value, fill = model)) +
    geom_boxplot(
      position     = position_dodge2(width = 0.8, preserve = "single"),
      width        = 0.7,
      outlier.size = 0.8
    ) +
    facet_grid(~ noisy + V, scales = "free_y") +
    scale_fill_manual(values = MODEL_COLORS, name = "Model") +
    labs(x = "Number of Observations (N)", y = "Noise Error") +
    bsfa_boxplot_theme()
}

# ---- Usage -------------------------------------------------------------------

results <- readRDS(here::here("tests", "results_noisy_variables.rds"))

samples_nv <- readRDS(here::here("tests", "sim_noisy_variables_samples.rds"))
results_ne <- add_noise_error_noisy_variables(results, samples_nv)
saveRDS(results_ne, here::here("tests", "results_noisy_variables_noise_error.rds"))

# Single combined plot
p_combined <- plot_noisy_variables(results)
p_combined
save_bsfa_plot(p_combined, "noisy_variables_boxplot.png", width = 10, height = 8)

# Two-panel version matching thesis layout (requires patchwork + cowplot).
# Uses noise_error (|true - estimated| noise ratio) in place of the raw
# estimated noise ratio, still labeled "Noise Ratio" in the plot.
p_patchwork <- plot_noisy_variables_patchwork_ne(results_ne)
p_patchwork
save_bsfa_plot(p_patchwork, "noisy_variables_boxplot_patchwork.png", width = 10, height = 6)

p_ne <- plot_noisy_variables_noise_error(results_ne)
p_ne
save_bsfa_plot(p_ne, "noisy_variables_noise_error_boxplot.png", width = 10, height = 6)

