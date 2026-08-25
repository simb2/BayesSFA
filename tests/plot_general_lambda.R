library(dplyr)
library(tidyr)
library(ggplot2)
source(here::here("tests", "sim_helpers.R"))

# Pivot results_general_lambda_noise_error to long format with model and metric
# columns. Each (V, model) combination now has N_REP rows (one per replicate).
# Uses noise *error* (|true_ratio - estimated noise ratio|) rather than the raw
# estimated noise ratio for the third panel.
pivot_general_lambda_results <- function(results_ne) {
  results_ne |>
    tidyr::pivot_longer(
      cols         = c(uglt_mse, uglt_crps, uglt_noise_error, splt_mse, splt_crps, splt_noise_error),
      names_to     = c("model", "metric"),
      names_pattern = "^(uglt|splt)_(.*)$",
      values_to    = "value"
    ) |>
    dplyr::mutate(
      model = dplyr::case_when(
        model == "splt" ~ "Sparse PLT",
        model == "uglt" ~ "UGLT"
      ),
      metric = dplyr::case_when(
        metric == "mse"         ~ "Goodness of Fit",
        metric == "crps"        ~ "CRPS",
        metric == "noise_error" ~ "Noise Error"
      ),
      model  = factor(model,  levels = c("Sparse PLT", "UGLT")),
      metric = factor(metric, levels = c("CRPS", "Goodness of Fit", "Noise Error")),
      V      = factor(V, levels = sort(unique(V)))
    )
}

# Boxplot across replicates, faceted by metric, split by V within each panel.
plot_general_lambda <- function(results_ne) {
  long <- pivot_general_lambda_results(results_ne)

  ggplot(long, aes(x = V, y = value, fill = model)) +
    geom_boxplot(
      position    = position_dodge2(width = 0.8, preserve = "single"),
      width       = 0.7,
      outlier.size = 0.8
    ) +
    facet_wrap(~metric, scales = "free_y") +
    scale_fill_manual(values = MODEL_COLORS[c("Sparse PLT", "UGLT")], name = "Model") +
    labs(x = "Number of Variables (V)", y = NULL) +
    bsfa_boxplot_theme()
}

# ---- Noise Error (true vs. estimated signal ratio) --------------------------
# Noise Error = |tr(Lambda Lambda^T)/tr(Omega) - tr(Lambda_est Lambda_est^T)/tr(cov(y))|
# See true_noise_ratio() in sim_helpers.R for the caveat on the estimated term.

add_noise_error_general_lambda <- function(results, samples) {
  true_ratio <- purrr::map_dbl(samples, ~ true_noise_ratio(.x$Lambda, .x$Sigma))
  results |>
    dplyr::mutate(
      true_ratio        = true_ratio,
      uglt_noise_error  = abs(true_ratio - uglt_noise),
      splt_noise_error  = abs(true_ratio - splt_noise)
    )
}

plot_general_lambda_noise_error <- function(results_ne) {
  long <- results_ne |>
    tidyr::pivot_longer(
      cols      = c(uglt_noise_error, splt_noise_error),
      names_to  = "model",
      values_to = "value"
    ) |>
    dplyr::mutate(
      model = dplyr::case_when(
        model == "splt_noise_error" ~ "Sparse PLT",
        model == "uglt_noise_error" ~ "UGLT"
      ),
      model = factor(model, levels = c("Sparse PLT", "UGLT")),
      V     = factor(V, levels = sort(unique(V)))
    )

  ggplot(long, aes(x = V, y = value, fill = model)) +
    geom_boxplot(
      position     = position_dodge2(width = 0.8, preserve = "single"),
      width        = 0.7,
      outlier.size = 0.8
    ) +
    scale_fill_manual(values = MODEL_COLORS[c("Sparse PLT", "UGLT")], name = "Model") +
    labs(x = "Number of Variables (V)", y = "Noise Error") +
    bsfa_boxplot_theme()
}

# ---- Usage -------------------------------------------------------------------

results    <- readRDS(here::here("tests", "results_general_lambda.rds"))
samples_gl <- readRDS(here::here("tests", "general_lambda_simstudy_samples.rds"))
results_ne <- add_noise_error_general_lambda(results, samples_gl)
saveRDS(results_ne, here::here("tests", "results_general_lambda_noise_error.rds"))

p <- plot_general_lambda(results_ne)
p
save_bsfa_plot(p, "general_lambda_boxplot.png")

p_ne <- plot_general_lambda_noise_error(results_ne)
p_ne
save_bsfa_plot(p_ne, "general_lambda_noise_error_boxplot.png")
