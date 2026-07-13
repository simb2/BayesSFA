library(dplyr)
library(tidyr)
library(ggplot2)
source(here::here("tests", "sim_helpers.R"))

# Pivot results_general_lambda to long format with model and metric columns.
# Each (V, model) combination now has N_REP rows (one per replicate).
pivot_general_lambda_results <- function(results) {
  results |>
    tidyr::pivot_longer(
      cols      = c(uglt_mse, uglt_crps, uglt_noise, splt_mse, splt_crps, splt_noise),
      names_to  = c("model", "metric"),
      names_sep = "_",
      values_to = "value"
    ) |>
    dplyr::mutate(
      model = dplyr::case_when(
        model == "splt" ~ "Sparse PLT",
        model == "uglt" ~ "UGLT"
      ),
      metric = dplyr::case_when(
        metric == "mse"   ~ "Goodness of Fit",
        metric == "crps"  ~ "CRPS",
        metric == "noise" ~ "Noise Ratio"
      ),
      model  = factor(model,  levels = c("Sparse PLT", "UGLT")),
      metric = factor(metric, levels = c("CRPS", "Goodness of Fit", "Noise Ratio")),
      V      = factor(V, levels = sort(unique(V)))
    )
}

# Boxplot across replicates, faceted by metric, split by V within each panel.
plot_general_lambda <- function(results) {
  long <- pivot_general_lambda_results(results)

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

# ---- Usage -------------------------------------------------------------------

results <- readRDS(here::here("tests", "results_general_lambda.rds"))
p <- plot_general_lambda(results)
p
save_bsfa_plot(p, "general_lambda_boxplot.png")
