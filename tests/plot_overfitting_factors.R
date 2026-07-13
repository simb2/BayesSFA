library(dplyr)
library(ggplot2)
source(here::here("tests", "sim_helpers.R"))

# Boxplot of the recovered factor count (posterior mean r) across replicates,
# grouped by fitting model and faceted by the true data-generating structure.
plot_overfitting_factors <- function(results, q_true) {
  results <- results |>
    dplyr::mutate(
      fit_model = factor(fit_model, levels = c("Sparse PLT", "UGLT")),
      truth     = factor(truth, levels = c("UGLT data", "PLT data"))
    )

  ggplot(results, aes(x = fit_model, y = mean_r, fill = fit_model)) +
    geom_boxplot(width = 0.6, outlier.size = 0.8) +
    geom_hline(yintercept = q_true, linetype = "dashed", color = "grey40") +
    facet_wrap(~truth) +
    scale_fill_manual(values = MODEL_COLORS[c("Sparse PLT", "UGLT")], guide = "none") +
    labs(x = NULL, y = "Posterior mean number of factors (r)") +
    bsfa_boxplot_theme()
}

# ---- Usage -------------------------------------------------------------------

overfit_r_results <- readRDS(here::here("tests", "sim_overfitting_factors_results.rds"))
p <- plot_overfitting_factors(overfit_r_results, q_true = 10)
p
save_bsfa_plot(p, "overfitting_factors_boxplot.png")
