#' Proportional Posterior Predictive Check for Tumor Volume or Confluency Data
#'
#' This function generates a posterior predictive check plot for tumor volume or confluency data.
#' It visualizes model uncertainty using a violin plot of posterior predictive samples,
#' scaled proportionally to the number of samples, overlaid with observed values.
#'
#' @param model A list containing the following elements:
#'   - `condition_v`: A vector of treatment group identifiers.
#'   - `fit_V`: A fitted Bayesian model object from which posterior predictive samples can be extracted.
#'   - `V`: Observed data (tumor volume or confluency).
#'   - `n_sam`: Number of posterior predictive samples to draw per treatment group.
#' @param q_up Upper quantile for uncertainty visualization. Default is 0.975.
#' @param q_lo Lower quantile for uncertainty visualization. Default is 0.025.
#'
#' @return A `ggplot` object showing a violin plot of posterior predictive uncertainty per treatment group,
#'         overlaid with observed data, scaled proportionally.
#'
#' @examples
#' \dontrun{
#' # Example not run: requires input list with components like 'condition_v', 'fit_V', etc.
#' model <- list(condition_v = ..., fit_V = ..., V = ..., n_sam = ...)
#' PPC_V(model)
#' }
#'
#' @import ggplot2
#' @export
PPC_V <- function(model, q_up = 0.975, q_lo = 0.025) {
  # Extract required elements from the model
  condition_v <- model$condition_v
  n_sam <- model$n_sam
  V <- model$V
  fit_V <- model$fit_V
  yrep_V <- rstan::extract(fit_V, "ypred")

  if(is.null(model$VT0)){
    uncertainty_V <- yrep_V$ypred
  }else{
    uncertainty_V = yrep_V$ypred/matrix(rep(model$VT0, each = nrow(yrep_V$ypred)),nrow(yrep_V$ypred))
  }

  # Prepare treatment-wise data
  group_V <- condition_v
# Prepare treatment-wise data
group_V <- condition_v
rep_per_treat <-table(factor(group_V, levels = unique(group_V)))

# Define sequences for each treatment group
start <- 1
end <- cumsum(rep_per_treat)
seq_list <- vector("list", length(rep_per_treat))
for (i in seq_along(rep_per_treat)) {
  seq_list[[i]] <- start:end[i]
  start <- end[i] + 1
}

# Extract and sample uncertainty per treatment
unc_V_treat <- unc_V_treat_ppc <- vector("list", length(rep_per_treat))
for (i in seq_along(seq_list)) {
  unc_V_treat[[i]] <- as.vector(uncertainty_V[, seq_list[[i]]])
  unc_V_treat_ppc[[i]] <- as.vector(yrep_V$ypred[, seq_list[[i]]])
}
sampled_elements <- lapply(unc_V_treat, function(x) sample(x, n_sam))
sampled_elements_ppc <- lapply(unc_V_treat_ppc, function(x) sample(x, n_sam))

# Prepare data for plotting
gr_lev <- rep(factor(unique(condition_v), levels = unique(condition_v)), each = n_sam)
approx_data <- data.frame(y = unlist(sampled_elements_ppc), gr_lev = gr_lev)
approx_data_Nor_V <- data.frame(y = unlist(sampled_elements), gr_lev = gr_lev)
data_v <- data.frame(V = V, treat = condition_v)
data_Nor_V <- data.frame(V = V/model$VT0, treat = condition_v)

# Determine plot label based on the data type
y_label <- if (all(V >= 0 & V <= 1)) {
  "Confluency & uncertainty"
} else if (all(V > 0)) {
  "Tumor volume & uncertainty"
} else {
  stop("Error: Data values must either be between 0 and 1 for confluency or strictly positive for tumor volume.")
}


# Determine plot label based on the data type
y_label_Nor <- if (all(V >= 0 & V <= 1)) {
  "Normalized confluency & uncertainty"
} else if (all(V > 0)) {
  "Normalized tumor volume & uncertainty"
} else {
  stop("Error: Data values must either be between 0 and 1 for confluency or strictly positive for tumor volume.")
}
# Generate the PPC plot with proportional scaling
ppc_plot <- ggplot() +
  geom_violin(data = approx_data, aes(x = gr_lev, y = y), fill = "deepskyblue", alpha = 0.7, scale = "area", trim = FALSE) +
  geom_point(data = data_v, aes(x = treat, y = V), color = "red", size = 2) +
  labs(x = "Treatment", y = y_label) +
  scale_fill_manual(values = "deepskyblue", guide = FALSE) +  # Updated color
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),  # Adds a black box around the plot
    axis.text.x = element_text(face = "bold", size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5, vjust = 0.5),
    strip.text = element_text(face = "bold", size = 14),
    plot.margin = margin(20, 20, 20, 20)
  )+
  scale_y_log10()

ppc_plot_Nor_V <- ggplot() +
  geom_violin(data = approx_data_Nor_V, aes(x = gr_lev, y = y), fill = "deepskyblue", alpha = 0.7, scale = "area", trim = FALSE) +
  geom_point(data = data_Nor_V, aes(x = treat, y = V), color = "red", size = 2) +
  labs(x = "Treatment", y = y_label_Nor) +
  scale_fill_manual(values = "deepskyblue", guide = FALSE) +  # Updated color
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),  # Adds a black box around the plot
    axis.text.x = element_text(face = "bold", size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.title.y = element_text(face = "bold", size = 16),
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5, vjust = 0.5),
    strip.text = element_text(face = "bold", size = 14),
    plot.margin = margin(20, 20, 20, 20)
  )+
  scale_y_log10()

  return(list(ppc_plot = ppc_plot,ppc_plot_Nor_V = ppc_plot_Nor_V,approx_data=approx_data, sampled_elements = sampled_elements))
}
