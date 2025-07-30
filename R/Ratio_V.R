#' Calculate and Visualize Ratio of Tumor Volumes or Confluency
#'
#' This function calculates the ratio of tumor volumes or confluency for different treatment groups
#' relative to a control group. The ratios are log-transformed and visualized using violin plots
#' to highlight the distribution across conditions.
#'
#' The `sampled_elements` input is derived using a Bayesian approach, where tumor volumes are
#' modeled with a lognormal likelihood and confluency with a beta likelihood. These ratios
#' provide insights into relative differences in treatment effects.
#'
#' @param model A list containing:
#'   \itemize{
#'     \item `condition_v`: A vector of condition labels for each sample.
#'     \item `K`: The number of unique conditions.
#'     \item `n_sam`: The number of samples per condition.
#'   }
#'   The `model` object is generated from the GPI function.
#' @param sampled_elements A named list of sampled tumor volume or confluency data.
#'   Names should correspond to the unique values in `condition_v`, and this data is typically
#'   derived from the `PPC_V` function.
#' @param q_up Numeric. The upper quantile threshold for filtering extreme values (default = 0.975).
#' @param q_lo Numeric. The lower quantile threshold for filtering extreme values (default = 0.025).
#' @param control_group Character. The name of the control group. Defaults to the first unique
#'   condition in `condition_v`.
#'
#' @return A list containing:
#'   - `plot_ratio_v`: A ggplot object visualizing the log-transformed ratios of tumor volumes
#'     or confluency across treatment groups.
#'   - `li_sam_ratio_V`: A list of computed ratios of tumor volumes or confluency.
#'
#' @examples
#' # Example usage
#' model <- list(condition_v = c("Control", "Treatment1", "Treatment2"), K = 3, n_sam = 100)
#' sampled_elements <- list(
#'   Control = rnorm(100, mean = 10),
#'   Treatment1 = rnorm(100, mean = 15),
#'   Treatment2 = rnorm(100, mean = 20)
#' )
#' result <- Ratio_v(model, sampled_elements, control_group = "Control")
#' print(result$plot_ratio_v)
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom stats ecdf quantile
#' @export
Ratio_v <- function(model, sampled_elements, q_up = 0.975, q_lo = 0.025, control_group = NULL) {
  # === Validations and Defaults ===
  if (is.null(sampled_elements)) {
    stop("`sampled_elements` is required but missing. Provide sampled tumor volume data.")
  }

  # Extract model parameters
  condition_v <- model$condition_v
  K <- model$K
  n_sam <- model$n_sam
  in_vivo = model$in_vivo
  names(sampled_elements) <- unique(condition_v)
  threshold <-0
  # Validate or set default control group
  if (is.null(control_group)) {
    control_group <- unique(condition_v)[1]
    message("Control group not specified. Defaulting to the first condition: ", control_group)
  } else if (!(control_group %in% condition_v)) {
    stop("The specified control_group must be present in condition_v.")
  }


  # Calculate ratios
  jj <- 1
  li_sam_ratio_V <- list()
  for (k in 1:K) {
    if (any(names(sampled_elements[k]) %in% control_group)) {
      expmus_ref <- sampled_elements[[k]]
    } else {
      li_sam_ratio_V[[jj]] <- sampled_elements[[k]] / expmus_ref
      jj <- jj + 1
    }
  }
  names(li_sam_ratio_V) <- names(sampled_elements)[-1]

  treatt_ratio_V <- rep(names(li_sam_ratio_V), each = n_sam)
  data_ratio_V <- data.frame(p_ratio = log(unlist(li_sam_ratio_V)), treat = treatt_ratio_V)

  # CDF calculation
  get_cdf_value <- function(x, threshold) {
    ecdf_x <- ecdf(x)
    ecdf_x(threshold)
  }
  cdf_data <- data_ratio_V %>%
    group_by(treat) %>%
    summarise(CDF = get_cdf_value(p_ratio, 0)) %>%
    arrange(CDF)

  # Filter data by quantiles
  filtered_data1 <- data_ratio_V %>%
    group_by(treat) %>%
    mutate(a = ifelse(p_ratio <= quantile(p_ratio, q_up) & p_ratio >= quantile(p_ratio, q_lo), p_ratio, NA)) %>%
    na.omit() %>%
    ungroup()

  # Merge filtered data with CDF
  merged_data <- left_join(filtered_data1, cdf_data, by = "treat")

  plot_ratio_v <-  ggplot(merged_data, aes(x = factor(treat, levels = unique(treat)), y = p_ratio)) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), fill = "deepskyblue", alpha = 0.7, trim = FALSE) +  # Updated color
    geom_hline(yintercept = threshold, linetype = "dashed", color = "blue", size = 2) +
    labs(
      x = "Treatment",
      y = ifelse(in_vivo, "Ratio of tumor volume (log scale)", "Ratio of confluency (log scale)")
    ) +
    scale_fill_manual(values = "deepskyblue", guide = FALSE) +  # Updated color
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2),  # Adds a black box around the plot
      axis.text.x = element_text(face = "bold", size = 18, angle = 45, hjust = 1),
      axis.text.y = element_text(face = "bold", size = 18),
      axis.title.x = element_text(face = "bold", size = 18),
      axis.title.y = element_text(face = "bold", size = 18),
      strip.text = element_text(face = "bold", size = 18),
      plot.title = element_text(size = 18, hjust = 0.5, vjust = 0.5),
      plot.margin = margin(20, 20, 40, 40)
    )
  


  return(list(plot_ratio_v=plot_ratio_v,merged_data=merged_data,li_sam_ratio_V=li_sam_ratio_V))
}

