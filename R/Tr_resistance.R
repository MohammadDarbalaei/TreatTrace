#' Treatment Resistance Analysis
#'
#' This function calculates treatment resistance from posterior samples, generates violin plots,
#' and summarizes resistance using mean values and CDF statistics.
#'
#' @param model A list containing model information:
#'   - `names_cell`: A vector of cell line names.
#'   - `n_sam`: The number of posterior samples.
#'   - `L`: The number of treatments.
#' @param li_sam_ratio_relative A named list of matrices representing relative ratios for each treatment.
#' @param li_sam_ratio_V A named list of vectors representing treatment-specific values.
#' @param q_up Upper quantile threshold for filtering (default: 0.975).
#' @param q_lo Lower quantile threshold for filtering (default: 0.025).
#' @return A list containing:
#'   - `treatment_resistance`: A ggplot object showing treatment resistance violin plots.
#'   - `cdf_data`: A dataframe summarizing cumulative distribution function (CDF) values.
#'   - `mean`: A vector of mean treatment resistance values.
#' @import ggplot2
#' @import dplyr
#' @export
#' @examples
#' \dontrun{
#'   model <- list(names_cell = c("Cell1", "Cell2"),
#'                 n_sam = 100,
#'                 L = 3)
#'   li_sam_ratio_relative <- list(Treatment1 = matrix(runif(300), nrow = 100),
#'                                 Treatment2 = matrix(runif(300), nrow = 100))
#'   li_sam_ratio_V <- list(Treatment1 = runif(100), Treatment2 = runif(100))
#'   result <- Ratio_Fraction(model, li_sam_ratio_relative, li_sam_ratio_V)
#'   print(result$treatment_resistance)
#' }

Tr_resistance <- function(model, li_sam_ratio_relative, li_sam_ratio_V, q_up = 0.975, q_lo = 0.025) {

  # Error checking
  if (is.null(li_sam_ratio_relative)) {
    stop("li_sam_ratio_relative is not available. Run Ratio_Fraction() first.")
  }
  if (is.null(li_sam_ratio_V)) {
    stop("li_sam_ratio_V is not available. Run Ratio_v() first.")
  }

  # Extract required model elements
  cell <- unique(model$names_cell)
  n_sam <- model$n_sam
  L <- model$L
  listm <- li_sam_ratio_relative
  listv <- li_sam_ratio_V

  result <- list()
  result_qu <- list()

  # Calculate treatment resistance metrics
  for (name in names(listm)) {
    matrix_val <- listm[[name]]
    vector_val <- listv[[name]]
    if (is.matrix(matrix_val) && is.vector(vector_val)) {
      result[[name]] <- matrix_val * matrix(vector_val, nrow = nrow(matrix_val), ncol = ncol(matrix_val))
      result_qu[[name]] <- matrix_val * quantile(vector_val, probs = 0.5)
    }
  }

  result <- Filter(Negate(is.null), result)
  result_qu <- Filter(Negate(is.null), result_qu)
  name_treat_comparison <- intersect(names(listm), names(listv))

  # Combine results for plotting
  all_eff <- unlist(result)
  qu_eff <- unlist(result_qu)
  treatt <- rep(name_treat_comparison, each = L * n_sam)
  celll <- rep(rep(cell, each = n_sam), length(name_treat_comparison))


  data_g <- data.frame(all_eff = log(all_eff), qu_eff = log(qu_eff), treat = treatt, cell = celll)
  data_g$cell <- factor(data_g$cell, levels = unique(data_g$cell))
  data_g$treat <- factor(data_g$treat, levels = unique(data_g$treat))
  # CDF calculation
  get_cdf_value <- function(x, threshold) {
    ecdf_x <- ecdf(x)
    ecdf_x(threshold)
  }
  cdf_data <- data_g %>%
    group_by(cell, treat) %>%
    summarise(CDF = get_cdf_value(all_eff, 0), .groups = "drop")
  cdf_data$CDF <- round(cdf_data$CDF, 2)

  # Filtered data for plotting
  filtered_data1 <- data_g %>%
    group_by(cell, treat) %>%
    mutate(a = ifelse(all_eff <= quantile(all_eff, q_up) & all_eff >= quantile(all_eff, q_lo), all_eff, NA)) %>%
    na.omit() %>%
    ungroup()

  # Mean values per treatment
  mean_data1 <- filtered_data1 %>%
    group_by(cell, treat) %>%
    summarize(mean_value = mean(all_eff), .groups = "drop")

  # Violin plot for treatment resistance
  treatment_resistance <- ggplot(filtered_data1, aes(x = factor(cell, levels = unique(cell)), y = all_eff, fill = factor(cell, levels = unique(cell)))) +
    geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), trim = FALSE) +
    facet_wrap(~ factor(treat, levels = unique(treat))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "blue", size = 2) +  # Example threshold
    labs(x = "Cell line", y = "Treatment resistance (log scale)", fill = "Cell line") +
    theme(axis.text.x = element_text(face = "bold", size = 15, angle = 45, hjust = 1),
          axis.text.y = element_text(face = "bold", size = 15),
          axis.title.x = element_text(face = "bold", size = 18),
          axis.title.y = element_text(face = "bold", size = 18),
          strip.text = element_text(face = "bold", size = 18),
          plot.title = element_text(size = 18, hjust = 0.5, vjust = 0.5),
          plot.margin = margin(20, 20, 40, 40)) +
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 15))

  return(list(treatment_resistance = treatment_resistance, cdf_data = cdf_data, mean = mean_data1$mean_value,
              filtered_data1=filtered_data1))
}
