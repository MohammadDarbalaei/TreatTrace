#' Treatment Resistance with Bubble Heatmap Plot
#'
#' This function calculates directional differences and mean differences within treatments,
#' and generates bubble heatmaps for visualizing resistance patterns.
#'
#' @param model A list containing model parameters, including `condition_count`.
#' @param cdf_data A data frame with columns: `cell` (cell lines), `CDF` (cumulative distribution function values),
#' and `treat` (treatment groups).
#' @param mean A numeric vector of mean values corresponding to the data.
#' @param q_up Upper quantile threshold for optional filtering (default = 0.975).
#' @param q_lo Lower quantile threshold for optional filtering (default = 0.025).
#' @param min_w Minimum bubble size for within-treatment heatmap (default = 1).
#' @param max_w Maximum bubble size for within-treatment heatmap (default = 8).
#' @param min_b Minimum bubble size for treatment-level heatmap (default = 1).
#' @param max_b Maximum bubble size for treatment-level heatmap (default = 8).
#'
#' @return A list containing two ggplot objects:
#' \itemize{
#'   \item `bubble_heatmap_w`: Heatmap showing directional differences between cell lines within treatments.
#'   \item `bubble_heatmap_b`: Heatmap summarizing treatment-level comparisons.
#' }
#'
#' @import ggplot2 dplyr
#' @export
#'
#' @examples
#' # Example usage:
#' model <- list(condition_count = c("Treatment1", "Treatment2"))
#' cdf_data <- data.frame(cell = c("A", "B"), CDF = c(0.1, 0.9), treat = c("Treatment1", "Treatment1"))
#' mean <- c(1.5, 2.5)
#' Tr_resistance_bubble(model, cdf_data, mean)

Tr_resistance_bubble <- function(model, cdf_data, mean, q_up = 0.975, q_lo = 0.025,
                                       min_w = 1, max_w = 8, min_b = 1, max_b = 8) {



  # Extract condition counts from the model
  condition_count <- model$condition_count

  # Combine cdf_data with the mean values
  cdf_mean <- cbind(cdf_data, mean)

  # Initialize an empty data frame for results
  result1 <- data.frame(
    new = character(0),
    new1 = character(0),
    Direction = numeric(0),
    treat = character(0),
    Mean = numeric(0)
  )

  # Iterate through unique treatments to compute pairwise differences
  for (treatment in unique(cdf_mean$treat)) {
    subset_data <- cdf_mean[cdf_mean$treat == treatment, ]

    for (i in 1:nrow(subset_data)) {
      for (j in 1:nrow(subset_data)) {
        if (i == j) {
          # Diagonal elements
          result1 <- rbind(
            result1,
            data.frame(
              new = subset_data$cell[i],
              new1 = subset_data$cell[j],
              Direction = 2 * ((1 - subset_data$CDF[i]) - 0.5),
              treat = treatment,
              Mean = abs(subset_data$mean[i])
            )
          )
        } else{
          # Off-diagonal elements
          result1 <- rbind(
            result1,
            data.frame(
              new = subset_data$cell[i],
              new1 = subset_data$cell[j],
              Direction = 0.5 * ((1 - subset_data$CDF[i]) - (1 - subset_data$CDF[j]) -
                                   (subset_data$CDF[i] - subset_data$CDF[j])),
              treat = treatment,
              Mean = abs(subset_data$mean[i] - subset_data$mean[j])
            )
          )
        }
      }
    }
  }

  # Format the result data frame
  result1$new <- factor(result1$new, levels = unique(result1$new))
  result1$new1 <- factor(result1$new1, levels = unique(result1$new1))
  result1$treat <- factor(result1$treat, levels = unique(condition_count))

  # Sort data frame
  result1 <- result1 %>%
    arrange(treat, new, new1, Direction, Mean)

  # Generate wide-format bubble heatmap
  bubble_heatmap_w <- ggplot(result1, aes(
    x = factor(new, levels = unique(new)),
    y = factor(new1, levels = unique(new1)),
    size = Mean, color = Direction
  )) +
    geom_point() +
    facet_wrap(~ treat) +
    geom_point(shape = 21, stroke = 1.2, color = "black") + 
    scale_radius(range = c(min_w,max_w)) +
    scale_color_gradientn(colors = c(
      "#00FF00",  # bright green
      "#00E57A",  # green-turquoise
      "#00D4A5",  # turquoise
      "#00BFD1",  # bright teal
      "#00AAE6",  # light blue-teal
      "#0094F0",  # lighter sky blue
      "#007FFF",  # sky blue
      "#FFFFFF",  # white (neutral midpoint)
      "#FFFFB3",  # pale yellow
      "#FFF176",  # banana yellow
      "#FFEB3B",  # vivid yellow
      "#FFC107",  # amber
      "#FF9800",  # orange
      "#FF5722",  # red-orange
      "#FF0000"   # red
    ), limits = c(-1, 1)) +
    labs(x = "Cell line", y = "Cell line") +
    theme(axis.text.x = element_text(face = "bold", size = 15, angle = 45, hjust = 1),
          axis.text.y = element_text(face = "bold", size = 15, angle = 45),
          axis.title.x = element_text(face = "bold", size = 18),
          axis.title.y = element_text(face = "bold", size = 18),
          strip.text = element_text(face = "bold", size = 18),
          plot.title = element_text(size = 18, hjust = 0.5, vjust = 0.5),
          plot.margin = margin(20, 20, 40, 40)) +
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 15))

  # Generate treatment-level bubble heatmap
  result_home <- result1[result1$new == result1$new1, ]

  bubble_heatmap_b <-  ggplot(result_home, aes(
    x = treat,
    y = new1,
    size = Mean,
    color = Direction
  )) +
    geom_point() +
    geom_point(shape = 21, stroke = 1.2, color = "black") + 
    scale_radius(range = c(min_b, max_b)) +  # Adjust bubble size
    scale_color_gradientn(colors = c(
      "#00FF00",  # bright green
      "#00E57A",  # green-turquoise
      "#00D4A5",  # turquoise
      "#00BFD1",  # bright teal
      "#00AAE6",  # light blue-teal
      "#0094F0",  # lighter sky blue
      "#007FFF",  # sky blue
      "#FFFFFF",  # white (neutral midpoint)
      "#FFFFB3",  # pale yellow
      "#FFF176",  # banana yellow
      "#FFEB3B",  # vivid yellow
      "#FFC107",  # amber
      "#FF9800",  # orange
      "#FF5722",  # red-orange
      "#FF0000"   # red
    ), limits = c(-1, 1)) +
    labs(x = "Treatment", y = "Cell line") +
    theme(
      axis.text.x = element_text(face = "bold", size = 15, angle = 45, hjust = 1),
      axis.text.y = element_text(face = "bold", size = 15, angle = 45),
      axis.title.x = element_text(face = "bold", size = 18),
      axis.title.y = element_text(face = "bold", size = 18),
      strip.text = element_text(face = "bold", size = 18),
      plot.title = element_text(size = 18, hjust = 0.5, vjust = 0.5),
      plot.margin = margin(20, 20, 40, 40),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 2)  # Adds a bold black box around the plot
    ) +
    theme(
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 15)
    )
  

  # Return the generated heatmaps as a list
  return(list(
    bubble_heatmap_w = bubble_heatmap_w,
    bubble_heatmap_b = bubble_heatmap_b
  ))
}
