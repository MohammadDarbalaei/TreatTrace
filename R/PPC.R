#' Posterior Predictive Check (PPC) with Visualization for Dirichlet-Multinomial Models
#'
#' This function performs posterior predictive checks (PPC) for count compositional data
#' modeled using a Dirichlet-multinomial distribution. It generates violin plots to
#' compare observed and predicted values, showcasing model uncertainty and fit.
#'
#' @param model A list containing the results from the `GPI` function, including:
#'   \itemize{
#'     \item \code{fit_dir_mult}: The fitted `rstan` model for Dirichlet-multinomial data.
#'     \item \code{group}: Grouping information for the data.
#'     \item \code{condition_count}: Treatment condition labels for each group.
#'     \item \code{chains}: Number of MCMC chains used in the model.
#'     \item \code{L}: Number of unique cell types in the data.
#'     \item \code{names_cell}: Cell line or type names corresponding to observations.
#'     \item \code{data_count}: Observed count data (matrix or data frame).
#'     \item \code{K}: Number of treatment groups.
#'   }
#' @param q_up Upper quantile for uncertainty visualization (default: 0.975).
#' @param q_lo Lower quantile for uncertainty visualization (default: 0.025).
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{ppc_plot}: A `ggplot2` object displaying the PPC visualization.
#'     \item \code{sampled_fraction}: A list of posterior predictive samples for compositional fractions.
#'   }
#'
#' @examples
#' # Assuming `model` is the result of the `GPI` function:
#' # result <- PPC_count(model)
#' # print(result$ppc_plot)
#'
#' @import ggplot2
#' @export
PPC_count <- function(model, q_up = 0.975, q_lo = 0.025) {
  # Helper function to extract and process posterior samples
  extract_ppc_data <- function(model) {
    # Extract raw data and model dimensions
    group <- model$group
    names_group <- model$condition_count
    n_sam <- model$n_sam
    chains <- model$chains
    fit_dir_mult <- model$fit_dir_mult
    L <- model$L
    names_cell <- as.character(model$names_cell)
    L_cell <- length(names_cell)
    y <- model$data_count
    n <- dim(y)[1]
    K <- model$K



    # Extract posterior predictions and other parameters
    yrep_ab <- rstan::extract(fit_dir_mult, "y_rep")$y_rep
    softmax_fraction <- rstan::extract(fit_dir_mult, pars = "theta")$theta
    kappa <- rstan::extract(fit_dir_mult, pars = "S")$S

    # Normalize predictions using the precision parameter
    for (j in 1:nrow(y)) {
      softmax_fraction[, j, ] <- t(apply(softmax_fraction[, j, ] * kappa, 1, function(row) row / sum(row)))
    }

  #  for (j in 1:nrow(y)) {
  #    softmax_fraction[, j, ] <- t(apply(softmax_fraction[, j, ] * kappa[, j, ], 1, function(row) row / sum(row)))
  #  }
    # Generate sample data for PPC
    final_matrix <- sampled_elements <- sampled_fraction <- fraction_matrix <- list()
    seq_list <- split(1:nrow(y), group)

    for (k in 1:K) {
      unc_volume_treat <- yrep_ab[, seq_list[[k]], ]
      vectors_list <- lapply(unc_volume_treat, as.vector)
      final_matrix[[k]] <- matrix(
        unlist(vectors_list),
        nrow = dim(unc_volume_treat)[1] * length(seq_list[[k]]),
        ncol = L_cell,
        byrow = FALSE
      )
      sampled_elements[[k]] <- c(final_matrix[[k]][sample(1:dim(unc_volume_treat)[1] * length(seq_list[[k]]), n_sam), ])
      fraction <- softmax_fraction[, seq_list[[k]], ]
      fraction <- lapply(fraction, as.vector)
      fraction_matrix[[k]] <- matrix(
        unlist(fraction),
        nrow = dim(unc_volume_treat)[1]* length(seq_list[[k]]),
        ncol = L_cell,
        byrow = FALSE
      )
      sampled_fraction[[k]] <- fraction_matrix[[k]][sample(1:dim(unc_volume_treat)[1]* length(seq_list[[k]]), n_sam), ]
    }

    # Save sampled_fraction in the global environment
    # Prepare data for visualization
    cell_ID <- ave(as.character(model$names_cell), model$names_cell, FUN = function(x) if (length(x) > 1) paste0(x, seq_along(x)) else x)
    cell_bar <- rep(rep(factor(cell_ID,levels =cell_ID ), each = n_sam), length(unique(names_group)))
    treatments <- rep(unique(names_group), each = n_sam * length(cell_ID))

    data_real <- data.frame(
      count = as.vector(t(y)) + 1,
      cell = rep(factor(cell_ID,levels =cell_ID ), n),
      treat = rep(names_group, each = L_cell)	
    )
    data_ppc <- data.frame(
      pred = unlist(sampled_elements) + 1,
      cell = cell_bar,
      treat = treatments
    )

    return(list(data_real = data_real, data_ppc = data_ppc, sampled_fraction = sampled_fraction))
  }

  # Extract PPC data
  ppc_data <- extract_ppc_data(model)
  data_real <- ppc_data$data_real
  data_ppc <- ppc_data$data_ppc
  sampled_fraction <- ppc_data$sampled_fraction


  # Generate the PPC plot
  ppc_plot <- ggplot(data_ppc, aes(x = cell, y = pred, fill = cell)) +
    geom_violin(scale = "width", trim = FALSE) +
    facet_wrap(~ factor(treat, levels = unique(treat))) +
    geom_point(data = data_real, aes(x = cell, y = count)) +
    labs(
      x = "Cell",
      y = "Absolute abundance & uncertainty",
      fill = "Cell"
    ) +
    theme(axis.text.x = element_text(face = "bold", size = 18, angle = 45, hjust = 1),
          axis.text.y = element_text(face = "bold", size = 18),
          axis.title.x = element_text(face = "bold", size = 18),
          axis.title.y = element_text(face = "bold", size = 18),
          strip.text = element_text(face = "bold", size = 18),
          plot.title = element_text(size = 20, hjust = 0.5, vjust = 0.5),
          plot.margin = margin(20, 20, 40, 40),
          legend.title = element_text(size = 18),
          legend.text = element_text(size = 15)) +
    guides(fill = guide_legend(ncol = 1)) +
    scale_y_log10()

  return(list(ppc_plot=ppc_plot,sampled_fraction=sampled_fraction))
}
