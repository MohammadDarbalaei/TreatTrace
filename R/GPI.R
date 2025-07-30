#' Generalized Probability Inference (GPI)
#'
#' Conducts Bayesian inference on barcoding data and tumor volume (in vivo) or confluency (in vitro) measurements.
#' Utilizes a Dirichlet-multinomial model for count data, a log-normal likelihood for tumor volume in vivo,
#' and a beta likelihood for confluency in vitro. Allows model selection based on the `dirr` parameter.
#'
#' @param data A list containing the following components:
#'   \itemize{
#'     \item \code{data_count}: A matrix or data frame of observed counts. Rows represent treatment groups, and columns represent cell line replicates.
#'     \item \code{condition_count}: A factor vector specifying treatment groups for \code{data_count}.
#'     \item \code{V}: A numeric vector of either tumor volumes (in vivo) or confluency values (in vitro).
#'     \item \code{condition_v}: A factor vector specifying treatment groups for \code{V}.
#'     \item \code{cell_line}: A factor vector of cell lines corresponding to the replicates in \code{data_count}.
#'     \item \code{cell_rep}: A character vector of replicate identifiers for each column in \code{data_count}.
#'     \item \code{VT0}: Optional numeric vector of baseline tumor volumes or confluency values for normalization.
#'   }
#' @param control_group Control group for comparative analysis.
#' @param dirr Logical flag indicating the choice of model for the Dirichlet-multinomial likelihood:
#'   \itemize{
#'     \item \code{TRUE}: Use the model defined in \code{dirr.stan}.
#'     \item \code{FALSE}: Use the model defined in \code{dir_new.stan}.
#'   }
#'  @param dispersion A list of parameters to control dispersion in the Dirichlet-multinomial model:
#'   \itemize{
#'     \item \code{psi_mean}: Mean for the log-normal prior on the precision parameter \code{psi} (default: 0).
#'     \item \code{psi_sd}: Standard deviation for the log-normal prior on \code{psi} (default: 0.5).
#'     \item \code{varphi_mean}: Mean for the log-normal prior on the variance parameter \code{varphi} (default: 0).
#'     \item \code{varphi_sd}: Standard deviation for the log-normal prior on \code{varphi} (default: 0.01).
#'   }
#' @param control List of control parameters:
#'   \itemize{
#'     \item \code{chains}: Number of MCMC chains (default: 3).
#'     \item \code{iter_count}: Total number of MCMC iterations for barcoding data (default: 1000).
#'     \item \code{iter_V}: Total number of MCMC iterations for tumor volume or confluency data (default: 10000).
#'     \item \code{cores}: Number of CPU cores to use (default: 3).
#'     \item \code{warm}: Logical, indicating whether warm-up iterations should be included (default: FALSE).
#'   }
#' @param q_up Upper quantile threshold (default: 0.975).
#' @param q_lo Lower quantile threshold (default: 0.025).
#' @param n_sam Number of posterior samples (defaults to \code{iter / 2}).
#' @param in_vivo Logical flag indicating whether the data are from an in vivo experiment (default: \code{TRUE}).
#'   If \code{TRUE}, a log-normal likelihood is used for tumor volume. If \code{FALSE}, a beta likelihood is used for confluency.
#' @param ... Other parameters for model fitting.
#' @return A list containing fitted Stan models and relevant metadata:
#'   \itemize{
#'     \item \code{fit_dir_mult}: Fitted Dirichlet-multinomial model.
#'     \item \code{fit_V}: Fitted tumor volume or confluency model.
#'     \item \code{L}: Number of unique cell lines.
#'     \item \code{iter_count}, \code{iter_V}, \code{chains}: MCMC parameters used.
#'     \item \code{data_count}, \code{V}: Input data.
#'     \item \code{group}, \code{group_V}: Numeric group identifiers.
#'     \item \code{condition_count}, \code{condition_v}: Group conditions for count and volume data.
#'     \item \code{control_group}: Control group identifier.
#'     \item \code{n_sam}: Number of posterior samples.
#'   }
#' @export

GPI <- function(data,control_group = NULL, dirr = TRUE,  dispersion = list(psi_mean = 0, psi_sd = 0.5, varphi_mean = 0, varphi_sd = 0.01),
                control = list(chains = 3, iter_count = 1000, iter_V = 10000, cores = 3, warm = FALSE),
                q_up = 0.975, q_lo = 0.025, n_sam = NULL, in_vivo = TRUE, ...) {
  if (is.null(n_sam)) {
    n_sam <- min(control$iter_count, control$iter_V) / 2
  } else {
    if (n_sam > min(control$iter_count, control$iter_V)) {
      n_sam <- min(control$iter_count, control$iter_V)
    }
  }

  if (is.null(control$chains)) control$chains <- 3
  else {
    if (length(control$chains) != 1 || !is.numeric(control$chains) || !is.finite(control$chains) || as.integer(control$chains) <= 0) {
      stop("The mcmc.chains parameter requires a positive integer value exceeding 0")
    }
  }

  if (is.null(control$iter_count)) control$iter_count <- 10000
  if (is.null(control$iter_V)) control$iter_V <- 10000
  if (is.null(control$cores)) control$cores <- FALSE
  else {
    if (length(control$cores) != 1 || !is.numeric(control$cores) || !is.finite(control$cores) || as.integer(control$cores) <= 0) {
      stop("The mcmc.cores parameter requires a positive integer value exceeding 0")
    }
  }

  # Validate dispersion parameters
  dispersion <- modifyList(list(psi_mean = 0, psi_sd = 0.5, varphi_mean = 0, varphi_sd = 0.01), dispersion)


  if (!is.null(data$pre_composition)) {
    # Perform the operation using data_count and pre_composition
    data$data_count <- t(apply(data$data_count, 1, function(row_x, y_mat) {
      i <- which(row_x == row_x[1])
      row_y <- y_mat[i, ]
      rounded <- round(row_x * (row_x / sum(row_x)) / row_y)
      ifelse(rounded == 0, 1, rounded)  # Replace zeros with 1
    }, y_mat = data$pre_composition))
  } else {
    # Print the message if pre_composition is not provided
    message("Each barcode distributed identically among conditions.")
  }
  factor_group <- factor(data$condition_count, levels = unique(data$condition_count))
  group <- as.numeric(factor_group)
  data_count <- data$data_count
  cell_line <- data$cell_line
  K <- length(unique(group))
  L <- length(unique(cell_line))
  n <- dim(data_count)[1]
  VT0 <- data$VT0



  dlist <- list(
    y = data_count,
    group = group,
    K = K,
    C = length(cell_line),
    n = n,
    psi_mean = dispersion$psi_mean,
    psi_sd = dispersion$psi_sd,
    varphi_mean = dispersion$varphi_mean,
    varphi_sd = dispersion$varphi_sd
  )

  # Fit Dirichlet-Multinomial Model Based on dirr Flag
  if (dirr) {
    fit <- rstan::stan(
      file = system.file("stan/dirr.stan", package = "GPIpackage"),
      data = dlist,
      chains = control$chains,
      iter = control$iter_count,
      cores = control$cores,
      control = list(adapt_delta = 0.999, stepsize = 0.5, max_treedepth = 12)
    )
  } else {
    fit <- rstan::stan(
      file = system.file("stan/dir_new2.stan", package = "GPIpackage"),
      data = dlist,
      chains = control$chains,
      iter = control$iter_count,
      cores = control$cores,
      control = list(adapt_delta = 0.999, stepsize = 0.5, max_treedepth = 18)
    )
  }

  # Verify treatment consistency
  if (!any(unique(data$condition_count) %in% unique(data$condition_v))) {
    stop("Error: Treatments in volume data do not match treatments in barcoding data.")
  }

  # Prepare tumor volume or confluency data
  V <- as.vector(data$V[data$condition_v %in% unique(data$condition_count)])
  group_V <- as.numeric(factor(data$condition_v, levels = unique(data$condition_v)))

  dlist_volume <- list(
    V = V,
    group = group_V,
    K = length(unique(data$condition_v)),
    n = length(group_V)
  )

  # Fit Appropriate Model Based on in_vivo Flag
  if (in_vivo) {
    # Error Check: Tumor volume values must be strictly positive
    if (any(V <= 0)) {
      stop("Error: Tumor volume values must be strictly positive for the log-normal likelihood.")
    }
    # Log-normal likelihood for tumor volume (in vivo)
    fit_V <- rstan::stan(
      file = system.file("stan/tumor_volume_lognormal1.stan", package = "GPIpackage"),
      data = dlist_volume,
      chains = control$chains,
      iter = control$iter_V,
      cores = control$cores,
      control = list(max_treedepth = 18)
    )
  } else {
    # Error Check: Confluency values must be between 0 and 1
    if (any(V <= 0 | V >= 1)) {
      stop("Error: Confluency values must be between 0 and 1 for the beta likelihood.")
    }
    # Beta likelihood for confluency (in vitro)
    fit_V <- rstan::stan(
      file = system.file("stan/confluency.stan", package = "GPIpackage"),
      data = dlist_volume,
      chains = control$chains,
      iter = control$iter_V,
      cores = control$cores,
      control = list(adapt_delta = 0.999, stepsize = 0.5, max_treedepth = 12)
    )
  }

  # Return results
  return(list(
    fit_dir_mult = fit,
    fit_V = fit_V,
    L = L,
    iter_count = control$iter_count,
    iter_V = control$iter_V,
    chains = control$chains,
    data_count = data$data_count,
    V = V,
    VT0 = VT0,
    K = K,
    names_cell = cell_line,
    group = group,
    group_V = group_V,
    condition_count = data$condition_count,
    condition_v = data$condition_v,
    control_group = control_group,
    n_sam = n_sam,
    in_vivo = in_vivo
  ))
}
