// -----------------------------------------------------
// Dirichlet-Multinomial Model in Stan
// -----------------------------------------------------

functions {
  // Log-probability mass function for the Dirichlet-Multinomial distribution
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha); // Sum of concentration parameters
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
           - lgamma(alpha_plus + sum(y)) - sum(lgamma(alpha));
  }

  // Random number generator for the Dirichlet-Multinomial distribution
  int[] dirichlet_multinomial_rng(vector alpha, int N) {
    return multinomial_rng(dirichlet_rng(alpha), N);
  }
}

data {
  int<lower=0> n;                     // Number of observations (e.g., pools or tumors)
  int<lower=0> C;                     // Number of categories (cell lines x replicates)
  int<lower=0> L;                     // Number of cell lines
  int y[n, C];                        // Observed cell counts (matrix)
  int<lower=1> K;                     // Number of groups (e.g., treatments)
  int<lower=1, upper=K> group[n];     // Group assignment for each observation
}

parameters {
  // Treatment scores and precision parameters
  matrix[K, C] beta;                  // Treatment scores for each group and category
  matrix<lower=0>[K, C] S;            // Precision parameters for Dirichlet-Multinomial

  // Hyperparameters for the Dirichlet precision
  matrix[K, C] psi;                   // Mean of precision parameters
  matrix<lower=0>[K, C] varphi;       // Std dev of precision parameters

  // Priors for category means
  vector[C] w;                        // Hyperprior mean for category effects
  matrix[K, C] mu;                    // Category mean for each group
  vector<lower=0>[C] tau;             // Variability of category means

  // Priors for treatment effects
  matrix<lower=0>[K, C] sigma;        // Variability of treatment scores
  matrix<lower=0>[K, C] nu;           // Degrees of freedom for Student-t noise
}

transformed parameters {
  // Derived parameters
  simplex[C] theta[n];                // Probabilities for each category
  vector[C] x_beta[n];                // Linear predictor for each observation
  vector[C] S_map[n];                 // Precision mapped to each observation

  for (j in 1:n) {
    x_beta[j] = to_vector(beta[group[j], ]);  // Group-specific linear predictor
    theta[j] = softmax(x_beta[j]);            // Convert linear predictor to probabilities
    S_map[j] = to_vector(S[group[j], ]);      // Map precision for each observation
  }
}

model {
  // Priors for hyperparameters
  w ~ normal(0, 1);                      // Prior for mean of category effects
  tau ~ gumbel(0.5, 0.1);                // Prior for variability of category means

  // Priors for treatment and precision parameters
  for (k in 1:K) {
    psi[k, ] ~ normal(0, 0.1);           // Prior for mean of precision
    varphi[k, ] ~ lognormal(0, 0.01);     // Prior for std dev of precision
    mu[k, ] ~ normal(w, tau);            // Prior for group-specific category means
    nu[k, ] ~ exponential(1);            // Prior for degrees of freedom (Student-t noise)
    sigma[k, ] ~ gumbel(0.5, 0.1);       // Prior for variability of treatment scores

    // Priors for treatment scores and precision
    beta[k, ] ~ student_t(nu[k, ], mu[k, ], sigma[k, ]);
    S[k, ] ~ lognormal(psi[k, ], varphi[k, ]);
  }

  // Likelihood: Dirichlet-Multinomial distribution for observed counts
  for (j in 1:n) {
    y[j] ~ dirichlet_multinomial(S_map[j] .* theta[j]);  // Scaled probabilities
  }
}

generated quantities {
  int y_rep[n, C];                // Simulated replicated counts
  real log_lik[n];                // Log-likelihood for each observation

  for (j in 1:n) {
    y_rep[j] = dirichlet_multinomial_rng(S_map[j] .* theta[j], sum(y[j]));  // Simulate data
    log_lik[j] = dirichlet_multinomial_lpmf(y[j] | S_map[j] .* theta[j]);   // Log-likelihood
  }
}
