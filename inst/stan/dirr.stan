// -----------------------------------------------------
// Define Dirichlet-Multinomial for likelihood estimation
// -----------------------------------------------------
functions {
  // Log-probability mass function for the Dirichlet-Multinomial distribution
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha); // Sum of concentration parameters
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y))) 
           - lgamma(alpha_plus + sum(y)) - sum(lgamma(alpha));
  }

  // Random number generator for Dirichlet-Multinomial distribution
  int[] dirichlet_multinomial_rng(vector alpha, int N) {
    return multinomial_rng(dirichlet_rng(alpha), N);
  }
}
// -----------------------------------------------------
// -----------------------------------------------------

data {
  int<lower=0> n;                     // Number of observations (e.g., pools or tumors)
  int<lower=0> C;                     // Number of cell lines multiplied by replicates per cell line
  int<lower=0> L;                     // Number of cell lines
  int y[n, C];                        // Observed cell counts
  int<lower=1> K;                     // Number of groups (treatments)
  int<lower=1, upper=K> group[n];     // Group assignment for each observation
}

parameters {
  matrix[K, C] beta;                  // Treatment scores for each group and category
  matrix<lower=0>[K, C] S;            // Precision parameter for the Dirichlet-Multinomial (one for each group and category)
  vector[C] psi;                      // Mean for the precision (S)
  vector<lower=0>[C] varphi;          // Standard deviation for the precision (S)
  vector[C] w;                        // Mean of category effects (hyperparameter for mu)
  matrix[C, K] mu;                    // Matrix of category means
  real<lower=0> tau;                  // Variability of category means
  vector<lower=0>[K] sigma;           // Variability of treatment scores
  vector<lower=0>[K] nu;              // Degrees of freedom for Student-t noise
}

transformed parameters {
  simplex[C] theta[n];                // Probabilities (softmax) for each category
  vector[C] x_beta[n];                // Linear combination for each observation

  // Calculate linear combination and convert to probabilities (softmax) for each observation
  for (j in 1:n) {
    x_beta[j] = to_vector(beta[group[j],]);  // Extract group-specific beta coefficients
    theta[j] = softmax(x_beta[j]);           // Apply softmax to convert to probabilities
  }
}

model {
  // Priors for precision-related parameters
  psi ~ normal(0, 0.1);
  varphi ~ lognormal(0, 0.1);
  for (k in 1:K) {
    S[k, ] ~ lognormal(psi, varphi); // Group-specific precisions for each category
  }

  // Priors for mean parameters
  w ~ normal(0, 1);
  tau ~ gumbel(0.5, 0.1);
  sigma ~ gumbel(0.5, 0.1);
  nu ~ exponential(1);

  // Prior distributions for category means and treatment effects
  for (k in 1:K) {
    for (l in 1:C) {
      mu[l, k] ~ normal(w[l], tau);                     // Prior for each category mean
    }
    beta[k,] ~ student_t(nu[k], mu[, k], sigma[k]);     // Prior for treatment effects with Student-t noise
  }

  // Likelihood function: Dirichlet-Multinomial distribution for observed counts
  for (j in 1:n) {
    y[j] ~ dirichlet_multinomial(to_vector(S[group[j], ]) .* theta[j]);
  }
}

generated quantities {
  int y_rep[n, C];                // Replicated (simulated) counts
  real log_lik[n];                // Log-likelihood for each observation

  // Simulate replicated data and compute log-likelihood
  for (j in 1:n) {
    y_rep[j] = dirichlet_multinomial_rng(to_vector(S[group[j], ]) .* theta[j], sum(y[j]));  // Generate replicated data
    log_lik[j] = dirichlet_multinomial_lpmf(y[j] | to_vector(S[group[j], ]) .* theta[j]);   // Compute log-likelihood
  }
}
