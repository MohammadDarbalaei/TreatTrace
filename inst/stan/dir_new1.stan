// -----------------------------------------------------
// define dirichlet_multinomial for likelihood estimation
// -----------------------------------------------------
functions {

  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
                - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }
  
  int[] dirichlet_multinomial_rng(vector alpha, int N) {
    return multinomial_rng(dirichlet_rng(alpha), N);
  }

}
// -----------------------------------------------------
// -----------------------------------------------------

data {
  int<lower=0> n; // Number of observations (pools or tomur)
  int<lower=0> C; // cell lines * (replicate of each cell line) >> for example 8*3
  int<lower=0> L; // number of cell lines
  int y [n,C]; // observed cell counts
  int<lower=1> K;                   // Number of groups (number of treatments)
  int<lower=1, upper=K> group[n];    // Group assignment for each observation
}

parameters {
  matrix[K,C] beta;   // treatment scores
  matrix<lower = 0>[K,C] S;  // precision (modified to a matrix)
  vector<lower = 0>[C] psi; // mean of S (precision)
  vector<lower = 0>[C] varphi; // standard deviation of S (precision)
  vector[C] w; // mean of mu (hyperparameter)
  matrix[C,K] mu; // matrix of observation means
  real<lower=0> tau; // between categories
  matrix<lower=0>[K,C] sigma; // within categories (modified to a matrix)
  vector<lower=0>[K] nu; // for Student-t noise
}

transformed parameters {
  simplex[C] theta[n]; 
  vector[C] x_beta[n];

  for (j in 1:n) {
     x_beta[j] = to_vector(beta[group[j],]);
    theta[j] = softmax(x_beta[j]);
  }
}

model {
  // priors
  psi ~ lognormal(0, 0.1);
  varphi ~ lognormal(0, 0.1);

  for (k in 1:K) {
    S[k,] ~ lognormal(psi, varphi);
  }

  w ~ normal(0, 1);
  tau ~ gumbel(0.5, 0.1);
  for (k in 1:K) {
    sigma[k,] ~ gumbel(0.5, 0.1);
  }
  nu ~ exponential(1);

  for (k in 1:K) {
    for (l in 1:C) {
      mu[l,k] ~ normal(w[l], tau);
    }
    beta[k,] ~ student_t(nu[k], mu[,k], sigma[k,]);
  }

  for (j in 1:n) {
    vector[C] scaled_theta = to_vector(S[group[j],]) .* to_vector(theta[j]);
    y[j] ~ dirichlet_multinomial(scaled_theta);
  }
}

generated quantities {
  int y_rep[n, C];
  real log_lik[n];

  // simulations
  for (j in 1:n) {
    vector[C] scaled_theta = to_vector(S[group[j],]) .* to_vector(theta[j]);
    y_rep[j] = dirichlet_multinomial_rng(scaled_theta, sum(y[j]));
    log_lik[j] = dirichlet_multinomial_lpmf(y[j] | scaled_theta);
  }
}
