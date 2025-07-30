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
  int<lower=0> C; // cell lines * (replecate of each cell line) >> for example 8*3
  int<lower=0> L; // number of cell lines
  int y [n,C]; // observed cell counts
  int<lower=1> K;                   // Number of groups (number of treatments)
  int<lower=1, upper=K> group[n];    // Group assignment for each observation
}

parameters {
  matrix[K,C] beta;   // treatment scores
  vector<lower = 0>[C] S;  // precision
  vector<lower = 0>[C] psi; // mean of S (precision)
  vector<lower = 0>[C] varphi;// standard deviation of S (precision)
  //vector[K] mu; // matrix of category means
  vector[C] w; // mean of mu (hyperparameter)
  matrix[C,K] mu; // matrix of observation means
 //matrix[C,J]  beta; // matrix of treatment effects
  real<lower=0>  tau; // between categories
  vector<lower=0>[K] sigma; // within categories
  vector<lower=0>[K] nu; // for Student-t noise
}

transformed parameters {
  simplex[C] theta[n]; 
  vector[C] x_beta[n];
  //  XB, for when there is a covariate as a predictor in the model. in our study, we only considered the group (treatment) effect."
  for (j in 1:n) {
//     x_beta[n] = to_vector(beta_0[group[n],] + beta_1[group[n],]*log_budget[n]);
     x_beta[j] = to_vector(beta[group[j],]);
    theta[j] = softmax(x_beta[j]);
    }
}


model {
// priors
     psi ~ lognormal(0, 0.1);
     varphi ~ lognormal(0,0.1);
     S ~ lognormal(psi,varphi);
     w ~ normal(0, 1);
     tau ~ gumbel(0.5, 0.1);
     sigma ~ gumbel(0.5, 0.1);
     nu ~ exponential(1);
  
  for(k in 1:K){

  for (l in 1:C) {
    mu[l,k] ~ normal( w[l], tau);
  }
  //beta[,k] ~ normal(muobs[,k], sigma_rawC[k]); 
 // beta_0[k,] ~ normal(muobs[,k], sigma_rawC[k]); 
   beta[k,] ~ student_t(nu[k], mu[,k], sigma[k]);
  }


for (j in 1:n) {
      y[j] ~  dirichlet_multinomial(S .*theta[j]);
}
}


generated quantities {
  int y_rep[n, C];
  real log_lik[n];
  
  
  // simulations
  for (j in 1:n) {
    y_rep[j] = dirichlet_multinomial_rng(S .*theta[j], sum(y[j]));
    log_lik[j] = dirichlet_multinomial_lpmf(y[j]| S .*theta[j]);
  }
}

