data {
  int<lower=1> N; // number of total data points
  int<lower=1> J; // number of people
  vector[N] y; // outcome
  vector[N] x; // level 1 predictor
  // vector[J] x2; // level 2 predictor
  array[N] int<lower=1, upper=J> id; // person id
}

parameters {
  vector[2] beta; // intercept and slope
  cholesky_factor_corr[2] L_tau; // cholesky for ranefs
  matrix[2, J] z_theta;  // std. normals for ncp??
  vector<lower=0>[2] tau; // RE sd
  real<lower=0> sigma;
}

transformed parameters{
  matrix[2, J] theta; // matrix to store individual random effects
  matrix[2, J] beta_j; // matrix to store individual coefficients
  
  
  for (j in 1:J) {
    // tau[1, i] = exp(beta[3] + beta[4] * x2[J]);
    // tau[2, i] = exp(beta[5] + beta[6] * x2[J]);

    theta[, j] = to_vector(diag_pre_multiply(tau, L_tau) * z_theta[, j]);
    beta_j[, j] = beta + theta[, j];
  }
  
}

model {
  real mu;
  
  // priors
  L_tau ~ lkj_corr_cholesky(1.0);
  to_vector(z_theta) ~ normal(0, 1);
  beta ~ normal(250, 10);
  tau ~ student_t(3, 0, 1);

  
  sigma ~ student_t(3, 25, 5);
  
  // likelihood
  
  for (i in 1:N) {
    mu = beta_j[1, id[i]] + beta_j[2, id[i]]*x[i];
    // sigma = exp(beta[2] + theta[2, id[i]]);
    
    y[i] ~ normal(mu, sigma);
  }
}

generated quantities {
    matrix[2, 2] R;
    
    R = L_tau * L_tau';
  }
  