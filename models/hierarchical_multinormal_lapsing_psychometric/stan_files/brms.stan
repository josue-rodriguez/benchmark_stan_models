functions {
  /* compute correlated group-level effects
   * Args:
   *   z: matrix of unscaled group-level effects
   *   SD: vector of standard deviation parameters
   *   L: cholesky factor correlation matrix
   * Returns:
   *   matrix of scaled group-level effects
   */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
}
data {
  // brms
  int<lower=1> N; // total number of observations
  vector[N] Y; // response variable
  int<lower=1> K; // number of population-level effects
  matrix[N, K] X; // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1; // number of grouping levels (i.e., number of unique units)
  int<lower=1> M_1; // number of coefficients per level
  array[N] int<lower=1> J_1; // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1; // predictors for random effects
  vector[N] Z_1_2;
  int<lower=1> NC_1; // number of group-level correlations
  int prior_only; // should the likelihood be ignored?
  
  // psychometric model
  // int<lower=1> N; // Total number of total trials
  // int<lower=1> N_1; // Number of unique people
  int<lower=2> n_levels; // Number of levels in factor
  array[N] real intensity; // Intensity of trial stimulus
  // int<lower=0> subject[N]; // Trial-level person id
  array[N] int<lower=0> level; // Trial-level level id
  array[N] int<lower=0, upper=1> correct; // Whether trial response was correct
  real<lower=0, upper=1> chance_performance; // Chance performance for experiment (e.g., 1/n for n-AFC)
}
transformed data {
  int Kc = K - 1;
  
  matrix[N, Kc] Xc; // centered version of X without an intercept
  
  vector[Kc] means_X; // column means of X before centering
  
  //2:K because discounting intercept
  for (i in 2 : K) {
    means_X[i - 1] = mean(X[ : , i]);
    Xc[ : , i - 1] = X[ : , i] - means_X[i - 1];
  }
  
  // Psychometric Model
  real<lower=0, upper=1> mean_beta; // mean of beta prior for individual lapse rate
  real<lower=0> betaEta; //Precision parameter of beta prior for individual lapse rate
  
  int df_levels_sA; // Number of degrees of freedom in interaction
  int n_levels_sA; // Number of levels in the interaction
  
  real lapse_alpha;
  real lapse_beta;
  
  df_levels_sA = (n_levels - 1) * (N_1 - 1);
  n_levels_sA = n_levels * N_1;
  
  mean_beta = 0.01;
  betaEta = 100;
  
  lapse_alpha = mean_beta * betaEta;
  lapse_beta = (1 - mean_beta) * betaEta;
}
parameters {
  vector[K] b; // population-level effects
  real Intercept; // temp. intercept for centered predictors
  real<lower=0> sigma; // dispersion parameter
  vector<lower=0>[M_1] sd_1; // group-level standard deviations
  matrix[M_1, N_1] z_1; // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;
  
  // Psychometric Model
  vector<lower=0, upper=1>[N_1] lapse; //Observer's lapse rate
  real mum;
  real<lower=0> muw;
  vector[n_levels - 1] fA;
  vector[N_1 - 1] sA;
  vector[N_1 - 1] sB;
  vector[df_levels_sA] fsA;
}
transformed parameters {
  matrix[N_1, M_1] r_1; // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  
  real lprior = 0; // prior contribution to the log posterior
  
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_1 = r_1[ : , 1];
  r_1_2 = r_1[ : , 2];
  // prior for beta coefs can be added here
  lprior += student_t_lpdf(Intercept | 3, 628, 314.4);
  lprior += student_t_lpdf(sigma | 3, 0, 314.4)
            - 1 * student_t_lccdf(0 | 3, 0, 314.4);
  lprior += student_t_lpdf(sd_1 | 3, 0, 314.4)
            - 2 * student_t_lccdf(0 | 3, 0, 314.4); // 2 corresponds to M_1, or number of random effects for the given grouping level
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1); // noninformative LKJ prior (i.e., eta = 1)
  
  // -- Psychometric Model
  vector[n_levels] factor_alpha = append_row(fA, -sum(fA));
  vector[N_1] subject_alpha = append_row(sA, -sum(sA));
  vector[N_1] subject_beta = append_row(sB, -sum(sB));
  
  // Uncomment below to allow for interaction between subject and factor
  //  (random slopes model)
  matrix[N_1, n_levels] interaction_alpha;
  interaction_alpha = rep_matrix(0, N_1, n_levels);
  for (sj in 1 : (N_1 - 1)) {
    for (l in 1 : (n_levels - 1)) {
      interaction_alpha[sj, l] = fsA[(sj - 1) * (n_levels - 1) + l];
    }
    interaction_alpha[sj, n_levels] = -1 * sum(interaction_alpha[sj,  : ]);
  }
  for (l in 1 : n_levels) {
    interaction_alpha[N_1, l] = -1 * sum(interaction_alpha[ : , l]);
  }
  
  // Psychometric Model
  array[N_1, n_levels] real pm_m;
  array[N_1, n_levels] real pm_w;
  
  vector[N] threshold;
  array[N] real width;
  
  vector[N] psi;
  for (sj in 1 : N_1) {
    for (l in 1 : n_levels) {
      //m[sj,l] = mum + factor_alpha[l] + subject_alpha[sj];
      // Comment above and uncomment below for subject/factor interaction
      pm_m[sj, l] = mum + factor_alpha[l] + subject_alpha[sj]
                    + interaction_alpha[sj, l];
      pm_w[sj, l] = exp(muw + subject_beta[sj]);
    }
  }
  
  for (tr in 1 : N) {
    threshold[tr] = pm_m[J_1[tr], level[tr]];
    //width[tr] = muw;
    width[tr] = pm_w[J_1[tr], level[tr]];
    psi[tr] = chance_performance
              + (1 - lapse[J_1[tr]] - chance_performance)
                * inv_logit((intensity[tr] - threshold[tr]) / width[tr]);
  }
  
  lprior += normal_lpdf(fA | 0, inv(sqrt(1 - inv(n_levels))));
  lprior += normal_lpdf(sA | 0, inv(sqrt(1 - inv(N_1))));
  lprior += normal_lpdf(sB | 0, inv(sqrt(1 - inv(N_1))));
  lprior += normal_lpdf(fsA | 0, inv(sqrt(1 - inv(n_levels_sA))));
  lprior += normal_lpdf(mum | 0, 100);
  lprior += gamma_lpdf(muw | 2, 0.5);
  lprior += beta_lpdf(lapse | lapse_alpha, lapse_beta);
}
model {
  // Append threshold values here for feeding into MLM
  vector[N] thresh_c = threshold - mean(threshold);
  matrix[N, K] X2;
  
  X2 = append_col(Xc, thresh_c);
  
  target += bernoulli_lpmf(correct | psi);
  
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = Intercept + rep_vector(0.0, N);
    for (n in 1 : N) {
      //add random effect terms to the linear predictor
      mu[n] += (r_1_1[J_1[n]] * Z_1_1[n]) + (r_1_2[J_1[n]] * Z_1_2[n]);
    }
    // normal likelihood for regression w/ with design matrix Xc, intercepts mu, 
    //  coefs b, and residual variance sigma
    target += normal_id_glm_lpdf(Y | X2, mu, b, sigma);
  }
  // add priors to target
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
}
generated quantities {
  vector[K] means_X2;
  means_X2 = append_row(means_X, mean(threshold));
  
  // recover actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X2, b);
  
  // compute unit-level correlations
  vector<lower=-1, upper=1>[NC_1] cor_1;
  
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  // extract upper diag of correlation matrix
  for (k in 1 : M_1) {
    for (j in 1 : (k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}


