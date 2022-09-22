data {  // Everything that must be input by the user
  int<lower=4> n_total;                 // Total number of trials to analyze
  int<lower=2> n_subjects;              // Number of unique observers
  int<lower=2> n_levels;                // Number of levels of Factor
  real intensity[n_total];              // Intensity of the stimulus on each trial
  int<lower=0> subject[n_total];        // Observer on each trial
  int<lower=0> level[n_total];          // Level of Factor on each trial
  int<lower=0,upper=1> correct[n_total];  // Whether the response was correct (1) on each trial
  real<lower=0,upper=1> chance_performance;  // Chance performance for experiment (e.g., 1/n for n-AFC)
}
transformed data {
  real<lower=0, upper=1> mean_beta;     // Mean of beta prior for individual lapse rate
  real<lower=0> betaEta;               // Precision parameter of beta prior for individual lapse rate
  
  // Uncomment below to allow for interaction between subject and factor
  //  (random slopes model)

  int df_levels_sA;       // Number of degrees of freedom in the interaction
  int n_levels_sA;        // Number of levels in the interaction
  
  df_levels_sA = (n_levels - 1) * (n_subjects - 1);
  n_levels_sA = n_levels * n_subjects;

  mean_beta  = 0.01;
  betaEta = 100;
}
parameters {
  vector<lower=0,upper=1>[n_subjects] lapse;        // Observer's lapse rate
  real mum;
  real<lower=0> muw;
  vector[n_levels-1] fA;
  vector[n_subjects-1] sA;
  vector[n_subjects-1] sB;
  vector[(n_subjects - 1) * (n_levels-1)] fsA;
}

transformed parameters {
  vector[n_levels] factor_alpha = append_row(fA, -sum(fA));
  vector[n_subjects] subject_alpha = append_row(sA, -sum(sA));
  vector[n_subjects] subject_beta = append_row(sB, -sum(sB));

  // Uncomment below to allow for interaction between subject and factor
  //  (random slopes model)
  matrix[n_subjects,n_levels] interaction_alpha;
  interaction_alpha = rep_matrix(0, n_subjects, n_levels);
  for(sj in 1:(n_subjects - 1)) {
    for(l in 1:(n_levels - 1)) {
      interaction_alpha[sj, l] = fsA[(sj - 1) * (n_levels - 1) + l];
    }
    interaction_alpha[sj, n_levels] = -1 * sum(interaction_alpha[sj,]);
  }
  for (l in 1:n_levels) {
    interaction_alpha[n_subjects,l] = -1 * sum(interaction_alpha[,l]);
  }
}

model {
  real m[n_subjects,n_levels];
  real w[n_subjects,n_levels];

  real threshold[n_total];
  real width[n_total];
  real lapse_alpha;
  real lapse_beta;
  vector [n_total] psi;
  for (sj in 1:n_subjects) {
    for (l in 1:n_levels) {
      //m[sj,l] = mum + factor_alpha[l] + subject_alpha[sj];
      // Comment above and uncomment below for subject/factor interaction
      m[sj,l] = mum + factor_alpha[l] + subject_alpha[sj] + interaction_alpha[sj, l];
      w[sj,l] = exp(muw + subject_beta[sj]);
    }
  } 

    for (tr in 1:n_total) {
      threshold[tr] = m[subject[tr],level[tr]];
      //width[tr] = muw;
      width[tr] = w[subject[tr],level[tr]];
      psi[tr] = chance_performance
                    + (1 - lapse[subject[tr]] - chance_performance)
                    * inv_logit((intensity[tr]-threshold[tr])/width[tr]);
    }

    lapse_alpha = mean_beta * betaEta;
    lapse_beta = (1-mean_beta) * betaEta ;


  //mean_beta ~ beta(1,60);
  //betaEta ~ gamma(1,0.01);
  fA ~ normal(0, inv(sqrt(1 - inv(n_levels))));
  sA ~ normal(0, inv(sqrt(1 - inv(n_subjects)))); // Produces a standard normal
  
  fsA ~ normal(0, inv(sqrt(1 - inv(n_levels_sA)))); // Produces a standard normal on on interaction_alpha
  sB ~ normal(0, inv(sqrt(1 - inv(n_subjects)))); // Produces a standard normal on on interaction_alpha
  
  mum ~ normal(0, 100);
  muw ~ gamma(2,.5);
  //mum ~ gamma(16,32);
  //muw ~ gamma(16,32);
  lapse ~ beta(lapse_alpha,lapse_beta);
  correct ~ bernoulli(psi);
}
