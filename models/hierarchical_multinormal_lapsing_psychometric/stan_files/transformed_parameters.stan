transformed parameters {
    matrix[N_1, M_1] r_1; // actual group-level effects
    // using vectors speeds up indexing in loops
    vector[N_1] r_1_1;
    vector[N_1] r_1_2;

    real lprior = 0; // prior contribution to the log posterior

    // compute actual group-level effects
    r_1 = scale_r_cor(z_1, sd_1, L_1);
    r_1_1 = r_1[, 1];
    r_1_2 = r_1[, 2];
    // prior for beta coefs can be added here
    lprior += student_t_lpdf(Intercept | 3, 628, 314.4);
    lprior += student_t_lpdf(sigma | 3, 0, 314.4) 
      - 1 * student_t_lccdf(0 | 3, 0, 314.4);
    lprior += student_t_lpdf(sd_1 | 3, 0, 314.4) 
      - 2 * student_t_lccdf(0 | 3, 0, 314.4); // 2 corresponds to M_1, or number of random effects for the given grouping level
    lprior += lkj_corr_cholesky_lpdf(L_1| 1); // noninformative LKJ prior (i.e., eta = 1)

    // -- Psychometric Model
    vector[n_levels] factor_alpha = append_row(fA, -sum(fA));
    vector[N_1] subject_alpha = append_row(sA, -sum(sA));
    vector[N_1] subject_beta = append_row(sB, -sum(sB));
    
    

    // Uncomment below to allow for interaction between subject and factor
    //  (random slopes model)
    matrix[N_1,n_levels] interaction_alpha;
    interaction_alpha = rep_matrix(0, N_1, n_levels);
    for(sj in 1:(N_1 - 1)) {
      for(l in 1:(n_levels - 1)) {
        interaction_alpha[sj, l] = fsA[(sj - 1) * (n_levels - 1) + l];
      }
      interaction_alpha[sj, n_levels] = -1 * sum(interaction_alpha[sj,]);
    }
    for (l in 1:n_levels) {
      interaction_alpha[N_1,l] = -1 * sum(interaction_alpha[,l]);
    }

    // Psychometric Model
    real pm_m[N_1,n_levels];
    real pm_w[N_1,n_levels];

    vector[N] threshold;
    real width[N];


    
    vector [N] psi;
    for (sj in 1:N_1) {
        for (l in 1:n_levels) {
        //m[sj,l] = mum + factor_alpha[l] + subject_alpha[sj];
        // Comment above and uncomment below for subject/factor interaction
        pm_m[sj,l] = mum + factor_alpha[l] + subject_alpha[sj] + interaction_alpha[sj, l];
        pm_w[sj,l] = exp(muw + subject_beta[sj]);
        }
    } 

    for (tr in 1:N) {
        threshold[tr] = pm_m[J_1[tr],level[tr]];
        //width[tr] = muw;
        width[tr] = pm_w[J_1[tr],level[tr]];
        psi[tr] = chance_performance
                        + (1 - lapse[J_1[tr]] - chance_performance)
                        * inv_logit((intensity[tr] - threshold[tr])/width[tr]);
        }

    lprior += normal_lpdf(fA | 0, inv(sqrt(1 - inv(n_levels))));
    lprior += normal_lpdf(sA | 0, inv(sqrt(1 - inv(N_1))));
    lprior += normal_lpdf(sB | 0, inv(sqrt(1 - inv(N_1))));
    lprior += normal_lpdf(fsA | 0, inv(sqrt(1 - inv(n_levels_sA))));
    lprior += normal_lpdf(mum | 0, 100);
    lprior += gamma_lpdf(muw | 2, 0.5);
    lprior += beta_lpdf(lapse | lapse_alpha, lapse_beta);

}