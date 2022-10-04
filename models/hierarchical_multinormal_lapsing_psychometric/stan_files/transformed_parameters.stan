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

    lprior += normal_lpdf(fA | 0, inv(sqrt(1 - inv(n_levels))));
    lprior += normal_lpdf(sA | 0, inv(sqrt(1 - inv(n_subjects))));
    lprior += normal_lpdf(sB | 0, inv(sqrt(1 - inv(n_subjects))));
    lprior += normal_lpdf(fsA | 0, inv(sqrt(1 - inv(n_levels_sA))));
    lprior += normal_lpdf(mum | 0, 100);
    lprior += gamma_lpdf(muw | 2, 0.5);
    lprior += beta_lpdf(lapse_alpha, lapse_beta);

}