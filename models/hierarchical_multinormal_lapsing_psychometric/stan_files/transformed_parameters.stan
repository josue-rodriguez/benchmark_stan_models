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
}