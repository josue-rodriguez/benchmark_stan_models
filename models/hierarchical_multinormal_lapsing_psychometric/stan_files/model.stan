model {
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

    matrix[N, K + 1] X2;
    X2 = append_col(Xc, threshold);
    target += bernoulli_lpmf(correct | psi);




    

    // likelihood including constants
    if (!prior_only) {
        // initialize linear predictor term
        vector[N] mu = Intercept + rep_vector(0.0, N);
        for (n in 1:N) {
            //add random effect terms to the linear predictor
            mu[n] += (r_1_1[J_1[n]] * Z_1_1[n]) + 
              (r_1_2[J_1[n]] * Z_1_2[n]);
        }
        // normal likelihood for regression w/ with design matrix Xc, intercepts mu, 
        //  coefs b, and residual variance sigma
        target += normal_id_glm_lpdf(Y | X2, mu, b, sigma);
    }
    // add priors to target
    target += lprior;
    target += std_normal_lpdf(to_vector(z_1));

    
}