model {
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
        target += normal_id_glm_lpdf(Y | Xc, mu, b, sigma);
    }
    // add priors to target
    target += lprior;
    target += std_normal_lpdf(to_vector(z_1));

    // Psychometric Model

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
        lapse_beta = (1-mean_beta) * betaEta;
        target += bernoulli_lpdf(correct | psi);
}