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