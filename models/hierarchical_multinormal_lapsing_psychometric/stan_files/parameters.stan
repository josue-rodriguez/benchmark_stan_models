parameters {
    vector[Kc] b; // population-level effects
    real Intercept; // temp. intercept for centered predictors
    real<lower=0> sigma; // dispersion parameter
    vector<lower=0>[M_1] sd_1; // group-level standard deviations
    matrix[M_1, N_1] z_1; // standardized group-level effects
    cholesky_factor_corr[M_1] L_1;
}