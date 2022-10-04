parameters {
    vector[Kc] b; // population-level effects
    real Intercept; // temp. intercept for centered predictors
    real<lower=0> sigma; // dispersion parameter
    vector<lower=0>[M_1] sd_1; // group-level standard deviations
    matrix[M_1, N_1] z_1; // standardized group-level effects
    cholesky_factor_corr[M_1] L_1;

    // Psychometric Model
    vector<lower=0, upper=1>[n_subjects] lapse; //Observer's lapse rate
    real mum;
    real<lower=0> muw;
    vector[n_levels - 1] fA;
    vector[n_subjects - 1] sA;
    vector[n_subjects - 1] sB;
    vector[df_levels_sA] fsA;

}