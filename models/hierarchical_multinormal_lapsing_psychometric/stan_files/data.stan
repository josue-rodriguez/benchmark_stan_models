data{
    // brms
    int<lower=1> N;  // total number of observations
    vector[N] Y;  // response variable
    int<lower=1> K;  // number of population-level effects
    matrix[N, K] X;  // population-level design matrix
    // data for group-level effects of ID 1
    int<lower=1> N_1;  // number of grouping levels
    int<lower=1> M_1;  // number of coefficients per level
    int<lower=1> J_1[N];  // grouping indicator per observation
    // group-level predictor values
    vector[N] Z_1_1;
    vector[N] Z_1_2;
    int<lower=1> NC_1;  // number of group-level correlations
    int prior_only;  // should the likelihood be ignored?

    // psychometric
}

