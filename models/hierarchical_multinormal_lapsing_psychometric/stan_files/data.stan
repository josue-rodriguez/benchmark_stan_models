data{
    // brms
    int<lower=1> N;  // total number of observations
    vector[N] Y;  // response variable
    int<lower=1> K;  // number of population-level effects
    matrix[N, K] X;  // population-level design matrix
    // data for group-level effects of ID 1
    int<lower=1> N_1;  // number of grouping levels (i.e., number of unique units)
    int<lower=1> M_1;  // number of coefficients per level
    int<lower=1> J_1[N];  // grouping indicator per observation
    // group-level predictor values
    vector[N] Z_1_1; // predictors for random effects
    vector[N] Z_1_2;
    int<lower=1> NC_1;  // number of group-level correlations
    int prior_only;  // should the likelihood be ignored?

    // psychometric model
    int<lower=1> n_total; // Total number of total trials
    int<lower=1> n_subject; // Number of unique people
    int<lower=2> n_levels; // Number of levels in factor
    real intensity[n_total]; // Intensity of trial stimulus
    int<lower=0> subject[n_total]; // Trial-level person id
    int<lower=0> level[n_total]; // Trial-level level id
    int<lower=0, upper=1> correct[n_total]; // Whether trial response was correct
    real<lower=0, upper=1> chance_performance; // Chance performance for experiment (e.g., 1/n for n-AFC)
}

