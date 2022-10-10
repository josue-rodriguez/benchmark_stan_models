transformed data {
    int Kc = K - 1;

    matrix[N, Kc] Xc; // centered version of X without an intercept
    
    vector[Kc] means_X; // column means of X before centering

    //2:K because discounting intercept
    for (i in 2:K) {
        means_X[i - 1] = mean(X[, i]);
        Xc[, i - 1] = X[, i] - means_X[i - 1];
    }


    // Psychometric Model
    real<lower=0, upper=1> mean_beta; // mean of beta prior for individual lapse rate
    real<lower=0> betaEta; //Precision parameter of beta prior for individual lapse rate

    int df_levels_sA; // Number of degrees of freedom in interaction
    int n_levels_sA; // Number of levels in the interaction
    
    real lapse_alpha;
    real lapse_beta;

    df_levels_sA = (n_levels - 1) * (N_1 - 1);
    n_levels_sA = n_levels * N_1;

    mean_beta = 0.01;
    betaEta = 100;

    lapse_alpha = mean_beta * betaEta;
    lapse_beta = (1-mean_beta) * betaEta;
}