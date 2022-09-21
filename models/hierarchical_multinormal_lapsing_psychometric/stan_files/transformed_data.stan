transformed data {
    int Kc = K - 1;
    matrix[N, Kc] Xc; // centered version of X without an intercept
    vector[Kc] means_X; // column means of X before centering
    
    //2:K because discounting intercept
    for (i in 2:K) {
        means_X[i - 1] = mean(X[, i]);
        Xc[, i - 1] = X[, i] - means_X[i - 1];
    }
}