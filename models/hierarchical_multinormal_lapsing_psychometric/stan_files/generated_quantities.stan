generated quantities {
    // recover actual population-level intercept
    real b_Intercept = Intercept - dot_product(means_X, b);

    // compute unit-level correlations
    vector<lower=-1, upper=1>[NC_1] cor_1;

    corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1); 
    // extract upper diag of correlation matrix
    for (k in 1:M_1) {
      for (j in 1:(k - 1)) {
          cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}