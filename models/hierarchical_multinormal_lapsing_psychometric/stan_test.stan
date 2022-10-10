data {
    int<lower=0> N;

    int<lower=0> p1;
    matrix[N, p1] X1;
    

    int<lower=0> p2;
    matrix[N, p2] X2;
    vector[N] y2;
}

// transformed data{
//     matrix[N, p2+1] X2new;
//     vector[N] y1;
    
//     y1 = X1 * Beta;

//     X2new = append_col(X2, y1);
// }

parameters {
    vector[p1] Beta;
    // real<lower=0> sigma;

    vector[p2 + 1] Beta2;
    real<lower=0> sigma2;
}

// transformed parameters {
    
// }

model {
    vector[N] y1;
    y1 = X1 * Beta;

    vector[N] b2vec = rep_vector(Beta[2], N);
    matrix[N, p2+1] X2new;
    X2new = append_col(X2, b2vec);
    
    Beta ~ normal(0, 10);
    // sigma ~ student_t(3, 0, 10);
    
    Beta2 ~ normal(0, 10);


    vector[N] yhat;
    yhat = X2new * Beta2;

    y2 ~ normal(yhat, sigma2);
}
