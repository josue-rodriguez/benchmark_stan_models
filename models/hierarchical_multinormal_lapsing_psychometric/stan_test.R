library(dplyr)
library(cmdstanr)


library(dplyr)

n <- 100


# first model
x1 <- rnorm(n)
x2 <- rnorm(n)

X1 <- cbind(1, x1)
B1 <- rbind(1, 1.5)

y1 <- X1 %*% B1 + rnorm(n)

# second model
X2 <- cbind(1, x2, B1)
B2 <- rbind(1, 2, 3)

y2 <-  X2 %*% B2 + rnorm(n)

stan_data <- list(
    N = n,
    p1 = ncol(X1),
    p2 = ncol(X2[, -3]),
    y1 = c(y1),
    y2 = c(y2),
    X1 = X1,
    X2 = X2[, -3]
)

mod <- cmdstan_model("stan_test.stan")

fit <- mod$sample(
    data = stan_data,
    iter_sampling = 5000,
    parallel_chains = 2
)

summ <- fit$summary()
summ |>
    filter(rhat > 1) |>
    print(n = 108)
