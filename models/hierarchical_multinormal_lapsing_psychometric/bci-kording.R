# -------------
# Estimating the probability of a common cause
# -------------

n <- 10000

p_common <- .27
p_indep <- 1 - p_common

C <- rbinom(n, size=1, prob = p_common)


sigma_p <- 1
sigma_v <- 1.5
sigma_a <- 0.5

x_v <- c()
x_a <- c()

for (i in 1:n) {
  # common cause
  if (i == 1) {
    s <- rnorm(1, 0, sigma_p)
    s_v <- s_a <- s
  } else {
    s_v <- rnorm(1, 0, sigma_p)
    s_a <- rnorm(1, 0, sigma_p)
  }
  x_v[i] <- rnorm(1, s_v, sigma_v)
  x_a[i] <- rnorm(1, s_a, sigma_a)
}



ll <- function(c, xv, xa, sv, sa, sp) {
    s2v <- sv^2
    s2a <- sa^2
    s2p <- sp^2
  
  if (c == 1) {
    var_sums <- (s2v * s2a) + (s2v * s2p) + (s2a * s2p)
    left <- 1 / (2 * pi * sqrt(var_sums))
    
    right <- -0.5 * ((((xv - xa)^2 * s2p) + ((xv - 0)^2 * s2a) + ((xa - 0)^2 * s2v)) / var_sums)
  }
  else {
   v_plus_p <- sv^2 + sp^2
   a_plus_p <- sa^2 + sp^2
   left <- 1 / (2 * pi * sqrt(v_plus_p + a_plus_p))
   right <- -0.5 * (((xv - 0)^2 / v_plus_p) + ((xa - 0)^2 / a_plus_p))
  }
    
  return(left * exp(right))
}


ll(c = 2, xv = x_v[i], xa = x_a[i], sv = sigma_v, sa = sigma_a, sp = sigma_p)

p_c1 <- c()


for (i in 1:n) {
  
  ll_common <- ll(c = 1, xv = x_v[i], xa = x_a[i], sv = sigma_v, sa = sigma_a, sp = sigma_p)
  
  ll_indep <- ll(c = 2, xv = x_v[i], xa = x_a[i], sv = sigma_v, sa = sigma_a, sp = sigma_p)
  
  p_c1[i] <- (ll_common * p_common) / ((ll_common * p_common) + (ll_indep * p_indep))
}


hist(p_c1)
mean(p_c1)