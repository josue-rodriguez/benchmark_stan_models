library(lme4)
library(brms)
library(rstan)
library(dplyr)
library(cmdstanr)

# for restarting sessions within vscodeß
radian_restart <- function(args = NULL) {
  getOption("rchitect.py_tools")$attach()
  os <- import("os")
  sys <- import("sys")
  os$execv(sys$executable, c("sys$executable", "-m", "radian", args))
}

dat_raw <- read.csv("ace_thresholds_data.csv")


# We explored four separate performance outcome measures: 
# performance efficiency (accuracy / mean response time (RT)), 
# consistency (1 - coefficient of variation (CV) where CV is standard deviation of RT/ mean RT),
#  average response window (RW),
#  and psychometric threshold. 

dat <-
    dat_raw |>
    mutate(pid = as.numeric(as.factor(pid)),
        correct = ifelse(correct_button == "correct", 1, 0)) |>
    group_by(pid) |>
    mutate(acc = mean(correct),
        pe = acc / mean(rt), 
        cv = sd(rt) / mean(rt),
        consistency = 1 - cv,
        avg_rw = mean(rw),
        thresh = mean(norm)) |>
        slice_sample(prop = 0.10) |>
    ungroup() |>
    select(pid, rt, timepoint, grade, acc:thresh)



N <- nrow(dat)
J <- length(unique(dat$pid))
y <- dat$rt
X <- model.matrix(~ avg_rw + pe + consistency + thresh, data = dat)

p <- ncol(X)

# bprior <- c(prior(normal(1,2), class = b, coef = pe))
stancode <- make_stancode(
    rt ~ avg_rw + pe + consistency + thresh + timepoint + (timepoint|pid),
    data = dat
)
stancode

standata <- make_standata(
    rt ~ avg_rw + pe + consistency + thresh + timepoint + (pe|pid),
    data = dat
)
str(standata)


path <- 
    c("~/GitHub/benchmark_stan_models/models/hierarchical_multinormal_lapsing_psychometric/stan_files/")
# mod <- stan("brms.stan", data = standata, chains = 4, iter = 5000)
 

# ~/.cmdstan/cmdstan-2.30.1/bin/stanc stan_files/brms.stan --include-paths ./
cmd_mod <- cmdstan_model("stan_files/brms.stan", force_recompile = TRUE, include_paths = "./stan_files")

cmd_fit <- fit <- cmd_mod$sample(
  data = standata, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  refresh = 500 # print update every 500 iters
)

# mod2 <- brm(rt ~ avg_rw + pe + consistency + thresh + timepoint + (pe|pid), data = dat, chains = 4, iter = 5000)
# cmd_mod$include_paths()
