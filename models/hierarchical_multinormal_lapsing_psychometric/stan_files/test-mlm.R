library(lme4)
library(brms)
# library(rstan)
library(dplyr)
library(cmdstanr)


## --- BRMS DATA ----
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
        # slice_sample(prop = 0.10) |>
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


## --- Psychometric DATA ----
dat <- read.csv("ace_thresholds_data.csv")
dat <- dat |> filter(timepoint == 1)

dat$condition <- factor(dat$condition)
dat$pid <- factor(dat$pid)

global_window_rescale <- 300

nTotal <- dim(dat)[1]
nCond <- length(unique(dat$condition))
nSubj <- length(unique(dat$pid))
intensity <- dat$rw / global_window_rescale
condition <- as.numeric(dat$condition)
pid <- as.numeric(dat$pid)
correct <- dat$correct_button == "correct"
chance_performance <- 0.5
stanData <- list(
  n_total = nTotal,
  n_levels = nCond,
  n_subjects = nSubj,
  subject = pid,
  intensity = intensity,
  level = condition,
  correct = correct,
  chance_performance = chance_performance
)

stan_data <- c(standata, stanData)


 
# system("~/.cmdstan/cmdstan-2.30.1/bin/stanc stan_files/brms.stan --include-paths ./stan_files")


# cmd_mod <- cmdstan_model(exe_file = "stan_files/brms", force_recompile = TRUE, include_paths = "./stan_files")
cmd_mod <- cmdstan_model("stan_files/brms.stan", force_recompile = TRUE, include_paths = "./stan_files")

cmd_fit <- cmd_mod$sample(
  data = stan_data, 
  seed = 123, 
  chains = 4, 
  parallel_chains = 4,
  iter_sampling = 100,
  refresh = 100 # print update every 500 iters
)

summ <- cmd_fit$summary()
summ
# mod2 <- brm(rt ~ avg_rw + pe + consistency + thresh + timepoint + (pe|pid), data = dat, chains = 4, iter = 5000)
# cmd_mod$include_paths()
