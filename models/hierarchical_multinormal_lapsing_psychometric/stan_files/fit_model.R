library(lme4)
library(brms)
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
        thresh = mean(norm))


N <- nrow(dat)
J <- length(unique(dat$pid))
y <- dat$rt
X <- model.matrix(~ avg_rw + pe + consistency + thresh, data = dat)

p <- ncol(X)

standata <- make_standata(
    rt ~ avg_rw + pe + consistency + thresh + timepoint + (pe|pid),
    data = dat
)

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


cmd_mod <- cmdstan_model("brms.stan", 
                         force_recompile = TRUE,
                         compile = TRUE)

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