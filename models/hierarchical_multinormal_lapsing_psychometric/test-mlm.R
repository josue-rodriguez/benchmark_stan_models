library(lme4)
library(cmdstanr)
library(dplyr)

dat_raw <- read.csv('ace_thresholds_data.csv')


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
    ungroup() |>       
    select(pid, avg_rw, pe, consistency, thresh, rt)  


N <- nrow(dat)
J <- length(unique(dat$pid))
y <- dat$rt
X <- model.matrix(~ avg_rw + pe + consistency + thresh, data = dat)