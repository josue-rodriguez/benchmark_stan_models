library(rstan)
library(psych)

dat <- read.csv("ace_thresholds_data.csv")
dat <- dat%>%filter(timepoint==1)

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
chance_performance <- 1/2
stanData <- list(n_total=nTotal, 
                 n_levels=nCond, 
                 n_subjects = nSubj, 
                 subject = pid, 
                 intensity=intensity, 
                 level=condition, 
                 correct=correct, 
                 chance_performance=chance_performance)

fit.norm_subj.cond.t1 <- stan(file="thresholds.1.newest.stan", 
                              data=stanData, cores = 4, chains = 4)


#posterior predictive
fit.samples <- extract(fit.norm_subj.cond.t1, 
                       c("mum", "muw", "factor_alpha",
                         "subject_alpha", "interaction_alpha", "subject_beta", "lapse"), permute=TRUE)


grand_mean <- mean(fit.samples$mum)

#condition
cond <- as.data.frame(fit.samples$factor_alpha)
cond_mean <- apply(cond, 2, mean)

#subject
subj <- as.data.frame(fit.samples$subject_alpha)
subj_mean <- apply(subj, 2, mean)

#subject by condition
subj_cond <- as.data.frame(fit.samples$interaction_alpha) #looks like these are the 14 conditions * the 17 subjects.
subj_cond_mean <- apply(subj_cond, 2, mean)

#width
width_mean <- mean(fit.samples$muw)

#width by condition
subj_width <- as.data.frame(fit.samples$subject_beta)
subj_width_mean <- apply(subj_width, 2, mean)

#multiply by number of conditions and subjects, second column is number of subjects
df <- data.frame(grand_mean = rep(grand_mean, 252), cond_mean = rep(cond_mean, each = 18), subj_mean = rep(subj_mean, 14), subj_cond_mean = subj_cond_mean, width_mean, subj_width_mean = rep(subj_width_mean, 14))
df <- df %>% 
  mutate(overall_midpoint = (grand_mean + cond_mean + subj_cond_mean),
         overall_width = width_mean + subj_width_mean,
         threshold = qlogis(.7)*overall_width + overall_midpoint,
         threshold = threshold * 300)
#qlogis(.7, df$overall_midpoint, df$overall_width)

#bringing PIDs and conditions to threshold dataset
pids <- dat%>%
  distinct(pid)%>%unlist()
conditions <- dat%>%
  distinct(condition)%>%unlist()
vals <- data.frame(pid = rep(pids, 14), conditions = rep(conditions, each = 18))%>%
  arrange(conditions)

final <- cbind(vals, df) #double check these values by pulling a couple from the dataframe.
#why am I getting values for pids that don't have values?

#this then should be ready to go into analysis, minus the negative and low values.

n_plot_samps <-100
plot_samps <- sample(x=length(fit.samples$mum), 
                     size=n_plot_samps, replace = FALSE)
xx <- seq(min(stanData$intensity), max(stanData$intensity), length.out=500)

psi <- matrix(NA, n_plot_samps, length(xx))
psi_cond <- array(NA, c(n_plot_samps, nCond, length(xx)))
psi_subj <- array(NA, c(n_plot_samps, nSubj, length(xx)))
for (s in 1:n_plot_samps) { 
  midpoint <- with(fit.samples, mum[plot_samps[s]])
  width <- with(fit.samples, muw[plot_samps[s]] + subject_beta[plot_samps[s]])
  psi[s,] <- chance_performance + (1 - chance_performance) * 
    logistic((xx-midpoint)/width);
  for (cn in 1:nCond) { 
    midpoint_cn <- midpoint + fit.samples$factor_alpha[plot_samps[s],cn]
    psi_cond[s,cn,] <- chance_performance + (1 - chance_performance) * 
      logistic((xx-midpoint_cn)/width);
  }
  for (sn in 1:nSubj) { 
    midpoint_sn <- midpoint + fit.samples$subject_alpha[plot_samps[s],sn]
    lapse_sn <- fit.samples$lapse[plot_samps[s],sn]
    psi_subj[s,sn,] <- chance_performance + 
      (1 - lapse_sn - chance_performance) * 
      logistic((xx-midpoint_sn)/width);
  }
  
}



#################################
## Group psychometric function ##
#################################
matplot(xx*global_window_rescale, t(psi), type='l', lty=1, 
        xlim=c(0,4000), ylim=c(.5,1), col=grey(0, alpha=.1),
        xlab="Response Window", ylab="P(Correct)", 
        main="Posterior Psychometric Functions")
abline(h=.5, lty=2)
for (cn in 1:nCond) { 
  matlines(xx*global_window_rescale, t(psi_cond[,cn,]), lty=1, 
           col=rainbow(nCond,alpha=.1)[cn])
}
legend("bottomright", legend=c("Group", levels(dat$condition)),
       lty=1, col=c("#000000", rainbow(nCond)))


####################################
## Subjects psychometric function ##
####################################
plot_subj <- c(7, 10, 17)
matplot(xx*global_window_rescale, t(psi), type='l', lty=1, 
        xlim=c(0,4000), ylim=c(.5,1), col=grey(0, alpha=.1),
        xlab="Response Window", ylab="P(Correct)", 
        main="Posterior Psychometric Functions")
abline(h=.5, lty=2)
for (sn in plot_subj) { 
  matlines(xx*global_window_rescale, t(psi_subj[,sn,]), lty=1, 
           col=rainbow(nSubj,alpha=.1)[sn])
}
legend("bottomright", legend=c("Group", levels(dat$pid)[plot_subj]),
       lty=1, col=c("#000000", rainbow(nSubj)[plot_subj]))
