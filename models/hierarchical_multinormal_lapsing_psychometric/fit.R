#install.packages("install_if_missing")

# Glossary ----

# n.b. the following employs a mix of snake_case and camelCase that is sure to
#  vex some, but represents the author"s best attempt to balance to the competing
#  aims of clarity & brevity.

# Y: observed outcome
# nY: number of observed outcomes
# X: predictor/contrast matrix
# nX: number of predictors (columns in the contrast matrix)
# rX: number of rows in the contrast matrix X
# (i)ndividual: a unit of observation within which correlated measurements may take place
# (c)ondition: a labelled set of observations within an individual that share some feature/predictor or conjunction of features/predictors
# Xc: condition-level contrast matrix
# nXc: number of predictors in the condition-level contrast matrix
# rXc: number of rows in the condition-level contrast matrix
# yXc: for each observation in y, an index indicating the associated row in Xc corresponding to that observation"s individual/condition combo
# (g)roup: a collection of individuals that share some feature/predictor
# Xg: group-level contrast matrix
# nXg: number of predictors in the group-level contrast matrix
# rXg: number of rows in the group-level contrast matrix
# Z: matrix of coefficient row-vectors to be dot-product"d with a contrast matrix
# iZc: matrix of coefficient row-vectors associated with each individual

# Preamble (options, installs, imports & custom functions) ----

options(warn = 1) #really should be default in R
`%!in%` = Negate(`%in%`) #should be in base R!

# specify the packages used:
required_packages <- c(
	"github.com/stan-dev/cmdstanr", #for Stan stuff
    "github.com/mike-lawrence/aria/aria", # for aria
    "tidyverse" #for all that is good and holy
)

# load the helper functions:
for (file in fs::dir_ls("r")) {
	cat("Loading function: ", fs::path_ext_remove(fs::path_file(file)), "()\n", sep="")
	source(file)
}

#install any required packages not already present
install_if_missing(required_packages)

# load tidyverse & aria
library(tidyverse)
library(magrittr)
library(aria)
library(cmdstanr)
#renv::install_github("mike-lawrence/aria/aria")
#remotes::install_github("mike-lawrence/aria/aria")
#install.packages("~/Desktop/benchmark_stan_models-main/aria_0.1.0.tar", repos = NULL, type = "source")

dat <-
(
	read_csv("ace_thresholds_data.csv")  |>
	separate(
		condition
		, into = c("task", "condition", "setsize")
		, fill = "right"
		, sep = "_"
	) |>
    mutate(
		pid = as.numeric(factor(pid)),
		correct_button = factor(correct_button,
								levels = c("incorrect", "correct", "no_response"))
	) |>
	rename(individual = pid)
)

#check separation of tasks/conditions
(
	dat |>
	group_by(task, condition, setsize) |>
	count()
)

# viz
(
	dat |>
	  unite(
		col = "task_condition_setsize", task, condition, setsize
	) |>
	ggplot() +
	  geom_point(aes(x = rw, y = individual), alpha = 0.5) +
	  facet_grid(task_condition_setsize ~ correct_button)
)

#treat misses as inaccurate & scale RW within each task
dat <-
(
	dat
	%>% group_by(task)
	%>% mutate(
			acc01 = case_when(
				correct_button == "correct" ~ 1,
				TRUE ~ 0
			),
			rw_scaled = scale(rw, center = FALSE)[, 1]
	)
)

# viz again
(
	dat  |>
	unite(
		col = "task_condition_setsize",
		task, condition, setsize
	)
	%>% ggplot() +
	  geom_point(aes(x = rw_scaled, y = individual), alpha = 0.5) +
	  facet_grid(task_condition_setsize ~ acc01)
)

# Prepare inputs to Stan ----
cmdstan_version()

# ungroup & sort by individual
dat <-
(
	dat
	# ungroup
	%>% ungroup()
	# ensure individual is a sequential numeric
	%>% mutate(
		individual = as.numeric(factor(individual))
	)
	# arrange rows by individual
	%>% arrange(individual)
	%>% mutate(dat_row = 1:n())
)

#compute group contrasts Xg from distinct combinations of groups
age_scaled <-
(
	dat |>
	#select down to any G vars
	  select(age) |>
	#collapse to distinct set of rows (combinations of grouping variables)
	  distinct() |>
	#arrange (not really necessary, but why not)
	  arrange(age) |>
	#first scale age
	  mutate(
		age_scaled = (age - median(age))/diff(range(age))
	)
)

(
age_scaled$contrasts <- get_contrast_matrix_rows_as_list(
							data = age_scaled,
							formula = ~ age_scaled,
							# half-sum contrasts are nice for 2-level variables bc they yield parameters whose value
							# is the difference between conditions
							contrast_kind = halfsum_contrasts)
)
Xg_with_vars <- age_scaled


#show contrasts
(
	Xg_with_vars
	%>% unnest(contrasts)
)

#join Xg with dat to label individuals with corresponding row from Xg
(
	Xg_with_vars
	# add row identifier
	%>% mutate(Xg_row=1:n())
	#join with dat, collapsed to 1-row per individual with their group info
	%>% right_join((#right-join to apply the row order from dat
		dat
		%>% select(individual, age)
		%>% distinct()
	))
	# grab the Xg row identifier (remember not to re-order dat from here on!)
	%>% pull(Xg_row)
) ->
	iXg

# compute Xc from distinct combinations of individuals & conditions
#   n.b. tweaked relative to parent hmg code to compute separate contrast matrices for each task (inc. more complexity in BOXED task)
(
	dat
	#select individual & any condition-defining columns
	%>% select(individual,task,condition,setsize)
	#collapse down to distinct rows (1 per individual/conditions combo)
	%>% distinct()
	%>% group_by(task)
	%>% group_split()
	# %>% pluck(1) -> x)
	%>% purrr::map_dfr(
		.f = function(x){
			task = x$task[1]
			if(task=="BOXED"){
				contrasts_formula = ~condition*setsize
				x = select(x,individual,condition,setsize)
			}else{
				contrasts_formula = ~condition
				x = select(x,individual,condition)
			}
			(
				x
				%>% as.list()
				%>% map(unique)
				%>% cross_df()
				# arrange (not really necessary, but why not)
				%>% arrange()
				# add the contrast matrix columns
				%>% mutate(
					contrasts = get_contrast_matrix(
						data = .
						, formula = contrasts_formula
						, contrast_kind = halfsum_contrasts
					)
					, task = task
				)
			) ->
				to_return
			names(to_return)[names(to_return)=="contrasts"] = paste0(task,"_contrasts")
			return(to_return)
		}
	)
	%>% {function(x){
		(
			x
			%>% select(contains("contrasts"))
			%>% as.matrix()
			%>% replace_na(0)
			%>% as_tibble()
			%>% bind_cols(
				select(x,!contains("contrasts"))
				,.
			)
		)
	}}()
) ->
	complete_Xc_with_vars

# show the unique contrasts
# (
# 	complete_Xc_with_vars
# 	%>% select(task,setsize,condition,contains("contrasts"))
# 	%>% distinct()
# 	%>% View()
# )

#subset down to just those individual-condition combos actually present in the data
#  it"s ok if there"s no missing data and nrow(complete_Xc_with_vars)==nrow(Xc_with_vars)
(
	complete_Xc_with_vars
	%>% semi_join(dat)
	%>% arrange()
	%>% arrange(individual)
) ->
	Xc_with_vars

#join Xc with dat to label observations with corresponding row from Xc
(
	Xc_with_vars
	# add row identifier
	%>% mutate(
		Xc_row=1:n()
	)
	# right-join with dat to preserve dat"s row order
	%>% left_join(
		(
			dat
			# %>% mutate(
			# 	setsize = replace_na(0)
			# )
		)
	)
	%>% arrange(dat_row)
	#pull the Xc row identifier
	%>% pull(Xc_row)
) ->
	yXc


# package for stan & sample ----

data_for_stan <- lst( #lst permits later entries to refer to earlier entries

	####
	# Entries we need to specify ourselves
	####

	# Xg: group-level predictor matrix
	Xg = (
		Xg_with_vars
		%>% select(contrasts)
		%>% unnest(contrasts)
		%>% as.matrix()
	)

	# iXg: which group each individual is associated with
	, iXg = iXg

	# Xc: condition-level predictor matrix
	, Xc = (
		Xc_with_vars
		%>% select(contains("contrasts"))
		%>% as.matrix()
		%>% replace_na(0)
	)


	# iXc: which individual is associated with each row in Xc
	, iXc = as.numeric(factor(Xc_with_vars$individual))

	# Y: observations
	, Y = dat$acc01

	# yXc: which row in Xc is associated with each observation in Y
	, yXc = yXc

	# intensity: intensity covariate for each observation
	, intensity = dat$rw_scaled

	# p_chance: probability of success @ chance for each observation
	, p_chance = rep(.5, times = length(Y))

	####
	# Entries computable from the above
	####

	# nXg: number of cols in the group-level predictor matrix
	, nXg = ncol(Xg)

	# rXg: number of rows in the group-level predictor matrix
	, rXg = nrow(Xg)

	# nI: number of individuals
	, nI = max(iXc)

	# nXc: number of cols in the condition-level predictor matrix
	, nXc = ncol(Xc)

	# rXc: number of rows in the condition-level predictor matrix
	, rXc = nrow(Xc)

	# nY: num entries in the observation vectors
	, nY = length(Y)

)

# double-check:
glimpse(data_for_stan)

#set the model path
mod_path <- "stan/hierarchical_multinormal_lapsing_psychometric.stan"

#set the model centered/non-centeredness
#  generally, if *either* nI_per_group *or* num_Y_per_q is small, non-centered will sample better than centered
data_for_stan$centered = FALSE
#conversion to 1/0 for stan
data_for_stan$centered = as.numeric(data_for_stan$centered)

mymodel_obj <- cmdstan_model(mod_path, cpp_options = list(stan_threads = TRUE))

# fill in NAs for Xc
data_for_stan$Xc[is.na(data_for_stan$Xc)] <- 0

#not used by cmdstan set_num_threads(8)
#---------------------------------------------
# Now fore sampling/optimizing
# ---   The command stan call to solve with MCMC:
post <- mymodel_obj$sample(data_for_stan, chains = 4, refresh = 1,
						   threads_per_chain = 8,
						   output_dir = ".", validate_csv = FALSE,
						   adapt_delta = 0.8)

stanfit <- rstan::read_stan_csv(post$output_files())

stanfit <- post$summary()
# View(summary)
# View(stanfit)
#set the posterior path (automated but you could do your own if you had multiple models)
# (
# 	mod_path
# 	%>% fs::path_file()
# 	%>% fs::path_ext_remove()
# 	%>% paste0(
# 		ifelse(data_for_stan$centered,"_c","_nc")
# 	)
# 	%>% fs::path(
# 		"posteriors"
# 		, .
# 		, ext = "netcdf4"
# 	)
# ) -> post_path

# ensure model is compiled
# aria:::check_and_compile(mod_path, block = TRUE)

# compose
# aria::compose(
# 	data = data_for_stan,
# 	code_path = mod_path,
# 	out_path = post_path,
# 	overwrite = T,
# 	block = T
# )

# check posterior diagnostics ----
# post = aria::coda(post_path)

# Check treedepth, divergences, & rebfmi
# (
# 	post$draws(group="sample_stats")
# 	%>% posterior::as_draws_df()
# 	%>% group_by(.chain)
# 	%>% summarise(
# 		max_treedepth = max(treedepth)
# 		, num_divergent = sum(divergent)
# 		, rebfmi = var(energy)/(sum(diff(energy)^2)/n()) #n.b. reciprocal of typical EBFMI, so bigger=bad, like rhat
# 	)
# )

# View(post$summary()) #view all of summary

# View(post$summary()%>%
# 	 	posterior::as_draws_df())
#
# post$summary()


# post$summary()%>%
# 	posterior::as_draws_df()%>%
# 	group_by(.chain)%>%
# 	summarise(
# 	max_treedepth = max(treedepth),
# 	num_divergent = sum(divergent),
# 	rebfmi = var(energy)/(sum(diff(energy)^2)/n()) #n.b. reciprocal of typical EBFMI, so bigger=bad, like rhat
# )

par_summary <- post$summary()

library(posterior)
draws_array <- post$draws()
str(draws_array)
draws_df <- as_draws_df(draws_array) # as_draws_matrix() for matrix
head(draws_df)



str(post$sampler_diagnostics())
diagnostics_df <- as_draws_df(post$sampler_diagnostics())

diagnostics_df |>
	summarise(
			max_treedepth = max(treedepth__),
			num_divergent = sum(divergent__),
			rebfmi = var(energy__)/(sum(diff(energy__)^2)/n()) #n.b. reciprocal of typical EBFMI, so bigger=bad, like rhat
		)





# # gather summary for core parameters (inc. r̂ & ess)
# (
# 	post$draws(group="parameters")
# 	%>% posterior::summarise_draws(.cores=parallel::detectCores())
# ) ->
# 	par_summary

# show the ranges of r̂/ess"s
(
	par_summary %>%
	select(rhat,contains("ess"))
	%>% summary()
)

#View those with suspect r̂
(
	par_summary
	%>% filter(rhat>1.01)
	%>% (function(suspects){
		if(nrow(suspects)>=1){
			View(suspects)
		}
		return(paste("# suspect parameters:",nrow(suspects)))
	})()
)

# Viz recovery of (some) non-correlation parameters ----


# str(post$draws())

post$draws(variables=c("Z","iZc_sd") )%>%
	posterior::as_draws_df() %>%
	select(-.draw)%>%
	pivot_longer(
	cols = -c(.chain,.iteration),
	names_to = "variable")%>%
	group_by(variable)%>%
	arrange(variable,.chain,.iteration)%>%
	summarise(
		rhat = 1.01<posterior::rhat(matrix(value,ncol=length(unique(.chain)))),
		ess_bulk = 100>posterior::ess_bulk(matrix(value,ncol=length(unique(.chain)))),
		ess_tail = 100>posterior::ess_tail(matrix(value,ncol=length(unique(.chain)))),
		as_tibble(t(posterior::quantile2(value,c(.05,.25,.5,.75,.95)))))

#should this be iZq or iZc
View(str(post$draws()))
(
	post$draws(variables=c("Z","iZc_sd"))
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	%>% pivot_longer(
		cols = -c(.chain,.iteration)
		, names_to = "variable"
	)
	%>% group_by(variable)
	%>% arrange(variable,.chain,.iteration)
	%>% summarise(
		rhat = 1.01<posterior::rhat(matrix(value,ncol=length(unique(.chain))))
		, ess_bulk = 100>posterior::ess_bulk(matrix(value,ncol=length(unique(.chain))))
		, ess_tail = 100>posterior::ess_tail(matrix(value,ncol=length(unique(.chain))))
		, as_tibble(t(posterior::quantile2(value,c(.05,.25,.5,.75,.95))))
	)
	# %>% mutate(variable = factor_1d(variable))
	%>% ggplot()
	+ geom_hline(yintercept = 0)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q5
			, ymax = q95
			, colour = ess_tail
		)
		, alpha = .5
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q25
			, ymax = q75
			, colour = ess_bulk
		)
		, size = 3
		, alpha = .5
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = q50
			, fill = rhat
		)
		, shape = 21
		, size = 3
	)
	+ coord_flip()
	+ scale_color_manual(
		values = lst(`TRUE`="red",`FALSE`="black")
		, labels = lst(`TRUE`="<100",`FALSE`=">=100")
	)
	+ scale_fill_manual(
		values = lst(`TRUE`="red",`FALSE`="white")
		, labels = lst(`TRUE`=">1.01",`FALSE`="<=1.01")
	)
	+ labs(
		y = "Posterior Value"
		, x = "Variable"
		, colour = "ESS"
		, fill = "Rhat"
	)
)


# Viz recovery of correlations ----
(
	post$draws(variables="iZq_r_vec")
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	%>% pivot_longer(
		cols = -c(.chain,.iteration)
		, names_to = "variable"
	)
	%>% group_by(variable)
	%>% arrange(variable,.chain,.iteration)
	%>% summarise(
		rhat = 1.01<posterior::rhat(matrix(value,ncol=length(unique(.chain))))
		, ess_bulk = 100>posterior::ess_bulk(matrix(value,ncol=length(unique(.chain))))
		, ess_tail = 100>posterior::ess_tail(matrix(value,ncol=length(unique(.chain))))
		, as_tibble(t(posterior::quantile2(value,c(.05,.25,.5,.75,.95))))
	)
	%>% mutate(variable = factor_1d(variable))
	%>% ggplot()
	+ geom_hline(yintercept = 0)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q5
			, ymax = q95
			, colour = ess_tail
		)
		, alpha = .5
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q25
			, ymax = q75
			, colour = ess_bulk
		)
		, size = 3
		, alpha = .5
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = q50
			, fill = rhat
		)
		, shape = 21
		, size = 3
	)
	+ coord_flip()
	+ scale_color_manual(
		values = lst(`TRUE`="red",`FALSE`="black")
		, labels = lst(`TRUE`="<100",`FALSE`=">=100")
	)
	+ scale_fill_manual(
		values = lst(`TRUE`="red",`FALSE`="white")
		, labels = lst(`TRUE`=">1.01",`FALSE`="<=1.01")
	)
	+ labs(
		y = "True & Posterior Value"
		, x = "Variable"
		, colour = "ESS"
		, fill = "Rhat"
	)
)

# Viz recovery of (some) non-correlation parameters ----
(
	post$draws(variables="iZq_")
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	%>% pivot_longer(
		cols = -c(.chain,.iteration)
		, names_to = "variable"
	)
	%>% group_by(variable)
	%>% arrange(variable,.chain,.iteration)
	%>% summarise(
		rhat = 1.01<posterior::rhat(matrix(value,ncol=length(unique(.chain))))
		, ess_bulk = 100>posterior::ess_bulk(matrix(value,ncol=length(unique(.chain))))
		, ess_tail = 100>posterior::ess_tail(matrix(value,ncol=length(unique(.chain))))
		, as_tibble(t(posterior::quantile2(value,c(.05,.25,.5,.75,.95))))
	)
	# %>% mutate(variable = factor_1d(variable))
	%>% ggplot()
	+ geom_hline(yintercept = 0)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q5
			, ymax = q95
			, colour = ess_tail
		)
		, alpha = .5
	)
	+ geom_linerange(
		mapping = aes(
			x = variable
			, ymin = q25
			, ymax = q75
			, colour = ess_bulk
		)
		, size = 3
		, alpha = .5
	)
	+ geom_point(
		mapping = aes(
			x = variable
			, y = q50
			, fill = rhat
		)
		, shape = 21
		, size = 3
	)
	+ coord_flip()
	+ scale_color_manual(
		values = lst(`TRUE`="red",`FALSE`="black")
		, labels = lst(`TRUE`="<100",`FALSE`=">=100")
	)
	+ scale_fill_manual(
		values = lst(`TRUE`="red",`FALSE`="white")
		, labels = lst(`TRUE`=">1.01",`FALSE`="<=1.01")
	)
	+ labs(
		y = "Posterior Value"
		, x = "Variable"
		, colour = "ESS"
		, fill = "Rhat"
	)
)


#getting "threshold_for_subj_by_cond"

# Viz individuals" functions ----
(
	post$draws("threshold_for_subj_by_cond")
	%>% posterior::as_draws_df()
	%>% select(-.draw)
	%>% pivot_longer(
		cols = c(-.chain,-.iteration)
		, names_prefix = fixed("threshold_for_subj_by_cond")
		, values_to = "threshold"
	)
	%>% left_join(
		(
			Xc_with_vars
			%>% select(-contrasts)
			%>% mutate(
				name = paste0("[",1:n(),"]")
			)
		)
		, by = "name"
	)
	%>% select(-name)
)


str(post$draws())
View(stanfit)
stanfit@model_pars

View(post$draws()%>%
	 	posterior::as_draws_df()%>%
	 	select(-.draw)%>%
	 	pivot_longer(
	 	cols = c(-.chain,-.iteration), names_prefix = fixed("threshold_for_subj_by_cond"),
	 	values_to = "threshold")%>%
	 	left_join(
	 	(
	 		Xc_with_vars%>% select(-contrasts)%>%
	 			mutate(
	 			name = paste0("[",1:n(),"]")
	 		)
	 	)
	 	, by = "name"
	 )%>%
	 	select(-name)
)
#so here I do get the code.


post$summary()

View(Xc_with_vars)
select(-contrasts)%>%
	mutate(
		name = paste0("[",1:n(),"]"))


View(post$draws()%>%
	 	posterior::as_draws_df()%>%
	 	select(-.draw)%>%
	 	pivot_longer(
	 		cols = c(-.chain,-.iteration), names_prefix = fixed("threshold_for_subj_by_cond"),
	 		values_to = "threshold"))
