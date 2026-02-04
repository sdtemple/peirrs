# Measles in Hagelloch analysis
# Seth Temple, sethtem@umich.edu

input_file = '../data/measles_hagelloch.csv'
p = 1
q = 1

# Load in packages --------------------------------------------------------

{
  
  library(peirrs)
  
}

# SEIR without classes ----------------------------------------------------

measles <- read.csv(input_file)
population_size <- 185

{
  # final form for package functions
  lag <- 10
  infections <- measles$infectious_time
  removals <- measles$removal_time
  removal_classes <- as.numeric(measles$class)
  infection_classes <- as.numeric(measles$class)
  infection_class_sizes <- c(90, 30, 65)
  epidemic_size <- length(removals)
}

# frequentist estimates
result_freq <- peirr_tau(removals,
                    infections,
                    population_size,
                    lag=lag)

result_freq_multitype <- peirr_tau_multitype(removals,
                                                infections,
                                                rep(1, population_size),
                                                infection_classes,
                                                infection_class_sizes,
                                                lag=lag
)

# bayesian estimates
result_bayes_multitype <- peirr_bayes_multitype(removals,
                                                infections,
                                                rep(1, population_size),
                                                infection_classes,
                                                infection_class_sizes,
                                                lag=lag,
                                                num_iter=10000,
                                                num_print=1e6,
                                                gamma_init=0.1,
                                                gamma_shape = 1,
                                                beta_init = rep(1,
                                                                length(unique(infection_classes))),
                                                beta_shape = rep(1,
                                                                 length(unique(infection_classes)))
)

result_bayes <- peirr_bayes(
  removals,
  infections,
  population_size,
  num_iter=10000,
  lag=lag,
  num_print=1e6,
  update_gamma=TRUE
)


# Credible intervals ------------------------------------------------------

# credible intervals for singular infection rate
quantile(result_bayes$infection_rate,c(0.025,0.975))

# credible intervals for multitype infection rates
# none of these overlap
apply(result_bayes_multitype$infection_rate, 1, quantile, c(0.025,0.975))

# Point estimates ---------------------------------------------------------

# results below are fundamentally all the same

# singular infection rates
result_freq$infection_rate
result_freq$removal_rate
result_freq$effective_number

mean(result_bayes$infection_rate)
mean(result_bayes$removal_rate)
mean(result_bayes$infection_rate/result_bayes$removal_rate)

# multitype infection rates
result_freq_multitype$infection_rate
apply(result_bayes_multitype$infection_rate,1,mean)
