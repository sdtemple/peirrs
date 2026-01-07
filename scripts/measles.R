# Measles in Hagelloch analysis
# Seth Temple, sethtem@umich.edu


# Load in packages --------------------------------------------------------

{

library(pblas)
library(peirrs)

library(outbreaks)
library(lubridate)

}

# Load and format the dataset ---------------------------------------------

{
  measles <- measles_hagelloch_1861
  measles <- measles[measles$age < 14,]
  population_size <- dim(measles)[1]
}
# do not analyze the 3 children 14 and 15
# not infected in the 1847 epidemic
# analysis likely robust to this exclusion

# convert the first symptoms date to numeric
{
  x <- as.numeric(ymd(measles$date_of_prodrome))
  x <- x + rnorm(length(x),sd=0.1)
  minx <- min(x)
  x <- x - minx
  measles$numeric_prodrome <- x
}

# convert the rash date to numeric
{
  x <- as.numeric(ymd(measles$date_of_rash))
  x <- x + rnorm(length(x),sd=0.1)
  x <- x - minx
  measles$numeric_rash <- x
}

# convert the death date to numeric
{
  x <- as.numeric(ymd(measles$date_of_death))
  x <- x + rnorm(length(x),sd=0.1)
  x <- x - minx
  measles$numeric_death <- x
}

{
  # transform to infectious and removal times
  measles$infectious_time <- measles$numeric_prodrome - 1
  measles$removal_time <- measles$numeric_rash + 3
  measles$removal_time <- pmin(measles$removal_time, measles$numeric_death, na.rm=TRUE)
  # removal time must end if there is death
}

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

# SEIR without classes ----------------------------------------------------


# frequentist

{
  num_boot <- 200
  num_runs <- 5
  ps <- c()
  qs <- c()
  betas <- c()
  gammas <- c()
  r0s <- c()
  betas_center <- c()
  gammas_center <- c()
  r0s_center <- c()
  betas_lower <- c()
  betas_upper <- c()
  gammas_lower <- c()
  gammas_upper <- c()
  r0s_lower <- c()
  r0s_upper <- c()
  betas_lower_center <- c()
  betas_upper_center <- c()
  gammas_lower_center <- c()
  gammas_upper_center <- c()
  r0s_lower_center <- c()
  r0s_upper_center <- c()
}

{

  start <- Sys.time()

  # complete data
  result <- peirr_tau(removals,
                      infections,
                      population_size,
                      lag=lag)
  beta <- result$infection_rate
  gamma <- result$removal_rate
  betas <- c(betas, beta)
  gammas <- c(gammas, gamma)
  r0s <- c(r0s, beta/gamma)
  ps <- c(ps, 1)
  qs <- c(qs, NA)

  boot_result <- peirr_bootstrap(num_boot,
                                 beta,
                                 gamma,
                                 population_size,
                                 epidemic_size,
                                 prop_complete = 1,
                                 prop_infection_missing = 1,
                                 num_renewals = 1,
                                 lag = lag,
                                 within= 0.05,
                                 peirr = peirr_tau
  )

  # extract bootstraps
  boot_beta <- boot_result$infection_rate[2:num_boot]
  boot_gamma <- boot_result$removal_rate[2:num_boot]
  boot_r0 <- boot_beta / boot_gamma

  # bias
  boot_bias_r0 <- mean(boot_r0 - beta / gamma, na.rm=TRUE)
  boot_bias_beta <- mean(boot_beta - beta, na.rm=TRUE)
  boot_bias_gamma <- mean(boot_gamma - gamma, na.rm=TRUE)

  # confidence bounds
  lower_cb_r0 <- quantile(boot_r0, 0.025, na.rm=TRUE)
  upper_cb_r0 <- quantile(boot_r0, 0.975, na.rm=TRUE)
  lower_cb_beta <- quantile(boot_beta, 0.025, na.rm=TRUE)
  upper_cb_beta <- quantile(boot_beta, 0.975, na.rm=TRUE)
  lower_cb_gamma <- quantile(boot_gamma, 0.025, na.rm=TRUE)
  upper_cb_gamma <- quantile(boot_gamma, 0.975, na.rm=TRUE)

  # r0 confidence intervals
  lower_reverse_r0 <- 2 * beta/gamma - upper_cb_r0
  lower_reverse_r0_center <- lower_reverse_r0 - boot_bias_r0
  upper_reverse_r0 <- 2 * beta/gamma - lower_cb_r0
  upper_reverse_r0_center <- upper_reverse_r0 - boot_bias_r0
  r0_center <- beta/gamma - boot_bias_r0

  # beta confidence intervals
  lower_reverse_beta <- 2 * beta - upper_cb_beta
  lower_reverse_beta_center <- lower_reverse_beta - boot_bias_beta
  upper_reverse_beta <- 2 * beta- lower_cb_beta
  upper_reverse_beta_center <- upper_reverse_beta - boot_bias_beta
  beta_center <- beta - boot_bias_beta

  # gamma confidence intervals
  lower_reverse_gamma <- 2 * gamma - upper_cb_gamma
  lower_reverse_gamma_center <- lower_reverse_gamma - boot_bias_gamma
  upper_reverse_gamma <- 2 * gamma- lower_cb_gamma
  upper_reverse_gamma_center <- upper_reverse_gamma - boot_bias_gamma
  gamma_center <- gamma - boot_bias_gamma

  # not centered recording
  betas_lower <- c(betas_lower, lower_reverse_beta)
  betas_upper <- c(betas_upper, upper_reverse_beta)
  gammas_lower <- c(gammas_lower, lower_reverse_gamma)
  gammas_upper <- c(gammas_upper, upper_reverse_gamma)
  r0s_lower <- c(r0s_lower, lower_reverse_r0)
  r0s_upper <- c(r0s_upper, upper_reverse_r0)

  # centered recording
  betas_center <- c(betas_center, beta_center)
  gammas_center <- c(gammas_center, gamma_center)
  r0s_center <- c(r0s_center, r0_center)
  betas_lower_center <- c(betas_lower_center, lower_reverse_beta_center)
  betas_upper_center <- c(betas_upper_center, upper_reverse_beta_center)
  gammas_lower_center <- c(gammas_lower_center, lower_reverse_gamma_center)
  gammas_upper_center <- c(gammas_upper_center, upper_reverse_gamma_center)
  r0s_lower_center <- c(r0s_lower_center, lower_reverse_r0_center)
  r0s_upper_center <- c(r0s_upper_center, upper_reverse_r0_center)

  end <- Sys.time()
  print(end - start)

}


matrix_time <- cbind(infections, removals)
for (run in 1:num_runs) {
  print(paste0("Run ", run))

  for (p in seq(0.4, 0.8, 0.2)) {

    for (q in seq(0.4, 0.8, 0.2)) {

      print(c(p,q))

      start <- Sys.time()

      # decomplete
      epidemic <- filter_sem(matrix_time)
      epidemic <- decomplete_sem(epidemic, p, q)
      epidemic <- sort_sem(epidemic)
      removals_d <- epidemic[,2]
      infections_d <- epidemic[,1]

      # estimate
      result_d <- peirr_tau(removals_d,
                            infections_d,
                            population_size,
                            lag=lag)
      beta_d <- result_d$infection_rate
      gamma_d <- result_d$removal_rate

      betas <- c(betas, beta_d)
      gammas <- c(gammas, gamma_d)
      r0s <- c(r0s, beta_d/gamma_d)
      ps <- c(ps, p)
      qs <- c(qs, q)

      # bootstrap routine
      boot_result_d <- peirr_bootstrap(num_boot,
                                       beta_d,
                                       gamma_d,
                                       population_size,
                                       epidemic_size,
                                       prop_complete = p,
                                       prop_infection_missing = q,
                                       num_renewals = 1,
                                       lag = lag,
                                       within= 0.05,
                                       peirr = peirr_tau
      )

      # extract bootstraps
      boot_beta <- boot_result_d$infection_rate[2:num_boot]
      boot_gamma <- boot_result_d$removal_rate[2:num_boot]
      boot_r0 <- boot_beta / boot_gamma

      # bias
      boot_bias_r0 <- mean(boot_r0 - beta_d / gamma_d, na.rm=TRUE)
      boot_bias_beta <- mean(boot_beta - beta_d, na.rm=TRUE)
      boot_bias_gamma <- mean(boot_gamma - gamma_d, na.rm=TRUE)

      # confidence bounds
      lower_cb_r0 <- quantile(boot_r0, 0.025, na.rm=TRUE)
      upper_cb_r0 <- quantile(boot_r0, 0.975, na.rm=TRUE)
      lower_cb_beta <- quantile(boot_beta, 0.025, na.rm=TRUE)
      upper_cb_beta <- quantile(boot_beta, 0.975, na.rm=TRUE)
      lower_cb_gamma <- quantile(boot_gamma, 0.025, na.rm=TRUE)
      upper_cb_gamma <- quantile(boot_gamma, 0.975, na.rm=TRUE)

      # r0 confidence intervals
      lower_reverse_r0 <- 2 * beta_d/gamma_d - upper_cb_r0
      lower_reverse_r0_center <- lower_reverse_r0 - boot_bias_r0
      upper_reverse_r0 <- 2 * beta_d/gamma_d - lower_cb_r0
      upper_reverse_r0_center <- upper_reverse_r0 - boot_bias_r0
      r0_center <- beta_d/gamma_d - boot_bias_r0

      # beta confidence intervals
      lower_reverse_beta <- 2 * beta_d - upper_cb_beta
      lower_reverse_beta_center <- lower_reverse_beta - boot_bias_beta
      upper_reverse_beta <- 2 * beta_d- lower_cb_beta
      upper_reverse_beta_center <- upper_reverse_beta - boot_bias_beta
      beta_center <- beta_d - boot_bias_beta

      # gamma confidence intervals
      lower_reverse_gamma <- 2 * gamma_d - upper_cb_gamma
      lower_reverse_gamma_center <- lower_reverse_gamma - boot_bias_gamma
      upper_reverse_gamma <- 2 * gamma_d- lower_cb_gamma
      upper_reverse_gamma_center <- upper_reverse_gamma - boot_bias_gamma
      gamma_center <- gamma_d - boot_bias_gamma

      # not centered recording
      betas_lower <- c(betas_lower, lower_reverse_beta)
      betas_upper <- c(betas_upper, upper_reverse_beta)
      gammas_lower <- c(gammas_lower, lower_reverse_gamma)
      gammas_upper <- c(gammas_upper, upper_reverse_gamma)
      r0s_lower <- c(r0s_lower, lower_reverse_r0)
      r0s_upper <- c(r0s_upper, upper_reverse_r0)

      # centered recording
      betas_center <- c(betas_center, beta_center)
      gammas_center <- c(gammas_center, gamma_center)
      r0s_center <- c(r0s_center, r0_center)
      betas_lower_center <- c(betas_lower_center, lower_reverse_beta_center)
      betas_upper_center <- c(betas_upper_center, upper_reverse_beta_center)
      gammas_lower_center <- c(gammas_lower_center, lower_reverse_gamma_center)
      gammas_upper_center <- c(gammas_upper_center, upper_reverse_gamma_center)
      r0s_lower_center <- c(r0s_lower_center, lower_reverse_r0_center)
      r0s_upper_center <- c(r0s_upper_center, upper_reverse_r0_center)

      end <- Sys.time()
      print(paste0(end - start))

    }
  }
}

{
  # saving the data
  table <- data.frame(
    cbind(ps,
          qs,
          betas,
          unname(betas_lower),
          unname(betas_upper),
          betas_center,
          unname(betas_lower_center),
          unname(betas_upper_center),
          gammas,
          unname(gammas_lower),
          unname(gammas_upper),
          gammas_center,
          unname(gammas_lower_center),
          unname(gammas_upper_center),
          r0s,
          unname(r0s_lower),
          unname(r0s_upper),
          r0s_center,
          unname(r0s_lower_center),
          unname(r0s_upper_center)
          )
  )

  colnames(table) <- c(
    "p",
    "q",
    "beta",
    "beta_lower",
    "beta_upper",
    "beta_center",
    "beta_lower_center",
    "beta_upper_center",
    "gamma",
    "gamma_lower",
    "gamma_upper",
    "gamma_center",
    "gamma_lower_center",
    "gamma_upper_center",
    "r0",
    "r0_lower",
    "r0_upper",
    "r0_center",
    "r0_lower_center",
    "r0_upper_center"
  )

  View(table)

  write.table(
    table,
    "measles_simple.tsv",
    sep="\t",
    row.names=FALSE
  )

}

# bayesian

result_bayes <- peirr_bayes(
  removals,
  infections,
  population_size,
  num_iter=10000,
  lag=lag,
  num_print=1e6
)

bayes_beta <- result_bayes$infection_rate
bayes_gamma <- result_bayes$removal_rate
bayes_r0 <- bayes_beta/bayes_gamma
bayes_beta_mean <- mean(bayes_beta)
bayes_gamma_mean <- mean(bayes_gamma)
bayes_r0_mean <- mean(bayes_r0)
quantile(bayes_beta, c(0.025,0.975))
quantile(bayes_gamma, c(0.025,0.975))
quantile(bayes_r0, c(0.025,0.975))

{
  start <- Sys.time()
  result_bayes <- peirr_bayes(
    removals_d,
    infections_d,
    population_size,
    num_iter=1000,
    lag=lag,
    num_print=1,
    gamma_init=0.1,
    update_gamma = TRUE
  )
  end <- Sys.time()
  end - start
}


# Pair-based likelihood approximation

peirr_pbla_infection_rate(removals,
                          infections,
                          population_size=population_size,
                          lag=lag)

peirr_pbla_both_rates(removals,
                      population_size = population_size,
                      lag = lag
                      )


# SEIR with classes -------------------------------------------------------

# not statistically significant to have different removal classes
peirr_tau_multitype(removals,
                    infections,
                    removal_classes,
                    infection_classes,
                    infection_class_sizes,
                    lag=lag
                    )

# all have the same removal class
peirr_tau_multitype(removals,
                    infections,
                    rep(1, epidemic_size),
                    infection_classes,
                    infection_class_sizes,
                    lag=lag
)

# effectively the same result
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

# may want to analyze the multitype with some partial data


