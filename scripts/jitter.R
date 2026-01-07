# How does rounding and jittering affect results --------------------------

library(peirrs)
library(pblas)

# Simulation setup --------------------------------------------------------

{
num_runs <- 30
not_jittered_betas <- c()
not_jittered_gammas <- c()
not_jittered_r0s <- c()
jittered_betas <- c()
jittered_gammas <- c()
jittered_r0s <- c()
sizes <- c()
betas <- c()
gammas <- c()
sigma <- 0.33
N <- 200
lag <- 10
}

# Simulations -------------------------------------------------------------

for (run in 1:num_runs) {
  print(run)
  for (beta in seq(0.2, 0.5, 0.1)) {
    for (gamma in c(0.05, 0.1, 0.2)) {
      simulated <- simulator(beta, 
                             gamma, 
                             N, 
                             lag=lag,
                             prop_complete=1,
                             prop_infection_missing=1
                             )
      r <- simulated$matrix_time[,2]
      i <- simulated$matrix_time[,1]
      i_rounded <- round(i)
      r_rounded <- round(r)
      i_jittered <- i_rounded + rnorm(length(i), 0, sigma)
      r_jittered <- r_rounded + rnorm(length(r), 0, sigma)
      result <- peirr_tau(r,
                          i,
                          N,
                          lag=lag
                          )
      result_jittered <- peirr_tau(
        r_jittered,
        i_jittered,
        N,
        lag=lag
      )
      jittered_beta <- result_jittered$infection_rate
      jittered_gamma <- result_jittered$removal_rate
      jittered_r0 <- jittered_beta / jittered_gamma
      not_jittered_beta <- result$infection_rate
      not_jittered_gamma <- result$removal_rate
      not_jittered_r0 <- not_jittered_beta / not_jittered_gamma
      betas <- c(betas, beta)
      gammas <- c(gammas, gamma)
      sizes <- c(sizes, length(i))
      not_jittered_betas <- c(not_jittered_betas, 
                              not_jittered_beta)
      not_jittered_gammas <- c(not_jittered_gammas, 
                              not_jittered_gamma)
      not_jittered_r0s <- c(not_jittered_r0s, 
                              not_jittered_r0)
      jittered_betas <- c(jittered_betas, 
                              jittered_beta)
      jittered_gammas <- c(jittered_gammas, 
                          jittered_gamma)
      jittered_r0s <- c(jittered_r0s, 
                          jittered_r0)
    }
  }
}


# Saving ------------------------------------------------------------------

table <- data.frame(cbind(betas, 
                       gammas,
                       sizes,
                       not_jittered_betas,
                       not_jittered_gammas,
                       not_jittered_r0s,
                       jittered_betas,
                       jittered_gammas,
                       jittered_r0s))
colnames(table) <- c("beta",
                  "gamma",
                  "sample_size",
                  "not_beta",
                  "not_gamma",
                  "not_r0",
                  "jitter_beta",
                  "jitter_gamma",
                  "jitter_r0"
                  )
View(table)

write.table(table,
            "jittered_results.tsv",
            sep='\t',
            row.names=FALSE
            )
