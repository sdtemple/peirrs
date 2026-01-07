# Analyzing the Abakaliki dataset -----------------------------------------

library(pblas)
library(peirrs)

# Setup -------------------------------------------------------------------

{
  num_iter <- 500
  num_boot <- 50
  population_size <- 15
  probs_complete <- c(0.4, 0.6, 0.8, 1.)
  probs_infection_missing <- c(0, 0.2, 0.4, 0.6, 0.8, 1) 
  num_datasets <- 100
  
  dataset_column <- c()
  ps <- c()
  qs <- c()
  bs <- c()
  bs2 <- c()
  gs <- c()
  r0s <- c()
  r0s2 <- c()
  
  bayes_bs <- c()
  bayes_bs2 <- c()
  bayes_gs <- c()
  bayes_gs2 <- c()
  bayes_r0s <- c()
  bayes_r0s2 <- c()
}


# Computation loop --------------------------------------------------------

for (d in 1:num_datasets) {
  
  print(paste0("Analyzing dataset ", d))
  
  filename <- paste0(
    "abakaliki-modified-",
    d,
    ".csv"
  )
  dataset <- read.csv(filename,
                      sep=',',
  )
  lag <- mean(dataset$FeverDay - dataset$ExposedDay)
  
  removals <- dataset$RecoveryDay
  infections <- dataset$FeverDay
  mn <- min(infections)
  removals <- removals - mn
  infections <- infections - mn
  matrix_time <- cbind(infections, removals)
  
  for (p in probs_complete) {
    for (q in probs_infection_missing) {
      # print(c(p,q))
      
      tryCatch ({
        matrix_decompleted <- decomplete_sem(matrix_time, 
                                             p,
                                             q
        )
        removals_decompleted <- matrix_decompleted[,2]
        infections_decompleted <- matrix_decompleted[,1]
        
        output <- peirr_tau(removals_decompleted,
                            infections_decompleted,
                            population_size,
                            lag
        )
        beta <- output$infection_rate
        gamma <- output$removal_rate
        r0 <- beta/gamma
        boot <- peirr_bootstrap(num_boot, 
                                beta, 
                                gamma, 
                                population_size, 
                                length(removals), 
                                p, 
                                q,
                                within=0.01
        )
        
        boot_beta <- boot$infection_rate[2:num_boot+1]
        boot_gamma <- boot$removal_rate[2:num_boot+1]
        boot_r0 <- boot_beta / boot_gamma
        g2 <- gamma + mean(gamma - boot_gamma, na.rm=TRUE)
        beta2 <- beta + mean(beta - boot_beta, na.rm=TRUE)
        r02 <- r0 + mean(r0 - boot_r0, na.rm=TRUE)
        
        bs <- c(bs, beta)
        bs2 <- c(bs2, beta2)
        r0s <- c(r0s, r0)
        r0s2 <- c(r0s2, r02)
        gs <- c(gs, gamma)
        ps <- c(ps, p)
        qs <- c(qs, q)
        dataset_column <- c(dataset_column, d)
        
        output_bayes <- peirr_bayes(removals_decompleted,
                                    infections_decompleted,
                                    population_size,
                                    beta_init=beta,
                                    gamma_init=gamma,
                                    num_iter=num_iter,
                                    num_print=1000,
                                    lag=lag
        )
        bayes_beta <- output_bayes$infection_rate[101:num_iter]
        bayes_gamma <- output_bayes$removal_rate[101:num_iter]
        bayes_r0 <- bayes_beta / bayes_gamma
        
        bayes_bs <- c(bayes_bs, mean(bayes_beta))
        bayes_bs2 <- c(bayes_bs2, quantile(bayes_beta, 0.5))
        
        bayes_gs <- c(bayes_gs, mean(bayes_gamma))
        bayes_gs2 <- c(bayes_gs2, quantile(bayes_gamma, 0.5))
        
        bayes_r0s <- c(bayes_r0s, mean(bayes_r0))
        bayes_r0s2 <- c(bayes_r0s2, quantile(bayes_r0, 0.5))
        
        
      },
      error = function(e) {
        message("An error occurred:", e$message)
        bs <- c(bs, NA)
        bs2 <- c(bs2, NA)
        r0s <- c(r0s, NA)
        r0s2 <- c(r0s2, NA)
        gs <- c(gs, NA)
        ps <- c(ps, p)
        qs <- c(qs, q)
        dataset_column <- c(dataset_column, d)
        bayes_bs <- c(bayes_bs, NA)
        bayes_bs2 <- c(bayes_bs2, NA)
        bayes_gs <- c(bayes_gs, NA)
        bayes_gs2 <- c(bayes_gs2, NA)
        bayes_r0s <- c(bayes_r0s, NA)
        bayes_r0s2 <- c(bayes_r0s2, NA)
        
      }
      )
    }
  }
}


# Saving ------------------------------------------------------------------

{
results <- cbind(dataset_column, 
                 ps, 
                 qs,
                 gs, 
                 bs, 
                 bs2, 
                 r0s, 
                 r0s2,
                 bayes_gs, 
                 bayes_gs2,
                 bayes_bs, 
                 bayes_bs2,
                 bayes_r0s, 
                 bayes_r0s2
                 )
results <- data.frame(results)
colnames(results) <- c("dataset",
                       "p",
                       "q",
                       "gamma",
                       "beta",
                       "beta_corrected",
                       "R0",
                       "R0_corrected",
                       "gamma_bayes_mean",
                       "gamma_bayes_median",
                       "beta_bayes_mean",
                       "beta_bayes_median",
                       "R0_bayes_mean",
                       "R0_bayes_median"
                       )

write.table(results,
            "abakaliki-results.tsv",
            sep='\t',
            row.names=FALSE
            )
}
