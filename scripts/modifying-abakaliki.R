# Modifying the Abakaliki dataset -----------------------------------------

dataset <- read.csv("abakaliki.csv")
dataset$JitteredDay <- dataset$Day + rnorm(13, 0, 0.5)

# Compute shape and rate of Gamma periods ---------------------------------

# mu and sigma of gamma periods
# from stockdale thesis 2019
{
  mu_fever_to_rash <- 2.49
  sigma_fever_to_rash <- 0.88
  mu_rash_to_recovery <- 16.0
  sigma_rash_to_recovery <- 2.83
}

{
  shape_func <- function(mu, sigma) {
    (mu / sigma) ** 2
  }
  rate_func <- function(mu, sigma) {
    mu / (sigma ** 2)
  }

  shape_fever_to_rash <- shape_func(
    mu_fever_to_rash,
    sigma_fever_to_rash
  )
  rate_fever_to_rash <- rate_func(
    mu_fever_to_rash,
    sigma_fever_to_rash
  )
  shape_rash_to_recovery <- shape_func(
    mu_rash_to_recovery,
    sigma_rash_to_recovery
  )
  rate_rash_to_recovery <- rate_func(
    mu_rash_to_recovery,
    sigma_rash_to_recovery
  )
}


# Augmenting fever and recovery dates -------------------------------------

# indicates data consistent with epidemic
check_if_epidemic = function(removals, infections, lag) {
  epidemic_size = length(removals)
  ind_matrix = matrix(0, nrow = epidemic_size, ncol = epidemic_size)
  for (j in 1:epidemic_size) {
    ind_matrix[j, ] = (infections[1:epidemic_size] < (infections[j] - lag)) * (removals > (infections[j] - lag))
  }
  chi_matrix = apply(ind_matrix, 1, sum)
  chi_matrix = chi_matrix[chi_matrix == 0]
  if (length(chi_matrix) > 1) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

num_datasets <- 100
for (d in 1:num_datasets) {
  rash_time <- dataset$JitteredDay
  removals <- rash_time + rgamma(13, 
                                 shape=shape_rash_to_recovery,
                                 rate=rate_rash_to_recovery
  )
  infections <- rash_time - rgamma(13,
                                   shape=shape_fever_to_rash,
                                   rate=rate_fever_to_rash
  )
  lag <- 11.6
  while (!check_if_epidemic(removals, infections, lag)) {
    removals <- rash_time + rgamma(13, 
                                   shape=shape_rash_to_recovery,
                                   rate=rate_rash_to_recovery
    )
    infections <- rash_time - rgamma(13,
                                     shape=shape_fever_to_rash,
                                     rate=rate_fever_to_rash
    )
  }
  
  # add augmented dates to dataset
  dataset$RecoveryDay <- removals
  dataset$FeverDay <- infections
  dataset$ExposedDay <- infections - lag
  
  
  # Saving ------------------------------------------------------------------
  
  write.table(dataset,
              paste0("abakaliki-modified-",d,".csv"),
              row.names=FALSE,
              sep=','
  )
  
}
