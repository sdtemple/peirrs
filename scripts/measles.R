# Measles in Hagelloch analysis
# Seth Temple, sethtem@umich.edu

# load in packages

{

# library(pblas)
# library(peirrs)

library(outbreaks)
library(lubridate)
  
}

# load in dataset
{
  measles <- measles_hagelloch_1861
  N <- dim(measles)[1] # Is this population size correct?
}

# convert the symptoms date to numeric
{
  x <- as.numeric(ymd(measles$date_of_prodrome))
  x <- x + rnorm(length(x),sd=0.1)
  minx <- min(x)
  x <- x - minx
  measles$numeric_procurement <- x
}

# convert the rash date to numeric
{
  x <- as.numeric(ymd(measles$date_of_rash))
  x <- x + rnorm(length(x),sd=0.1)
  x <- x - minx
  measles$numeric_rash <- x
}

# convert the date death to numeric
{
  x <- as.numeric(ymd(measles$date_of_death))
  x <- x + rnorm(length(x),sd=0.1)
  x <- x - minx
  measles$numeric_death <- x
}

# set the time of infection
{
  fixed.time.to.symptoms <- 10
  measles$numeric_infection <- measles$numeric_procurement - fixed.time.to.symptoms
}

# convert time of infectious and removal
{
  infectious.sd <- 0.5
  infectious.mean <- 4
  measles$numeric_infectious <- measles$numeric_rash - 
    rnorm(length(x),
          mean=infectious.mean,
          sd=infectious.sd)
  measles$numeric_removal <- measles$numeric_rash + 
    rnorm(length(x),
          mean=infectious.mean,
          sd=infectious.sd)
}

# get the final data
{
  scalar <- 1
  i <- measles$numeric_infectious
  mi <- min(i)
  i <- i - mi
  r <- measles$numeric_removal
  r <- r - mi
  i <- i/scalar
  r <- r/scalar
}

# analysis with complete data
{
  complete.mle.estimates <- mle.complete(r,i,N)
  complete.mle.estimates  
}


# analysis with pbla
{
  pbla.estimates <- peirr_pbla(r,N)
  pbla.estimates
}

# analysis with peirr pbla
{
  peirr.pbla.estimates <- peirr_likelihood(r,i,N)
  peirr.pbla.estimates  
}


# analysis with peirr pbla (50% missingness)
{
  p <- 0.5
  q <- 1.
  epi <- cbind(i,r)
  iepi <- decomplete_gsem(epi, p, q)
  iepi <- sort_gsem(iepi)
  r2 <- iepi[,2]
  i2 <- iepi[,1]
  
  peirr.pbla.estimates <- peirr_likelihood(r2,i2,N)
  peirr.pbla.estimates
}

# analysis with peirr pbla (90% missingness)
{
  p <- 0.1
  q <- 1.
  epi <- cbind(i,r)
  iepi <- decomplete_gsem(epi, p, q)
  iepi <- sort_gsem(iepi)
  r2 <- iepi[,2]
  i2 <- iepi[,1]
  
  peirr.pbla.estimates <- peirr_likelihood(r2,i2,N)
  peirr.pbla.estimates
}

# analysis with peirr tau (50% missingness)

{
  p <- 0.5
  q <- 0.5
  epi <- cbind(i,r)
  iepi <- decomplete_gsem(epi, p, q)
  iepi <- sort_gsem(iepi)
  r2 <- iepi[,2]
  i2 <- iepi[,1]
  peirr.tau.estimates <- peirr_tau(r2,i2,N)
  peirr.tau.estimates
}

# analysis with peirr tau (90% missingness)
{
  p <- 0.1
  q <- 0.5
  epi <- cbind(i,r)
  iepi <- decomplete_gsem(epi, p, q)
  iepi <- sort_gsem(iepi)
  r2 <- iepi[,2]
  i2 <- iepi[,1]
  peirr.tau.estimates <- peirr_tau(r2,i2,N)
  peirr.tau.estimates
}

# analysis with peirr tau (90% missingness, all infectious times)
{
  p <- 0.1
  q <- 1.
  epi <- cbind(i,r)
  iepi <- decomplete_gsem(epi, p, q)
  iepi <- sort_gsem(iepi)
  r2 <- iepi[,2]
  i2 <- iepi[,1]
  peirr.tau.estimates <- peirr_tau(r2,i2,N)
  peirr.tau.estimates
}

# analysis with peirr bayes


# there are different classes to consider


