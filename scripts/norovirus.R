# Norovirus in English schools
# Seth Temple, sethtem@umich.edu

# library(pblas)
# library(peirrs)

library(outbreaks)

# prepare data
{
  noro <- norovirus_derbyshire_2001_school 
  noro.sick <- noro[noro$start_illness>0,]
  i <- noro.sick$start_illness + rnorm(dim(noro.sick)[1],0,sd=0.1)
  r <- noro.sick$end_illness + rnorm(dim(noro.sick)[1],0,sd=0.1)
  N <- dim(noro)[1]  
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
