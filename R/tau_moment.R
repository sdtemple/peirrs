#' Expected value for tau with different missingness
#'
#' Compute expectation of tau with (in)complete observations.
#'
#' @param rk numeric: removal time rk
#' @param rj numeric: removal time rj
#' @param ik numeric: infection time ik
#' @param ij numeric: infection time ij
#' @param lambdak numeric: removal rate of k
#' @param lambdaj numeric: removal rate of j
#' @param med bool: use median imputation if true
#'
#' @return numeric: expectation of some pair-wise tau
#'
#' @export
tau_moment <- function(rk,rj,ik,ij,lambdak,lambdaj, med=T){

  # Hypoexponential CDF with two different rates
  hypo2.cdf <- function(x, lambdak, lambdaj){
    if(lambdak==lambdaj){
      return(pgamma(x,shape=2,rate=lambdaj))
    }
    first <- lambdaj / (lambdaj - lambdak) * exp(-lambdak*x)
    second <- lambdak / (lambdak-lambdaj) * exp(-lambdaj*x)
    return(1 - first - second)
  }

  # Hypoexponential survival with two different rates
  hypo2.sf <- function(x, lambdak, lambdaj){
    return(1-hypo2.cdf(x,lambdak,lambdaj))
  }

  # Expected tau: rk, rj, ik, ij observed
  E.tau.rk.rj.ik.ij <- function(rk,rj,ik,ij,lambdak,lambdaj, med=T){
    # med does nothing
    val <- 0
    if(ij<ik){
      val <- 0
    } else if(ij > rk){
      val <- rk - ik
    } else{
      val <- ij - ik
    }
    return(val)
  }

  # Expected tau: rk, ik, ij observed
  E.tau.rk.ik.ij <- function(rk,rj,ik,ij,lambdak,lambdaj, med=T){
    # med does nothing
    val <- 0
    if(ij<ik){
      val <- 0
    } else if(ij > rk){
      val <- rk - ik
    } else{
      val <- ij - ik
    }
    return(val)
  }

  # Expected tau: rj, ik, ij observed
  E.tau.rj.ik.ij <- function(rk,rj,ik,ij,lambdak,lambdaj,med=T){
    # med does nothing
    val <- 0
    if(ij<ik){
      val <- 0
    } else if(ij > ik){
      H4 <- pexp(ij-ik,rate=lambdak,lower.tail=FALSE)
      H5 <- (exp(-lambdak*(ij-ik))*(lambdak*(ik-ij)-1)+1)/lambdak
      val <- H4 + H5
    }
    return(val)
  }

  # Expected tau: rk, rj, ik observed
  E.tau.rk.rj.ik <- function(rk,rj,ik,ij,lambdak,lambdaj, med=T){
    # med does nothing
    val <- 0
    if(rj < ik){
      val <- 0
    } else if(rj < rk){
      H6 <- (exp(-lambdaj*(rj-ik))-lambdaj*(ik-rj)-1)/lambdaj
      val <- H6
    } else{
      H8 <- (rk-ik)*pexp(rj-rk,rate=lambdaj,lower.tail=TRUE)
      H7 <- (exp(-lambdaj*rj)*(exp(lambdaj*ik)+exp(lambdaj*rk)*(lambdaj*(rk-ik)-1)))/lambdaj
      val <- H7+H8
    }
    return(val)
  }

  # Expected tau: rk, rj, ij observed
  E.tau.rk.rj.ij <- function(rk,rj,ik,ij,lambdak,lambdaj, med=T){
    median.scalar = 1
    if(med){median.scalar=log(2)}
    val <- 0
    if(ij < rk){
      val <- (exp(lambdak*ij)-lambdak*ij - 1)*exp(-lambdak*rk)/lambdak
    } else{
      val <- 1 / lambdak * median.scalar
    }
    return(val)
  }

  # Expected tau: rk, rj observed
  E.tau.rk.rj <- function(rk,rj,ik,ij,lambdak,lambdaj,med=T){
    median.scalar=1
    if(med){median.scalar=log(2)}
    val <- 0
    if(rj < rk){
      H1 <- pexp(rk-rj,rate=lambdak,lower.tail=FALSE)*lambdaj/lambdak/(lambdak+lambdaj)
      val <- H1
    } else{
      H2 <- pexp(rj-rk,rate=lambdaj,lower.tail=FALSE)*lambdaj/lambdak/(lambdak+lambdaj)
      H3 <- pexp(rj-rk,rate=lambdaj,lower.tail=TRUE)/lambdak
      val <- H2 + H3
    }
    return(val*median.scalar)
  }

  # Expected tau: rk, ij observed
  # rj is not useful if ij is known
  E.tau.rk.ij <- function(rk,rj,ik,ij,lambdak,lambdaj,med=T){
    return(E.tau.rk.rj.ij(rk,rj,ik,ij,lambdak,lambdaj,med))
  }

  # Expected tau: ik, ij observed
  # rj is not useful if ij is known
  E.tau.ik.ij <- function(rk,rj,ik,ij,lambdak,lambdaj,med=T){
    return(E.tau.rj.ik.ij(rk,rj,ik,ij,lambdak,lambdaj,med))
  }


  # integral from a to b
  integral.upp.low <- function(func,low,upp,rate1,rate2,rate3){
    return(func(upp,rate1,rate2,rate3)-func(low,rate1,rate2,rate3))
  }

  # Expected tau: rj, ik observed
  # This integral problem is extremely hard
  E.tau.rj.ik <- function(rk,rj,ik,ij,lambdak,lambdaj,med=T){
    # med does nothing
    val <- 0
    func1 <- function(x,rate1,rate2,rate3){
      return((1-(rate2-rate1)*x)*exp((rate2-rate1)*x))
    }
    func2 <- function(x,rate1,rate2,rate3){
      return(exp((rate2-rate1)*x)/(rate2-rate1))
    }
    if(rj<ik){
      val <- 0
    } else {
      # ij > rk case
      if(lambdak==lambdaj){
        H14 <- lambdak*(rj^2-ik^2)/2
        H13 <- rj - ik
        H27 <- rj - ik
      } else{
        H15 <- lambdak/((lambdaj-lambdak)^2)*integral.upp.low(func1,ik,rj,lambdak,lambdaj,0)
        H14 <- H15
        H13 <- integral.upp.low(func2,rj,ik,lambdak,lambdaj,0)
        H27 <- (exp((lambdaj-lambdak)*rj)-exp((lambdaj-lambdak)*ik))/(lambdaj-lambdak)
      }
      H12 <- H13 + H14
      H11 <- (1+lambdak*ik)*exp(-lambdak*ik)*(exp(lambdaj*rj)-exp(lambdaj*ik))/lambdaj
      H10 <- (H11-H12)/(lambdak^2)
      H9 <- ((exp(lambdaj*rj)-exp(lambdaj*ik))*exp(-lambdak*ik)/lambdaj - H27)*ik/lambdak
      E1 <- (H10-H9)*lambdak*lambdaj*exp(-lambdaj*rj)*exp(lambdak*ik)

      # ik < ij < rj < rk case
      H18 <- exp(-lambdak*rj)*integral.upp.low(func1,rj,ik,0,lambdaj,0)/lambdak/(lambdaj^2)
      H19 <- ik/lambdak/lambdaj*exp(-lambdak*rj)*(exp(lambdaj*rj)-exp(lambdaj*ik))
      H16 <- lambdak*lambdaj*exp(-lambdaj*rj)*exp(lambdak*ik)*(H18-H19)

      # ik < ij < rk < rj case
      if(lambdaj==lambdak){
        H22 <- (rj^2-ik^2)/2
        H25 <- rj-ik
      } else{
        H24 <- integral.upp.low(func1,ik,rj,lambdak,lambdaj,0)/((lambdaj-lambdak)^2)
        H22 <- H24
        H25 <- (exp((lambdaj-lambdak)*rj)-exp((lambdaj-lambdak)*ik))/(lambdaj-lambdak)
      }
      H26 <- (exp(lambdaj*rj)-exp(lambdaj*ik))*exp(-lambdak*rj)/lambdaj
      H21 <- ik*(H25-H26)/lambdak
      H23 <- integral.upp.low(func1,rj,ik,0,lambdaj,0)/(lambdaj^2)*exp(-lambdak*rj)
      H20 <- (H22-H23)/lambdak
      H17 <- lambdak*lambdaj*exp(-lambdaj*rj)*exp(lambdak*ik)*(H20-H21)

      # combine as ik < ij < rk
      E2 <- H16 + H17

      # final result
      val <- E1 + E2
    }

    return(val)
  }

  itr <- 0

  # should not both be NA
  # means individual not infected
  if(is.na(rk)){
    if(is.na(ik)){
      itr <- 2
      val <- NULL
    }
  }

  # should not both be NA
  # means individual not infected
  if(is.na(rj)){
    if(is.na(ij)){
      itr <- 2
      val <- NULL
    }
  }

  # case 1
  if(is.na(ij)){
    if(is.na(ik)){
      itr <- itr + 1
      val <- E.tau.rk.rj(rk, rj, ik, ij, lambdak, lambdaj, med)
    }
  }
  # case 2
  if(is.na(ij)){
    if(is.na(rk)){
      itr <- itr + 1
      val <- E.tau.rj.ik(rk, rj, ik, ij, lambdak, lambdaj, med)
    }
  }
  # case 3
  if(is.na(ik)){
    if(is.na(rj)){
      itr <- itr + 1
      val <- E.tau.rk.ij(rk, rj, ik, ij, lambdak, lambdaj, med)
    }
  }
  # case 4
  if(is.na(rk)){
    if(is.na(rj)){
      itr <- itr + 1
      val <- E.tau.ik.ij(rk, rj, ik, ij, lambdak, lambdaj, med)
    }
  }
  # case 5
  if(is.na(rk)&itr==0){
    itr <- itr + 1
    val <- E.tau.rj.ik.ij(rk, rj, ik, ij, lambdak, lambdaj, med)
  }
  # case 6
  if(is.na(rj)&itr==0){
    itr <- itr + 1
    val <- E.tau.rk.ik.ij(rk, rj, ik, ij, lambdak, lambdaj, med)
  }
  # case 7
  if(is.na(ik)&itr==0){
    itr <- itr + 1
    val <- E.tau.rk.rj.ij(rk, rj, ik, ij, lambdak, lambdaj, med)
  }
  # case 8
  if(is.na(ij)&itr==0){
    itr <- itr + 1
    val <- E.tau.rk.rj.ik(rk, rj, ik, ij, lambdak, lambdaj, med)
  }

  # case 9
  if(!is.na(rk+rj+ik+ij)){
    itr <- itr + 1
    val <- E.tau.rk.rj.ik.ij(rk, rj, ik, ij, lambdak, lambdaj, med)
  }

  result <- tryCatch({
    if(is.null(val)){
      stop("Non-case data is input")
    }
    return(val)
  }, warning=function(w){
    message("This is a warning that non-case data is input")
    print(c(rk,rj,ik,ij))
    return(NULL)
  }, error=function(e){
    message("This is an error that non-case data is input")
    print(c(rk,rj,ik,ij))
    stop("Non-case data is input")
  }, finally={
  })

  result <- tryCatch({
    if(val<0){
      stop("tau is negative")
    }
    return(val)
  }, warning=function(w){
    message("This is a warning that tau is negative")
    return(NULL)
  }, error=function(e){
    message("This is an error that tau is negative")
    stop("tau is negative")
  }, finally={
  })

  result <- tryCatch({
    if(itr!=1){
      stop("Too many cases triggered")
    }
    return(val)
  }, warning=function(w){
    message("This is a warning that too many cases triggered")
    return(NULL)
  }, error=function(e){
    message("This is an error that too many cases triggered")
    stop("Too many cases triggered")
  }, finally={
  })

  return(val)

}
