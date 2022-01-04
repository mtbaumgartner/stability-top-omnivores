runODE <- function(mod, times){
  
  Bs <- mod$Bs
  
  #---- build ODE system ----
  ODEmodel <- function(t, Bs, pars){
    
    if(any(is.nan(Bs))){ Bs[is.nan(Bs)] <- 0 }
    if(any(Bs < 1e-30)){ Bs[Bs < 1e-30] <- 0 }
    
    dBs <- Bs
    
    # Functional response
    Fji <- matrix(0, nrow = pars$S, ncol = pars$S)
    for(j in 1:pars$S){
      Fji[j,] <- (pars$wi * Bs[j]^pars$q) /
        (pars$B0^pars$q + Bs * pars$B0^pars$q +
           sum(pars$wi[-j] * Bs[-j]^pars$q))
    }
    # for(j in 1:pars$S){
    #   for(i in 1:pars$S){
    #     Fji[j,i] <- (pars$wi[i] * Bs[j]^pars$q) /
    #       (pars$B0^pars$q + Bs[i] * pars$B0^pars$q +
    #          sum(pars$wi[i] * Bs[-j]^pars$q))
    #   }
    # }

    # ODE solution
    for(i in 1:pars$S){

      if(pars$TL[i] == 1){ # producers
        #dBs[i] <- pars$ei[i] * pars$xi[i] * Bs[i] * exp((1-sum(Bs))/pars$K) -
        dBs[i] <- pars$ei[i] * pars$xi[i] * Bs[i] * (1 - (sum(BS) / pars$K)) -
          sum((pars$y/pars$ei[i]) * Bs * Fji[i,]) -
          (1 - pars$ei[i]) * pars$xi[i] * Bs[i]
      }

      if(pars$TL[i] > 1){ # consumers
        dBs[i] <- sum(pars$xi[i] * pars$y * Bs[i] * Fji[i,]) -
          sum(pars$xi * (pars$y/pars$ei) * Bs * Fji[i,]) -
          pars$xi[i] * Bs[i]
      }
    }
    
    dBs[is.nan(Bs)] <- 0
    dBs[Bs < 1e-30] <- 0
    
    return(list(as.vector(dBs)))
  }
  
  # run ODE model
  output <- ode(Bs, times = seq(0, times, by = 0.5), ODEmodel, parms = mod,
                method = "radau", atol = 1e-4, rtol = 1e-4)
  
  mod$Bs <- output[nrow(output),-1]
  alive <- mod$Bs > 1e-30
  
  # keep time series resumed
  output <- output[floor(seq(1, nrow(output), length.out = 100)),]
  
  out <- c(mod, list(ODEsteps = output, alive = alive))
  
  return(out)
}