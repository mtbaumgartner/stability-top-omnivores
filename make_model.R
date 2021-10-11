#' Make community model for simulations
#' S = number of species
#' C = connectance (L/S^2)
#' propPP = minimum and maximum proportion of primary producers
#' maxit = maximum number of trials

make_model <- function(S, C, propPP = c(0.1, 0.3), maxit = 10){
  
  comm <- T
  nit <- 1
  while(comm){ # loop until proportion of producers is achieved
    
    # niche model based on Williams and Martinez (2000) - Nature
    NM <- niche_model(S = S, C = C)
    web <- NM$web
    sp_niche <- NM$sp.niche
    #web <- niche_model(S = S, C = C)
    
    # prey-averaged trophic level following Levine (1980) - J Theor Biol
    suppressWarnings(
      TL <- get_TL(web)
    )
    
    # check if minimum proportion of producers is achieved
    if(all(!is.na(TL)) & mean(TL == 1) >= 0.1 & mean(TL == 1) <= 0.3){
      comm <- F
    }
    
    # check if maximum TL < S
    if(max(TL) > S){
      comm <- F
    }
    
    nit <- nit + 1
    if(nit == maxit){
      stop("Maxit reached")
    }
    
  }
  
  #---- Model Parameters ----
  Mi = ai = ei = rep(NA, S)
  for(i in 1:S){
    
    # Depending on trophic level:
    # body masses (Mi)
    # specific allometric constant (ai)
    # assimilation efficiency (ei)
    if(TL[i] == 1){
      Mi[i] = 1
      ai[i] = 1
      ei[i] = 0.83
    }
    if(TL[i] == 2){
      #Mi[i] = 10 ^ (runif(1, 2, 10) * sp_niche[i])
      Mi[i] = rlnorm(1, 0.65, 1.52)
      ai[i] = 0.314
      ei[i] = 0.45
    }
    if(TL[i] > 2){
      #Mi[i] = 10 ^ (runif(1, 2, 10) * sp_niche[i])
      ndraw = round(TL[i]) - 1
      Mi[i] = prod(rlnorm(ndraw, 2.73, 1.60))
      ai[i] = 0.88
      ei[i] = 0.85
    }
    
  }
  
  # consumer preference (wi)
  wi <- 1/colSums(web)
  wi[is.infinite(wi)] <- 0
  
  # metabolic rate (xi)
  xi <- ai * Mi^-0.25
  
  # starting biomass
  Bs <- rep(1, S)
  
  pars <- list(S = S, C = C, web = web, TL = TL,
               Bi = 10, xi = xi, y = 10, ei = ei,
               K = 10^6, wi = wi, B0 = 1500, ai = ai,
               q = runif(1, 1, 2), Mi = Mi, Bs = Bs)
  
  return(pars)
}
