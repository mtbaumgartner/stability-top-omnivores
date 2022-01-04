#' Main script with simulations for 'The sensitivity of complex food webs to the loss
#' of top omnivorous species' by Matheus T. Baumgartner and Lucas del Bianco Faria
#' 
#' This script was coded by MTB and runs simulations of complex food webs 
#' assuming allometric constraints under the Allometric Trophic Network (ATN) 
#' framework.

# Load required packages and functions
rm(list = ls())
load("inits.RData")
load("stable_communities.RData")

library(deSolve)
library(rootSolve)
library(numDeriv)
library(foreach)
library(doParallel)
library(doSNOW)
source("identify_omnivores.R")
source("runODE.R")
source("keep_living.R")
source("ODEstruc.R")

# Keep only living species
Nwebs <- length(stable_comm)
stable_comm <- lapply(stable_comm, keep_living)

#---- Remove random species ----
# Ensure that they are not omnivores
# Select species that resembles the trophic level of removed omnivore
load("without_omnivores.RData")

random_removed <- data.frame(ID = seq(Nwebs),
                             SPremoved = 0, TLremoved = 0,
                             TLdelta = 0)

without_random <- list()
for(i in 1:Nwebs){
  foodWeb <- stable_comm[[i]] # pass community
  
  omniv <- identify_omnivores(foodWeb$web)
  if(nrow(omniv) > 0){ # at least one omnivore in web
    not_omniv <- seq(foodWeb$S)[-omniv$sp] # identify non-omnivores
    
    delta_TL <- abs(foodWeb$TL[not_omniv] - omnivores_removed$TLremoved[i])
    sel_species <- not_omniv[which.min(delta_TL)]
    
    # Record information
    random_removed$SPremoved[i] <- sel_species # removed species
    random_removed$TLremoved[i] <- foodWeb$TL[sel_species] # trophic level
    random_removed$TLdelta[i] <- delta_TL[which.min(delta_TL)] # delta in trophic level
    
    # Remove from community
    foodWeb$alive[sel_species] <- FALSE
    foodWeb <- keep_living(foodWeb)
    
    # Pass to main list
    without_random[[i]] <- foodWeb
    cat("Random species removed from community", i, "\r")
    
  }else{
    without_random[[i]] <- NULL
  }
}
#save.image("without_random.RData")


#---- Set up restarting communities
rm(list = ls())
load("without_random.RData")
rm(stable_comm) # remove to save memory
rm(without_omnivores)

# Keep webs with at least one omnivore
without_random <- without_random[sapply(without_random, function(x){!is.null(x)})]

# Update number of webs
Nwebs <- length(without_random)

for(i in 1:Nwebs){ # initial set up for function 'runODE'
  without_random[[i]] <- without_random[[i]][1:15]
}


#---- Re-iterate each community without random species for 5000 time steps ----
# parallel iterations
numCores = 3
registerDoParallel(numCores)
cl <- makeSOCKcluster(numCores)
registerDoSNOW(cl)

# progress bar
pb <- txtProgressBar(min = 1, max = Nwebs, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

without_random_stable <- foreach(i = 1:Nwebs,
                                 .packages = "deSolve",
                                 .options.snow = opts) %dopar%
  {
    runODE(mod = without_random[[i]], times = 5000)
  }
#save(without_random_stable, file = "without_random_stable.RData")


#---- Calculate community matrices for stable communities with random species removed ----
# keep only living species
without_random_stable <- lapply(without_random_stable, keep_living)


# Calculate Jacobian matrix
jacobians_random <- foreach(i = 1:Nwebs,
                            .packages = "numDeriv",
                            .options.snow = opts) %dopar%
  {
    jacobian(ODEstruc, without_random_stable[[i]]$Bs,
             pars = without_random_stable[[i]], method = "simple")
  }
#save(jacobians_random, file = "jacobians_random.RData")


# Calculate stability metrics
metrics_without_random <- data.frame(ID = omnivores_removed$webID[omnivores_removed$Nomniv > 0])
metrics_without_random$S <- sapply(without_random_stable, function(x){x$S}) # Species richness
metrics_without_random$C <- sapply(without_random_stable, function(x){x$C}) # Connectance
metrics_without_random$q <- sapply(without_random_stable, function(x){x$q}) # Hill's exponent
metrics_without_random$rand_TL <- random_removed$TLremoved[omnivores_removed$Nomniv > 0] # TL of removed species
metrics_without_random$delta_TL <- random_removed$TLdelta[omnivores_removed$Nomniv > 0] # delta TL between random and omnivore
metrics_without_random$maxeig_Re <- sapply(jacobians_random, function(x){Re(eigen(x)$value[1])})
metrics_without_random$maxeig_Im <- sapply(jacobians_random, function(x){Im(eigen(x)$value[1])})


#save(metrics_without_random, file = "metrics_without_random.RData")
