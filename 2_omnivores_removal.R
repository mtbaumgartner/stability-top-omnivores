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
library(pbapply)
library(foreach)
library(doParallel)
library(doSNOW)
source("identify_omnivores.R")
source("runODE.R")
source("keep_living.R")
source("ODEstruc.R")

# Keep only living species
Nwebs <- length(stable_comm)
stable_comm <- pblapply(stable_comm, keep_living)

#---- Identify omnivores, remove and record their trophic level ----
omnivores_removed <- data.frame(webID = seq(Nwebs),
                                Nomniv = 0, SPremoved = 0,
                                TLremoved = 0, TYremoved = "NA")

without_omnivores <- list()
for(i in 1:Nwebs){
  foodWeb <- stable_comm[[i]] # pass community
  
  omniv <- identify_omnivores(foodWeb$web)
  if(nrow(omniv) > 0){ # if there are omnivores in the web
    sel_species <- omniv[sample(nrow(omniv), 1),] # select species to remove
    
    # Record information
    omnivores_removed$Nomniv[i] <- nrow(omniv) # number of omnivores in the web
    omnivores_removed$SPremoved[i] <- sel_species$sp # identity of removed species
    omnivores_removed$TLremoved[i] <- foodWeb$TL[sel_species$sp] # trophic level
    omnivores_removed$TYremoved[i] <- ifelse(sel_species$CL == 1, "CL", "MP") # type
    
    # Remove from community
    foodWeb$alive[names(foodWeb$alive) == sel_species$sp] <- FALSE
    foodWeb <- keep_living(foodWeb)
    
    # Pass to main list
    without_omnivores[[i]] <- foodWeb
    cat("Omnivore removed from community", i, "\r")
    
  }else{
    without_omnivores[[i]] <- NULL
  }
}
#save.image("without_omnivores.RData")


#---- Set up communities to run simulations ----
rm(list = ls())
load("without_omnivores.RData")
rm(stable_comm) # remove to save memory space

# Keep webs with at least one omnivore
without_omnivores <- without_omnivores[sapply(without_omnivores, function(x){!is.null(x)})]

# Update number of webs
Nwebs <- length(without_omnivores)

for(i in 1:Nwebs){ # initial set up for function 'runODE'
  without_omnivores[[i]] <- without_omnivores[[i]][1:15]
}


#---- Re-iterate each community without omnivores for 5000 time steps ----
# parallel iterations
numCores = 3
registerDoParallel(numCores)
cl <- makeSOCKcluster(numCores)
registerDoSNOW(cl)

# progress bar
pb <- txtProgressBar(min = 1, max = Nwebs, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

without_omnivores_stable <- foreach(i = 1:Nwebs,
                                    .packages = "deSolve",
                                    .options.snow = opts) %dopar%
  {
    runODE(mod = without_omnivores[[i]], times = 5000)
  }
#save(without_omnivores_stable, file = "without_omnivores_stable.RData")


#---- Calculate community matrices for stable communities with omnivores removed ----
# keep only living species
without_omnivores_stable <- lapply(without_omnivores_stable, keep_living)


# Calculate Jacobian matrix
jacobians_omnivores <- foreach(i = 1:Nwebs,
                               .packages = "numDeriv",
                               .options.snow = opts) %dopar%
  {
    jacobian(ODEstruc, without_omnivores_stable[[i]]$Bs,
             pars = without_omnivores_stable[[i]], method = "simple")
  }
#save(jacobians_omnivores, file = "jacobians_omnivores.RData")


# Calculate stability metrics
metrics_without_omnivores <- data.frame(ID = omnivores_removed$webID[omnivores_removed$Nomniv > 0])
metrics_without_omnivores$S <- sapply(without_omnivores_stable, function(x){x$S}) # Species richness
metrics_without_omnivores$C <- sapply(without_omnivores_stable, function(x){x$C}) # Connectance
metrics_without_omnivores$omniv_N <- omnivores_removed$Nomniv[omnivores_removed$Nomniv > 0] # number of omnivores
metrics_without_omnivores$omniv_TL <- omnivores_removed$TLremoved[omnivores_removed$Nomniv > 0] # TL of removed omnivores
metrics_without_omnivores$omniv_TY <- omnivores_removed$TYremoved[omnivores_removed$Nomniv > 0] # type of omnivory
metrics_without_omnivores$q <- sapply(without_omnivores_stable, function(x){x$q}) # Hill's exponent
metrics_without_omnivores$maxeig_Re <- sapply(jacobians_omnivores, function(x){Re(eigen(x)$value[1])})
metrics_without_omnivores$maxeig_Im <- sapply(jacobians_omnivores, function(x){Im(eigen(x)$value[1])})

#save(metrics_without_omnivores, file = "metrics_without_omnivores.RData")
