#' Main script with simulations for 'The sensitivity of complex food webs to the loss
#' of top omnivorous species' by Matheus T. Baumgartner and Lucas del Bianco Faria
#' 
#' This script was coded by MTB and runs simulations of complex food webs 
#' assuming allometric constraints under the Allometric Trophic Network (ATN) 
#' framework.

# Load required packages and functions
rm(list = ls())

library(magrittr)
library(igraph)
library(deSolve)
library(rootSolve)
library(numDeriv)
library(pbapply)
library(foreach)
library(doParallel)
library(doSNOW)
source("get_TL.R")
source("niche_model.R")
source("make_model.R")
source("runODE.R")
source("ODEstruc.R")
source("keep_living.R")

#---- Initial setups ----
Nwebs = 1e4 # Approx. 30 webs for each S*C combination (change to 3e4!!!)
inits <- data.frame(S = sample(15:50, Nwebs, replace = T),
                    C = runif(Nwebs, min = 0.1, max = 0.3))
plot(inits)
table(inits$S)
#save(inits, file = "inits.RData")


#---- Generate starting communities ----
#load("inits.RData")

start_comm <- list()

t.begin <- Sys.time()
for(i in 1:Nwebs){
  run <- NULL
  
  while(is.null(run)){
    run <- tryCatch(
      {
        make_model(S = inits$S[i], C = inits$C[i])
      },
      error = function(x){ NULL }
    )
  }
  
  start_comm[[i]] <- run
  cat("Starting community", i, "generated. Elapsed time:", Sys.time() - t.begin, "\n")
}
#save(start_comm, file = "starting_communities.RData")


#---- Iterate each starting community for 5000 time steps ----
#load("starting_communities.RData")

# Parallel runs
numCores = 3
registerDoParallel(numCores)
cl <- makeSOCKcluster(numCores)
registerDoSNOW(cl)

# Progress bar
pb <- txtProgressBar(min = 1, max = Nwebs, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Run simulations
stable_comm <- foreach(i = 1:Nwebs,
                       .packages = "deSolve",
                       .options.snow = opts) %dopar%
  {
    runODE(mod = start_comm[[i]], times = 300) # Change to 5000
  }
#save(stable_comm, file = "stable_communities.RData")

# Check the distribution of the proportion of living species
hist(sapply(stable_comm, function(x){ mean(x$alive) }), xlim = c(0,1))


#---- Calculate community matrices for stable communities ----
#load("stable_communities.RData")

# Keep only living species
stable_comm <- pblapply(stable_comm, keep_living)

# Plot initial S x stable S
plot(sapply(stable_comm, function(x){x$S}) ~ inits$S,
     xlab = "Starting", ylab = "Stable", pch = 16)
abline(a = 0, b = 1, lwd = 4, col = "red")

# Calculate Jacobian matrices for stable states
jacobians_stable <- foreach(i = 1:Nwebs,
                            .packages = "numDeriv",
                            .options.snow = opts) %dopar%
  {
    jacobian(ODEstruc, stable_comm[[i]]$Bs,
             pars = stable_comm[[i]], method = "simple")
  }
#save(jacobians_stable, file = "jacobians_stable.RData")

# Calculate stability metrics
metrics_stable <- data.frame(ID = seq(Nwebs))
metrics_stable$S <- sapply(stable_comm, function(x){x$S}) # Species richness
metrics_stable$C <- sapply(stable_comm, function(x){x$C}) # Connectance
metrics_stable$q <- sapply(stable_comm, function(x){x$q}) # Hill's exponent
metrics_stable$prop_alive <- metrics_stable$S / inits$S
metrics_stable$maxeig_Re <- sapply(jacobians_stable, function(x){Re(eigen(x)$values[1])})
metrics_stable$maxeig_Im <- sapply(jacobians_stable, function(x){Im(eigen(x)$values[1])})

hist(metrics_stable$maxeig_Re)
plot(maxeig_Re ~ S, metrics_stable)
hist(metrics_stable$prop_alive)

#save(metrics_stable, file = "metrics_stable.RData")
