#' Main script with simulations for 'The sensitivity of complex food webs to the loss
#' of top omnivorous species' by Matheus T. Baumgartner and Lucas del Bianco Faria
#' 
#' This script was coded by MTB and runs simulations of complex food webs 
#' assuming allometric constraints under the Allometric Trophic Network (ATN) 
#' framework. Any question should be addressed to matheustbs@gmail.com

# Load required packages and functions
rm(list = ls())

library(magrittr)
library(igraph)
library(deSolve)
library(rootSolve)
library(foreach)
library(doParallel)
library(doSNOW)
source("body_masses.R")
source("get_growth.R")
source("get_TL.R")
source("niche_model.R")
source("make_model.R")
source("runODE.R")
source("keep_living.R")
source("ODEmodel.R")
source("stability.R")

#---- Initial setups ----
Nwebs = 300
inits <- data.frame(S = sample(15:50, Nwebs, replace = T),
                    C = runif(Nwebs, min = 0.1, max = 0.3))
plot(inits)
table(inits$S)
#save.image("inits.RData")


#---- Generate the 1000 starting communities ----
load("inits.RData")

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
#save.image("starting_communities.RData")


#---- Iterate each starting community for 1000 time steps ----
load("starting_communities.RData")

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
    runODE(mod = start_comm[[i]], times = 100)
  }
#save(stable_comm, file = "stable_communities.RData")

# Check the distribution of the proportion of living species
hist(sapply(stable_comm, function(x){ mean(x$alive) }), xlim = c(0,1))


#---- Calculate community matrices for stable communities ----
# Keep only living species
stable_comm_alive <- list()
for(i in 1:Nwebs){
  stable_comm_alive[[i]] <- keep_living(stable_comm[[i]])
}
#save(stable_comm_alive, file = "stable_comm_alive.RData")

# Plot initial S x stable S
plot(sapply(stable_comm_alive, function(x){x$S}) ~ inits$S,
     xlab = "Stable", ylab = "Starting", pch = 16)
abline(a = 0, b = 1, lwd = 4, col = "red")

# Calculate Jacobian matrix
jacobians_stable <- list()
for(i in 1:Nwebs){
  jacobians_stable[[i]] <- jacobian.full(stable_comm_alive[[i]]$Bs, ODEmodel,
                                         parms = stable_comm_alive[[i]], pert = 1e-2)
  cat("Jacobian matrix calculated for stable community", i, "\n")
}
#save(jacobians_stable, file = "jacobians_stable.RData")


# Calculate stability metrics
metrics_stable <- data.frame(ID = seq(Nwebs))
metrics_stable$S <- sapply(stable_comm_alive, function(x){x$S}) # Species richness
metrics_stable$C <- sapply(stable_comm_alive, function(x){x$C}) # Connectance
metrics_stable$q <- sapply(stable_comm_alive, function(x){x$q}) # Hill's exponent
metrics_stable <- data.frame(metrics_stable, stability(jacobians_stable[[1]]))
for(i in 1:Nwebs){
  metrics_stable[i,5:10] <- stability(jacobians_stable[[i]])
  cat("Stability metrics for stable community", i, "calculated\n")
}
#save(metrics_stable, file = "metrics_stable.RData")
