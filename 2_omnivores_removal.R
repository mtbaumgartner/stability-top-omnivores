rm(list = ls())
load("inits.RData")
load("starting_communities.RData")
load("stable_comm_alive.RData")
load("jacobians_stable.RData")
load("metrics_stable.RData")

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
source("identify_omnivores.R")
source("niche_model.R")
source("make_model.R")
source("runODE.R")
source("keep_living.R")
source("ODEmodel.R")
source("stability.R")


#---- Identify omnivores, remove and record trophic level ----
omnivores_removed <- data.frame(webID = seq(Nwebs),
                                Nomniv = 0, SPremoved = 0,
                                TLremoved = 0)

without_omnivores <- list()
for(i in 1:Nwebs){
  without_omnivores[[i]] <- stable_comm_alive[[i]] # transfer community
  omniv <- identify_omnivores(without_omnivores[[i]]$web) # identify omnivores
  extinct <- sample(omniv$sp, 1) # select species to extinct
  
  omnivores_removed$Nomniv[i] <- nrow(omniv) # Record the number of omnivores
  omnivores_removed$SPremoved[i] <- extinct
  omnivores_removed$TLremoved[i] <- without_omnivores[[i]]$TL[extinct]
  
  # Remove from community
  without_omnivores[[i]]$alive[names(without_omnivores[[i]]$alive) == extinct] <- FALSE
  without_omnivores[[i]] <- keep_living(without_omnivores[[i]])
  cat("Omnivore removed from community", i, "\n")
}
#save.image("without_omnivores.RData")


#---- Set up restarting communities ----
rm(list = ls())
load("without_omnivores.RData")

for(i in 1:Nwebs){
  without_omnivores[[i]] <- without_omnivores[[i]][names(without_omnivores[[i]]) %in%
                                                     names(start_comm[[i]])]
}

#---- Re-iterate each community without omnivores for 1000 time steps ----
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
    runODE(mod = without_omnivores[[i]], times = 100)
  }
#save(without_omnivores_stable, file = "without_omnivores_stable.RData")


#---- Calculate community matrices for stable communities with omnivores removed ----
# keep only living species
for(i in 1:Nwebs){
  without_omnivores_stable[[i]] <- keep_living(without_omnivores_stable[[i]])
}

# Calculate Jacobian matrix
jacobians_omnivores <- list()
for(i in 1:Nwebs){
  jacobians_omnivores[[i]] <- jacobian.full(without_omnivores_stable[[i]]$Bs, ODEmodel,
                                            parms = without_omnivores_stable[[i]], pert = 1e-2)
  cat("Jacobian matrix calculated for community without omnivores", i, "\n")
}
#save(jacobians_omnivores, file = "jacobians_omnivores.RData")


# Calculate stability metrics
metrics_without_omnivores <- data.frame(ID = seq(Nwebs))
metrics_without_omnivores$S <- sapply(without_omnivores_stable, function(x){x$S}) # Species richness
metrics_without_omnivores$C <- sapply(without_omnivores_stable, function(x){x$C}) # Connectance
metrics_without_omnivores$q <- sapply(without_omnivores_stable, function(x){x$q}) # Hill's exponent
metrics_without_omnivores <- data.frame(metrics_without_omnivores, omnivores_removed[,-1]) # Removed omnivore
metrics_without_omnivores <- data.frame(metrics_without_omnivores, stability(jacobians_omnivores[[1]]))
for(i in 1:Nwebs){
  metrics_without_omnivores[i,8:13] <- stability(jacobians_omnivores[[i]])
  cat("Stability metrics for community without omnivores", i, "calculated\n")
}
#save(metrics_without_omnivores, file = "metrics_without_omnivores.RData")
