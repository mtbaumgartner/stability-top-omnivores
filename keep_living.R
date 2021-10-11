keep_living <- function(stab){
  
  # Function to take the simulated model from 'runODE' and
  # remove vanishing species and their parameters
  
  stab$web <- stab$web[stab$alive,stab$alive]
  stab$S <- ncol(stab$web)
  stab$C <- sum(stab$web > 0) / stab$S^2
  stab$TL <- stab$TL[stab$alive]
  stab$xi <- stab$xi[stab$alive]
  stab$ei <- stab$ei[stab$alive]
  stab$wi <- stab$wi[stab$alive]
  stab$ai <- stab$ai[stab$alive]
  stab$Mi <- stab$Mi[stab$alive]
  stab$ODEsteps <- data.frame(stab$ODEsteps[,1], stab$ODEsteps[,stab$alive])
  stab$Bs <- stab$Bs[stab$alive]
  stab$alive <- stab$alive[stab$alive]
  
  return(stab)
}
