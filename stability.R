stability <- function(comm){
  
  # Ives et al. (2003) - Ecological Monographs
  detB <- det(comm) # determinant
  lamb_max <- Re(eigen(comm)$values)[1] # lambda max
  lamb_min <- rev(Re(eigen(comm)$values))[1] # lambda min
  abs_lamb_max <- abs(lamb_max) # absolute lambda max
  #maxBkB <- max(abs(eigen(kronecker(comm, comm))$values)) # lambda max of Kronecker product
  
  # Arnoldi et al. (2016) - J Theor Biol
  Rinf <- -Re(eigen(comm)$values)[1]
  R0 <- -0.5 * eigen(comm + t(comm))$values[1]
  
  out <- data.frame(detB, lamb_max, lamb_min, abs_lamb_max, Rinf, R0)
  return(out)
}

