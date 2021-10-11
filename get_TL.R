get_TL <- function(web){
  # Based on Owen Petchey's function
  # Calculates the prey-averaged trophic level based on
  # the matrix inversion method by Levine (1980 - J Theor Biol)
  # takes predation matrix with consumers in columns
  
  diag(web)<-0 # Remove cannibal links
  
  tweb <- t(web)
  sp<- length(tweb[,1])
  rs <- rowSums(tweb)
  for(i in 1:sp) tweb[i,tweb[i,]==1] = 1/rs[i]
  nb.TL <- try(solve(diag(sp) - tweb), T)
  if(class(nb.TL)=="try-error")return(rep(NA, sp))
  if(class(nb.TL)!="try-error")return(rowSums(nb.TL))
}