#' This function takes an interaction matrix and identifies omnivores as:
#' Classic, Multi-resource, and Mutual predation (as in McLeod and Leroux 2020)
#' IMPORTANT: The interaction matrix must be binary with columns 'eating' rows

identify_omnivores <- function(i_matrix){
  
  # Outputs
  # omn <- data.frame(sp = seq(ncol(i_matrix)),
  #                   CL = 0, MP = 0, MR = 0)
  omn <- data.frame(sp = seq(ncol(i_matrix)),
                    CL = 0, MP = 0)
  
  
  # Count motifs by following triplets
  for(r in 1:nrow(i_matrix)){
    for(c in 1:ncol(i_matrix)){
      if(i_matrix[r,c] == 1){ # if "c" eats "r"
        for(x in 1:nrow(i_matrix)){
          if(i_matrix[x,c] == 1 & r != x){ # if "c" eats "x"
            if(i_matrix[x,r] == 1 & r != x){ # if "r" eats "x"
              if(i_matrix[c,r] == 1){ # if "r" eats "c"
                # Mutual predation
                omn[c,"MP"] <- 1
                omn[r,"MP"] <- 1
              }
              else{
                # Classic
                omn[c,"CL"] <- 1
              }
            }
            # else { # maybe multi-resource
            #   for(y in 1:nrow(i_matrix)){
            #     if(i_matrix[y,r] == 1 & x != y){ # if "r" eats "y"
            #       # Multi-resource
            #       omn[c,"MR"] <- 1
            #     }
            #   }
            # }
          }
        }
      }
    }
  }
  
  omn <- omn[rowSums(omn[,-1]) > 0,]
  return(omn)
}
