
# -------------------------------------------------------------------------
# Idea:
# Write one function that obtains centrality measures from posterior arrays
# The one function that plots these centrality measures
# currently supports instrength, outstrength, density
# Things to think about: does this also work with the Stan samples?

get_centrality <- function(fitobj,
                           burnin = 500){
  
  # TODO add check that the diagonal is zero
  # TODO allow that we don't take absolute values?
  # TODO ignore autoregressive effects???
  
  
  #--- BGGM
  # Obtain samples
  beta_samps <- abs(fitobj$fit$beta[,, -c(1:burnin)])
  pcor_samps <- abs(fitobj$fit$pcors[,, -c(1:burnin)])
  
  cnames <- colnames(fitobj$Y)
  
  #--- Centrality measures
  # In-strength
  instrength <- t(apply(beta_samps, MARGIN = 3, FUN = rowSums))
  colnames(instrength) <- cnames
  
  # Out-strength
  outstrength <- t(apply(beta_samps, MARGIN = 3, FUN = colSums))
  colnames(outstrength) <- cnames
  
  # Contemporaneous strength
  strength <- t(apply(pcor_samps, MARGIN = 3, FUN = rowSums))
  colnames(strength) <- cnames
  
  # Density
  density_beta <- apply(beta_samps, MARGIN = 3, FUN = sum)
  density_pcor <- apply(pcor_samps, MARGIN = 3, FUN = sum)
  
  #--- Return
  return(list(instrength = instrength,
              outstrength = outstrength,
              strength = strength,
              density_beta = density_beta,
              density_pcor = density_pcor))
  
}



