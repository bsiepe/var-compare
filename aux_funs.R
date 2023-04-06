# -------------------------------------------------------------------------
# Data Generation ---------------------------------------------------------
# -------------------------------------------------------------------------
# Temporary Server Imports Matrixcalc -------------------------------------
# package matrixcalc was not installed on the server
# so I copied source code
is.square.matrix <- function( x )
{
  ###
  ### determines if the given matrix is a square matrix
  ###
  ### arguments
  ### x = a matrix object
  ###
  if ( !is.matrix( x ) )
    stop( "argument x is not a matrix" )
  return( nrow(x) == ncol(x) )
}

is.symmetric.matrix <- function( x )
{
  ###
  ### this function determines if the matrix is symmetric
  ###
  ### argument
  ### x = a numeric matrix object
  ###
  if ( !is.matrix( x ) ) {
    stop( "argument x is not a matrix" )
  }
  if ( !is.numeric( x ) ) {
    stop( "argument x is not a numeric matrix" )
  }    
  if ( !is.square.matrix( x ) )
    stop( "argument x is not a square numeric matrix" )
  return( sum( x == t(x) ) == ( nrow(x) ^ 2 ) )
}

is.positive.semi.definite <- function( x, tol=1e-8 )
{
  ###
  ### this function determines if the given real symmetric matrix is positive semi definite
  ### parameters
  ### x = a square numeric matrix object
  ### tol = tolerance level for zero
  ###
  if ( !is.square.matrix( x ) )
    stop( "argument x is not a square matrix" )
  if ( !is.symmetric.matrix( x ) )
    stop( "argument x is not a symmetric matrix" )
  if ( !is.numeric( x ) )
    stop( "argument x is not a numeric matrix" )
  eigenvalues <- eigen(x, only.values = TRUE)$values
  n <- nrow( x )
  for ( i in 1: n ) {
    if ( abs( eigenvalues[i] ) < tol ) {
      eigenvalues[i] <- 0
    }
  }    
  if ( any( eigenvalues < 0 ) ) {
    return( FALSE )
  }
  return( TRUE )
}




# Change Graph ------------------------------------------------------------
# Function to change "true graph" according to specified changes of max value
# and noise
# Uses uniform noise for now


#' Change Graph
#' Changes graph structure (beta and kappa matrix) according to prescpecified changes. 
#' @param truegraph The true graph (should contain a beta and a kappa matrix)
#' @param changemax Vector of factors with which the largest matrix element
#' should be multiplied.
#' @param noise Vector of uniform noise that is added to the matrixes.
#'
#' @return List containing truegraph and all changed graphs as elements. 
#' @export
#'
change_graphs <- function(truegraph = NULL, 
                          changemax = NULL,
                          noise,
                          permute_index = NULL,
                          permute_active = FALSE, # if matrix should be permuted
                          seed = 2022){
  
  set.seed(seed)
  
  ### Create storage
  l_out <- list()
  l_out[[1]] <- truegraph
  l_out[[1]]$args <- "truegraph"
  names(l_out)[[1]] <- "truegraph"
  
  ### Obtain matrix characteristics
  ## Beta
  m_beta <- truegraph$beta
  b_i <- nrow(m_beta)
  b_j <- ncol(m_beta)
  
  ## Kappa
  m_kappa <- truegraph$kappa
  k_i <- nrow(m_kappa)
  k_j <- ncol(m_kappa)
  
  ### Find max
  ## Beta
  b_max <- which.max(m_beta)
  
  ## Kappa - only search for off-diagonal
  # Initialize the largest element to be zero
  largest_element <- 0
  max_row <- NA
  max_col <- NA
  # Loop through kappa and find the largest off-diagonal element
  for (i in 1:k_i) {
    for (j in 1:k_j) {
      # Skip the diagonal elements (i.e. the elements where i == j)
      if (i != j) {
        # Update the largest element if we find a larger one
        if (abs(m_kappa[i,j]) > abs(largest_element)) {
          largest_element <- m_kappa[i,j]
          max_row <- i
          max_col <- j
        }
      }
    }
  }
  
  
  
  
  ### Change max
  if(!is.null(changemax)){
    l_out_change <- list()
    ## Beta
    for(c in seq(changemax)){
      # create storage for each changed graph
      l_out_change[[c]] <- list()
      tmp_beta <- m_beta
      tmp_beta[b_max] <- m_beta[b_max]*changemax[c]
      l_out_change[[c]]$beta <- tmp_beta
    }
    
    ## Kappa
    # change the element in the matrix on both sides of diagonal
    for(c in seq(changemax)){
      # Check for positive-semidefiniteness
      psd_check <- FALSE
      check_count_max <- 0
      while(psd_check == FALSE & check_count_max <= 100){
        check_count_max <- check_count_max + 1
        tmp_kappa <- m_kappa
        tmp_kappa[max_row,max_col] <- m_kappa[max_row,max_col]*changemax[c]
        tmp_kappa[max_col,max_row] <- m_kappa[max_col,max_row]*changemax[c]

        
        psd <- matrixcalc::is.positive.semi.definite(tmp_kappa)
        
        if(psd){
          psd_check <- TRUE
        }
        else{
          print(paste0("Changemax leads to matrix that is not positive definite for ", changemax[c], " attempt ", check_count_max))
        }
        
      }
      l_out_change[[c]]$kappa <- tmp_kappa
      l_out_change[[c]]$args <- paste0("change", changemax[c])
      names(l_out_change)[[c]] <- paste0("change", changemax[c])
    }
  #   ## Checking
  #   # values outside unit circle
  #   # TODO adapt to not setting kappa diagonal to 1
  #   l_out_change <- lapply(l_out_change, function(x){
  #     lapply(x, function(mat) {
  #       mat[mat > 1] <- 1
  #       mat[mat < -1] <- -1
  #       mat
  #     })
  #   })
  }
  
  
  
  
  ### Add noise
  if(!is.null(noise)){
    l_out_noise <- list()
    ## Beta
    for(n in seq_along(noise)){
      l_out_noise[[n]] <- list()
      noise_beta <- matrix(runif(n = b_i*b_j, min = -noise[n], max = noise[n]), 
                           nrow = b_i, ncol = b_j)
      tmp_beta <- m_beta + noise_beta
      l_out_noise[[n]]$beta <- tmp_beta
      l_out_noise[[n]]$noisebeta <- noise_beta
    }
    
    ## Kappa
    noise_mat <- list()
    for(n in seq(noise)){
      # create change matrix for kappa, which is then scaled w.r.t. the sqrt of  diagonal elements
      noise_mat <- matrix(data = runif(n = k_i*k_j, min = -noise[n], max = noise[n]),
                          nrow = k_i, ncol = k_j)
      # scaling loop
      for(i in seq(k_i)){      # loop over rows
        for(j in seq(k_j)){    # loop over cols
          noise_mat[i,j] <- noise_mat[i,j]*sqrt(m_kappa[i,i]*m_kappa[j,j])
        }     
      }
      noise_mat <- as.matrix(Matrix::forceSymmetric(noise_mat))
      

      psd_check <- FALSE
      check_count_noise <- 0
      while(psd_check == FALSE & check_count_noise <= 100){
        check_count_noise <- check_count_noise + 1
        tmp_kappa <- m_kappa + noise_mat
        # tmp_kappa <- as.matrix(Matrix::forceSymmetric(tmp_kappa))
        psd <- matrixcalc::is.positive.semi.definite(tmp_kappa)  
        if(psd){
          psd_check <- TRUE
        } else{
          print(paste0("Adding noise to Kappa leads to matrix that is not positive semi-definite for ", noise[n], " attempt ", check_count_noise))
          # Redraw matrix
          noise_mat <- matrix(data = runif(n = k_i*k_j, min = -noise[n], max = noise[n]),
                              nrow = k_i, ncol = k_j)
          # scaling loop
          for(i in seq(k_i)){      # loop over rows
            for(j in seq(k_j)){    # loop over cols
              noise_mat[i,j] <- noise_mat[i,j]*sqrt(m_kappa[i,i]*m_kappa[j,j])
            }     
          }
          noise_mat <- as.matrix(Matrix::forceSymmetric(noise_mat))
          
        } # end else
      } # end while
      
      
      l_out_noise[[n]]$kappa <- tmp_kappa
      l_out_noise[[n]]$args <- paste0("noise", noise[n])
      l_out_noise[[n]]$noisekappa <- noise_mat
      names(l_out_noise)[[n]] <- paste0("noise", noise[n])
    }
  } # end noise
  
  
  ### Permute matrix
  if(isTRUE(permute_active)){
    l_out_perm <- list()
    l_out_perm[[1]] <- list()
    l_out_perm[[1]]$beta <- permute_mat_col(m_beta, permute_index)
    l_out_perm[[1]]$kappa <- permute_mat_col(m_kappa, permute_index)  
    l_out_perm[[1]]$args <- "perm"
    names(l_out_perm) <- "perm"
    l_out <- c(l_out, l_out_change, l_out_noise, l_out_perm)
    
  }

  
  
  ### Output
  # Combine truegraph, maxchange, and added noise
  l_out <- c(l_out, l_out_change, l_out_noise)
  return(l_out)
  
}



# Permute matrix columns --------------------------------------------------
# This permutes matrix columns while keeping diagonal elements on the diagonal
permute_mat_col <- function(mat,
                            symmetric = FALSE,
                            permute_index = NULL,
                            remove_names = TRUE){    # names removed for identical check
  # number of matrix columns
  p <- ncol(mat)
  
  if(remove_names){
    colnames(mat) <- NULL
    rownames(mat) <- NULL
  }
  
  # if no index is provided
  if(is.null(permute_index)){
    permute_index <- sample(p, p, replace = FALSE)
  }
  
  # Create a matrix
  diag_m <- diag(mat)
  
  # Permutation index for columns
  perm_ind <- permute_index
  
  # Permute the matrix
  perm_mat <- mat[, perm_ind]
  
  # Give new index
  new_ind <- 1:p
  
  # Calculate difference between old and new index
  diff_ind <- perm_ind - new_ind
  
  # Setup new matrix
  perm_mat_c <- perm_mat
  
  # Change position of diagonal elements
  for(i in 1:ncol(mat)){
    # value that should be on the diagonal
    diag_val <- diag_m[perm_ind[i]]
    
    # value that currently is on diagonal 
    swap_val <- diag(perm_mat)[i]
    
    
    # Set correct diagonal value
    perm_mat_c[i, i] <- diag_val
    perm_mat_c[i+diff_ind[i], i] <- swap_val
    
  }
  
  # Double check
  id_check <- identical(sort(diag(mat)), sort(diag(perm_mat_c)))
  if(!id_check){
    stop("Something went wrong. Diagonal not correct.")
  }
  
  # For symmetric matrix, copy upper diagonal
  if(isTRUE(symmetric)){
    perm_mat_c <- as.matrix(Matrix::forceSymmetric(perm_mat_c))
  }
  
  
  return(perm_mat_c)

}


# Simulate raw data -------------------------------------------------------
# 
#' Simulate raw data
#' Simulates raw data from fixed generating process to obtain $n$ new individuals. 
#' Data are standardized by default. 
#'
#' @param dgp Graph containing beta and kappa matrix. 
#' @param n Number of individuals to simulate
#' @param tp Timepoints per time series.
#' @param means Mean vector. Defaults to zero. 
#' @param standardize Should data be z-standardized? Defaults to true. 
#'
#' @return List with data and function arguments. 
#' @export
#'
#' 
sim_raw_parallel <- function(dgp,
                             n,
                             tp,
                             means = 0,
                             seed,
                             standardize = TRUE){
  
  # save function arguments in output
  args <- as.list(environment())
  args$dgp <- dgp$args
  

  # Reproducible loops - not used anymore, as we need different seeds
  # across conditions. 
  # registerDoRNG(seed)
  
  data <- foreach(i = seq(n), .packages = "graphicalVAR") %dopar% {
    raw_data <- list()
    raw_data$data <- as.data.frame(graphicalVAR::graphicalVARsim(nTime = tp,
                                                                 beta = dgp$beta,
                                                                 kappa = dgp$kappa,
                                                                 mean = means))
    if(standardize == TRUE){
      # Standardize data
      raw_data$data <- apply(raw_data$data, 2, scale)
    }
    # return name of data-generating process
    raw_data$args <- args
    raw_data 
    
  }
  return(data)

}





# Fit VAR parallel merged -------------------------------------------------
# This is a merge of fit_var_parallel and fit_var_parallel_post into one function

fit_var_parallel_merged <- function(data, 
                                     n,         # number of individuals
                                     nds,       # number of datasets
                                     rho_prior, 
                                     beta_prior,
                                     seed,
                                     iterations,
                                     get_kappa = FALSE,
                                     posteriorsamples = FALSE,
                                     multigroup = FALSE,
                                     pruneresults = FALSE, 
                                     summarize_post = TRUE, # summarize posterior samples into intervals 
                                     cred_int = 0.95,       # credible interval for posterior summary
                                     select = FALSE,       # apply selection based on CI, currently only supported for non single dataset fitting
                                     save_files = FALSE,
                                     dgp_name = NULL){     # option to return as .rds
  
  if(n != length(data)){
    warning("The n provided does not match the number of available data frames")
  }
  require(doParallel)

  # reproducible parallelization
  doRNG::registerDoRNG(seed)
  
  
  # Input checks
  if(isTRUE(select) & isTRUE(pruneresults)){
    stop("Using select and pruneresults at the same time is not supported.")
  }
  

  if(isFALSE(posteriorsamples) & isFALSE(multigroup)){
    print("Fitting to raw data")
    fit <- foreach(i = seq(nds), .packages = "BGGM", .export = "sim_select") %dopar% {

      fit_ind <- list()
      if(is.list(data[[i]]$data) | is.numeric(data[[i]]$data)){
        fit_ind <- tryCatch({BGGM::var_estimate(as.data.frame(data[[i]]$data),
                                                rho_sd = rho_prior,
                                                beta_sd = beta_prior,
                                                iter = iterations,
                                                progress = FALSE,
                                                seed = seed)}, error = function(e) NULL)
        
        
        # Delete irrelevant matrices
        if(is.list(fit_ind)){
          fit_ind$fit$fisher_z <- NULL
        }
        
        
        if(isTRUE(get_kappa)){
          # check if fitting worked
          if(is.list(fit_ind)){
            # Invert covariance matrix of residuals to obtain precision matrix
            fit_ind$fit$kappa <- array(apply(fit_ind$fit$Sigma, 3, solve), 
                                       dim = dim(fit_ind$fit$Sigma))
            
            # Calculate mean of kappa
            fit_ind$kappa_mu <- apply(fit_ind$fit$kappa, c(1,2), mean)
            
          }

        }
        
        # prune results for comparison purposes
        if(isTRUE(pruneresults) & is.list(fit_ind)){
          # beta <- fit_ind$fit$beta
          # kappa <- fit_ind$fit$kappa
          beta_mu <- fit_ind$beta_mu 
          pcor_mu <- fit_ind$pcor_mu
          # kappa_mu <- fit_ind$kappa_mu
          args <- data[[i]]$args
          fit_ind <- list()
          fit_ind$beta_mu <- beta_mu
          fit_ind$pcor_mu <- pcor_mu
          fit_ind$kappa_mu <- NA
          fit_ind$args <- args
          # fit_ind$fit$beta <- beta
          # fit_ind$fit$kappa <- kappa

        } # end if isTRUE 
        
        # Select option
        if(isTRUE(select) & is.list(fit_ind)){
          sel <- sim_select(fit_ind)
          if(isTRUE(summarize_post)){
            cred_interval <- summarize_post(fit_ind, cred = cred_int)
          }
          fit_ind <- list()
          fit_ind <- stats::setNames(sel, names(sel))
          fit_ind$args <- data[[i]]$args
          if(isTRUE(summarize_post)){
            fit_ind$cred_interval <- cred_interval
          }
          
          
        } # end isTRUE select
        
        if(isFALSE(select) & isTRUE(summarize_post) & is.list(fit_ind)){
          beta_mu <- fit_ind$beta_mu 
          pcor_mu <- fit_ind$pcor_mu
          # kappa_mu <- fit_ind$kappa_mu
          args <- data[[i]]$args
          cred_interval <- summarize_post(fit_ind, cred = cred_int)
          fit_ind <- list()
          fit_ind$beta_mu <- beta_mu
          fit_ind$pcor_mu <- pcor_mu
          fit_ind$kappa_mu <- NA
          fit_ind$cred_interval <- cred_interval
          fit_ind$args <- args
          
        }
        

        
      } # end if is.list
      # else fit_ind <- NA
      
      return(fit_ind)  
    } # end foreach
    # Cut away nonconverged attempts
    # only for non-select
    # if(isFALSE(select)){
    #   lapply(fit, function(x){
    #     if(length(x$fit) == 0){
    #       warning("Some models did not converge!")}
    #   })
    #   fit <- fit[!sapply(fit, function(x) length(x$fit) == 0)]
    #   
    #   
    #   # Return list with desired length
    #   fit <- fit[c(1:n)]
    # }


    
    if(isFALSE(save_files)){
      return(fit)
      
    }
    # If outputs should be saved as RDS
    if(isTRUE(save_files)){
      saveRDS(fit, file = here::here(paste0("data/compare_sim_fits/fit_",dgp_name, "_", data[[1]]$args$dgp, "_", data[[1]]$args$tp,".RDS")))
      
    }
    
    
  } # end isFALSE(posteriorsamples)
  
  
  # if we fit the data on samples generated from the posterior
  if(isTRUE(posteriorsamples)){
    print("Fitting to posterior samples")
    
    # counter for converged models 
    m <- list()
    # counter for attempted models
    c <- list()
    
    
    # loop across individuals
    foreach(i = seq(n), .packages = "BGGM") %dopar% {
      # Counter for converged models
      m[[i]] <- 0
      # Counter for estimated models
      c[[i]] <- 0
      fit <- list()
      # loop across datasets
      for(d in seq(nds)) {
        if(m[[i]] >= 100){
          break 
        }
        
        fit_ind <- list()
        if(is.list(data[[i]][[d]])){
          fit_ind <- tryCatch({BGGM::var_estimate(data[[i]][[d]],
                                                  rho_sd = rho_prior,
                                                  beta_sd = beta_prior,
                                                  iter = iterations,
                                                  progress = FALSE,
                                                  seed = seed)}, error = function(e) NULL)
          # Add to counter for estimated models
          c[[i]] <- c[[i]]+1
          
          # check if fitting worked
          if(is.list(fit_ind)){
            m[[i]] <- m[[i]]+1
            # print(m[[i]])
            
            
            if(isTRUE(get_kappa)){
              # Invert covariance matrix of residuals to obtain precision matrix
              fit_ind$fit$kappa <- array(apply(fit_ind$fit$Sigma, 3, solve), 
                                         dim = dim(fit_ind$fit$Sigma))
              # Calculate mean of kappa
              fit_ind$kappa_mu <- apply(fit_ind$fit$kappa, c(1,2), mean)
              
            }
            
            
            # prune results for comparison purposes
            # When fitting to posterior data, we no longer need 
            # raw kappa/beta matrizes
            if(isTRUE(pruneresults)){
              # kappa <- fit_ind$fit$kappa
              # beta <- fit_ind$fit$beta
              beta_mu <- fit_ind$beta_mu 
              kappa_mu <- fit_ind$kappa_mu
              pcor_mu <- fit_ind$pcor_mu
              n_attempts <- c[[i]]
              fit_ind <- list()
              # fit_ind$fit <- list()
              fit_ind$beta_mu <- beta_mu
              fit_ind$kappa_mu <- kappa_mu
              fit_ind$pcor_mu <- pcor_mu
              fit_ind$n_attempts <- n_attempts
              # fit_ind$fit$beta <- beta
              # fit_ind$fit$kappa <- kappa
              
            } 
            
            fit[[d]] <- fit_ind  
          }    # end is.list(fit_ind)
          
          
        } # end is.list(data[[i]])  
        
        
      } # end for loop 
      
      return(fit)
      
    }
    
    
  } # end isTRUE posteriorsamples
  if(isTRUE(multigroup)){
    print("Fitting to multigroup posterior samples")
    
    # counter for converged models 
    m <- list()
    # counter for attempted models
    c <- list()
    
    fit <- list()
    
    # Progressbar
    pb <- utils::txtProgressBar(0, nds, style = 3)
    
    # loop across samples
    fit <- foreach(i = seq(nds), .packages = "BGGM") %dopar% {
      utils:: setTxtProgressBar(pb, i)
      
      # TODO reimplement!
      # # Counter for converged models
      # m[[i]] <- 0
      # # Counter for estimated models
      # c[[i]] <- 0
      fit_ind <- list()
      # loop across groups
      for(d in seq(n)) {
        # if(m[[i]] >= 100){
        #   break 
        # }
        
        if(is.list(data[[i]])){
          if(is.list(data[[i]][[d]])){
            fit_ind[[d]] <- tryCatch({BGGM::var_estimate(as.data.frame(data[[i]][[d]]$data),
                                                         rho_sd = rho_prior,
                                                         beta_sd = beta_prior,
                                                         iter = iterations,
                                                         progress = FALSE,
                                                         seed = seed)}, error = function(e) NA)
            # # Add to counter for estimated models
            # c[[i]] <- c[[i]]+1
            
            # check if fitting worked (length condition needed when first d iteration
            # fails and list is empty)
            if(length(fit_ind) > 0){
              if(is.list(fit_ind[[d]])){
                # m[[i]] <- m[[i]]+1
                # print(m[[i]])
                
                
                if(isTRUE(get_kappa)){
                  # Invert covariance matrix of residuals to obtain precision matrix
                  fit_ind[[d]]$fit$kappa <- array(apply(fit_ind[[d]]$fit$Sigma, 3, solve), 
                                                  dim = dim(fit_ind[[d]]$fit$Sigma))
                  # Calculate mean of kappa
                  fit_ind[[d]]$kappa_mu <- apply(fit_ind[[d]]$fit$kappa, c(1,2), mean)
                  
                }
                
                
                # prune results for comparison purposes
                # When fitting to posterior data, we no longer need 
                # raw kappa/beta matrizes
                if(isTRUE(pruneresults)){
                  # kappa <- fit_ind[[d]]$fit$kappa
                  # beta <- fit_ind[[d]]$fit$beta
                  beta_mu <- fit_ind[[d]]$beta_mu 
                  kappa_mu <- fit_ind[[d]]$kappa_mu
                  pcor_mu <- fit_ind[[d]]$pcor_mu
                  # n_attempts <- c[[i]]
                  fit_ind[[d]] <- list()
                  # fit_ind[[d]]$fit <- list()
                  fit_ind[[d]]$beta_mu <- beta_mu
                  fit_ind[[d]]$kappa_mu <- kappa_mu
                  fit_ind[[d]]$pcor_mu <- pcor_mu
                  # fit_ind[[d]]$n_attempts <- n_attempts
                  # fit_ind[[d]]$fit$beta <- beta
                  # fit_ind[[d]]$fit$kappa <- kappa
                  
                } # end isTRUE(pruneresults)
                
                
              }    # end is.list(fit_ind)
            }
            
            
            
          } # end is.list(data[[i]])  
        }
        
        
        
      } # end loop across groups 
      
      return(fit_ind)
      
    } # end foreach 
    return(fit)
    
  } # end isTRUE multigroup
  
  

  
}  













# Sim from posterior ------------------------------------------------------
# TODO save function arguments

#' Simulate from Posterior Samples
#' This function simulates a specified number of datasets from the posterior
#' of a given model. Uses the graphicalVARsim function from the graphicalVAR
#' package to generate data. 
#' @param fitobj BGGM fit object containing all posterior samples
#' @param n_datasets Number of datasets to create
#' @param n Number of individuals
#' @param tp Number of timepoints
#' @param iterations Number of iterations used in BGGM sampling
#' @param means Mean vector
#' @param convert_bggm DEPRECATED: Should results be converted to BGGM List format?
#'
#' @return List of datasets in dataframe format. 
#' @export

sim_from_post_parallel <- function(fitobj, 
                                   n_datasets, 
                                   n,
                                   tp,
                                   iterations,
                                   seed,  
                                   means = 0,
                                   convert_bggm = FALSE){
  # Extract parameters from fitobject
  # delete first 50 samples
  l_params <- list()
  for(i in seq(n)){
    l_params[[i]] <- list()
    l_params[[i]]$beta <- fitobj[[i]]$fit$beta[,,51:iterations]
    l_params[[i]]$kappa <- fitobj[[i]]$fit$kappa[,,51:iterations]
  }
  
  # reproducible parallelization
  registerDoRNG(seed)
  
  # Loop to create new datasets from posterior samples
  post_data <- foreach(i = seq(n), .packages = "graphicalVAR") %dopar% {
    dat <- list()
    # Loop over number of datasets to create from posterior
    for(j in seq(n_datasets)){
      # get random posterior sample
      # Needs transposing of beta matrix!
      smp <- sample(iterations, size = 1)
      dat[[j]] <- try(as.data.frame(graphicalVAR::graphicalVARsim(nTime = tp,
                                                                  beta = t(l_params[[i]]$beta[,,smp]),
                                                                  kappa = l_params[[i]]$kappa[,,smp],
                                                                  mean = means))) 
      
    }
    
    dat
    
    
    
  } # end parallel
  
  # # No longer in use: Should data be converted to bggm format?
  # if(isTRUE(convert_bggm)){
  #   post_data <- lapply(post_data, format_bggm_list)
  #   
  # }
  post_data
  
}



# Summarize posterior -----------------------------------------------------
# summarize_post <- function(res,
#                            cred = 0.95){
#   
#   # Lower and upper bound
#   lb <- (1-cred)/2
#   ub <- 1-lb
#   
#   # for beta
#   beta_lb <- apply(res$fit$beta, c(1,2), stats::quantile, lb)
#   beta_ub <- apply(res$fit$beta, c(1,2), stats::quantile, ub)
#   
#   # for pcor
#   pcor_lb <- apply(res$fit$pcors, c(1,2), stats::quantile, lb)
#   pcor_ub <- apply(res$fit$pcors, c(1,2), stats::quantile, ub)  
#   
#   # Output
#   out <- list(beta_lb = beta_lb,
#               beta_ub = beta_ub, 
#               pcor_lb = pcor_lb, 
#               pcor_ub = pcor_ub)
# 
#   return(out)
# }


summarize_post <- function(res, cred = c(0.95)) {
  
  # Check if the "cred" argument is a numeric vector
  if (!is.numeric(cred) | any(cred < 0 | cred > 1)) {
    stop("The 'cred' argument must be a numeric vector between 0 and 1.")
  }
  
  # Convert "cred" to a sorted vector of unique values
  cred <- sort(unique(cred))
  
  # Initialize output lists
  beta_lb_list <- vector("list", length(cred))
  beta_ub_list <- vector("list", length(cred))
  pcor_lb_list <- vector("list", length(cred))
  pcor_ub_list <- vector("list", length(cred))
  
  # Compute bounds for each "cred" value
  for (i in seq_along(cred)) {
    
    # Lower and upper bound
    lb <- (1 - cred[i])/2
    ub <- 1 - lb
    
    # for beta
    beta_lb_list[[i]] <- apply(res$fit$beta, c(1, 2), stats::quantile, lb)
    beta_ub_list[[i]] <- apply(res$fit$beta, c(1, 2), stats::quantile, ub)
    
    # for pcor
    pcor_lb_list[[i]] <- apply(res$fit$pcors, c(1, 2), stats::quantile, lb)
    pcor_ub_list[[i]] <- apply(res$fit$pcors, c(1, 2), stats::quantile, ub)  
    
  }
  names(beta_lb_list) <- paste0("lb_", cred)
  names(beta_ub_list) <- paste0("ub_", cred)
  names(pcor_lb_list) <- paste0("lb_", cred)
  names(pcor_ub_list) <- paste0("ub_", cred)
  
  # Combine output into a list
  out <- list(beta_lb = beta_lb_list,
              beta_ub = beta_ub_list, 
              pcor_lb = pcor_lb_list, 
              pcor_ub = pcor_ub_list)
  
  return(out)
}





# -------------------------------------------------------------------------
# Comparison Functions ----------------------------------------------------
# -------------------------------------------------------------------------
# Distance between empirical and posterior --------------------------------
# THIS IS A NEW BETA VERSION! Not yet in use
postemp_distance_new <- function(post,
                             emp,
                             comp,
                             mod){
  
  # storage
  dist_out <- list()
  
  
  # define the distance function based on comp
  distance_fn_beta <- switch(comp,
                             frob =   {function(x) norm(emp[[mod]]$beta_mu-x$beta_mu, type = "F")},
                             maxdiff = {function(x) max(abs((emp[[mod]]$beta_mu-x$beta_mu)))},
                             l1 = {function(x) sum(abs((emp[[mod]]$beta_mu-x$beta_mu)))}
  )
  distance_fn_pcor <- switch(comp,
                             frob = {function(x) norm(emp[[mod]]$pcor_mu-x$pcor_mu, type = "F")},
                             maxdiff = {function(x) max(abs((emp[[mod]]$pcor_mu-x$pcor_mu)))},
                             l1 = {function(x) sum(abs((emp[[mod]]$pcor_mu-x$pcor_mu)))}
  )
  
  
## Check if estimation worked
# Should be unneccessary if non-converged attempts were deleted
if(!is.list(post[[mod]]$fit) | !is.list(post[[mod]]$fit)){
    beta_distance <- NA
    pcor_distance <- NA
    stop("Input not a list, probably estimation did not converge.")
    
} 

  # if both elements are lists
  else{
    dist_out[["beta"]] <- unlist(lapply(post[[mod]]$fit, distance_fn_beta(x)))
    dist_out[["pcor"]] <- unlist(lapply(post[[mod]]$fit, distance_fn_pcor(x)))
    
  }   

  
  
  return(dist_out)
}


# TODO output structure
postemp_distance <- function(post,
                             emp,
                             comp,
                             mod){

  # storage
  dist_out <- list()



  if(comp == "frob"){
    normtype = "F"
    frob_beta <- unlist(lapply(post[[mod]]$fit, function(x){
      if(length(x) == 0 | !is.list(x))
        {NA}
      else
        {norm(emp[[mod]]$beta_mu-x$beta_mu, type = normtype)}}))

    frob_pcor <- unlist(lapply(post[[mod]]$fit, function(x){
      if(length(x) == 0 | !is.list(x))
        {NA}
      else
        {norm(emp[[mod]]$pcor_mu-x$pcor_mu, type = normtype)}}))

    dist_out[["beta"]] <- frob_beta
    dist_out[["pcor"]] <- frob_pcor

    }  # end frob




  if(comp == "maxdiff"){
    maxdiff_beta <- unlist(lapply(post[[mod]]$fit, function(x){
      if(length(x) == 0 | !is.list(x))
        {NA}
      else
        {max(abs(emp[[mod]]$beta_mu-x$beta_mu))}}))

    maxdiff_pcor <- unlist(lapply(post[[mod]]$fit, function(x){
      if(length(x) == 0 | !is.list(x))
        {NA}
      else
        {max(abs(emp[[mod]]$pcor_mu-x$pcor_mu))}}))

    dist_out[["beta"]] <- maxdiff_beta
    dist_out[["pcor"]] <- maxdiff_pcor

    }  # end maxdiff



  if(comp == "l1"){
    l1_beta <-  unlist(lapply(post[[mod]]$fit, function(x){
      if(length(x) == 0 | !is.list(x))
        {NA}
      else
        {sum(abs(emp[[mod]]$beta_mu-x$beta_mu))}}))

    l1_pcor <-  unlist(lapply(post[[mod]]$fit, function(x){
      if(length(x) == 0 | !is.list(x))
      {NA}
      else
      {sum(abs(emp[[mod]]$pcor_mu-x$pcor_mu))}}))

      dist_out[["beta"]] <- l1_beta
      dist_out[["pcor"]] <- l1_pcor


    }  # end l1

  if(comp == "pcmd"){
    # only suitable for (partial) correlations
    pcmd_pcor <-  unlist(lapply(post[[mod]]$fit, function(x){
      if(length(x) == 0 | !is.list(x))
      {NA}
      else
      {1 - (sum(diag(emp[[mod]]$pcor_mu %*% x$pcor_mu)) / (norm(emp[[mod]]$pcor_mu, type = "F") * norm(x$pcor_mu, type = "F")))}}))
    # dist_out[["beta"]] <- NULL
    dist_out[["pcor"]] <- pcmd_pcor


  }


  return(dist_out)
}


# Distance between posterior samples --------------------------------------
postpost_distance <- function(post_a,
                             post_b,
                             comp,
                             mod_a,
                             mod_b){
  
  # storage
  dist_out <- list()
  
  
  
  if(comp == "frob"){
    normtype = "F"
    frob_beta <- tryCatch(mapply(function(x,y)
    {
      if(!is.list(x) | !is.list(y))
      {NA} 
      else 
        norm(x$beta_mu-y$beta_mu, type = normtype)}, post_a[[mod_a]], post_b[[mod_b]]), 
    error = function(e) {NA})
    
    frob_pcor <- tryCatch(mapply(function(x,y)
    {
      if(!is.list(x) | !is.list(y))
      {NA} 
      else 
        norm(x$pcor_mu-y$pcor_mu, type = normtype)}, post_a[[mod_a]], post_b[[mod_b]]), 
    error = function(e) {NA})
    
    dist_out[["beta"]] <- frob_beta
    dist_out[["pcor"]] <- frob_pcor
    
  }  # end frob
  
  
  
  
  if(comp == "maxdiff"){
    maxdiff_beta <- tryCatch(mapply(function(x,y)
    {
      if(!is.list(x) | !is.list(y))
      {NA} 
      else 
        max(abs((x$beta_mu-y$beta_mu)))}, post_a[[mod_a]], post_b[[mod_b]]), 
    error = function(e) {NA})
    
    maxdiff_pcor <- tryCatch(mapply(function(x,y)
    {
      if(!is.list(x) | !is.list(y))
      {NA} 
      else 
        max(abs((x$pcor_mu-y$pcor_mu)))}, post_a[[mod_a]], post_b[[mod_b]]), 
    error = function(e) {NA})
    
    dist_out[["beta"]] <- maxdiff_beta
    dist_out[["pcor"]] <- maxdiff_pcor
    
  }  # end maxdiff
  
  
  
  if(comp == "l1"){
    l1_beta <-  tryCatch(mapply(function(x,y)
    {
      if(!is.list(x) | !is.list(y))
      {NA} 
      else 
        sum(abs((x$beta_mu-y$beta_mu)))}, post_a[[mod_a]], post_b[[mod_b]]), 
    error = function(e) {NA})
    
    l1_pcor <-  tryCatch(mapply(function(x,y)
    {
      if(!is.list(x) | !is.list(y))
      {NA} 
      else 
        sum(abs((x$beta_mu-y$beta_mu)))}, post_a[[mod_a]], post_b[[mod_b]]), 
    error = function(e) {NA})
    
    dist_out[["beta"]] <- l1_beta
    dist_out[["pcor"]] <- l1_pcor  
    
    
  }  # end l1
  
  return(dist_out)
}



# Distance within posterior predictive samples ---------------------------------
#' Calculates distances between pairs of fitted models using the posterior samples or posterior predictive draws
#'
#' This function computes distances between a specified number of pairs of fitted models, which can be obtained either from posterior samples or posterior predictive draws. The distance between two models can be calculated based on three options: Frobenius norm, maximum difference, or L1 norm. The function allows for comparison of posterior samples or posterior predictive draws, and beta coefficients or partial correlations can be used as inputs.
#'
#' @param post An object of class \code{posterior}, which contains either posterior samples or posterior predictive draws.
#' @param comp A character string indicating the type of distance between models that should be calculated. The options include: "frob" (Frobenius norm), "maxdiff" (maximum difference), or "l1" (L1 norm).
#' @param pred A logical indicating whether the input is posterior predictive draws (TRUE) or posterior samples (FALSE).
#' @param draws An integer specifying the number of random pairs of models that should be compared.
#'
#' @return A list of distances between the specified pairs of fitted models. The list has length equal to the specified number of random pairs. Each list element contains two distance values, one for beta coefficients and one for partial correlations.
#'
#'
#' @export post_distance_within
post_distance_within <- function(post, 
                                 comp,
                                 pred,         # posterior predictive?
                                 draws = 1000){
  
  # storage
  dist_out <- list()
  
  
  # for posterior predictive approach
  if(isTRUE(pred)){
    # define the distance function based on comp
    distance_fn_beta <- switch(comp,
                               frob =   {function(x, y, mod_one, mod_two) norm(x$fit[[mod_one]]$beta_mu-y$fit[[mod_two]]$beta_mu, type = "F")},
                               maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit[[mod_one]]$beta_mu-y$fit[[mod_two]]$beta_mu)))},
                               l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit[[mod_one]]$beta_mu-y$fit[[mod_two]]$beta_mu)))}
    )
    distance_fn_pcor <- switch(comp,
                               frob = {function(x, y, mod_one, mod_two) norm(x$fit[[mod_one]]$pcor_mu-y$fit[[mod_two]]$pcor_mu, type = "F")},
                               maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit[[mod_one]]$pcor_mu-y$fit[[mod_two]]$pcor_mu)))},
                               l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit[[mod_one]]$pcor_mu-y$fit[[mod_two]]$pcor_mu)))}
    )
    
    
    
    
    # Obtain number of models
    n_mod <- length(post$fit)
    
  }
  
  
  # for posteriors of empirical models
  if(isFALSE(pred)){
    # define the distance function based on comp
    # draw from all posterior samples
    

    distance_fn_beta <- switch(comp,
                               frob =   {function(x, y, mod_one, mod_two) norm(x$fit$beta[,,mod_one]-y$fit$beta[,,mod_two], type = "F")},
                               maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit$beta[,,mod_one]-y$fit$beta[,,mod_two])))},
                               l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit$beta[,,mod_one]-y$fit$beta[,,mod_two])))}
    )
    distance_fn_pcor <- switch(comp,
                               frob = {function(x, y, mod_one, mod_two) norm(x$fit$pcors[,,mod_one]-y$fit$pcors[,,mod_two], type = "F")},
                               maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit$pcors[,,mod_one]-y$fit$pcors[,,mod_two])))},
                               l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit$pcors[,,mod_one]-y$fit$pcors[,,mod_two])))}
    )
    
    # Obtain number of posterior samples
    n_mod <- dim(post$fit$beta)[3]
    
  }
  
  
  
  ## Draw two random models
  # delete burn-in iterations (to 50)
  # n_mod <- n_mod[-c(1:50)]
  
  # TODO: actually, "samples" would be more fitting here than "models"
  # "model" is still a residue from posterior predictive approach
  # Draw models spaced apart so that we don't have autocorrelation from sampling
  mod_pairs <- array(NA, dim = c(2, draws))
  # draw from first half of samples
  mod_pairs[1,1:(draws)] <- seq(51, draws+50, by = 1)
  
  # draw from second half of samples
  mod_pairs[2,1:(draws)] <- seq((n_mod/2)+51, (n_mod/2)+50+(draws), by = 1)
  
  # mod_pairs <- replicate(draws, sample(1:n_mod, size = 2, replace = TRUE))
  
  for(i in seq(draws)){
    # storage
    dist_out[[i]] <- list()
    mod_one <- mod_pairs[1,i]
    mod_two <- mod_pairs[2,i]
    
    # if mod_one and mod_two are equal, redraw
    if(mod_one == mod_two){
      mod_two <- sample(1:n_mod, size = 1)
    }
    
    ## Check if estimation worked
    # Should be unneccessary if non-converged attempts were deleted
    if(isTRUE(pred)){
      if(!is.list(post$fit[[mod_one]]) | !is.list(post$fit[[mod_two]])){
        beta_distance <- NA
        pcor_distance <- NA
        stop("Not a list.")
        
        
      } 
      # if both elements are lists
      else{
        beta_distance <- distance_fn_beta(post, post, mod_one, mod_two)
        pcor_distance <- distance_fn_pcor(post, post, mod_one, mod_two)
        
      }  
    }
    
    if(isFALSE(pred)){
      if(!is.list(post) | !is.list(post)){
        beta_distance <- NA
        pcor_distance <- NA
        
        stop("Not a list.")
        
      } 
      # if both elements are lists
      else{
        beta_distance <- distance_fn_beta(post, post, mod_one, mod_two)
        pcor_distance <- distance_fn_pcor(post, post, mod_one, mod_two)
        
      }  
    }
    
    
    
    
    dist_out[[i]]$comp <- comp
    dist_out[[i]]$mod_one <- mod_one
    dist_out[[i]]$mod_two <- mod_two
    dist_out[[i]]$beta <- beta_distance
    dist_out[[i]]$pcor <- pcor_distance  
    
    
  } # end for loop  
  out <- do.call(rbind, dist_out)
  out <- as.data.frame(out)
  
  
  return(out)
}



# # Cross-compare all posterior samples -------------------------------------

cross_compare <- function(
    postpost = FALSE,        # compute differences between posteriors
    fitpost_a = l_postres,
    fitpost_b = l_postres,
    fitemp_a = l_res,
    fitemp_b = l_res,
    mod_a = 1, 
    mod_b = 2,
    n_datasets = 100,
    comparison = "frob",
    ...){
  if(!is.numeric(mod_a) | !is.numeric(mod_b)){
    stop("Error: Model needs to have numerical index")
  }
  
  if(comparison == "frob"){
    normtype = "F"
  }
    # Distance within posterior samples
    null_a <- postemp_distance(post = fitpost_a, emp = fitemp_a, 
                               comp = comparison, mod = mod_a)
    null_b <- postemp_distance(post = fitpost_b, emp = fitemp_b, 
                               comp = comparison, mod = mod_b)
    
    # Compute empirical distance as test statistic
    if(isFALSE(postpost)){
      
      if(comparison == "frob"){
        # Compute Distance of empirical betas between a and b
        emp_beta <- tryCatch(norm(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu, type = normtype), error = function(e) {NA})
        
        # Compute Distance of empirical pcors between a and b
        emp_pcor <- tryCatch(norm(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu, type = normtype), error = function(e) {NA})
        
      }

      if(comparison == "maxdiff"){
        # Compute maxdiff of empirical betas between a and b
        emp_beta <- tryCatch(max(abs(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu)), error = function(e) {NA})
        
        # Compute maxdiff of empirical pcors between a and b
        emp_pcor <- tryCatch(max(abs(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu)), error = function(e) {NA})
        
      }
      
      if(comparison == "l1"){
        # Compute l1 of empirical betas between a and b
        emp_beta <- tryCatch(sum(abs(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu)), error = function(e) {NA})
        
        # Compute l1 of empirical pcors between a and b
        emp_pcor <- tryCatch(sum(abs(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu)), error = function(e) {NA})
        
        
      }
      
      # Save results
      cc_res_beta <- data.frame(model_ind = c(rep(mod_a, n_datasets), rep(mod_b, n_datasets)),
                                null = c(unlist(null_a[["beta"]]), unlist(null_b[["beta"]])),
                                emp = rep(emp_beta, n_datasets*2),
                                comp = rep(comparison, n_datasets*2),
                                type = rep("postemp", n_datasets*2))
      
      
      cc_res_pcor <- data.frame(model_ind = c(rep(mod_a, n_datasets), rep(mod_b, n_datasets)),
                                null = c(unlist(null_a[["pcor"]]), unlist(null_b[["pcor"]])),
                                emp = rep(emp_pcor, n_datasets*2),
                                comp = rep(comparison, n_datasets*2),
                                type = rep("postemp", n_datasets*2))
      
      
      
      
    } # end isFALSE(postpost)
    
    
    # Compute posterior distances as test distribution
    if(isTRUE(postpost)){
      
      # Distance posterior estimates between A and B
      post <- postpost_distance(post_a = fitpost_a, post_b = fitpost_b, 
                                comp = comparison,mod_a = mod_a, mod_b = mod_b)
      
      # Save results
      cc_res_beta <- data.frame(model_ind = c(rep(mod_a, n_datasets), rep(mod_b, n_datasets)),
                                null = c(null_a[["beta"]], null_b[["beta"]]),
                                emp = rep(post[["beta"]], 2),    
                                comp = rep(comparison, n_datasets*2),
                                type = rep("postpost", n_datasets*2))
      
      
      cc_res_pcor <- data.frame(model_ind = c(rep(mod_a, n_datasets), rep(mod_b, n_datasets)),
                                null = c(null_a[["pcor"]], null_b[["pcor"]]),
                                emp = rep(post[["pcor"]], 2),
                                comp = rep(comparison, n_datasets*2),
                                type = rep("postpost", n_datasets*2))
      
      
      
      
      
    } # end isTRUE(postpost)


  
  l_cc_res <- list()
  l_cc_res[["beta"]] <- cc_res_beta
  l_cc_res[["pcor"]] <- cc_res_pcor
  
  return(l_cc_res)
} # end function



# Within posterior comparison ---------------------------------------------
# This function uses bootstrapping to generate a difference distribution
# within the posterior

# Kappa not used for now
within_compare <- function(
    fitpost_a = l_postres,
    fitpost_b = l_postres,
    fitemp_a = l_res,
    fitemp_b = l_res,
    mod_a = 1, 
    mod_b = 2,
    n_datasets = 100,
    comparison = "frob",
    n_draws = 1000,
    postpred = FALSE,          # do we use the posterior predictive approach?
    # parallel = FALSE,
    # ncores = NULL, 
    ...){
  if(!is.numeric(mod_a) | !is.numeric(mod_b)){
    stop("Error: Model needs to have numerical index")
  }
  
  # If one of the models did not converge
  if(isFALSE(postpred) && length(fitpost_a[[mod_a]]$fit) == 0 |
     isFALSE(postpred) && length(fitpost_a[[mod_b]]$fit) == 0 ){
    cc_res <- data.frame(mat = c("beta", "pcor", "beta", "pcor"),
                         null = c(NA, NA, NA, NA),
                         model_ind = c(mod_a, mod_a, mod_b, mod_b),
                         mod = c("mod_a", "mod_a", "mod_b", "mod_b"),
                         emp = c(NA, NA, NA, NA),
                         comp = rep(comparison, 4))
    return(cc_res)
  }
  
  
  if(comparison == "frob"){
    normtype = "F"
  }
    # # Parallel or not?
    # if(isFALSE(parallel)){
      # Distance empirical to posterior sample estimates
      null_a <- post_distance_within(post = fitpost_a[[mod_a]], 
                                     comp = comparison, 
                                     draws = n_draws,
                                     pred = postpred)
      null_b <- post_distance_within(post = fitpost_b[[mod_b]], 
                                     comp = comparison, 
                                     draws = n_draws,
                                     pred = postpred)
      
    # }
    # 
    # if(isTRUE(parallel)){
    #   if(is.null(ncores)){
    #     stop("Error: Please specify >1 Core for parallelization")
    #   }
    #   null_a <- post_distance_within_par(post = fitpost_a[[mod_a]], 
    #                                      comp = comparison, 
    #                                      draws = n_draws,
    #                                      pred = postpred,
    #                                      ncores = ncores)
    #   null_b <- post_distance_within_par(post = fitpost_b[[mod_b]], 
    #                                      comp = comparison, 
    #                                      draws = n_draws,
    #                                      pred = postpred,
    #                                      ncores = ncores)
    # }

      
    


  
  # Compute empirical distance as test statistic
    if(comparison == "frob"){
      # Compute Distance of empirical betas between a and b
      emp_beta <- tryCatch(norm(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu, type = normtype), error = function(e) {NA})
      
      # Compute Distance of empirical pcors between a and b
      emp_pcor <- tryCatch(norm(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu, type = normtype), error = function(e) {NA})
      
      # Compute Distance of empirical kappas between a and b
      # emp_kappa <- tryCatch(norm(fitemp_a[[mod_a]]$kappa_mu - fitemp_b[[mod_b]]$kappa_mu, type = normtype), error = function(e) {NA})
      
    }
    
    if(comparison == "maxdiff"){
      # Compute maxdiff of empirical betas between a and b
      emp_beta <- tryCatch(max(abs(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu)), error = function(e) {NA})
      
      # Compute maxdiff of empirical pcors between a and b
      emp_pcor <- tryCatch(max(abs(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu)), error = function(e) {NA})
      
      # Compute maxdiff of empirical kappas between a and b
      # emp_kappa <- tryCatch(max(abs(fitemp_a[[mod_a]]$kappa_mu - fitemp_b[[mod_b]]$kappa_mu)), error = function(e) {NA})
    }
    
    if(comparison == "l1"){
      # Compute l1 of empirical betas between a and b
      emp_beta <- tryCatch(sum(abs(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu)), error = function(e) {NA})
      
      # Compute l1 of empirical pcors between a and b
      emp_pcor <- tryCatch(sum(abs(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu)), error = function(e) {NA})
      
      # Compute l1 of empirical kappas between a and b
      # emp_kappa <- tryCatch(sum(abs(fitemp_a[[mod_a]]$kappa_mu - fitemp_b[[mod_b]]$kappa_mu)), error = function(e) {NA})
    }
    
    # Save results
    cc_res_beta <- data.frame(null = c(unlist(null_a[["beta"]]), unlist(null_b[["beta"]])),
                              model_ind = c(rep(mod_a, n_draws), rep(mod_b, n_draws)),
                              mod = c(rep("mod_a", n_draws), rep("mod_b", n_draws)),
                              emp = rep(emp_beta, n_draws*2),
                              comp = rep(comparison, n_draws*2))
    
    
    cc_res_pcor <- data.frame(null = c(unlist(null_a[["pcor"]]), unlist(null_b[["pcor"]])),
                              model_ind = c(rep(mod_a, n_draws), rep(mod_b, n_draws)),
                              mod = c(rep("mod_a", n_draws), rep("mod_b", n_draws)),
                              emp = rep(emp_pcor, n_draws*2),
                              comp = rep(comparison, n_draws*2))
    
    # cc_res_kappa <- data.frame(null = c(unlist(null_a[["kappa"]]), unlist(null_b[["kappa"]])),
    #                            emp = rep(emp_kappa, n_draws*2),
    #                            comp = rep(comparison, n_draws*2))
    
    
 
  
  l_cc_res <- list()
  l_cc_res[["beta"]] <- cc_res_beta
  l_cc_res[["pcor"]] <- cc_res_pcor
  # l_cc_res[["kappa"]] <- cc_res_kappa
  
  cc_res <- dplyr::bind_rows(l_cc_res, .id = "mat")
  
  return(cc_res)
} # end function










# -------------------------------------------------------------------------
# Comparison Evaluation ---------------------------------------------------
# -------------------------------------------------------------------------
# Evaluate Cross-Comparison with Matrix Norm ------------------------------
#' Evaluate Cross-Comparison with Matrix Norm
#'
#' @param df_res Dataframe with posterior predictive and empirical distance measures. 
#' @param plot Should results be plotted with densities? Defaults to FALSE.
#'
#' @return List with number of posterior predictive differences that were smaller than empirical difference. 
#' @export
#'
#' 

cross_compare_eval <- function(l_res){
  
  ### Betas
  df_res_beta <- l_res[["beta"]]
  # Obtain model indexes
  model_ind_a <- unique(df_res_beta$model_ind)[1]
  model_ind_b <- unique(df_res_beta$model_ind)[2]
  
  # Check if model even converged
  # if(is.na(df_res_beta$null[df_res_beta$model_ind == model_ind_a])){
  #   teststat_a_beta <- NULL
  # }
  # if(is.na(df_res_beta$null[df_res_beta$model_ind == model_ind_b])){
  #   teststat_b_beta <- NULL
  # }
  
  # Replace infinite with NA
  df_res_beta$null[is.infinite(df_res_beta$null)] <- NA
  df_res_beta$emp[is.infinite(df_res_beta$emp)] <- NA
  
    # Number of posterior difference > empirical difference
    teststat_a_beta <- sum(df_res_beta$null[df_res_beta$model_ind == model_ind_a] > df_res_beta$emp[df_res_beta$model_ind == model_ind_a])
    teststat_b_beta <- sum(df_res_beta$null[df_res_beta$model_ind == model_ind_b] > df_res_beta$emp[df_res_beta$model_ind == model_ind_b])
    
    
    

  ### Pcor
  df_res_pcor <- l_res[["pcor"]]
  # Obtain model indexes
  model_ind_a <- unique(df_res_pcor$model_ind)[1]
  model_ind_b <- unique(df_res_pcor$model_ind)[2]
  
  # Check if model even converged
  # if(is.na(df_res_pcor$null[df_res_pcor$model_ind == model_ind_a])){
  #   teststat_a_pcor <- NULL
  # }
  # if(is.na(df_res_pcor$null[df_res_pcor$model_ind == model_ind_b])){
  #   teststat_b_pcor <- NULL
  # }
  #     
  
  # Replace infinite with NA
  df_res_pcor$null[is.infinite(df_res_pcor$null)] <- NA
  df_res_pcor$emp[is.infinite(df_res_pcor$emp)] <- NA
  
      # Number of posterior difference > empirical difference
      teststat_a_pcor <- sum(df_res_pcor$null[df_res_pcor$model_ind == model_ind_a] > df_res_pcor$emp[df_res_pcor$model_ind == model_ind_a])
      teststat_b_pcor <- sum(df_res_pcor$null[df_res_pcor$model_ind == model_ind_b] > df_res_pcor$emp[df_res_pcor$model_ind == model_ind_b])
    


  ### Store Output
  testres <- list(beta_a = teststat_a_beta,
                  beta_b = teststat_b_beta,
                  pcor_a = teststat_a_pcor,
                  pcor_b = teststat_b_pcor,
                  comp = df_res_beta$comp[[1]], # get type of comparison
                  dgp = l_res[["params"]]$dgp,
                  tp = l_res[["params"]]$tp,
                  comp_graph = l_res[["params"]]$comp_graph)  
  testres
}







# Evaluate Within-Comparison ----------------------------------------------
within_compare_eval <- function(l_res,
                                pcor = TRUE,
                                kappa = TRUE){
  
  # if input is not a list, i.e. did not converge
  if(!is.list(l_res)){
    wcompres <- list(beta_a = NA,
                     beta_b = NA,
                     pcor_a = NA,
                     pcor_b = NA,
                     mod_a = NA, 
                     mod_b = NA,
                     comp = NA)
    return(wcompres)
  }
  
  ### Betas
  df_res <- as.data.frame(l_res)
  df_res_beta <- subset(df_res, mat == "beta")
  
  # Obtain model indexes
  model_ind_a <- unique(df_res_beta$model_ind)[1]
  model_ind_b <- unique(df_res_beta$model_ind)[2]
  
  # if both indexes are the same, there is only one unique element
  if(length(unique(df_res_beta$model_ind)) == 1){
    model_ind_b <- unique(df_res_beta$model_ind)[1]
  }
  
  # Number of posterior difference > empirical difference
  teststat_a_beta <- sum(df_res_beta$null[df_res_beta$mod == "mod_a"] > df_res_beta$emp[df_res_beta$mod == "mod_a"], na.rm = TRUE)
  teststat_b_beta <- sum(df_res_beta$null[df_res_beta$mod == "mod_b"] > df_res_beta$emp[df_res_beta$mod == "mod_b"], na.rm = TRUE)
  
  
  if(isTRUE(pcor)){
    ### Pcor
    df_res_pcor <- subset(df_res, mat == "pcor")
    # # Obtain model indexes
    # "mod_a" <- unique(df_res_pcor$mod)[1]
    # "mod_b" <- unique(df_res_pcor$mod)[2]
    
    # Number of posterior difference > empirical difference
    teststat_a_pcor <- sum(df_res_pcor$null[df_res_pcor$mod == "mod_a"] > df_res_pcor$emp[df_res_pcor$mod == "mod_a"], na.rm = TRUE)
    teststat_b_pcor <- sum(df_res_pcor$null[df_res_pcor$mod == "mod_b"] > df_res_pcor$emp[df_res_pcor$mod == "mod_b"], na.rm = TRUE)
    
  }
  
  # if(isTRUE(kappa)){
  #   ### kappa
  #   df_res_kappa <- subset(df_res, mat == "kappa")
  #   # # Obtain model indexes
  #   # "mod_a" <- unique(df_res_kappa$mod)[1]
  #   # "mod_b" <- unique(df_res_kappa$mod)[2]
  #   
  #   # Number of posterior difference > empirical difference
  #   teststat_a_kappa <- sum(df_res_kappa$null[df_res_kappa$mod == "mod_a"] > df_res_kappa$emp[df_res_kappa$mod == "mod_a"], na.rm = TRUE)
  #   teststat_b_kappa <- sum(df_res_kappa$null[df_res_kappa$mod == "mod_b"] > df_res_kappa$emp[df_res_kappa$mod == "mod_b"], na.rm = TRUE)
  #   
  # }

  
  

  wcompres <- list(beta_a = teststat_a_beta,
                   beta_b = teststat_b_beta,
                   pcor_a = teststat_a_pcor,
                   pcor_b = teststat_b_pcor,
                   # kappa_a = teststat_a_kappa, 
                   # kappa_b = teststat_b_kappa,
                   mod_a = model_ind_a, 
                   mod_b = model_ind_b,
                   comp = df_res_beta$comp[[1]]) # get type of comparison
                   # dgp = l_res$params$dgp,
                   # tp = l_res$params$tp,
                   # comp_graph = l_res$params$comp_graph)  
  wcompres
}







# Compare to DGP ----------------------------------------------------------
# This function compares the BGGM VAR estimates to the true data-generating process
compare_dgp <- function(true, 
                        est_bggm,
                        comp_gvar = TRUE,
                        est_gvar,
                        n = n_ind){
  
  # replicate dgp 
  dgp <- list()
  for(i in 1:n){
    dgp[[i]] <- list()
    dgp[[i]]$beta <- true$beta
    # dgp[[i]]$kappa <- true$kappa
    dgp[[i]]$PCC <- true$PCC
  }
  
  # storage
  l_diff_beta <- list()
  l_diff_pcor <- list()
  out <- list()
  ## Differences
  # Look at difference to true values for bggm
  l_diff_beta[["true_bggm"]] <- map2(dgp, est_bggm, .f = function(x,y){abs(x$beta-t(y$beta_mu))})
  l_diff_pcor[["true_bggm"]] <- map2(dgp, est_bggm, .f = function(x,y){abs(x$PCC-y$pcor_mu)})
  
  if(isTRUE(comp_gvar)){
    # Look at difference to true values for gvar
    # delete intercepts
    l_diff_beta[["true_gvar"]] <- map2(dgp, est_gvar, .f = function(x,y){abs(x$beta-y$beta[,-1])})
    l_diff_pcor[["true_gvar"]] <- map2(dgp, est_gvar, .f = function(x,y){abs(x$PCC-y$PCC)})
    
    # Look at difference between bggm and gvar
    # delete intercepts
    l_diff_beta[["bggm_gvar"]] <- map2(est_bggm, est_gvar, .f = function(x,y){abs(t(x$beta_mu)-y$beta[,-1])})
    l_diff_pcor[["bggm_gvar"]] <- map2(est_bggm, est_gvar, .f = function(x,y){abs(x$pcor_mu-y$PCC)})
    
  }
  
  
  # Aggregate differences
  # TODO add standard deviations? or even look at raw differences across all individuals

  
  out[["diff_beta_mean"]] <- lapply(l_diff_beta, function(x){apply(simplify2array(x), 1:2, mean)})
  out[["diff_pcor_mean"]] <- lapply(l_diff_pcor, function(x){apply(simplify2array(x), 1:2, mean)})
  
  ## Correlations
  # Look at correlations between bggm and gvar
  # # delete intercepts
  # 
  # out[["cor"]] <- map2(est_bggm, est_gvar, .f = function(x,y){cor(t(x$beta_mu)-y$beta[,-1])})
  # 
  # # aggregate correlations
  
  out
  
}











# -------------------------------------------------------------------------
# BGGM-VAR-Simulation -----------------------------------------------------
# -------------------------------------------------------------------------
# Correlation with zero check ---------------------------------------------
# Inspired by Mansueto et al. 2021
cor_zero <- function(x,y){
  if(sd(x)==0 & sd(y) == 0){
    co <- 1
  } else if(sd(x)==0 | sd(y) == 0){
    co <- 0
  }  else {
    co <- cor(x, y)
  }
  co
}


# Compare BGGM GVAR -------------------------------------------------------
# Function to compare BGGM and GVAR fitting results per simulation 
compare_bggm_gvar <- function(fit_bggm,
                              fit_gvar){
  # Output
  l_out <- list()
  
  # Create vectors
  beta_bggm_vec <- c(fit_bggm$beta_mu)
  pcor_bggm_vec <- c(fit_bggm$pcor_mu[upper.tri(fit_bggm$pcor_mu,
                                                diag = FALSE)])
  
  beta_sel_bggm_vec <- c(fit_bggm$beta_weighted_adj)
  pcor_sel_bggm_vec <- c(fit_bggm$pcor_weighted_adj[upper.tri(fit_bggm$pcor_weighted_adj,
                                                              diag = FALSE)])
  
  beta_gvar_vec <- c(t(fit_gvar$beta[,-1]))
  pcor_gvar_vec <- c(fit_gvar$PCC[upper.tri(fit_gvar$PCC, diag = FALSE)])
  
  # Compute correlation without selection
  l_out$cor_beta <- cor_zero(beta_bggm_vec, beta_gvar_vec)
  l_out$cor_pcor <- cor_zero(pcor_bggm_vec, pcor_gvar_vec)
  
  # Compute correlation with selection
  l_out$cor_beta_sel <- cor_zero(beta_sel_bggm_vec, beta_gvar_vec)
  l_out$cor_pcor_sel <- cor_zero(pcor_sel_bggm_vec, pcor_gvar_vec)
  
  
  # Absolute differences
  l_out$diff_beta <- mean(abs(fit_bggm$beta_mu - t(fit_gvar$beta[,-1])))
  l_out$diff_pcor <- mean(abs(fit_bggm$pcor_mu - fit_gvar$PCC))
  
  return(l_out)
  
}


# Sim Select --------------------------------------------------------------
# Loop over all graphs, extract BGGM fitting object and apply the selection function 
# extend select function to also include the lower and upper bound 
# https://github.com/donaldRwilliams/BGGM/blob/master/R/select.VAR_estimate.R
sim_select <- function(simobj,
                       cred = 0.95,
                       include_ci = TRUE){
  
  # Input check 
  
  ## Use BGGM select
  sel <- BGGM::select(simobj, cred = 0.95, alternative = "two.sided")
  
  # Prune results
  pcor_adj <- sel$pcor_adj
  beta_adj <- sel$beta_adj
  pcor_weighted_adj <- sel$pcor_weighted_adj
  beta_weighted_adj <- sel$beta_weighted_adj
  pcor_mu <- sel$pcor_mu
  beta_mu <- sel$beta_mu
  
  
  ## Confidence Intervals
  # Extract samples
  pcors <- simobj$fit$pcors[,,51:(simobj$iter +50)]
  beta <- simobj$fit$beta[,,51:(simobj$iter +50)]
  
  # Define bounds
  lb <- (1-cred)/2
  ub <- 1-lb
  
  # Obtain CIs
  lb_pcor <- apply(pcors, 1:2, quantile, lb)
  ub_pcor <- apply(pcors, 1:2, quantile, ub)
  lb_beta <- apply(beta, 1:2, quantile, lb)
  ub_beta <- apply(beta, 1:2, quantile, ub)
  
  
  ## Return
  res <- list(
    pcor_adj = pcor_adj,
    beta_adj = beta_adj,
    pcor_weighted_adj = pcor_weighted_adj,
    beta_weighted_adj = beta_weighted_adj,
    pcor_mu = pcor_mu,
    beta_mu = beta_mu,
    lb_pcor = lb_pcor,
    ub_pcor = ub_pcor,
    lb_beta = lb_beta,
    ub_beta = ub_beta
  )
  return(res)
  
}


# Fit graphicalVAR parallel -----------------------------------------------

fit_graphicalvar_parallel <- function(data,
                                      n, 
                                      pruneresults = TRUE, 
                                      ...){     # other arguments passed to graphicalVAR
  # Save arguments
  
  require(doParallel)
  
  # reproducible parallelization
  doRNG::registerDoRNG(seed)
  
  # Foreach
  fit <- foreach(i = seq(n), .packages = "graphicalVAR") %dopar% {
    fit_ind <- tryCatch({graphicalVAR::graphicalVAR(data = data[[i]]$data,
                                                    ...)}, 
                        error = function(e) NULL)
    
    
    
    if(isTRUE(pruneresults) & is.list(fit_ind)){
      fit_ind[c("path", "allResults", "data", "labels")] <- NULL
      
    }
    fit_ind$args <- data[[i]]$args
    return(fit_ind)
  } # end foreach
  
  
  # Output
  return(fit)
  
}     




# Evaluate BGGM Simulation ------------------------------------------------
bias <- function(e,t){
  b <- mean(abs(e - t), na.rm = TRUE)
  return(b)
}

rmse <- function(e,t){
  r <- sqrt(mean((t - e)^2, na.rm = TRUE))
  return(r)
}


eval_bggm <- function(fit,
                      cred_int = c(0.9, 0.95, 0.99),    # different credible intervals
                      nds = 100,         # number of datasets per simulation condition
                      dgp_list = l_graphs){
  
  
  # Prepare output
  l_out <- list()
  
  
  # Save arguments
  args <- fit$args
  l_out$dgp_ind <- fit$sim_cond$dgp
  l_out$tp_ind <- fit$sim_cond$n_tp
  l_out$rho_prior <- fit$sim_cond$rho_prior
  l_out$beta_prior <- fit$sim_cond$beta_prior
  
  
  ## Find corresponding true graph
  # transposing bc we use graphicalVARsim
  true_graph <- dgp_list[[l_out$dgp_ind]]
  beta_true_vec <- t(true_graph$beta)
  kappa_true <- true_graph$kappa
  
  # Calculate PCOR 
  pcor_true_vec <- -1*stats::cov2cor(kappa_true)
  
  
  # Vectors for correlations
  beta_true_vec <- c(beta_true_vec)
  pcor_true_vec <- c(pcor_true_vec[upper.tri(pcor_true_vec, diag = FALSE)])
  
  
  #--- Nonselect Method ---#
  # Point estimates
  beta_est <- fit$beta_mu
  pcor_est <- fit$pcor_mu
  
  # Vectors for correlations
  beta_est_vec <- c(beta_est)
  pcor_est_vec <- pcor_est[upper.tri(pcor_est, diag = FALSE)]
  
  
  # Compute Bias
  l_out$bias_beta <- bias(beta_est_vec, beta_true_vec)
  l_out$bias_pcor <- bias(pcor_est_vec, pcor_true_vec)
  
  # RMSE
  l_out$rmse_beta <- rmse(beta_est_vec, beta_true_vec)
  l_out$rmse_pcor <- rmse(pcor_est_vec, pcor_true_vec)
  
  # Correlations
  # TODO should I do it like this? just ignore matrix structure?
  l_out$cor_beta <- cor_zero(beta_est_vec, beta_true_vec)
  l_out$cor_pcor <- cor_zero(pcor_est_vec, pcor_true_vec)
  
  
  
  #--- Select Method ---#
  # Obtain different credible intervals
  cred_interval <- fit$cred_interval
  
  
  # Obtain estimates with selection
  beta_est_sel_vec <- fit$beta_weighted_adj
  pcor_est_sel_vec <- fit$pcor_weighted_adj
  
  # Vectors for correlations
  beta_est_sel_vec <- c(beta_est_sel_vec)
  pcor_est_sel_vec <- c(pcor_est_sel_vec[upper.tri(pcor_est_sel_vec, diag = FALSE)])
  
  # Correlations
  l_out$cor_beta_sel <- cor_zero(beta_est_sel_vec, beta_true_vec)
  l_out$cor_pcor_sel <- cor_zero(pcor_est_sel_vec, pcor_true_vec)
  
  # Bias
  l_out$bias_beta_sel <- bias(beta_est_sel_vec, beta_true_vec)
  l_out$bias_pcor_sel <- bias(pcor_est_sel_vec, pcor_true_vec)
  
  # rmse
  l_out$rmse_beta_sel <- rmse(beta_est_sel_vec, beta_true_vec)
  l_out$rmse_pcor_sel <- rmse(pcor_est_sel_vec, pcor_true_vec)
  
  # Amount of zeros
  l_out$zeros_beta_sel <- sum(beta_est_sel_vec == 0)
  l_out$zeros_pcor_sel <- sum(pcor_est_sel_vec == 0)
  
  
  ## True/False Positive/Negative
  # TP
  l_out$true_pos_beta <- sum(beta_true_vec != 0 & beta_est_sel_vec != 0)
  l_out$true_pos_pcor <- sum(pcor_true_vec != 0 & pcor_est_sel_vec != 0)
  
  # FP
  l_out$fal_pos_beta <- sum(beta_true_vec == 0 & beta_est_sel_vec != 0)
  l_out$fal_pos_pcor <- sum(pcor_true_vec == 0 & pcor_est_sel_vec != 0)  
  
  # TN
  l_out$true_neg_beta <- sum(beta_true_vec == 0 & beta_est_sel_vec == 0)
  l_out$true_neg_pcor <- sum(pcor_true_vec == 0 & pcor_est_sel_vec == 0)
  
  # FN
  l_out$fal_neg_beta <- sum(beta_true_vec != 0 & beta_est_sel_vec == 0)
  l_out$fal_neg_pcor <- sum(pcor_true_vec != 0 & pcor_est_sel_vec == 0)
  
  ## Sensitivity
  l_out$sens_beta <- l_out$true_pos_beta / (l_out$true_pos_beta + l_out$fal_neg_beta)
  l_out$sens_pcor <- l_out$true_pos_pcor / (l_out$true_pos_pcor + l_out$fal_neg_pcor)
  
  ## Specificity
  l_out$spec_beta <- l_out$true_neg_beta / (l_out$true_neg_beta + l_out$fal_pos_beta)
  l_out$spec_pcor <- l_out$true_neg_pcor / (l_out$true_neg_pcor + l_out$fal_pos_pcor)
  
  
  
  ## ...
  
  
  ## Coverage
  # Loop across different cred. ints
  df_ci <- data.frame(ci = cred_int,
                      sum_cover_beta = rep(NA, length(cred_int)),
                      sum_cover_pcor = rep(NA, length(cred_int)),
                      width_beta = rep(NA, length(cred_int)),
                      width_pcor = rep(NA, length(cred_int)) )
  for(i in 1:length(cred_int)){
    lb_beta <- c(cred_interval$beta_lb[[i]])
    ub_beta <- c(cred_interval$beta_ub[[i]])
    lb_pcor <- cred_interval$pcor_lb[[i]]
    lb_pcor <- lb_pcor[upper.tri(lb_pcor)]
    ub_pcor <- cred_interval$pcor_ub[[i]]  
    ub_pcor <- ub_pcor[upper.tri(ub_pcor)]
    
    m_cover_beta <- beta_true_vec >= lb_beta & beta_true_vec <= ub_beta 
    m_cover_pcor <- pcor_true_vec >= lb_pcor & pcor_true_vec <= ub_pcor 
    
    # Only consider upper diagonal of pcor
    df_ci[i, "sum_cover_beta"] <- sum(m_cover_beta)
    df_ci[i, "sum_cover_pcor"] <- sum(m_cover_pcor)
    
    
  }
  
  ## Average width
  
  for(i in 1:length(cred_int)){
    lb_beta <- cred_interval$beta_lb[[i]]
    ub_beta <- cred_interval$beta_ub[[i]]
    lb_pcor <- cred_interval$pcor_lb[[i]]
    ub_pcor <- cred_interval$pcor_ub[[i]]
    
    df_ci[i,"width_beta"] <- mean(ub_beta - lb_beta)
    df_ci[i,"width_pcor"] <- mean(ub_pcor - lb_pcor)
    
  }
  
  l_out$ci <- df_ci
  
  
  
  #--- Output ---#
  
  return(l_out)
}



# Evaluate GVAR Simulation ------------------------------------------------
eval_gvar <- function(fit,
                      nds = 1000,         # number of datasets per simulation condition
                      dgp_list = l_graphs){
  
  
  # Prepare output
  l_out <- list()
  
  
  # Save arguments
  args <- fit$args
  l_out$dgp_ind <- fit$sim_cond$dgp
  l_out$tp_ind <- fit$sim_cond$n_tp
  l_out$ebic <- fit$sim_cond$gamma_ebic
  l_out$lambda <- fit$sim_cond$lambda
  
  
  ## Find corresponding true graph
  # transposing because we use graphicalVARsim
  true_graph <- dgp_list[[l_out$dgp_ind]]
  beta_true <- t(true_graph$beta)
  kappa_true <- true_graph$kappa
  
  # Calculate PCOR 
  pcor_true<- -1*stats::cov2cor(kappa_true)
  
  # Vectors for correlations
  beta_true_vec <- c(beta_true)
  pcor_true_vec <- c(pcor_true[upper.tri(pcor_true, diag = FALSE)])
  
  #--- Nonselect Method ---#
  # Point estimates
  beta_est <- t(fit$beta[,-1])
  pcor_est <- fit$PCC
  
  # Vectors for correlations
  beta_est_vec <- c(beta_est)
  pcor_est_vec <- c(pcor_est[upper.tri(pcor_est, diag = FALSE)])
  
  
  # Compute Bias
  l_out$bias_beta <- bias(beta_est_vec, beta_true_vec)
  l_out$bias_pcor <- bias(pcor_est_vec, pcor_true_vec)
  
  # Compute rmse
  l_out$rmse_beta <- rmse(beta_est_vec, beta_true_vec)
  l_out$rmse_pcor <- rmse(pcor_est_vec, pcor_true_vec)
  
  # Correlations
  # add conditions for very sparse matrices
  l_out$cor_beta <- cor_zero(beta_est_vec, beta_true_vec)
  l_out$cor_pcor <- cor_zero(pcor_est_vec, pcor_true_vec)
  
  # Sum of zeros
  l_out$zeros_beta <- sum(beta_est_vec == 0)
  l_out$zeros_pcor <- sum(pcor_est_vec == 0)
  
  
  ## True/False Positive/Negative
  # TP
  l_out$true_pos_beta <- sum(beta_true_vec != 0 & beta_est_vec != 0)
  l_out$true_pos_pcor <- sum(pcor_true_vec != 0 & pcor_est_vec != 0)
  
  # FP
  l_out$fal_pos_beta <- sum(beta_true_vec == 0 & beta_est_vec != 0)
  l_out$fal_pos_pcor <- sum(pcor_true_vec == 0 & pcor_est_vec != 0)  
  
  # TN
  l_out$true_neg_beta <- sum(beta_true_vec == 0 & beta_est_vec == 0)
  l_out$true_neg_pcor <- sum(pcor_true_vec == 0 & pcor_est_vec == 0)
  
  # FN
  l_out$fal_neg_beta <- sum(beta_true_vec != 0 & beta_est_vec == 0)
  l_out$fal_neg_pcor <- sum(pcor_true_vec != 0 & pcor_est_vec == 0)
  
  ## Sensitivity
  l_out$sens_beta <- l_out$true_pos_beta / (l_out$true_pos_beta + l_out$fal_neg_beta)
  l_out$sens_pcor <- l_out$true_pos_pcor / (l_out$true_pos_pcor + l_out$fal_neg_pcor)
  
  ## Specificity
  l_out$spec_beta <- l_out$true_neg_beta / (l_out$true_neg_beta + l_out$fal_pos_beta)
  l_out$spec_pcor <- l_out$true_neg_pcor / (l_out$true_neg_pcor + l_out$fal_pos_pcor)
  
  
  
  
  #--- Output ---#
  
  return(l_out)
}



# TODO add args
eval_across_gvar <- function(fit, 
                             true_graph,
                             n_rep = 1000){
  # output
  l_out <- list()
  
  # Convert data into matrices
  beta_est_vec <- unlist(lapply(fit, function(x) as.numeric(c(t(x$beta[,-1])))))
  beta_true_vec <- rep(as.numeric(c(t(true_graph$beta))), n_rep)
  beta_full_mat <- matrix(c(beta_est_vec, beta_true_vec), ncol = 2)
  
  pcor_est_vec <- unlist(lapply(fit, function(x) as.numeric(c(x$PCC[upper.tri(x$PCC, diag = FALSE)]))))
  pcor_true <- -1*stats::cov2cor(true_graph$kappa)
  pcor_true_vec <- rep(as.numeric(c(pcor_true[upper.tri(pcor_true, diag = FALSE)])), n_rep)
  pcor_full_mat <- matrix(c(pcor_est_vec, pcor_true_vec), ncol = 2)
  
  # Example calculation
  
  # Bias
  l_out$beta_bias<- mean(abs(beta_full_mat[,1]-beta_full_mat[,2]))
  l_out$pcor_bias<- mean(abs(pcor_full_mat[,1]-pcor_full_mat[,2]))
  
  # RMSE
  beta_squared_error <- (beta_full_mat[,1]-beta_full_mat[,2])^2
  l_out$beta_mse <- mean(beta_squared_error)
  l_out$beta_rmse <- sqrt(l_out$beta_mse)
  pcor_squared_error <- (pcor_full_mat[,1]-pcor_full_mat[,2])^2
  l_out$pcor_mse <- mean(pcor_squared_error)
  l_out$pcor_rmse <- sqrt(l_out$pcor_mse)
  
  # number of estimates
  nrow_beta <- nrow(beta_full_mat)
  nrow_pcor <- nrow(pcor_full_mat)
  # TP
  l_out$true_pos_beta <- sum(beta_full_mat[,2] != 0 & beta_full_mat[,1] != 0)/nrow_beta
  l_out$true_pos_pcor <- sum(pcor_full_mat[,2] != 0 & pcor_full_mat[,2] != 0)/nrow_pcor
  
  # FP
  l_out$fal_pos_beta <- sum(beta_full_mat[,2] == 0 & beta_full_mat[,1] != 0)/nrow_beta
  l_out$fal_pos_pcor <- sum(pcor_full_mat[,2] == 0 & pcor_full_mat[,2] != 0)/nrow_pcor  
  
  # TN
  l_out$true_neg_beta <- sum(beta_full_mat[,2] == 0 & beta_full_mat[,1] == 0)/nrow_beta
  l_out$true_neg_pcor <- sum(pcor_full_mat[,2] == 0 & pcor_full_mat[,2] == 0)/nrow_pcor
  
  # FN
  l_out$fal_neg_beta <- sum(beta_full_mat[,2] != 0 & beta_full_mat[,1] == 0)/nrow_beta
  l_out$fal_neg_pcor <- sum(pcor_full_mat[,2] != 0 & pcor_full_mat[,2] == 0)/nrow_pcor
  
  ## Sensitivity
  l_out$sens_beta <- l_out$true_pos_beta / (l_out$true_pos_beta + l_out$fal_neg_beta)
  l_out$sens_pcor <- l_out$true_pos_pcor / (l_out$true_pos_pcor + l_out$fal_neg_pcor)
  
  ## Specificity
  l_out$spec_beta <- l_out$true_neg_beta / (l_out$true_neg_beta + l_out$fal_pos_beta)
  l_out$spec_pcor <- l_out$true_neg_pcor / (l_out$true_neg_pcor + l_out$fal_pos_pcor)
  
  
  ## Correlations
  l_out$cor_beta <- cor_zero(beta_full_mat[,1], beta_full_mat[,2])
  l_out$cor_pcor <- cor_zero(pcor_full_mat[,1], pcor_full_mat[,2])
  
  
  # return(beta_full_mat)
  return(l_out)
  
}






# -------------------------------------------------------------------------
# Empirical Example -------------------------------------------------------
# -------------------------------------------------------------------------
# Effective Sample Size VAR -----------------------------------------------
#' Compute Effective Sample Sizes for MCMC Samples of Beta and Partial Correlation Coefficients
#'
#' This function computes the effective sample sizes (ESS) of MCMC samples of VAR and partial correlation coefficients (pcor) based on the provided MCMC fit object.
#'
#' @param fitobj A list containing a BGGM fit object.
#' @param burnin An integer indicating the number of burn-in iterations to discard. Default is 50.
#'
#' @return A list with two elements: ess_beta and ess_pcor. ess_beta contains the ESS of MCMC samples of VAR, and ess_pcor contains the ESS of MCMC samples of partial correlation coefficients.
#'

#'
#' @import coda
#' @export

var_ess <- function(fitobj,
                    burnin = 50){
  # Input Information
  it <- fitobj$iter
  p <- fitobj$p
  
  ## Get samples
  beta <- fitobj$fit$beta[,,(burnin+1):(iterations+burnin)]
  pcor <- fitobj$fit$pcors[,,(burnin+1):(iterations+burnin)]
  
  # Transform to mcmc objects
  mcmc_beta <- coda::as.mcmc(t(matrix(beta, p*p, iterations)))
  mcmc_pcor <- coda::as.mcmc(t(matrix(pcor, p*p, iterations)))
  
  # correct variable names
  # column after column 
  cnames <- colnames(fitobj$Y)
  cnames_lag <- paste0(colnames(fitobj$Y), ".l1")
  
  beta_names <- c(sapply(cnames, paste, cnames_lag, sep = "--"))
  pcor_names <- c(sapply(cnames, paste, cnames, sep = "--"))
  
  ## Calculate ESS
  ess_beta <- coda::effectiveSize(mcmc_beta)
  ess_pcor <- coda::effectiveSize(mcmc_pcor)
  
  names(ess_beta) <- beta_names
  names(ess_pcor) <- pcor_names
  
  ## Return
  l_out <- list(
    ess_beta = ess_beta,
    ess_pcor = ess_pcor
  )
  return(l_out)
  
}

# Compare VAR -------------------------------------------------------------
compare_var_old <- function(fit_a, 
                        fit_b, 
                        cutoff = 5,           # percentage level of test
                        dec_rule = "OR",
                        n_draws = 1000,
                        comp = "frob",
                        return_all = FALSE){  # return all distributions?
  
  require(magrittr)
  
  ## Helper function for computing distance metrics
  compute_metric <- function(a, b, metric) {
    tryCatch({
      if (metric == "frob") {
        norm(a - b, type = "F")
      } else if (metric == "maxdiff") {
        max(abs(a - b))
      } else if (metric == "l1") {
        sum(abs(a - b))
      }
    }, error = function(e) NA)
  }
  
  ## Create reference distributions for both models
  ref_a <- post_distance_within(fit_a, comp = comp, pred = FALSE, draws = n_draws)
  ref_b <- post_distance_within(fit_b, comp = comp, pred = FALSE, draws = n_draws)
  
  ## Empirical distance
  # Compute empirical distance as test statistic
  if(comp == "frob"){
    normtype = "F"
    # Compute Distance of empirical betas between a and b
    emp_beta <- tryCatch(norm(fit_a$beta_mu - fit_b$beta_mu, type = normtype), error = function(e) {NA})

    # Compute Distance of empirical pcors between a and b
    emp_pcor <- tryCatch(norm(fit_a$pcor_mu - fit_b$pcor_mu, type = normtype), error = function(e) {NA})

  }

  if(comp == "maxdiff"){
    # Compute maxdiff of empirical betas between a and b
    emp_beta <- tryCatch(max(abs(fit_a$beta_mu - fit_b$beta_mu)), error = function(e) {NA})

    # Compute maxdiff of empirical pcors between a and b
    emp_pcor <- tryCatch(max(abs(fit_a$pcor_mu - fit_b$pcor_mu)), error = function(e) {NA})

  }

  if(comp == "l1"){
    # Compute l1 of empirical betas between a and b
    emp_beta <- tryCatch(sum(abs(fit_a$beta_mu - fit_b$beta_mu)), error = function(e) {NA})

    # Compute l1 of empirical pcors between a and b
    emp_pcor <- tryCatch(sum(abs(fit_a$pcor_mu - fit_b$pcor_mu)), error = function(e) {NA})

  }
  emp_beta <- compute_metric(fit_a$beta_mu, fit_b$beta_mu, comp)
  emp_pcor <- compute_metric(fit_a$pcor_mu, fit_b$pcor_mu, comp)
  
  
  ## Combine results
  res_beta <- data.frame(null = c(unlist(ref_a[["beta"]]), unlist(ref_b[["beta"]])),
                         mod = c(rep("mod_a", n_draws), rep("mod_b", n_draws)),
                         emp = rep(emp_beta, n_draws*2),
                         comp = rep(comp, n_draws*2))
  
  
  res_pcor <- data.frame(null = c(unlist(ref_a[["pcor"]]), unlist(ref_b[["pcor"]])),
                         mod = c(rep("mod_a", n_draws), rep("mod_b", n_draws)),
                         emp = rep(emp_pcor, n_draws*2),
                         comp = rep(comp, n_draws*2))
  
  ## Implement decision rule "OR"
  if(dec_rule == "OR"){
    suppressWarnings(sig_beta <- res_beta %>% 
                       dplyr::group_by(mod) %>% 
                       dplyr::summarize(sum_larger = sum(null > emp)) %>% 
                       dplyr::summarize(sig = ifelse(sum_larger < cutoff * (n_draws/100), 1, 0)) %>% 
                       dplyr::summarize(sig_decision = sum(sig)) %>% 
                       dplyr::pull(sig_decision))
    
    suppressWarnings(larger_beta <- res_beta %>% 
                       dplyr::group_by(mod) %>% 
                       dplyr::summarize(sum_larger = sum(null > emp))) %>% 
      dplyr::pull(sum_larger)
    
    suppressWarnings(sig_pcor <- res_pcor %>% 
                       dplyr::group_by(mod) %>% 
                       dplyr::summarize(sum_larger = sum(null > emp)) %>% 
                       dplyr::summarize(sig = ifelse(sum_larger < cutoff * (n_draws/100), 1, 0)) %>% 
                       dplyr::summarize(sig_decision = sum(sig)) %>% 
                       dplyr::pull(sig_decision))
    
    suppressWarnings(larger_pcor<- res_pcor %>% 
                       dplyr::group_by(mod) %>% 
                       dplyr::summarize(sum_larger = sum(null > emp))) %>% 
      dplyr::pull(sum_larger)
    
    
  }
  
  # sig_beta <- as.numeric(sig_beta)
  # larger_beta <- as.numeric(larger_beta)
  # sig_pcor <- as.numeric(sig_pcor)
  # larger_pcor <- as.numeric(larger_pcor)
  
  
  if(!return_all){
    l_res <- list(sig_beta = sig_beta,
                  sig_pcor = sig_pcor,
                  # res_beta = res_beta,
                  # res_pcor = res_pcor,
                  emp_beta = emp_beta,
                  emp_pcor = emp_pcor,
                  larger_beta = larger_beta,
                  larger_pcor = larger_pcor)
    
  }
  if(isTRUE(return_all)){
    l_res <- list(sig_beta = sig_beta,
                  sig_pcor = sig_pcor,
                  res_beta = res_beta,
                  res_pcor = res_pcor,
                  emp_beta = emp_beta,
                  emp_pcor = emp_pcor,
                  larger_beta = larger_beta,
                  larger_pcor = larger_pcor)
    
  }
  
  
  
  return(l_res)
  
  
  
}

##### New compare_var

# TODO this roxygen is still incorrect, e.g. larger_beta

#' Compare Variance Components between Two Models
#'
#' Computes the empirical distance between two models based on the variance components
#' and compares them using reference distributions. Returns the p-value for the comparison
#' based on a decision rule specified by the user.
#'
#' @param fit_a Fitted model object for Model A
#' @param fit_b Fitted model object for Model B
#' @param cutoff The percentage level of the test (default: 5%)
#' @param dec_rule The decision rule to be used. Currently only supports default "OR".
#' @param n_draws The number of draws to use for reference distributions (default: 1000)
#' @param comp The distance metric to use. Should be one of "frob" (Frobenius norm), 
#' "maxdiff" (maximum  difference), or "l1" (L1 norm) (default: "frob")
#' @param return_all Logical indicating whether to return all distributions (default: FALSE)
#'
#' @return A list containing the results of the comparison. The list includes:
#'   \item{sig_beta}{The decision on whether there is a significant difference between the variance components for Model A and Model B (based on the beta parameter)}
#'   \item{sig_pcor}{The decision on whether there is a significant difference between the variance components for Model A and Model B (based on the partial correlation parameter)}
#'   \item{res_beta}{The null distribution for the variance components (based on the beta parameter) for both models}
#'   \item{res_pcor}{The null distribution for the variance components (based on the partial correlation parameter) for both models}
#'   \item{emp_beta}{The empirical distance between the two models (based on the beta parameter)}
#'   \item{emp_pcor}{The empirical distance between the two models (based on the partial correlation parameter)}
#'   \item{larger_beta}{The number of times the null hypothesis (based on the beta parameter) was rejected across all draws}
#'   \item{larger_pcor}{The number of times the null hypothesis (based on the partial correlation parameter) was rejected across all draws}
#'
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarize pull
#' @importFrom stats norm max abs sum
#'
#' @export

compare_var <- function(fit_a, 
                        fit_b, 
                        cutoff = 5,           # percentage level of test
                        dec_rule = "OR",
                        n_draws = 1000,
                        comp = "frob",
                        return_all = FALSE){  # return all distributions?
  
  require(magrittr)
  
  ## Helper function for computing distance metrics
  compute_metric <- function(a, b, metric) {
    tryCatch({
      if (metric == "frob") {
        norm(a - b, type = "F")
      } else if (metric == "maxdiff") {
        max(abs(a - b))
      } else if (metric == "l1") {
        sum(abs(a - b))
      }
    }, error = function(e) NA)
  }
  
  ## Create reference distributions for both models
  ref_a <- post_distance_within(fit_a, comp = comp, pred = FALSE, draws = n_draws)
  ref_b <- post_distance_within(fit_b, comp = comp, pred = FALSE, draws = n_draws)
  
  ## Empirical distance
  # Compute empirical distance as test statistic
  emp_beta <- compute_metric(fit_a$beta_mu, fit_b$beta_mu, comp)
  emp_pcor <- compute_metric(fit_a$pcor_mu, fit_b$pcor_mu, comp)
  
  
  ## Combine results
  res_beta <- data.frame(null = c(unlist(ref_a[["beta"]]), unlist(ref_b[["beta"]])),
                         mod = c(rep("mod_a", n_draws), rep("mod_b", n_draws)),
                         emp = rep(emp_beta, n_draws*2),
                         comp = rep(comp, n_draws*2))
  
  
  res_pcor <- data.frame(null = c(unlist(ref_a[["pcor"]]), unlist(ref_b[["pcor"]])),
                         mod = c(rep("mod_a", n_draws), rep("mod_b", n_draws)),
                         emp = rep(emp_pcor, n_draws*2),
                         comp = rep(comp, n_draws*2))
  
  ## Implement decision rule "OR"
  # Helper function
  compute_stats <- function(data, var, cutoff, n_draws) {
    sig_decision <- data %>%
      dplyr::group_by(mod) %>%
      dplyr::summarize(sum_larger = sum(null > emp)) %>%
      dplyr::summarize(sig_decision = sum(sum_larger < cutoff * (n_draws/100))) %>%
      dplyr::pull(sig_decision)
    
    sum_larger <- data %>%
      dplyr::group_by(mod) %>%
      dplyr::summarize(sum_larger = sum(null > emp))
    
    return(list(sig_decision = sig_decision, sum_larger = sum_larger))
  }
  
  if(dec_rule == "OR"){
    sig_beta <- compute_stats(res_beta, "null", cutoff, n_draws)$sig_decision
    larger_beta <- compute_stats(res_beta, "null", cutoff, n_draws)$sum_larger
    sig_pcor <- compute_stats(res_pcor, "null", cutoff, n_draws)$sig_decision
    larger_pcor <- compute_stats(res_pcor, "null", cutoff, n_draws)$sum_larger
    
  }

  if(!return_all){
    l_res <- list(sig_beta = sig_beta,
                  sig_pcor = sig_pcor,
                  emp_beta = emp_beta,
                  emp_pcor = emp_pcor,
                  larger_beta = larger_beta,
                  larger_pcor = larger_pcor)
    
  }
  if(isTRUE(return_all)){
    l_res <- list(sig_beta = sig_beta,
                  sig_pcor = sig_pcor,
                  res_beta = res_beta,
                  res_pcor = res_pcor,
                  emp_beta = emp_beta,
                  emp_pcor = emp_pcor,
                  larger_beta = larger_beta,
                  larger_pcor = larger_pcor)
    
  }
  
  
  
  return(l_res)
  
  
  
}








# Plotting method
# THIS IS ONLY TEMPORARY!
plot.compare_var <- function(compres,
                             ...){
  require(ggplot2)
  require(cowplot)
  # create df
  # dat <- rbind(compres$res_beta, compres$res_pcor)
  
  
  
  # Plotting
  plt_beta <- ggplot(compres$res_beta, 
                     aes(x = null, fill = mod))+
    geom_density(alpha = .7)+
    theme_classic()+
    ggokabeito::scale_fill_okabe_ito()+
    geom_vline(aes(xintercept = compres$emp_beta), 
               col = "red", lty = 1, linewidth = .75)+
    labs(title = "Temporal")
  
  plt_pcor <- ggplot(compres$res_pcor, 
                     aes(x = null, fill = mod))+
    geom_density(alpha = .7)+
    theme_classic()+
    ggokabeito::scale_fill_okabe_ito()+
    geom_vline(aes(xintercept = compres$emp_pcor), 
               col = "red", lty = 1, linewidth = .75)+
    labs(title = "Contemporaneous")
  
  leg <- get_legend(plt_beta)
  
  # Plot
  plt_tmp <- cowplot::plot_grid(plt_beta + theme(legend.position = "none"),
                                plt_pcor + theme(legend.position = "none"))
  
  # Add legend
  plt <- plot_grid(plt_tmp, leg, rel_widths = c(3, .4))
  plt
  
}






# Posterior Matrix Plot ---------------------------------------------------
# use inspiration from stat_wilke function here: 
# https://bookdown.org/content/8ba612b7-90f2-4ebc-b329-0159008e2340/metric-predicted-variable-with-multiple-metric-predictors.html#metric-predicted-variable-with-multiple-metric-predictors

#' posterior_plot
#'
#' Plots posterior distributions of betas and partial correlations (pcors).
#'
#' @param object An object of class 'bggm'.
#' @param mat A matrix to use for plotting. Default is beta.
#' @param cis A numeric vector of credible intervals to use for plotting. Default is c(0.8, 0.9, 0.95).
#'
#' @import ggdist
#' @import tidyr
#' @import dplyr
#' @importFrom BGGM posterior_samples
#' @importFrom ggplot aes facet_grid geom_vline labs scale_alpha scale_fill_brewer
#' @importFrom ggdist stat_pointinterval stat_slab
#' @importFrom tidyr separate_wider_delim pivot_longer
#' @export

posterior_plot <- function(object,
                           mat = beta,
                           cis = c(0.8, 0.9, 0.95)){   # credible intervals for plotting
  
  require(ggdist)   # visualize uncertainty
  
  
  # Obtain samples
  samps <- BGGM::posterior_samples(object)
  
  # throw a bug if one of the variable names contains an underscore
  if(length(grep("_", colnames(object$Y))) > 0) {
    stop("Column names must not contain an underscore. Please rename.")}
  
  # Split into betas and pcors
  beta_cols <- grep(".l1", colnames(samps), value = TRUE)
  pcor_cols <- grep("--", colnames(samps), value = TRUE)
  
  beta_samps <- as.data.frame(samps[,beta_cols])
  pcor_samps <- as.data.frame(samps[,pcor_cols])
  
  # Pivot longer
  beta <- beta_samps %>%
    as.data.frame() %>%
    dplyr::mutate(iteration = dplyr::row_number()) %>% 
    tidyr::pivot_longer(cols = !iteration, names_to = "edge", values_to = "value") %>% 
    # split edge description into nodes
    tidyr::separate_wider_delim(cols = edge, delim = "_",
                                names = c("dv" ,"iv"))
  
  pcor <- pcor_samps %>%
    as.data.frame() %>%
    dplyr::mutate(iteration = dplyr::row_number()) %>% 
    tidyr::pivot_longer(cols = !iteration,names_to = "edge", values_to = "value") %>% 
    # split edge description into nodes
    tidyr::separate_wider_delim(cols = edge, delim = "--",
                                names = c("dv" ,"iv"))
  
  
  
  # Create matrix layout
  
  if(mat == beta){
    # Start plotting
    # TODO add option for numerical value per plot
    beta_plot <- beta %>% 
      dplyr::group_by(dv, iv) %>% 
      dplyr::mutate(mean_value = mean(value, na.rm = TRUE)) %>% 
      dplyr::ungroup() %>% 
      ggplot(aes(x = value))+
      # ggdist::stat_halfeye(aes(fill = after_stat(level)), .width = cis)+
      ggdist::stat_slab(aes(fill = after_stat(level), alpha = abs(mean_value)), .width = c(cis, 1)) +
      ggdist::stat_pointinterval(aes(alpha = abs(mean_value)), size = 1) +
      scale_alpha(guide = "none")+
      facet_grid(iv~dv,
                 switch = "y")+
      ggdist::theme_ggdist()+
      geom_vline(xintercept = 0, linetype = "dashed")+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())+
      scale_fill_brewer() +
      labs(y = "",
           fill = "CI")+
      ylim(-0.1, 1)
    

    print(beta_plot)
  }
  
  if(mat == pcor){
    # has to be made symmetric, maybe reorder variable
    # create duplicate plot
    pcor_tmp1<- pcor %>% 
      dplyr::group_by(dv, iv) %>% 
      dplyr::mutate(mean_value = mean(value, na.rm = TRUE)) %>% 
      dplyr::ungroup()
    pcor_tmp2 <- pcor %>% 
      dplyr::group_by(dv, iv) %>% 
      dplyr::mutate(mean_value = mean(value, na.rm = TRUE)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(dv2 = dv, iv2 = iv) %>% 
      dplyr::mutate(dv = iv2, iv = dv2) %>% 
      dplyr::select(-c(iv2, dv2))
    pcor <- rbind(pcor_tmp1, pcor_tmp2)
    
    pcor_plot <- pcor %>% 
      dplyr::group_by(dv, iv) %>% 
      dplyr::mutate(mean_value = mean(value, na.rm = TRUE)) %>% 
      dplyr::ungroup() %>% 
      ggplot(aes(x = value))+
      # ggdist::stat_halfeye(aes(fill = after_stat(level)), .width = cis)+
      ggdist::stat_slab(aes(fill = after_stat(level), alpha = abs(mean_value)), .width = c(cis, 1)) +
      ggdist::stat_pointinterval(aes(alpha = abs(mean_value)), size = 1) +
      # scale_alpha_manual(guide = "none")+
      facet_grid(iv~dv,
                 switch = "y")+
      ggdist::theme_ggdist()+
      scale_alpha(guide = "none")+
      geom_vline(xintercept = 0, linetype = "dashed")+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())+
      scale_fill_brewer() +
      labs(y = "",
           fill = "CI",
           alpha = "")+
      ylim(-0.1, 1)
    
    print(pcor_plot)
  }


  # return(beta)
  # return(pcor)
  
}




# Plot test results -------------------------------------------------------
#' Plot posterior difference test results
#'
#' This function plots the posterior difference test results using ggplot2.
#'
#' @param comp_obj Comparison object for the posterior difference test.
#' @param modmat Model matrix used.
#' @param ref_dist Reference distribution used.
#' @param emp_diff Empirical difference used.
#' @param ind Index of the model used.
#' @param comp_type Type of comparison used.
#'
#' @return A ggplot2 object.
#'
#' @import ggplot2
#' @import ggokabeito
#' @import dplyr
#'
#' @export

plot_test <- function(comp_obj,
                      modmat,
                      ref_dist = null,
                      emp_diff = emp,
                      ind = model_ind,
                      comp_type = comp){
  
  # Store params
  pr <- comp_obj$params
  
  # get comp type
  ct <- comp_obj$res %>% distinct({{comp_type}})
  
  # Get matrix as character
  c_matrix <- deparse(substitute(modmat))
  
  
  comp_obj$res %>% 
    filter(mat == c_matrix) %>%  
    ggplot()+
    geom_histogram(aes(x = {{ref_dist}}, fill = as.factor({{ind}})), alpha = 0.65,  position = "identity", bins = 100)+
    geom_vline(aes(xintercept = max({{emp_diff}})))+
    theme_minimal()+
    labs(x = paste0(ct, " Norm Value"),
         fill = "Model",
         caption = paste0("DGP: ", pr$dgp,", TP: ", pr$tp, ", Comparison Graph: ", pr$comp_graph, ", Matrix: ", c_matrix))+
    ggokabeito::scale_fill_okabe_ito()
  
}








# -------------------------------------------------------------------------
# Miscellaneous -----------------------------------------------------------
# -------------------------------------------------------------------------
# Expand grid unique ------------------------------------------------------
# Adapted from GIMME
# https://rdrr.io/cran/gimme/src/R/expand.grid.unique.R

expand_grid_unique <- function(mod_a, mod_b){
  g <- function(i){
    z <- setdiff(mod_b, mod_a[seq_len(i - 1)])
    if(length(z)) cbind(mod_a[i], z, deparse.level = 0)
  }
  comb <- do.call(rbind, lapply(seq_along(mod_a), g))
  # delete where mod_a == mod_b 
  comb <- comb[comb[,1] != comb[,2],]
  comb
}


# Plotting Theme ----------------------------------------------------------
theme_compare <- function(){
  # add google font
  sysfonts::font_add_google("News Cycle", "news")
  # use showtext
  showtext::showtext_auto()
  # theme
  ggplot2::theme_minimal(base_family = "news") +
    ggplot2::theme(
      # remove minor grid
      panel.grid.minor = ggplot2::element_blank(),
      # Title and Axis Texts
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1.1), hjust = 0.5),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1)),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, b = 10)),
      
      # Faceting
      strip.text = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.1), hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      # Grid
      panel.grid = ggplot2::element_line(colour = "#F3F4F5"),
      # Legend
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "top",
      legend.justification = 1,
      # Panel/Facets
      panel.spacing.y = ggplot2::unit(1.5, "lines")
    )
}



