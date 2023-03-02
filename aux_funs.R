# Find ToDos --------------------------------------------------------------
# todor::todor(c("TODO"))



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
      tmp_beta <- m_beta + runif(n = b_i*b_j, min = -noise[n], max = noise[n])
      l_out_noise[[n]]$beta <- tmp_beta
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
      names(l_out_noise)[[n]] <- paste0("noise", noise[n])
    }
  } # end noise
  
  
  ### Permute matrix
  l_out_perm <- list()
  l_out_perm[[1]] <- list()
  l_out_perm[[1]]$beta <- permute_mat_col(m_beta, permute_index)
  l_out_perm[[1]]$kappa <- permute_mat_col(m_kappa, permute_index)  
  l_out_perm[[1]]$args <- "perm"
  names(l_out_perm) <- "perm"
  
  
  ### Output
  # Combine truegraph, maxchange, and added noise
  l_out <- c(l_out, l_out_change, l_out_noise, l_out_perm)
  return(l_out)
  
}



# Permute matrix columns --------------------------------------------------
# This permutes matrix columns while keeping diagonal elements on the diagonal
permute_mat_col <- function(mat,
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
  
  # ncores = parallel::detectCores() - 2
  # cl = makeCluster(ncores)
  # registerDoParallel(cl)
  
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
      # return name of data-generating process
      raw_data$args <- args
      raw_data 
      
    }

    
  }
  return(data)
  # stopCluster(cl)
}


# Fit VAR parallel --------------------------------------------------------
# TODO add saving of function arguments

fit_var_parallel <- function(data, 
                             n,
                             rho_prior, 
                             beta_prior,
                             seed,
                             iterations,
                             get_kappa = TRUE,
                             posteriorsamples = FALSE,
                             pruneresults = FALSE){
 
  if(n != length(data)){
    warning("The n provided does not match the number of available data frames")
  }
  require(doParallel)
  # require("BGGM", lib.loc = "C:/Users/Bjoern/R-dev")
  
  # reproducible parallelization
  registerDoRNG(seed)
  # if we take simulated data in a list object
  if(isFALSE(posteriorsamples)){
    fit <- foreach(i = seq(n), .packages = "BGGM") %dopar% {
      fit_ind <- list()
      if(is.list(data[[i]]$data) | is.numeric(data[[i]]$data)){
        fit_ind <- tryCatch({BGGM::var_estimate(data[[i]]$data,
                                          rho_sd = rho_prior,
                                          beta_sd = beta_prior,
                                          iter = iterations,
                                          progress = FALSE,
                                          seed = seed)}, error = function(e) NULL)
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
          beta <- fit_ind$fit$beta
          kappa <- fit_ind$fit$kappa
          beta_mu <- fit_ind$beta_mu 
          pcor_mu <- fit_ind$pcor_mu
          kappa_mu <- fit_ind$kappa_mu
          fit_ind <- list()
          fit_ind$fit <- list()
          fit_ind$beta_mu <- beta_mu
          fit_ind$pcor_mu <- pcor_mu
          fit_ind$kappa_mu <- kappa_mu
          fit_ind$fit$beta <- beta
          fit_ind$fit$kappa <- kappa
          
        } 
        
      }
      else fit_ind <- NA
      
      fit_ind  
    }
    
    
  }
  
  
  # if we fit the data on samples generated from the posterior
  if(isTRUE(posteriorsamples)){
    print("Fitting to posterior samples")

    fit <- foreach(i = seq(n), .packages = "BGGM") %dopar% {
      fit_ind <- list()
      if(is.list(data[[i]])){
        fit_ind <- tryCatch({BGGM::var_estimate(data[[i]],
                                          rho_sd = rho_prior,
                                          beta_sd = beta_prior,
                                          iter = iterations,
                                          progress = FALSE,
                                          seed = seed)}, error = function(e) NULL)
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
          beta_mu <- fit_ind$beta_mu 
          kappa_mu <- fit_ind$kappa
          pcor_mu <- fit_ind$pcor_mu
          fit_ind <- list()
          fit_ind$beta_mu <- beta_mu
          fit_ind$kappa_mu <- kappa_mu
          fit_ind$pcor_mu <- pcor_mu
          
        } 
        
        
      }
      else fit_ind <- NA
      
      fit_ind  
    }  

    
  } # end isTRUE posteriorsamples
  

  
    return(fit)

  
  
}



# Fit var parallel to posterior samples -----------------------------------
# TODO maybe merge this to the other fit function?
# TODO counter for attempts is only implemented for posteriorsamples = TRUE
# TODO Should maybe parallelize at the level of posterior datasets, not individuals
# this is especially relevant going forward, when this function should be used by researchers


fit_var_parallel_post <- function(data, 
                                  n,
                                  nds = n_postds, 
                                  rho_prior, 
                                  beta_prior,
                                  seed,
                                  iterations,
                                  get_kappa = TRUE,
                                  posteriorsamples = FALSE,
                                  pruneresults = FALSE){
  
  if(n != length(data)){
    warning("The n provided does not match the number of available data frames")
  }
  require(doParallel)
  # require("BGGM", lib.loc = "C:/Users/Bjoern/R-dev")
  
  # reproducible parallelization
  registerDoRNG(seed)
  
  # if we take simulated data in a list object
  if(isFALSE(posteriorsamples)){
    fit <- foreach(i = seq(n), .packages = "BGGM") %dopar% {
      fit_ind <- list()
      if(is.list(data[[i]]$data) | is.numeric(data[[i]]$data)){
        fit_ind <- try(BGGM::var_estimate(data[[i]]$data,
                                          rho_sd = rho_prior,
                                          beta_sd = beta_prior,
                                          iter = iterations,
                                          progress = FALSE,
                                          seed = seed), silent = TRUE)
        if(isTRUE(get_kappa)){
          # check if fitting worked
          if(is.list(fit_ind)){
            # Invert covariance matrix of residuals to obtain precision matrix
            fit_ind$fit$kappa <- array(apply(fit_ind$fit$Sigma, 3, solve), 
                                       dim = dim(fit_ind$fit$Sigma))
          }
          
        }
      }
      else fit_ind <- NA
      
      fit_ind  
    }
    
    
  }
  
  
  # if we fit the data on samples generated from the posterior
  if(isTRUE(posteriorsamples)){
    print("Fitting to posterior samples")
    
    # counter for converged models 
    m <- list()
    # counter for attempted models
    c <- list()
    
    
    # loop across individuals
  l_out <- foreach(i = seq(n), .packages = "BGGM") %dopar% {
      # Counter for converged models
      m[[i]] <- 0

      fit <- list()
      # loop across datasets
      for(d in seq(nds)) {
        # we want 100 converged models
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
            if(isTRUE(pruneresults)){
              # kappa <- fit_ind$fit$kappa
              # beta <- fit_ind$fit$beta
              beta_mu <- fit_ind$beta_mu 
              kappa_mu <- fit_ind$kappa_mu
              pcor_mu <- fit_ind$pcor_mu
              # n_attempts <- c[[i]]
              fit_ind <- list()
              # fit_ind$fit <- list()
              fit_ind$beta_mu <- beta_mu
              fit_ind$kappa_mu <- kappa_mu
              fit_ind$pcor_mu <- pcor_mu
              # fit_ind$n_attempts <- n_attempts
              # fit_ind$fit$beta <- beta
              # fit_ind$fit$kappa <- kappa
              
            } 
            
            fit[[d]] <- fit_ind  
          }    # end is.list(fit_ind)
          
          
        } # end is.list(data[[i]])  
        
        
      } # end for loop 
      res <- list()
      res$fit <- fit
      
      # store parameters
      res$params <- list()
      # number of failed attempts
      res$params$n_notconv <- length(res$fit)-100
      
      # Cut away failed attempts
      res$fit <- res$fit[!sapply(res$fit, is.null)]
      
      
      return(res)
      
    }  # end foreach
  # cut away failed attempts

  return(l_out)
  } # end isTRUE posteriorsamples

  
}  # end function




# Fit VAR parallel merged -------------------------------------------------
# This is a merge of fit_var_parallel and fit_var_parallel_post into one function

fit_var_parallel_merged <- function(data, 
                                     n,
                                     nds,       # number of datasets
                                     rho_prior, 
                                     beta_prior,
                                     seed,
                                     iterations,
                                     get_kappa = TRUE,
                                     posteriorsamples = FALSE,
                                     multigroup = FALSE,
                                     pruneresults = FALSE, 
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
        if(isTRUE(get_kappa)){
          # check if fitting worked
          if(is.list(fit_ind)){
            # Invert covariance matrix of residuals to obtain precision matrix
            fit_ind$fit$kappa <- array(apply(fit_ind$fit$Sigma, 3, solve), 
                                       dim = dim(fit_ind$fit$Sigma))
            
            # Calculate mean of kappa
            fit_ind$kappa_mu <- apply(fit_ind$fit$kappa, c(1,2), mean)
            
          }
          # else{
            # n_nonconv <- 1     # binary indicator of convergence now
          # }
          # fit_ind$n_nonconv <- n_nonconv
        }
        
        # prune results for comparison purposes
        if(isTRUE(pruneresults) & is.list(fit_ind)){
          # beta <- fit_ind$fit$beta
          # kappa <- fit_ind$fit$kappa
          beta_mu <- fit_ind$beta_mu 
          pcor_mu <- fit_ind$pcor_mu
          kappa_mu <- fit_ind$kappa_mu
          args <- data[[i]]$args
          # n_nonconv <- fit_ind$n_nonconv
          fit_ind <- list()
          # fit_ind$fit <- list()       - this is probably why everything got deleted below!
          fit_ind$beta_mu <- beta_mu
          fit_ind$pcor_mu <- pcor_mu
          fit_ind$kappa_mu <- kappa_mu
          fit_ind$args <- args
          # fit_ind$n_nonconv <- n_nonconv
          # fit_ind$fit$beta <- beta
          # fit_ind$fit$kappa <- kappa

        } # end if isTRUE 
        
        if(isTRUE(select) & is.list(fit_ind)){
          sel <- sim_select(fit_ind)
          fit_ind <- list()
          fit_ind <- stats::setNames(sel, names(sel))
          fit_ind$args <- data[[i]]$args
          
        } # end isTRUE select
        
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
      saveRDS(fit, file = here::here(paste0("data/compare_sim_fits/fit_",dgp_name, "_", fit[[1]]$args$dgp, "_", fit[[1]]$args$tp,".RDS")))
      
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







# # External data to BGGM format --------------------------------------------
# # External data needs to have two objects
# # Y: response data
# # X: lagged response data
# 
# format_bggm <- function(Y){
#   if(is.data.frame(Y)){
#     Y <- scale(na.omit(Y))
#     p <- ncol(Y)
#     n <- nrow(Y)
#     Y_lag <- rbind(NA, Y)
#     colnames(Y_lag) <- paste0(colnames(Y), ".l1")
#     Y_all <- na.omit(cbind.data.frame(rbind(Y, NA), Y_lag))
#     Y <- as.matrix(Y_all[, 1:p])
#     X <- as.matrix(Y_all[, (p + 1):(p * 2)])
#     out <- list()
#     out$Y <- Y
#     out$X <- X
#     
#     
#   }
#   else out <- NA
#   
#   out
# }
# 
# 
# 
# # External data from nested list to BGGM format ---------------------------
# # same as format_bggm, but has nested list as input
# format_bggm_list <- function(listname){
#   l_out <- lapply(listname, format_bggm)
#   l_out
# }






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
# PCMD not yet implemented in other functions
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


# # Distance within normal posterior samples --------------------------------
# # This function resamples from the original posterior distribution
# # so no posterior predictive stuff involved
# 
# post_distance_within <- function(post, 
#                                  comp, 
#                                  draws = 1000){
#   
#   # storage
#   dist_out <- list()
#   
#   # define the distance function based on comp
#   distance_fn_beta <- switch(comp,
#                              frob =   {function(x, y, mod_one, mod_two) norm(x$fit[[mod_one]]$beta_mu-y$fit[[mod_two]]$beta_mu, type = "F")},
#                              maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit[[mod_one]]$beta_mu-y$fit[[mod_two]]$beta_mu)))},
#                              l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit[[mod_one]]$beta_mu-y$fit[[mod_two]]$beta_mu)))}
#   )
#   distance_fn_pcor <- switch(comp,
#                              frob = {function(x, y, mod_one, mod_two) norm(x$fit[[mod_one]]$pcor_mu-y$fit[[mod_two]]$pcor_mu, type = "F")},
#                              maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit[[mod_one]]$pcor_mu-y$fit[[mod_two]]$pcor_mu)))},
#                              l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit[[mod_one]]$pcor_mu-y$fit[[mod_two]]$pcor_mu)))}
#   )
#   
#   
#   
#   
#   ## Draw two random models
#   # Obtain number of models
#   n_mod <- length(post$fit)
#   
#   # Draw pairs of models
#   mod_pairs <- replicate(draws, sample(1:n_mod, size = 2, replace = TRUE))
#   
#   for(i in seq(draws)){
#     # storage
#     dist_out[[i]] <- list()
#     mod_one <- mod_pairs[1,i]
#     mod_two <- mod_pairs[2,i]
#     
#     # if mod_one and mod_two are equal, redraw
#     if(mod_one == mod_two){
#       mod_two <- sample(1:n_mod, size = 1)
#     }
#     
#     ## Check if estimation worked
#     # Should be unneccessary if non-converged attempts were deleted
#     if(!is.list(post$fit[[mod_one]]) | !is.list(post$fit[[mod_two]])){
#       beta_distance <- NA
#       pcor_distance <- NA
#       stop("Not a list.")
#       
#       
#     } 
#     # if both elements are lists
#     else{
#       beta_distance <- distance_fn_beta(post, post, mod_one, mod_two)
#       pcor_distance <- distance_fn_pcor(post, post, mod_one, mod_two)
#       
#     }  
#     dist_out[[i]]$comp <- comp
#     dist_out[[i]]$mod_one <- mod_one
#     dist_out[[i]]$mod_two <- mod_two
#     dist_out[[i]]$beta <- beta_distance
#     dist_out[[i]]$pcor <- pcor_distance  
#     
#   } # end for loop  
#   out <- do.call(rbind, dist_out)
#   out <- as.data.frame(out)
#   
#   
#   return(out)
# }




# Distance within posterior predictive samples ---------------------------------
# Looks at differences between models sampled from the same
# "original" model, so similar to bootstrapping
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
    
    # Convert Kappas to Pcors
    
    
    
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




# 
# 
# ### Parallel beta
# post_distance_within_par <- function(post, 
#                                  comp,
#                                  pred,         # posterior predictive?
#                                  draws = 1000,
#                                  ncores){
#   
#   # storage
#   dist_out <- list()
#   
#   
#   # for posterior predictive approach
#   if(isTRUE(pred)){
#     # define the distance function based on comp
#     distance_fn_beta <- switch(comp,
#                                frob =   {function(x, y, mod_one, mod_two) norm(x$fit[[mod_one]]$beta_mu-y$fit[[mod_two]]$beta_mu, type = "F")},
#                                maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit[[mod_one]]$beta_mu-y$fit[[mod_two]]$beta_mu)))},
#                                l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit[[mod_one]]$beta_mu-y$fit[[mod_two]]$beta_mu)))}
#     )
#     distance_fn_pcor <- switch(comp,
#                                frob = {function(x, y, mod_one, mod_two) norm(x$fit[[mod_one]]$pcor_mu-y$fit[[mod_two]]$pcor_mu, type = "F")},
#                                maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit[[mod_one]]$pcor_mu-y$fit[[mod_two]]$pcor_mu)))},
#                                l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit[[mod_one]]$pcor_mu-y$fit[[mod_two]]$pcor_mu)))}
#     )
#     distance_fn_kappa <- switch(comp,
#                                 frob = {function(x, y, mod_one, mod_two) norm(x$fit[[mod_one]]$kappa_mu-y$fit[[mod_two]]$kappa_mu, type = "F")},
#                                 maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit[[mod_one]]$kappa_mu-y$fit[[mod_two]]$kappa_mu)))},
#                                 l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit[[mod_one]]$kappa_mu-y$fit[[mod_two]]$kappa_mu)))}
#     )
#     
#     
#     
#     # Obtain number of models
#     n_mod <- length(post$fit)
#     
#   }
#   
#   
#   # for posteriors of empirical models
#   if(isFALSE(pred)){
#     # define the distance function based on comp
#     # draw from all posterior samples
#     
#     # Convert Kappas to Pcors
#     post$fit$pcor<- array(apply(post$fit$kappa, 3, function(x){-1*cov2cor(x)}), dim = dim(post$fit$kappa))
#     
#     
#     
#     distance_fn_beta <- switch(comp,
#                                frob =   {function(x, y, mod_one, mod_two) norm(x$fit$beta[,,mod_one]-y$fit$beta[,,mod_two], type = "F")},
#                                maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit$beta[,,mod_one]-y$fit$beta[,,mod_two])))},
#                                l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit$beta[,,mod_one]-y$fit$beta[,,mod_two])))}
#     )
#     distance_fn_pcor <- switch(comp,
#                                frob = {function(x, y, mod_one, mod_two) norm(x$fit$pcor[,,mod_one]-y$fit$pcor[,,mod_two], type = "F")},
#                                maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit$pcor[,,mod_one]-y$fit$pcor[,,mod_two])))},
#                                l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit$pcor[,,mod_one]-y$fit$pcor[,,mod_two])))}
#     )
#     distance_fn_kappa <- switch(comp,
#                                 frob = {function(x, y, mod_one, mod_two) norm(x$fit$kappa[,,mod_one]-y$fit$kappa[,,mod_two], type = "F")},
#                                 maxdiff = {function(x, y, mod_one, mod_two) max(abs((x$fit$kappa[,,mod_one]-y$fit$kappa[,,mod_two])))},
#                                 l1 = {function(x, y, mod_one, mod_two) sum(abs((x$fit$kappa[,,mod_one]-y$fit$kappa[,,mod_two])))}
#     )
#     
#     # Obtain number of posterior samples
#     n_mod <- dim(post$fit$beta)[3]
#     
#   }
#   
#   
#   
#   ## Draw two random models
#   # TODO:
#   # Draw pairs of models, spaced far apart so we don't have autocorrelation
#   # delete burn-in iterations (to 50)
#   mod_pairs <- replicate(draws, sample(1:n_mod, size = 2, replace = TRUE))
#   
#   
#   require(doParallel)
#   cl = parallel::makeCluster(ncores)
#   registerDoParallel(cl)
#   
#   dist_out <- foreach(i = seq(draws)) %dopar% {
#     # storage
#     d_out <- list()
#     mod_one <- mod_pairs[1,i]
#     mod_two <- mod_pairs[2,i]
#     
#     # if mod_one and mod_two are equal, redraw
#     if(mod_one == mod_two){
#       mod_two <- sample(1:n_mod, size = 1)
#     }
#     
#     ## Check if estimation worked
#     # Should be unneccessary if non-converged attempts were deleted
#     if(isTRUE(pred)){
#       if(!is.list(post$fit[[mod_one]]) | !is.list(post$fit[[mod_two]])){
#         beta_distance <- NA
#         pcor_distance <- NA
#         kappa_distance <- NA
#         stop("Not a list.")
#         
#         
#       } 
#       # if both elements are lists
#       else{
#         beta_distance <- distance_fn_beta(post, post, mod_one, mod_two)
#         pcor_distance <- distance_fn_pcor(post, post, mod_one, mod_two)
#         kappa_distance <- distance_fn_kappa(post, post, mod_one, mod_two)
#       }  
#     }
#     
#     if(isFALSE(pred)){
#       if(!is.list(post) | !is.list(post)){
#         beta_distance <- NA
#         pcor_distance <- NA
#         kappa_distance <- NA
#         stop("Not a list.")
#         
#       } 
#       # if both elements are lists
#       else{
#         beta_distance <- distance_fn_beta(post, post, mod_one, mod_two)
#         pcor_distance <- distance_fn_pcor(post, post, mod_one, mod_two)
#         kappa_distance <- distance_fn_kappa(post, post, mod_one, mod_two)
#       }  
#     }
#     
#     
#     
#     
#     d_out$comp <- comp
#     d_out$mod_one <- mod_one
#     d_out$mod_two <- mod_two
#     d_out$beta <- beta_distance
#     d_out$pcor <- pcor_distance  
#     d_out$kappa <- kappa_distance
#     
#   } # end for loop  
#   stopCluster(cl)
#   out <- do.call(rbind, dist_out)
#   out <- as.data.frame(out)
#   
#   
#   return(out)
# }

 
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
# Change 07.02.:  deleted $res, no longer necessary. also changed model_ind == to mod == , which fixes errors when both model indicators are identical
within_compare_eval <- function(l_res,
                                pcor = TRUE,
                                kappa = TRUE){
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
  
  if(isTRUE(kappa)){
    ### kappa
    df_res_kappa <- subset(df_res, mat == "kappa")
    # # Obtain model indexes
    # "mod_a" <- unique(df_res_kappa$mod)[1]
    # "mod_b" <- unique(df_res_kappa$mod)[2]
    
    # Number of posterior difference > empirical difference
    teststat_a_kappa <- sum(df_res_kappa$null[df_res_kappa$mod == "mod_a"] > df_res_kappa$emp[df_res_kappa$mod == "mod_a"], na.rm = TRUE)
    teststat_b_kappa <- sum(df_res_kappa$null[df_res_kappa$mod == "mod_b"] > df_res_kappa$emp[df_res_kappa$mod == "mod_b"], na.rm = TRUE)
    
  }

  
  

  wcompres <- list(beta_a = teststat_a_beta,
                   beta_b = teststat_b_beta,
                   pcor_a = teststat_a_pcor,
                   pcor_b = teststat_b_pcor,
                   kappa_a = teststat_a_kappa, 
                   kappa_b = teststat_b_kappa,
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
    dgp[[i]]$kappa <- true$kappa
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





# Effective Sample Size VAR -----------------------------------------------
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








# -------------------------------------------------------------------------
# BGGM-VAR-Simulation -----------------------------------------------------
# -------------------------------------------------------------------------



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






















