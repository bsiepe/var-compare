# -------------------------------------------------------------------------
# Data Generation ---------------------------------------------------------
# -------------------------------------------------------------------------
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
                          const,
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

  }
  
  
  
  
  ### Add noise
  if(!is.null(noise)){
    l_out_noise <- list()
    ## Beta
    for(n in seq_along(noise)){
      l_out_noise[[n]] <- list()
      # Add eigenvalue check, similar to gVARsimulator function in graphicalVAR
      eigen_cond <- FALSE
      eigen_iter <- 0
      while(isFALSE(eigen_cond) & eigen_iter < 100){
        eigen_iter <- eigen_iter + 1
        noise_beta <- matrix(runif(n = b_i*b_j, min = -noise[n], max = noise[n]), 
                             nrow = b_i, ncol = b_j)
        tmp_beta <- m_beta + noise_beta
        eigen_beta <- eigen(tmp_beta)$values
        eigen_cond <- all(Re(eigen_beta)^2 + Im(eigen_beta)^2 <1)
      }
      l_out_noise[[n]]$eigen_iter <- eigen_iter
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
        eigen_kappa <- eigen(tmp_kappa)$values
        # tmp_kappa <- as.matrix(Matrix::forceSymmetric(tmp_kappa))
        psd <- matrixcalc::is.positive.semi.definite(tmp_kappa) & all(eigen_kappa > 0)  
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
  
  ### Add constant
  if(!is.null(const)){
    ## Beta
    l_out_const <- list()
    ## Beta
    for(n in seq_along(const)){
      
      # add Eigenvalue check
      eigen_cond <- FALSE
      eigen_count <- 0
      while(isFALSE(eigen_cond) & eigen_iter < 100){
        eigen_count <- eigen_count + 1
        l_out_const[[n]] <- list()
        const_beta <- matrix(sample(size = b_i*b_j, 
                                    x = c(-const[n], max = const[n]),
                                    replace = TRUE),
                                    nrow = b_i, ncol = b_j)
        tmp_beta <- m_beta + const_beta
        eigen_beta <- eigen(tmp_beta)$values
        eigen_cond <- all(Re(eigen_beta)^2 + Im(eigen_beta)^2 <1)
    
      } 
      l_out_const[[n]]$beta <- tmp_beta
      l_out_const[[n]]$constbeta <- const_beta
      l_out_const[[n]]$eigen_iter <- eigen_count
      
    }
    
    ## Kappa
    noise_mat <- list()
    for(n in seq(const)){
      # create change matrix for kappa, which is then scaled w.r.t. the sqrt of  diagonal elements
      const_mat <- matrix(sample(size = b_i*b_j, 
                                 x = c(-const[n], max = const[n]),
                                 replace = TRUE),
                          nrow = b_i, ncol = b_j)
      # scaling loop
      for(i in seq(k_i)){      # loop over rows
        for(j in seq(k_j)){    # loop over cols
          const_mat[i,j] <- const_mat[i,j]*sqrt(m_kappa[i,i]*m_kappa[j,j])
        }     
      }
      const_mat <- as.matrix(Matrix::forceSymmetric(const_mat))
      
      
      psd_check <- FALSE
      check_count_const <- 0
      while(psd_check == FALSE & check_count_const <= 100){
        check_count_const <- check_count_const + 1
        tmp_kappa <- m_kappa + const_mat
        eigen_kappa <- eigen(tmp_kappa)$values
        
        psd <- matrixcalc::is.positive.semi.definite(tmp_kappa) & all(eigen_kappa > 0)
        if(psd){
          psd_check <- TRUE
        } else{
          print(paste0("Adding const to Kappa leads to matrix that is not positive semi-definite for ", const[n], " attempt ", check_count_const))
          # Redraw matrix
          const_mat <- matrix(sample(size = b_i*b_j, 
                                     x = c(-const[n], max = const[n]),
                                     replace = TRUE),
                              nrow = b_i, ncol = b_j)
          # scaling loop
          for(i in seq(k_i)){      # loop over rows
            for(j in seq(k_j)){    # loop over cols
              const_mat[i,j] <- const_mat[i,j]*sqrt(m_kappa[i,i]*m_kappa[j,j])
            }     
          }
          const_mat <- as.matrix(Matrix::forceSymmetric(const_mat))
          
        } # end else
      } # end while
      
      
      l_out_const[[n]]$kappa <- tmp_kappa
      l_out_const[[n]]$args <- paste0("const", const[n])
      l_out_const[[n]]$constkappa <- const_mat
      names(l_out_const)[[n]] <- paste0("const", const[n])
    }
    
  } # end constant
  
  
  
  ### Permute matrix
  if(isTRUE(permute_active)){
    l_out_perm <- list()
    l_out_perm[[1]] <- list()
    l_out_perm[[1]]$beta <- permute_mat_col(m_beta, permute_index)
    l_out_perm[[1]]$kappa <- permute_mat_col(m_kappa, permute_index)  
    l_out_perm[[1]]$args <- "perm"
    names(l_out_perm) <- "perm"
    l_out <- c(l_out, l_out_change, l_out_noise, l_out_const, l_out_perm)
    
  }

  
  
  ### Output
  # Combine truegraph, maxchange, and added noise
  l_out <- c(l_out, l_out_change, l_out_noise, l_out_const)
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


# Fit graphicalVAR parallel -----------------------------------------------
fit_graphicalvar_parallel <- function(data,
                                      n, 
                                      pruneresults = TRUE, 
                                      seed,
                                      ...){     # other arguments passed to graphicalVAR
  # Save arguments
  
  require(doParallel)
  
  # reproducible parallelization
  doRNG::registerDoRNG(seed = seed)
  
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


# Fit VAR parallel merged -------------------------------------------------
# In simulation 1, we just use the default options. 

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





# Summarize posterior -----------------------------------------------------
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
  beta_mean_list <- vector("list", length(cred))
  pcor_lb_list <- vector("list", length(cred))
  pcor_ub_list <- vector("list", length(cred))
  pcor_mean_list <- vector("list", length(cred))
  
  # Compute bounds for each "cred" value
  for (i in seq_along(cred)) {
    
    # Lower and upper bound
    lb <- (1 - cred[i])/2
    ub <- 1 - lb
    
    # for beta
    beta_lb_list[[i]] <- apply(res$fit$beta, c(1, 2), stats::quantile, lb)
    beta_ub_list[[i]] <- apply(res$fit$beta, c(1, 2), stats::quantile, ub)
    beta_mean_list[[i]] <- apply(res$fit$beta, c(1, 2), mean)
    
    # for pcor
    pcor_lb_list[[i]] <- apply(res$fit$pcors, c(1, 2), stats::quantile, lb)
    pcor_ub_list[[i]] <- apply(res$fit$pcors, c(1, 2), stats::quantile, ub)
    pcor_mean_list[[i]] <- apply(res$fit$pcors, c(1, 2), mean)
    
  }
  names(beta_lb_list) <- paste0("lb_", cred)
  names(beta_ub_list) <- paste0("ub_", cred)
  names(beta_mean_list) <- paste0("mean_", cred)
  names(pcor_lb_list) <- paste0("lb_", cred)
  names(pcor_ub_list) <- paste0("ub_", cred)
  names(pcor_mean_list) <- paste0("mean_", cred)
  
  # Combine output into a list
  out <- list(beta_lb = beta_lb_list,
              beta_ub = beta_ub_list, 
              beta_mean = beta_mean_list,
              pcor_lb = pcor_lb_list, 
              pcor_ub = pcor_ub_list,
              pcor_mean = pcor_mean_list)
  
  return(out)
}





# -------------------------------------------------------------------------
# Comparison Functions ----------------------------------------------------
# -------------------------------------------------------------------------
# Within posterior comparison ---------------------------------------------
# This function uses bootstrapping to generate a difference distribution
# within the posterior

# will be replaced by the "compare_gvar" function in the package

within_compare <- function(
    fitpost_a = l_postres,
    fitpost_b = l_postres,
    fitemp_a = l_res,
    fitemp_b = l_res,
    mod_a = 1, 
    mod_b = 2,
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
  
  # Helper to select upper triangle elements of matrix
  ut <- function(x) {
    matrix(x[upper.tri(x, diag = FALSE)])
  }
  
  
  null_a <- post_distance_within(fitpost_a[[mod_a]], 
                                 comp = comparison, 
                                 draws = n_draws,
                                 pred = postpred,
                                 sampling_method = "random")
  null_b <- post_distance_within(fitpost_b[[mod_b]], 
                                 comp = comparison, 
                                 draws = n_draws,
                                 pred = postpred, 
                                 sampling_method = "random")
  
  
  
  # Compute empirical distance as test statistic
  if(comparison == "frob"){
    # Compute Distance of empirical betas between a and b
    emp_beta <- tryCatch(norm(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu, type = normtype), error = function(e) {NA})
    
    # Compute Distance of empirical pcors between a and b
    emp_pcor <- tryCatch(norm(ut(fitemp_a[[mod_a]]$pcor_mu) - ut(fitemp_b[[mod_b]]$pcor_mu), type = normtype), error = function(e) {NA})
    
    
  }
  
  if(comparison == "maxdiff"){
    # Compute maxdiff of empirical betas between a and b
    emp_beta <- tryCatch(max(abs(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu)), error = function(e) {NA})
    
    # Compute maxdiff of empirical pcors between a and b
    emp_pcor <- tryCatch(max(abs(ut(fitemp_a[[mod_a]]$pcor_mu) - ut(fitemp_b[[mod_b]]$pcor_mu))), error = function(e) {NA})
    
  }
  
  if(comparison == "l1"){
    # Compute l1 of empirical betas between a and b
    emp_beta <- tryCatch(sum(abs(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu)), error = function(e) {NA})
    
    # Compute l1 of empirical pcors between a and b
    emp_pcor <- tryCatch(sum(abs(ut(fitemp_a[[mod_a]]$pcor_mu) - ut(fitemp_b[[mod_b]]$pcor_mu))), error = function(e) {NA})
    
    
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
  
  
  
  
  
  l_cc_res <- list()
  l_cc_res[["beta"]] <- cc_res_beta
  l_cc_res[["pcor"]] <- cc_res_pcor
  
  cc_res <- dplyr::bind_rows(l_cc_res, .id = "mat")
  
  return(cc_res)
}


# -------------------------------------------------------------------------
# Between posterior comparison
# -------------------------------------------------------------------------
# Write this anew to be able to compare between models
# Why does this compare the beta_mu and pcor_mu? Why not the whole matrix?
# mod_a seemingly reflects a posterior sample here, not a model
post_distance_between <- function(fitobj_a,
                                  fitobj_b,
                                  comp,
                                  burnin = 50, 
                                  draws = 1000) {
  #--- Storage
  # browser()
  dist_out <- list()
  
  # Helper to select upper triangle elements of matrix
  ut <- function(x) {
    matrix(x[upper.tri(x, diag = FALSE)])
  }
  
  
  # for posteriors of empirical models
  # define the distance function based on comp
  # draw from all posterior samples
  
  distance_fn_beta <- switch(comp,
                             frob = {
                               function(x, y, mod_one, mod_two) norm(x$fit$beta[, , mod_one] - y$fit$beta[, , mod_two], type = "F")
                             },
                             maxdiff = {
                               function(x, y, mod_one, mod_two) max(abs((x$fit$beta[, , mod_one] - y$fit$beta[, , mod_two])))
                             },
                             l1 = {
                               function(x, y, mod_one, mod_two) sum(abs((x$fit$beta[, , mod_one] - y$fit$beta[, , mod_two])))
                             }
  )
  distance_fn_pcor <- switch(comp,
                             frob = {
                               function(x, y, mod_one, mod_two) norm(ut(x$fit$pcors[, , mod_one]) - ut(y$fit$pcors[, , mod_two]), type = "F")
                             },
                             maxdiff = {
                               function(x, y, mod_one, mod_two) max(abs((ut(x$fit$pcors[, , mod_one]) - ut(y$fit$pcors[, , mod_two]))))
                             },
                             l1 = {
                               function(x, y, mod_one, mod_two) sum(abs((ut(x$fit$pcors[, , mod_one]) - ut(y$fit$pcors[, , mod_two]))))
                             }
  )
  
  
  
  
  
  #--- Draw random samples from posterior
  # Obtain number of posterior samples
  # Assume that both models have the same number of samples
  n_post_samples <- dim(fitobj_a$fit$beta)[3]
  
  # Draw pairs of samples 
  sample_pairs <- replicate(draws, 
                            cbind(sample((burnin+1):n_post_samples, 2, 
                                         replace = TRUE)))
  
  # save the samples as an array
  dim(sample_pairs) <- c(2, draws)
  
  
  for (d in seq(draws)) {
    # storage
    dist_out[[d]] <- list()
    
    if (!is.list(fitobj_a) | !is.list(fitobj_b)) {
      beta_distance <- NA
      pcor_distance <- NA
      
      stop("Not a list.")
    }
    # if both elements are lists
    else {
      beta_distance <- distance_fn_beta(fitobj_a, fitobj_b, 
                                        sample_pairs[1, d], sample_pairs[2, d])
      pcor_distance <- distance_fn_pcor(fitobj_a, fitobj_b, 
                                        sample_pairs[1, d], sample_pairs[2, d])
    }
    
    
    # Store results
    dist_out[[d]]$comp <- comp
    dist_out[[d]]$beta <- beta_distance
    dist_out[[d]]$pcor <- pcor_distance
  } # end for loop
  out <- do.call(rbind.data.frame, dist_out)
  
  
  return(out)
}


# Between-compare with posterior distances ------------------------------------
# This function compares within-posterior distances to between-posterior distances
compare_gvar_between <- function(fit_a, 
                                 fit_b, 
                                 comp = "l1",
                                 n_draws = 1000,
                                 # Combine within-posterior uncertainty
                                 combine_post = TRUE,
                                 # sampling method for within-posterior comparison
                                 sampling_method = "random",
                                 burnin = 50
){
  
  #--- Obtain within-posterior uncertainty
  ref_a <- post_distance_within(fit_a,
                                comp = comp,
                                pred = FALSE,
                                draws = n_draws,
                                sampling_method = sampling_method
  )
  ref_b <- post_distance_within(fit_b,
                                comp = comp,
                                pred = FALSE,
                                draws = n_draws,
                                sampling_method = sampling_method
  )
  
  # Combine into single distribution
  if(isTRUE(combine_post)){
    ref <- data.frame(
      beta_ref = c(ref_a$beta[1:(n_draws/2)], ref_b$beta[1:(n_draws/2)]),
      pcor_ref = c(ref_a$pcor[1:(n_draws/2)], ref_b$pcor[1:(n_draws/2)])
    )
  } else
    ref <- data.frame(
      beta_ref_a = ref_a$beta,
      beta_ref_b =  ref_b$beta,
      pcor_ref_a = ref_a$pcor,
      pcor_ref_b = ref_b$pcor
    )
  
  
  
  #--- Obtain between posterior distances
  dist_between <- post_distance_between(fit_a, 
                                        fit_b, 
                                        burnin = burnin, 
                                        comp = comp, 
                                        draws = n_draws)
  
  #--- Return Results
  # Combine into single data frame
  out <- cbind(ref, dist_between)
  out <- out[,!colnames(out) %in% "comp"]
  
  # class(out) <- c("compare_gvar", class(out))
  
  return(out)
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


# TEMPORARY FROM TSNET
# ONLY USES UPPER.TRI for PCOR COMP
post_distance_within <- function(fitobj,
                                 comp,
                                 pred, # posterior predictive?
                                 draws = 1000,
                                 sampling_method = "random") {
  # storage
  dist_out <- list()
  
  # Helper to select upper triangle elements of matrix
  ut <- function(x) {
    matrix(x[upper.tri(x, diag = FALSE)])
  }
  
  # for posterior predictive approach
  if (isTRUE(pred)) {
    # define the distance function based on comp
    distance_fn_beta <- switch(comp,
                               frob = {
                                 function(x, y, mod_one, mod_two) norm(x$fit[[mod_one]]$beta_mu - y$fit[[mod_two]]$beta_mu, type = "F")
                               },
                               maxdiff = {
                                 function(x, y, mod_one, mod_two) max(abs((x$fit[[mod_one]]$beta_mu - y$fit[[mod_two]]$beta_mu)))
                               },
                               l1 = {
                                 function(x, y, mod_one, mod_two) sum(abs((x$fit[[mod_one]]$beta_mu - y$fit[[mod_two]]$beta_mu)))
                               }
    )
    distance_fn_pcor <- switch(comp,
                               frob = {
                                 function(x, y, mod_one, mod_two) norm(ut(x$fit[[mod_one]]$pcor_mu) - ut(y$fit[[mod_two]]$pcor_mu), type = "F")
                               },
                               maxdiff = {
                                 function(x, y, mod_one, mod_two) max(abs((ut(x$fit[[mod_one]]$pcor_mu) - ut(y$fit[[mod_two]]$pcor_mu))))
                               },
                               l1 = {
                                 function(x, y, mod_one, mod_two) sum(abs((ut(x$fit[[mod_one]]$pcor_mu) - ut(y$fit[[mod_two]]$pcor_mu))))
                               }
    )
    
    # Obtain number of models
    n_mod <- length(fitobj$fit)
  }
  
  
  # for posteriors of empirical models
  if (isFALSE(pred)) {
    # define the distance function based on comp
    # draw from all posterior samples
    
    distance_fn_beta <- switch(comp,
                               frob = {
                                 function(x, y, mod_one, mod_two) norm(x$fit$beta[, , mod_one] - y$fit$beta[, , mod_two], type = "F")
                               },
                               maxdiff = {
                                 function(x, y, mod_one, mod_two) max(abs((x$fit$beta[, , mod_one] - y$fit$beta[, , mod_two])))
                               },
                               l1 = {
                                 function(x, y, mod_one, mod_two) sum(abs((x$fit$beta[, , mod_one] - y$fit$beta[, , mod_two])))
                               }
    )
    distance_fn_pcor <- switch(comp,
                               frob = {
                                 function(x, y, mod_one, mod_two) norm(ut(x$fit$pcors[, , mod_one]) - ut(y$fit$pcors[, , mod_two]), type = "F")
                               },
                               maxdiff = {
                                 function(x, y, mod_one, mod_two) max(abs((ut(x$fit$pcors[, , mod_one]) - ut(y$fit$pcors[, , mod_two]))))
                               },
                               l1 = {
                                 function(x, y, mod_one, mod_two) sum(abs((ut(x$fit$pcors[, , mod_one]) - ut(y$fit$pcors[, , mod_two]))))
                               }
    )
    
    # Obtain number of posterior samples
    n_mod <- dim(fitobj$fit$beta)[3]
  }
  
  
  ## Draw two random models
  # delete burn-in iterations (to 50)
  # n_mod <- n_mod[-c(1:50)]
  
  # "samples" would be more fitting here than "models"
  # "model" is still a residue from posterior predictive approach
  # Draw models spaced apart so that we don't have autocorrelation from sampling
  mod_pairs <- array(NA, dim = c(2, draws))
  
  
  if(sampling_method == "sequential"){
    # draw from first half of samples
    mod_pairs[1, 1:(draws)] <- seq(51, draws + 50, by = 1)
    
    # draw from second half of samples
    mod_pairs[2, 1:(draws)] <- seq((n_mod / 2) + 51, (n_mod / 2) + 50 + (draws), by = 1)
    
  }
  if(sampling_method == "random"){
    # Determine the valid range of values for drawing pairs
    # Leave out middle 100 to keep certain distance between halves
    valid_range <- setdiff(1:n_mod, ((n_mod/2 - 50):(n_mod/2 + 49)))
    
    # Draw pairs randomly from the first half and second half
    mod_pairs[1, 1:draws] <- sample(valid_range[1:(length(valid_range)/2)], draws)
    mod_pairs[2, 1:draws] <- sample(valid_range[(length(valid_range)/2 + 1):length(valid_range)], draws)
    
    
  }
  
  
  
  for (i in seq(draws)) {
    # storage
    dist_out[[i]] <- list()
    mod_one <- mod_pairs[1, i]
    mod_two <- mod_pairs[2, i]
    
    # if mod_one and mod_two are equal, redraw
    if (mod_one == mod_two) {
      mod_two <- sample(1:n_mod, size = 1)
    }
    
    ## Check if estimation worked
    # Should be unneccessary if non-converged attempts were deleted
    if (isTRUE(pred)) {
      if (!is.list(fitobj$fit[[mod_one]]) | !is.list(fitobj$fit[[mod_two]])) {
        beta_distance <- NA
        pcor_distance <- NA
        stop("Not a list.")
      }
      # if both elements are lists
      else {
        beta_distance <- distance_fn_beta(fitobj, fitobj, mod_one, mod_two)
        pcor_distance <- distance_fn_pcor(fitobj, fitobj, mod_one, mod_two)
      }
    }
    
    if (isFALSE(pred)) {
      if (!is.list(fitobj) | !is.list(fitobj)) {
        beta_distance <- NA
        pcor_distance <- NA
        
        stop("Not a list.")
      }
      # if both elements are lists
      else {
        beta_distance <- distance_fn_beta(fitobj, fitobj, mod_one, mod_two)
        pcor_distance <- distance_fn_pcor(fitobj, fitobj, mod_one, mod_two)
      }
    }
    
    # Store results
    dist_out[[i]]$comp <- comp
    dist_out[[i]]$mod_one <- mod_one
    dist_out[[i]]$mod_two <- mod_two
    dist_out[[i]]$beta <- beta_distance
    dist_out[[i]]$pcor <- pcor_distance
  } # end for loop
  out <- do.call(rbind.data.frame, dist_out)
  
  
  return(out)
}



# -------------------------------------------------------------------------
# Comparison Evaluation ---------------------------------------------------
# -------------------------------------------------------------------------
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

    # Number of posterior difference > empirical difference
    teststat_a_pcor <- sum(df_res_pcor$null[df_res_pcor$mod == "mod_a"] > df_res_pcor$emp[df_res_pcor$mod == "mod_a"], na.rm = TRUE)
    teststat_b_pcor <- sum(df_res_pcor$null[df_res_pcor$mod == "mod_b"] > df_res_pcor$emp[df_res_pcor$mod == "mod_b"], na.rm = TRUE)
    
  }



  wcompres <- list(beta_a = teststat_a_beta,
                   beta_b = teststat_b_beta,
                   pcor_a = teststat_a_pcor,
                   pcor_b = teststat_b_pcor,
                   mod_a = model_ind_a, 
                   mod_b = model_ind_b,
                   comp = df_res_beta$comp[[1]]) # get type of comparison

  wcompres
}




# Within Compare Eval Revision --------------------------------------------
# Change: Allow combination of posteriors 
within_compare_eval_rev <- function(l_res,
                                pcor = TRUE){
  
  # if input is not a list, i.e. did not converge
  if(!is.list(l_res)){
    wcompres <- list(beta = NA,
                     pcor = NA,
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
  teststat_beta <- sum(df_res_beta$null > df_res_beta$emp, na.rm = TRUE)
  
  
  if(isTRUE(pcor)){
    ### Pcor
    df_res_pcor <- subset(df_res, mat == "pcor")
    
    # Number of posterior difference > empirical difference
    teststat_pcor <- sum(df_res_pcor$null > df_res_pcor$emp, na.rm = TRUE)
    
  }
  
  
  
  wcompres <- list(beta = teststat_beta,
                   pcor = teststat_pcor,
                   mod_a = model_ind_a, 
                   mod_b = model_ind_b,
                   comp = df_res_beta$comp[[1]]) # get type of comparison
  
  wcompres
}



# Compare to DGP ----------------------------------------------------------
# This function compares the BGGM VAR estimates to the true data-generating process
# Not used for manuscript, can just be helpful for some checks
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







# Evaluate BGGM Simulation ------------------------------------------------
bias <- function(e,t){
  b <- mean(abs(e - t), na.rm = TRUE)
  return(b)
}

# for calculation of mcse of bias
bias_sq <- function(e,t){
  b <- mean((e - t)^2, na.rm = TRUE)
  return(b)
}

rmse <- function(e,t){
  r <- sqrt(mean((t - e)^2, na.rm = TRUE))
  return(r)
}


eval_bggm <- function(fit,
                      cred_int = c(0.9, 0.95, 0.99),    # different credible intervals
                      nds = 1000,         # number of datasets per simulation condition
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
  l_out$bias_sq_beta <- bias_sq(beta_est_vec, beta_true_vec)
  l_out$bias_pcor <- bias(pcor_est_vec, pcor_true_vec)
  l_out$bias_sq_pcor <- bias_sq(pcor_est_vec, pcor_true_vec)
  
  # RMSE
  l_out$rmse_beta <- rmse(beta_est_vec, beta_true_vec)
  l_out$rmse_pcor <- rmse(pcor_est_vec, pcor_true_vec)
  
  # Correlations
  l_out$cor_beta <- cor_zero(beta_est_vec, beta_true_vec)
  l_out$cor_pcor <- cor_zero(pcor_est_vec, pcor_true_vec)
  
  
  
  #--- Select Method ---#
  # Obtain different credible intervals
  cred_interval <- fit$cred_interval
  
  
  # Obtain estimates with selection
  adj_beta <- fit$beta_weighted_adj
  adj_pcor <- fit$pcor_weighted_adj
  
  # Vectors for correlations
  adj_beta <- c(adj_beta)
  adj_pcor <- c(adj_pcor[upper.tri(adj_pcor, diag = FALSE)])
  
  # Correlations
  l_out$cor_sel_beta <- cor_zero(adj_beta, beta_true_vec)
  l_out$cor_sel_pcor <- cor_zero(adj_pcor, pcor_true_vec)
  
  # Bias
  l_out$bias_sel_beta <- bias(adj_beta, beta_true_vec)
  l_out$bias_sq_sel_beta <- bias_sq(adj_beta, beta_true_vec)
  l_out$bias_sel_pcor <- bias(adj_pcor, pcor_true_vec)
  l_out$bias_sq_sel_pcor <- bias_sq(adj_pcor, pcor_true_vec)
  
  # rmse
  l_out$rmse_sel_beta <- rmse(adj_beta, beta_true_vec)
  l_out$rmse_sel_pcor <- rmse(adj_pcor, pcor_true_vec)
  
  # Amount of zeros
  l_out$zeros_sel_beta <- sum(adj_beta == 0)
  l_out$zeros_sel_pcor <- sum(adj_pcor == 0)
  
  
  ## True/False Positive/Negative
  # TP
  l_out$true_pos_beta <- sum(beta_true_vec != 0 & adj_beta != 0)
  l_out$true_pos_pcor <- sum(pcor_true_vec != 0 & adj_pcor != 0)
  
  # FP
  l_out$fal_pos_beta <- sum(beta_true_vec == 0 & adj_beta != 0)
  l_out$fal_pos_pcor <- sum(pcor_true_vec == 0 & adj_pcor != 0)  
  
  # TN
  l_out$true_neg_beta <- sum(beta_true_vec == 0 & adj_beta == 0)
  l_out$true_neg_pcor <- sum(pcor_true_vec == 0 & adj_pcor == 0)
  
  # FN
  l_out$fal_neg_beta <- sum(beta_true_vec != 0 & adj_beta == 0)
  l_out$fal_neg_pcor <- sum(pcor_true_vec != 0 & adj_pcor == 0)
  
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
                      width_pcor = rep(NA, length(cred_int)),
                      true_pos_beta = rep(NA, length(cred_int)),
                      true_pos_pcor = rep(NA, length(cred_int)),
                      fal_pos_beta = rep(NA, length(cred_int)),
                      fal_pos_pcor = rep(NA, length(cred_int)),
                      true_neg_beta = rep(NA, length(cred_int)),
                      true_neg_pcor = rep(NA, length(cred_int)),
                      fal_neg_beta = rep(NA, length(cred_int)),
                      fal_neg_pcor = rep(NA, length(cred_int)),
                      sens_beta = rep(NA, length(cred_int)),
                      sens_pcor = rep(NA, length(cred_int)),
                      spec_beta = rep(NA, length(cred_int)),
                      spec_pcor = rep(NA, length(cred_int)))
  for(i in 1:length(cred_int)){
    lb_beta <- c(cred_interval$beta_lb[[i]])
    ub_beta <- c(cred_interval$beta_ub[[i]])
    lb_pcor <- cred_interval$pcor_lb[[i]]
    lb_pcor <- c(lb_pcor[upper.tri(lb_pcor)])
    ub_pcor <- cred_interval$pcor_ub[[i]]  
    ub_pcor <- c(ub_pcor[upper.tri(ub_pcor)])
    
    m_cover_beta <- beta_true_vec >= lb_beta & beta_true_vec <= ub_beta 
    m_cover_pcor <- pcor_true_vec >= lb_pcor & pcor_true_vec <= ub_pcor 
    
    # Only consider upper diagonal of pcor
    df_ci[i, "sum_cover_beta"] <- sum(m_cover_beta)
    df_ci[i, "sum_cover_pcor"] <- sum(m_cover_pcor)
    
    # Get adjacency matrix
    adj_beta <- lb_beta > 0 | ub_beta < 0
    adj_pcor <- lb_pcor > 0 | ub_pcor < 0
    
    # Compute stats
    ## True/False Positive/Negative
    # TP
    df_ci[i, "true_pos_beta"] <- sum(beta_true_vec != 0 & adj_beta != 0)
    df_ci[i, "true_pos_pcor"] <- sum(pcor_true_vec != 0 & adj_pcor != 0)
    
    # FP
    df_ci[i, "fal_pos_beta"] <- sum(beta_true_vec == 0 & adj_beta != 0)
    df_ci[i, "fal_pos_pcor"] <- sum(pcor_true_vec == 0 & adj_pcor != 0)  
    
    # TN
    df_ci[i, "true_neg_beta"] <- sum(beta_true_vec == 0 & adj_beta == 0)
    df_ci[i, "true_neg_pcor"] <- sum(pcor_true_vec == 0 & adj_pcor == 0)
    
    # FN
    df_ci[i, "fal_neg_beta"] <- sum(beta_true_vec != 0 & adj_beta == 0)
    df_ci[i, "fal_neg_pcor"] <- sum(pcor_true_vec != 0 & adj_pcor == 0)
    
    ## Sensitivity
    df_ci[i, "sens_beta"] <- df_ci[i, "true_pos_beta"] / (df_ci[i, "true_pos_beta"] + df_ci[i, "fal_neg_beta"])
    df_ci[i, "sens_pcor"] <- df_ci[i, "true_pos_pcor"] / (df_ci[i, "true_pos_pcor"] + df_ci[i, "fal_neg_pcor"])
    
    ## Specificity
    df_ci[i, "spec_beta"] <- df_ci[i, "true_neg_beta"] / (df_ci[i, "true_neg_beta"] + df_ci[i, "fal_pos_beta"])
    df_ci[i, "spec_pcor"] <- df_ci[i, "true_neg_pcor"] / (df_ci[i, "true_neg_pcor"] + df_ci[i, "fal_pos_pcor"])
    
    
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
  l_out$bias_sq_beta <- bias_sq(beta_est_vec, beta_true_vec)
  l_out$bias_pcor <- bias(pcor_est_vec, pcor_true_vec)
  l_out$bias_sq_pcor <- bias_sq(pcor_est_vec, pcor_true_vec)
  
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




# Or, from the tsnet package directly: 

#' Compute Effective Sample Sizes for MCMC Samples of BGGM GVAR models
#'
#' This function computes the effective sample size (ESS) of MCMC samples of
#' temporal and contemporaneous network parameters based on the provided
#' [BGGM::var_estimate()] fit object. It uses the default functionality of the
#' `coda` package to compute ESS. Effective sample size estimates for
#' [stan_gvar()] models can be found in their model output itself. Currently
#' used for internal development only.
#'
#' @param fitobj A [BGGM::var_estimate()] fit object.
#' @param burnin An integer indicating the number of burn-in iterations to
#'   discard. Default is 0.
#'
#' @return A list with two elements: ess_beta and ess_pcor. ess_beta contains
#'   the ESS of MCMC samples of the temporal coefficients, and ess_pcor contains
#'   the ESS of MCMC samples of partial correlation coefficients.
#'
#' @examples
#' \dontrun{
#' # Load data
#' data(ts_data)
#' example_data <- ts_data[1:100,]
#'
#' # Estimate a GVAR model
#' fit <- BGGM::var_estimate(example_data)
#' ess <- ess_gvar(fit)
#' }
#'
#' @importFrom coda as.mcmc effectiveSize
#' @noRd

ess_gvar <- function(fitobj,
                     burnin = 0) {


  # Input Information
  it <- fitobj$iter
  p <- fitobj$p

  ## Get samples
  beta <- fitobj$fit$beta[, , (burnin + 1):(it + burnin)]
  pcor <- fitobj$fit$pcors[, , (burnin + 1):(it + burnin)]

  # Transform to mcmc objects
  mcmc_beta <- coda::as.mcmc(t(matrix(beta, p * p, it)))
  mcmc_pcor <- coda::as.mcmc(t(matrix(pcor, p * p, it)))

  # correct variable names
  # column after column
  cnames <- colnames(fitobj$Y)
  cnames_lag <- paste0(colnames(fitobj$Y), ".l1")

  beta_names <- as.vector(vapply(cnames,
                                  function(x) paste(x, cnames_lag, sep = "--"),
                                  character(length(cnames))))
  pcor_names <- as.vector(vapply(cnames,
                                 function(x) paste(x, cnames, sep = "--"),
                                 character(length(cnames))))

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
compare_var <- function(fit_a, 
                        fit_b, 
                        cutoff = 5,           # percentage level of test
                        dec_rule = "OR",
                        n_draws = 1000,
                        comp = "frob",
                        return_all = FALSE){  # return all distributions?
  
  require(magrittr)
  
  # Store arguments
  args <- list(match.call)
  
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
                  larger_pcor = larger_pcor,
                  args = args)
    
  }
  if(isTRUE(return_all)){
    l_res <- list(sig_beta = sig_beta,
                  sig_pcor = sig_pcor,
                  res_beta = res_beta,
                  res_pcor = res_pcor,
                  emp_beta = emp_beta,
                  emp_pcor = emp_pcor,
                  larger_beta = larger_beta,
                  larger_pcor = larger_pcor,
                  args = args)
    
  }
  
  
  
  return(l_res)
  
  
  
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
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.15)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1.1)),
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, b = 10)),
      
      # Faceting
      strip.text = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.1), hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = NA, color = NA),
      # Grid
      panel.grid = ggplot2::element_line(colour = "#F3F4F5"),
      # Legend
      legend.title = ggplot2::element_text(face = "plain"),
      legend.position = "top",
      legend.justification = 1,
      # Panel/Facets
      panel.spacing.y = ggplot2::unit(1.5, "lines")
    )
}



