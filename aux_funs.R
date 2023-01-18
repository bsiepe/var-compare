# Find ToDos --------------------------------------------------------------
# todor::todor(c("TODO"))



# Change Kappa ------------------------------------------------------------
# Changes maximum value of kappa
changemax_kappa <- function(m_kappa,
                            change){
  # Find the dimensions of the m_kappa
  n <- nrow(m_kappa)
  m <- ncol(m_kappa)
  
  # Initialize the largest element to be zero
  largest_element <- 0
  row <- NA
  col <- NA
  
  # Loop through the m_kappa and find the largest off-diagonal element
  for (i in 1:n) {
    for (j in 1:m) {
      # Skip the diagonal elements (i.e. the elements where i == j)
      if (i != j) {
        # Update the largest element if we find a larger one
        if (abs(m_kappa[i,j]) > abs(largest_element)) {
          largest_element <- m_kappa[i,j]
          row <- i
          col <- j
        }
      }
    }
  }
  print(largest_element)
  # change the element in the matrix on both sides of diagonal
  m_kappa[row,col] <- m_kappa[row,col]*change
  m_kappa[col,row] <- m_kappa[col,row]*change
  
  
  
  # Multiply the largest element by a factor and return the result
  return(m_kappa)
}


# Change Graph ------------------------------------------------------------
# Function to change "true graph" according to specified changes of max value
# and noise
# Uses uniform noise for now
# TODO make sure that everything is stationary/below 1

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
change_graphs <- function(truegraph, 
                          changemax, 
                          noise){
  ### Create storage
  l_out <- list()
  l_out[[1]] <- truegraph
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
  l_out_change <- list()
  ## Beta
  for(c in seq_along(changemax)){
    # create storage for each changed graph
    l_out_change[[c]] <- list()
    tmp_beta <- m_beta
    tmp_beta[b_max] <- m_beta[b_max]*changemax[c]
    l_out_change[[c]]$beta <- tmp_beta
  }
  
  ## Kappa
  # change the element in the matrix on both sides of diagonal
  for(c in seq_along(changemax)){
    tmp_kappa <- m_kappa
    tmp_kappa[max_row,max_col] <- m_kappa[max_row,max_col]*changemax[c]
    tmp_kappa[max_col,max_row] <- m_kappa[max_col,max_row]*changemax[c]
    l_out_change[[c]]$kappa <- tmp_kappa
    names(l_out_change)[[c]] <- paste0("change_", changemax[c])
  }
  
  
  
  ### Add noise
  l_out_noise <- list()
  ## Beta
  for(n in seq_along(noise)){
    l_out_noise[[n]] <- list()
    tmp_beta <- m_beta + runif(n = b_i*b_j, min = -noise[n], max = noise[n])
    l_out_noise[[n]]$beta <- tmp_beta
  }
  
  ## Kappa
  for(n in seq_along(noise)){
    tmp_kappa <- m_kappa + runif(n = k_i*k_j, min = -noise[n], max = noise[n])
    tmp_kappa <- as.matrix(Matrix::forceSymmetric(tmp_kappa))
    l_out_noise[[n]]$kappa <- tmp_kappa
    names(l_out_noise)[[n]] <- paste0("noise_", noise[n])
  }
  
  
  ### Checking 
  ## Check if max change only changed one edge
  
  ## Check if stationarity holds
  
  ## Check if everything is <= 1
  l_out_change <- lapply(l_out_change, function(x){
    lapply(x, function(mat) {
      mat[mat > 1] <- 1
      mat
    })
  })


  l_out_noise <- lapply(l_out_noise, function(x){
    lapply(x, function(mat) {
      mat[mat > 1] <- 1
      mat
    })
  })
  
  
  
  ### Output
  # Combine truegraph, maxchange, and added noise
  l_out <- c(l_out, l_out_change, l_out_noise)
  return(l_out)
  
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
  args$dgp <- deparse(substitute(dgp))
  
  # ncores = parallel::detectCores() - 2
  # cl = makeCluster(ncores)
  # registerDoParallel(cl)
  
  # Reproducible loops
  registerDoRNG(seed)
  
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
# TODO needs different options for different provided data formats
# because $data does not work for posterior samples data
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
    # TODO does this call the correct BGGM version?
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
                                     pruneresults = FALSE){
  
  if(n != length(data)){
    warning("The n provided does not match the number of available data frames")
  }
  require(doParallel)

  # reproducible parallelization
  registerDoRNG(seed)
  

  # if we take simulated data in a list object
  # We sequence along number of posterior datasets and then cut away redundant datasets afterwards
  if(isFALSE(posteriorsamples)){
    fit <- foreach(i = seq(nds), .packages = "BGGM") %dopar% {
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
    # Cut away nonconverged attempts
    fit <- fit[!sapply(fit, is.null)]
    
    # Return list with desired length
    fit <- fit[c(1:n)]
    return(fit)
    
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

  
}  













# Sim from posterior ------------------------------------------------------
# TODO save function arguments
# TODO INSANELY IMPORTANT! DO I NEED TO TRANSPOSE? I THINK YES

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
                                   seed = seed,
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
      smp <- sample(iterations, size = 1)
      dat[[j]] <- try(as.data.frame(graphicalVAR::graphicalVARsim(nTime = tp,
                                                                  beta = t(l_params[[i]]$beta[,,smp]),
                                                                  kappa = l_params[[i]]$kappa[,,smp],
                                                                  mean = means))) 
      
    }
    
    dat
    
    
    
  } # end parallel
  
  # DEPRECATED: Should data be converted to bggm format?
  if(isTRUE(convert_bggm)){
    post_data <- lapply(post_data, format_bggm_list)
    
  }
  post_data
  
}
# Predict function for external data --------------------------------------
# Taken from https://github.com/donaldRwilliams/BGGM/blob/master/R/predict.estimate.R
# adapted for using external data for refitting of the model
# data needs to be of the same length, which is no issue when using simulated data

predict.var_estimate <- function(object,
                                 data,
                                 summary = TRUE,
                                 cred = 0.95,
                                 iter = NULL,
                                 progress = TRUE,
                                 ...){
  
  
  # lower bound
  lb <- (1 - cred) / 2
  
  # uppder bound
  ub <- 1 - lb
  
  
  # data matrix
  # B: changed by me to accomodate external object
  # needs to be in the proper (lagged) formatting
  X <- data$X
  n <- nrow(X)
  
  
  if(is.null(iter)){
    
    iter <- object$iter
    
  }
  
  p <- object$p
  
  post_names <- sapply(1:p, function(x) paste0(
    colnames(data$Y)[x], "_",  colnames(data$X))
  )
  
  post_samps <- BGGM::posterior_samples(object)
  
  if(isTRUE(progress)){
    pb <- utils::txtProgressBar(min = 0, max = p, style = 3)
  }
  
  yhats <- lapply(1:p, function(x){
    
    yhat_p <- post_samps[, post_names[,x]] %*% t(X)
    
    
    if(isTRUE(progress)){
      utils::setTxtProgressBar(pb, x)
    }
    
    yhat_p
    
  })
  
  if(isTRUE(summary)){
    
    
    
    fitted_array <- array(0, dim = c(n, 4, p))
    
    dimnames(fitted_array)[[2]] <- c("Post.mean",
                                     "Post.sd",
                                     "Cred.lb",
                                     "Cred.ub")
    
    dimnames(fitted_array)[[3]] <- colnames(data$Y)
    
    
    for(i in 1:p){
      
      fitted_array[,,i] <- cbind(colMeans(yhats[[i]]),
                                 apply(yhats[[i]], 2, sd),
                                 apply(yhats[[i]], 2, quantile, lb),
                                 apply(yhats[[i]], 2, quantile, ub)
      )
    }
    
    
  } else {
    
    fitted_array <- array(0, dim = c(iter, n, p))
    
    dimnames(fitted_array)[[3]] <- colnames(data$Y)
    
    for(i in 1:p){
      
      fitted_array[,,i] <- t(as.matrix(yhats[[i]]))
      
    }
    
  }
  
  return(fitted_array)
  
}





# Predict with posterior mean ---------------------------------------------
# IDEA: add option to manually provide a beta matrix with posterior means
# instead of building the mean of the posterior samples

# Change the function above to only use posterior mean for prediction
predict_pmu.var_estimate <- function(object,
                                 data,
                                 summary = TRUE,
                                 cred = 0.95,
                                 iter = NULL,
                                 progress = TRUE,
                                 ...){
  
  
  # lower bound
  lb <- (1 - cred) / 2
  
  # uppder bound
  ub <- 1 - lb
  
  
  # data matrix
  # B: changed by me to accomodate external object
  # needs to be in the proper (lagged) formatting
  X <- data$X
  n <- nrow(X)
  
  
  if(is.null(iter)){
    
    iter <- object$iter
    
  }
  
  p <- object$p
  
  post_names <- sapply(1:p, function(x) paste0(
    colnames(data$Y)[x], "_",  colnames(data$X))
  )
  
  post_samps <- BGGM::posterior_samples(object)
  
  # Only use relevant Beta samples
  post_samps_b <- post_samps[,post_names]
  colnames(post_samps_b) <- post_names
  post_samps_mu <- matrix(colMeans(post_samps_b), ncol = 6)
  # Very important: Do I need to do byrow TRUE or not?
  # I don't think so, because it is done variable-wise
  # so the first 6 values are the the coefficients on V1
  
  
  if(isTRUE(progress)){
    pb <- utils::txtProgressBar(min = 0, max = p, style = 3)
  }
  
  yhats <- lapply(1:p, function(x){
    
    yhat_p <- post_samps_mu[,x] %*% t(X)
    
    
    if(isTRUE(progress)){
      utils::setTxtProgressBar(pb, x)
    }
    
    yhat_p
    
  })
  
  if(isTRUE(summary)){
    
    
    
    fitted_array <- array(0, dim = c(n, 1, p))
    
    dimnames(fitted_array)[[2]] <- c("Post.mean")
    
    dimnames(fitted_array)[[3]] <- colnames(data$Y)
    
    
    for(i in 1:p){
      
      fitted_array[,,i] <- cbind(colMeans(yhats[[i]])
      )
    }
    
    
  } else {
    
    fitted_array <- array(0, dim = c(iter, n, p))
    
    dimnames(fitted_array)[[3]] <- colnames(data$Y)
    
    for(i in 1:p){
      
      fitted_array[,,i] <- t(as.matrix(yhats[[i]]))
      
    }
    
  }
  
  return(fitted_array)
  
}







# External data to BGGM format --------------------------------------------
# External data needs to have two objects
# Y: response data
# X: lagged response data

format_bggm <- function(Y){
  if(is.data.frame(Y)){
    Y <- scale(na.omit(Y))
    p <- ncol(Y)
    n <- nrow(Y)
    Y_lag <- rbind(NA, Y)
    colnames(Y_lag) <- paste0(colnames(Y), ".l1")
    Y_all <- na.omit(cbind.data.frame(rbind(Y, NA), Y_lag))
    Y <- as.matrix(Y_all[, 1:p])
    X <- as.matrix(Y_all[, (p + 1):(p * 2)])
    out <- list()
    out$Y <- Y
    out$X <- X
    
    
  }
  else out <- NA
  
  out
}



# External data from nested list to BGGM format ---------------------------
# same as format_bggm, but has nested list as input
format_bggm_list <- function(listname){
  l_out <- lapply(listname, format_bggm)
  l_out
}







# Plot error distribution -------------------------------------------------
plot_error_dist <- function(dat, errorcol = mse){
  ggplot(dat, aes(x = {{errorcol}}))+
  ggdist::stat_dist_halfeye(fill = ggokabeito::palette_okabe_ito()[2],
                            color = "black")+
  theme_minimal()
}






# Distance between empirical and posterior --------------------------------

# THIS IS A NEW BETA VERSION!
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


# Distance within posterior samples ---------------------------------------
# Looks at differences between models sampled from the same
# "original" model, so similar to bootstrapping
# TODO this is still in the works
# TODO make other postpost functions shorter as well
postpost_distance_within <- function(post_a, 
                      comp, 
                      draws = 1000){
  
  # storage
  dist_out <- list()
  
  # define the distance function based on comp
  distance_fn_beta <- switch(comp,
                             frob =   {function(x, y) norm(x$beta_mu-y$beta_mu, type = "F")},
                             maxdiff = {function(x, y) max(abs((x$beta_mu-y$beta_mu)))},
                             l1 = {function(x, y) sum(abs((x$beta_mu-y$beta_mu)))}
  )
  distance_fn_pcor <- switch(comp,
                             frob = {function(x, y) norm(x$pcor_mu-y$pcor_mu, type = "F")},
                             maxdiff = {function(x, y) max(abs((x$pcor_mu-y$pcor_mu)))},
                             l1 = {function(x, y) sum(abs((x$pcor_mu-y$pcor_mu)))}
  )
  
  
  
  
  ## Draw two random models
  # Obtain number of models
  n_mod <- length(post_a)
  
  # Draw pairs of models
  mod_pairs <- replicate(draws, sample(1:n_mod, size = 2, replace = TRUE))
  
for(i in seq(draws)){
  # storage
  dist_out[[i]] <- list()
  mod_a <- mod_pairs[1,i]
  mod_b <- mod_pairs[2,i]
  
# if mod_a and mod_b are equal, redraw
  if(mod_a == mod_b){
    mod_b <- sample(1:n_mod, size = 1)
}
  
  ## Check if estimation worked
  # Should be unneccessary if non-converged attempts were deleted
  if(!is.list(post_a[[mod_a]]) | !is.list(post_a[[mod_b]])){
    beta_distance <- NA
    pcor_distance <- NA
    stop("Not a list.")
    
    
  } 
  # if both elements are lists
  else{
    beta_distance <- distance_fn_beta(post_a[[mod_a]], post_a[[mod_b]])
    pcor_distance <- distance_fn_pcor(post_a[[mod_a]], post_a[[mod_b]])
    
  }  
  dist_out[[i]]$comp <- comp
  dist_out[[i]]$mod_a <- mod_a
  dist_out[[i]]$mod_b <- mod_b
  dist_out[[i]]$beta <- beta_distance
  dist_out[[i]]$pcor <- pcor_distance  
  
} # end for loop  
  out <- do.call(rbind, dist_out)
  
  
  
  return(out)
}



 
# # Cross-compare all posterior samples -------------------------------------
# TODO implement for all comparison types
# does not work yet, emp gives NA
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
    # Distance empirical to posterior sample estimates
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
                                null = c(null_a[["beta"]], null_b[["beta"]]),
                                emp = rep(emp_beta, n_datasets*2),
                                comp = rep(comparison, n_datasets*2),
                                type = rep("postemp", n_datasets*2))
      
      
      cc_res_pcor <- data.frame(model_ind = c(rep(mod_a, n_datasets), rep(mod_b, n_datasets)),
                                null = c(null_a[["pcor"]], null_b[["pcor"]]),
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
                                emp = rep(post[["beta"]], 2),    # TODO repeat the distribution?
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





# Evaluate Cross-Comparison with Matrix Norm ------------------------------
# TODO make this applicable to both postpost and postemp comparison types?

#' Evaluate Cross-Comparison with Matrix Norm
#'
#' @param df_res Dataframe with posterior predictive and empirical distance measures. 
#' @param plot Should results be plotted with densities? Defaults to FALSE.
#'
#' @return List with number of posterior predictive differences that were smaller than empirical difference. 
#' @export
#'
#' 

cross_compare_eval <- function(l_res, 
                               pcor = TRUE){
  
  ### Betas
  df_res_beta <- l_res[["beta"]]
  # Obtain model indexes
  model_ind_a <- unique(df_res_beta$model_ind)[1]
  model_ind_b <- unique(df_res_beta$model_ind)[2]
  
  # Number of posterior difference > empirical difference
  teststat_a_beta <- sum(df_res_beta$null[df_res_beta$model_ind == model_ind_a] > df_res_beta$emp[df_res_beta$model_ind == model_ind_a], na.rm = TRUE)
  teststat_b_beta <- sum(df_res_beta$null[df_res_beta$model_ind == model_ind_b] > df_res_beta$emp[df_res_beta$model_ind == model_ind_b], na.rm = TRUE)
  
    
    
  if(isTRUE(pcor)){
    ### Pcor
    df_res_pcor <- l_res[["pcor"]]
    # Obtain model indexes
    model_ind_a <- unique(df_res_pcor$model_ind)[1]
    model_ind_b <- unique(df_res_pcor$model_ind)[2]
    
    # Number of posterior difference > empirical difference
    teststat_a_pcor <- sum(df_res_pcor$null[df_res_pcor$model_ind == model_ind_a] > df_res_pcor$emp[df_res_pcor$model_ind == model_ind_a], na.rm = TRUE)
    teststat_b_pcor <- sum(df_res_pcor$null[df_res_pcor$model_ind == model_ind_b] > df_res_pcor$emp[df_res_pcor$model_ind == model_ind_b], na.rm = TRUE)
    
  }

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




























