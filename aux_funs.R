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



# Add noise kappa ---------------------------------------------------------
# # Add noise to kappa
# addnoise_kappa <- function(m_kappa,
#                             ){
#   # Find the dimensions of the m_kappa
#   n <- nrow(m_kappa)
#   m <- ncol(m_kappa)
#   
# 
#   # Loop through the m_kappa
#   for (i in 1:n) {
#     for (j in 1:m) {
#       # Skip the diagonal elements (i.e. the elements where i == j)
#       if (i != j) {
#         # Update the largest element if we find a larger one
#         if (abs(m_kappa[i,j]) > abs(largest_element)) {
#           largest_element <- m_kappa[i,j]
#           row <- i
#           col <- j
#         }
#       }
#     }
#   }
#   print(largest_element)
#   # change the element in the matrix on both sides of diagonal
#   m_kappa[row,col] <- m_kappa[row,col]*change
#   m_kappa[col,row] <- m_kappa[col,row]*change
#   
#   
#   
#   # Multiply the largest element by a factor and return the result
#   return(m_kappa)
# }



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
                             standardize = TRUE){
  # save function arguments in output
  args <- as.list(environment())
  args$dgp <- deparse(substitute(dgp))
  
  # ncores = parallel::detectCores() - 2
  # cl = makeCluster(ncores)
  # registerDoParallel(cl)
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
          }

        }
        # prune results for comparison purposes
        if(isTRUE(pruneresults) & is.list(fit_ind)){
          beta_mu <- fit_ind$beta_mu 
          kappa_mu <- fit_ind$kappa_mu
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
# TODO: Keep counter of failed attempts?

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
  require("BGGM", lib.loc = "C:/Users/Bjoern/R-dev")
  
  
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
    
    # storage for overall results
    
    
    # loop across individuals
    foreach(i = seq(n), .packages = "BGGM") %dopar% {
      m[[i]] <- 0
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
          
          # check if fitting worked
          if(is.list(fit_ind)){
            m[[i]] <- m[[i]]+1
            # print(m[[i]])
            
            
            if(isTRUE(get_kappa)){
              # Invert covariance matrix of residuals to obtain precision matrix
              fit_ind$fit$kappa <- array(apply(fit_ind$fit$Sigma, 3, solve), 
                                         dim = dim(fit_ind$fit$Sigma))
            }
            
            
            # prune results for comparison purposes
            if(isTRUE(pruneresults)){
              beta_mu <- fit_ind$beta_mu 
              kappa_mu <- fit_ind$kappa_mu
              pcor_mu <- fit_ind$pcor_mu
              fit_ind <- list()
              fit_ind$beta_mu <- beta_mu
              fit_ind$kappa_mu <- kappa_mu
              fit_ind$pcor_mu <- pcor_mu
              
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
sim_from_post_parallel <- function(fitobj, 
                                   n_datasets, 
                                   n,
                                   tp,
                                   iterations,
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








# Combine y and yhat ------------------------------------------------------
# Function to get yhat and actual data into the same dataframe
# data = data object
# refit = predict.var_estimate object
f_merge_pred <- function(data, refit){
  y_pred <- array(NA, dim = c(dim(refit)[1],dim(refit)[3]))
  # loop over number of variables in refit object
  for(p in 1:dim(refit)[3]){
    y_pred[,p] <- refit[,,p]
    dimnames(y_pred) <- list(NULL, paste0("V",1:6,"_pred"))
  }
  res <- cbind(data$Y, y_pred)
  return(res)
}




# Compute RMSE ------------------------------------------------------------
# everything should be standardized
f_rmse <- function(res){
  # loop through the merged results
  # 1st column with 7th, 2nd column with 8th etc.
  ncol <- dim(res)[2]/2
  res_rmse <- array(NA, dim = c(1, ncol))
  for(j in 1:ncol){
    res_rmse[,j] <- sqrt(sum((res[,j]-res[,j+ncol])^2))
  }
  colnames(res_rmse)[1:ncol] <- paste0(colnames(res)[1:ncol],"_rmse")
  return(res_rmse)
}





# Evaluate refit ----------------------------------------------------------
# function that combines everything
f_eval_refit <- function(data,
                         refit){
  comb <- f_merge_pred(data = data, refit = refit)
  fit_stat <- f_rmse(comb)
  mean_rmse <- rowMeans(fit_stat)
  
}





# Plot error distribution -------------------------------------------------
plot_error_dist <- function(dat, errorcol = mse){
  ggplot(dat, aes(x = {{errorcol}}))+
  ggdist::stat_dist_halfeye(fill = ggokabeito::palette_okabe_ito()[2],
                            color = "black")+
  theme_minimal()
}





# JSD between reference and other error distributions ---------------------
# ref_model = number of model for reference
# n = number of generated datasets
f_comp_jsd <- function(df = df_errors, ref_model = 1, n){
  # create storage
  l_err <- list()
  l_ecdf <- list()
  df_jsd <- data.frame(model = rep(NA, n),
                       jsd = rep(NA, n))
  
  
  # Obtain characteristics of reference model
  # obtain RMSEs
  tmp <- subset(df_errors, model == ref_model, select = rmse)
  rmse_refmod <- tmp$rmse
  
  # obtain ECDFs
  f_ecdf_ref <- stats::ecdf(rmse_refmod)
  ecdf_refmod <- f_ecdf_ref(rmse_refmod)
  
  # setup loop
  for(i in seq(n)){
    # if model has errors stored
    if(nrow(subset(df_errors, model == i, select = rmse)) > 0){
      # obtain RMSEs
      tmp <- subset(df_errors, model == i, select = rmse)
      rmse_mod <- tmp$rmse
      
      # obtain ECDFs
      f_ecdf <- stats::ecdf(rmse_mod)
      ecdf_mod <- f_ecdf(rmse_mod)
      
      # compute JSD to reference distribution
      v_ecdf <- rbind(ecdf_refmod, ecdf_mod)
      jsd <- philentropy::JSD(v_ecdf)
      
      # store values
      df_jsd[i,"model"] <- i
      df_jsd[i, "jsd"] <- jsd
    }

   
  }
  
  return(df_jsd)
}


# Compute likelihood ------------------------------------------------------
# https://github.com/cran/graphicalVAR/blob/master/R/graphicalVAR.R
# starting with line 158
# kappa = precision matrix of residuals
# sigma = unconstrained covariance matrix of residuals

logll <- function(kappa, sigma, n){
  lik1 <- determinant(kappa)$modulus[1]
  lik2 <- sum(diag(kappa %*% sigma))
  
  llk <- (n/2)*(lik1-lik2)
  llk
}





# Distance between empirical and posterior --------------------------------

# TODO output structure
postemp_distance <- function(post,
                             emp,
                             comp,
                             mod){
  
  # storage
  dist_out <- list()
    
  
  
  if(comp == "frob"){
    normtype = "F"
    frob_beta <- unlist(lapply(post[[mod]], function(x){
      if(length(x) == 0 | !is.list(x))
        {NA} 
      else
        {norm(emp[[mod]]$beta_mu-x$beta_mu, type = normtype)}}))
      
    frob_pcor <- unlist(lapply(post[[mod]], function(x){
      if(length(x) == 0 | !is.list(x))
        {NA} 
      else
        {norm(emp[[mod]]$pcor_mu-x$pcor_mu, type = normtype)}}))
    
    dist_out[["beta"]] <- frob_beta
    dist_out[["pcor"]] <- frob_pcor
    
    }  # end frob


  
  
  if(comp == "maxdiff"){
    maxdiff_beta <- unlist(lapply(post[[mod]], function(x){
      if(length(x) == 0 | !is.list(x))
        {NA}
      else
        {max(abs(emp[[mod]]$beta_mu-x$beta_mu))}}))
    
    maxdiff_pcor <- unlist(lapply(post[[mod]], function(x){
      if(length(x) == 0 | !is.list(x))
        {NA}
      else
        {max(abs(emp[[mod]]$pcor_mu-x$pcor_mu))}}))
    
    dist_out[["beta"]] <- maxdiff_beta
    dist_out[["pcor"]] <- maxdiff_pcor
    
    }  # end maxdiff
  
  
  
  if(comp == "l1"){
    l1_beta <-  unlist(lapply(post[[mod]], function(x){
      if(length(x) == 0 | !is.list(x))
        {NA}
      else
        {sum(abs(emp[[mod]]$beta_mu-x$beta_mu))}}))
    
    l1_pcor <-  unlist(lapply(post[[mod]], function(x){
      if(length(x) == 0 | !is.list(x))
      {NA}
      else
      {sum(abs(emp[[mod]]$pcor_mu-x$pcor_mu))}}))
    
    dist_out[["beta"]] <- l1_beta
    dist_out[["pcor"]] <- l1_pcor  
    
    
    }  # end l1
           
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






# Cross compare two models with posterior ---------------------------------
#' Cross-compare two models with their posterior
#'
#' @param mod_a Numerical indicator of model/person A. 
#' @param mod_b Numerical indicator of model/person B. 
#' @param n_datasets Number of datasets sampled from posterior. 
#' @param comparison Which comparison to use. "Frob" for Frobenius-Norm, "maxdiff" for maximum edge difference. 
#' @param ... Currently ignored. 
#' @param fitpost_a Posterior fit objects for model A. 
#' @param fitpost_b Posterior fit objects for model B.
#' @param fitemp_a Empirical fit object for model A.
#' @param fitemp_b Empirical fit object for model B. 
#'
#' @return Dataframe with null distributions for Frobenius norm of both models under the Null and empirical Frobenius norm between both models 
#' @export
#'
#' 
cross_compare_emp <- function(
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
  null_a <- postemp_distance(post = fitpost_a, emp = fitemp_a, comp = "frob", mod = mod_a)
  
  
  
  
  
  # Frobenius Norm Comparison 
  if(comparison == "frob"){
    normtype = "F"
    # Compute Distance of empirical estimates to posterior samples estimates
    frob_null_a <- postemp_distance(post = fitpost_a, emp = fitemp_a, comp = "frob", mod = mod_a)
    frob_null_b <- postemp_distance(post = fitpost_b, emp = fitemp_b, comp = "frob", mod = mod_b)
    
    # Compute Distance of empirical betas between a and b
    frob_emp_beta <- tryCatch(norm(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu), error = function(e) {NA})
    
    # Compute Distance of empirical pcors between a and b
    frob_emp_pcor <- tryCatch(norm(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu), error = function(e) {NA})
    
    cc_res_beta <- data.frame(model_ind = c(rep(mod_a, n_datasets), rep(mod_b, n_datasets)),
                         null = c(frob_null_a[["beta"]], frob_null_b[["beta"]]),
                         emp = rep(frob_emp_beta, n_datasets*2),
                         comp = rep("frob", n_datasets*2))
    
    
    cc_res_pcor <- data.frame(model_ind = c(rep(mod_a, n_datasets), rep(mod_b, n_datasets)),
                         null = c(frob_null_a[["pcor"]], frob_null_b[["pcor"]]),
                         emp = rep(frob_emp_pcor, n_datasets*2),
                         comp = rep("frob", n_datasets*2))
    
    
  }
  # Max Difference Comparison
  if(comparison == "maxdiff"){
    # Compute maximum distance of empirical estimates to posterior samples estimates
    maxdiff_null_a <- postemp_distance(post = fitpost_a, emp = fitemp_a, comp = "maxdiff", mod = mod_a)
    maxdiff_null_b <- postemp_distance(post = fitpost_b, emp = fitemp_b, comp = "maxdiff", mod = mod_b)
   
    # Compute maxdiff of empirical betas between a and b
    maxdiff_emp_beta <- tryCatch(max(abs(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu)), error = function(e) {NA})
    
    # Compute maxdiff of empirical pcors between a and b
    maxdiff_emp_pcor <- tryCatch(max(abs(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu)), error = function(e) {NA})
    
    cc_res_beta <- data.frame(model_ind = c(rep(mod_a, n_datasets), rep(mod_b, n_datasets)),
                         null = c(maxdiff_null_a[["beta"]], maxdiff_null_b[["beta"]]),
                         emp = rep(maxdiff_emp_beta, n_datasets*2),
                         comp = rep("maxdiff", n_datasets*2))
    

    
    cc_res_pcor <- data.frame(model_ind = c(rep(mod_a, n_datasets), rep(mod_b, n_datasets)),
                              null = c(maxdiff_null_a[["pcor"]], maxdiff_null_b[["pcor"]]),
                              emp = rep(maxdiff_emp_pcor, n_datasets*2),
                              comp = rep("maxdiff", n_datasets*2))
    
    
    
  } # end maxdiff
  
  
  if(comparison == "l1"){
    l1_null_a <- postemp_distance(post = fitpost_a, emp = fitemp_a, comp = "l1", mod = mod_a)
    l1_null_b <- postemp_distance(post = fitpost_b, emp = fitemp_b, comp = "l1", mod = mod_b)
    
    # Compute l1 of empirical betas between a and b
    l1_emp_beta <- tryCatch(sum(abs(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu)), error = function(e) {NA})
    
    # Compute l1 of empirical pcors between a and b
    l1_emp_pcor <- tryCatch(sum(abs(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu)), error = function(e) {NA})
    
    cc_res_beta <- data.frame(model_ind = c(rep(mod_a, n_datasets), rep(mod_b, n_datasets)),
                              null = c(l1_null_a[["beta"]], l1_null_b[["beta"]]),
                              emp = rep(l1_emp_beta, n_datasets*2),
                              comp = rep("l1", n_datasets*2))
    
    
    
    
    cc_res_pcor <- data.frame(model_ind = c(rep(mod_a, n_datasets), rep(mod_b, n_datasets)),
                              null = c(l1_null_a[["pcor"]], l1_null_b[["pcor"]]),
                              emp = rep(l1_emp_pcor, n_datasets*2),
                              comp = rep("l1", n_datasets*2))
    
  }
  
  l_cc_res <- list()
  l_cc_res[["beta"]] <- cc_res_beta
  l_cc_res[["pcor"]] <- cc_res_pcor
  
  return(l_cc_res)
} # end function








# 
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
    null_b <- postemp_distance(post = fitpost_a, emp = fitemp_a, 
                               comp = comparison, mod = mod_b)
    
    # Compute empirical distance as test statistic
    if(isFALSE(postpost)){
      
      if(comparison == "frob"){
        # Compute Distance of empirical betas between a and b
        emp_beta <- tryCatch(norm(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu), error = function(e) {NA})
        
        # Compute Distance of empirical pcors between a and b
        emp_pcor <- tryCatch(norm(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu), error = function(e) {NA})
        
      }

      if(comparison == "maxdiff"){
        # Compute maxdiff of empirical betas between a and b
        maxdiff_emp_beta <- tryCatch(max(abs(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu)), error = function(e) {NA})
        
        # Compute maxdiff of empirical pcors between a and b
        maxdiff_emp_pcor <- tryCatch(max(abs(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu)), error = function(e) {NA})
        
      }
      
      if(comparison == "l1"){
        # Compute l1 of empirical betas between a and b
        l1_emp_beta <- tryCatch(sum(abs(fitemp_a[[mod_a]]$beta_mu - fitemp_b[[mod_b]]$beta_mu)), error = function(e) {NA})
        
        # Compute l1 of empirical pcors between a and b
        l1_emp_pcor <- tryCatch(sum(abs(fitemp_a[[mod_a]]$pcor_mu - fitemp_b[[mod_b]]$pcor_mu)), error = function(e) {NA})
        
        
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
                  comp = df_res_beta$comp[[1]])  # get type of comparison
  testres
}






# Frobenius norm two datasets from each posterior sample ------------------
f_post_frob <- function(sample, seed = 2022,
                        rho_sd = 0.5, beta_sd = 1){
  set.seed = seed
  
  # Simulate
  d1 <- mlVAR::simulateVAR(pars = sample,
                           means = 0,
                           Nt = 200,
                           residuals = .1)
  d2 <- mlVAR::simulateVAR(pars = sample,
                           means = 0,
                           Nt = 200,
                           residuals = .1)
  # Estimate
  m1 <- try(BGGM::var_estimate(d1,
                               rho_sd = rho_sd,
                               beta_sd = beta_sd,
                               iter = n_iter,
                               progress = FALSE,
                               seed = seed))
  m2 <- try(BGGM::var_estimate(d2,
                               rho_sd = rho_sd,
                               beta_sd = beta_sd,
                               iter = n_iter,
                               progress = FALSE,
                               seed = seed))
  
  # Compute Frobenius Norm
  if(is.list(m1) & is.list(m2)){
    frob_norm <- norm(m1$beta_mu - m2$beta_mu, type = "F")
    return(frob_norm)
    
  }
  else return(NA)
  
  
}


# Compare to DGP ----------------------------------------------------------
# TODO: Do I need to takek absolute differences here or does it work this way?
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
  l_diff_beta[["true_bggm"]] <- map2(dgp, est_bggm, .f = function(x,y){x$beta-t(y$beta_mu)})
  l_diff_pcor[["true_bggm"]] <- map2(dgp, est_bggm, .f = function(x,y){x$PCC-y$pcor_mu})
  
  if(isTRUE(comp_gvar)){
    # Look at difference to true values for gvar
    # delete intercepts
    l_diff_beta[["true_gvar"]] <- map2(dgp, est_gvar, .f = function(x,y){x$beta-y$beta[,-1]})
    l_diff_pcor[["true_gvar"]] <- map2(dgp, est_gvar, .f = function(x,y){x$PCC-y$PCC})
    
    # Look at difference between bggm and gvar
    # delete intercepts
    l_diff_beta[["bggm_gvar"]] <- map2(est_bggm, est_gvar, .f = function(x,y){t(x$beta_mu)-y$beta[,-1]})
    l_diff_pcor[["bggm_gvar"]] <- map2(est_bggm, est_gvar, .f = function(x,y){x$pcor_mu-y$PCC})
    
  }
  
  
  # Aggregate differences
  out[["diff_beta"]] <- lapply(l_diff_beta, function(x){apply(simplify2array(x), 1:2, mean)})
  out[["diff_pcor"]] <- lapply(l_diff_pcor, function(x){apply(simplify2array(x), 1:2, mean)})
  
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














# Posterior sampmles covariance matrix ------------------------------------
# res = results of var_estimate
f_postcov <- function(res){
  beta_posterior <- res$fit$beta
  # delete warm-up samples
  beta_posterior <- beta_posterior[,,51:5050]

  # get mean and SD of posterior estimates
  beta_mu <- round(apply(beta_posterior,1:2,mean), digits = 3)
  beta_sd <- round(apply(beta_posterior,1:2,sd), digits = 3)

  # obtain the covariance matrix of estimates
  dimnames(beta_posterior)[[1]] <- c("V1.l1", "V2.l1", "V3.l1", "V4.l1", "V5.l1", "V6.l1")
  dimnames(beta_posterior)[[2]] <- c("V1", "V2", "V3", "V4", "V5", "V6")

  # convert array to list
  l_beta_posterior <- lapply(seq(dim(beta_posterior)[3]), function(x) beta_posterior[,,x])

  ldf_beta_posterior <- lapply(l_beta_posterior, function(x){reshape2::melt(as.matrix(x))})

  # keep sample index
  df_beta_posterior <- purrr::map_dfr(ldf_beta_posterior, .f = rbind, .id = "index")

  # pivot wider to obtain cov-matrix of predictors across posterior samples
  vcov_beta <- df_beta_posterior |>
    tidyr::pivot_wider(id_cols = index,
                names_from = c(Var1, Var2)) |>
    dplyr::select(!index) |>
    stats::cov()
  return(vcov_beta)
}


# Get Beta Variance -------------------------------------------------------
# Input: fit object res
# iter: number of iterations
f_betavar <- function(res){
  iter <- res$iter
  beta_var <- apply(res$fit$beta[,,51:(res$iter+50)], 1:2, var)   # delete first 50 samples
  beta_var
  
}
  
  


# Compute reliability -----------------------------------------------------
# Input: Vector of beta weights and posterior samples covariance matrix
f_rel <- function(beta, covmat){
  var_bw <- sd(beta)
  var_wi <- tr(covmat)
}


















