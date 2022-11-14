# Find ToDos --------------------------------------------------------------
todor::todor(c("TODO"))




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
# TODO: add option to manually provide a beta matrix with posterior means
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
  # TODO very important: Do I need to do byrow TRUE or not?
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
  out
}


# Combine y and yhat ------------------------------------------------------
# Function to get yhat and actual data into the same dataframe
# data = data object
# refit = predict.var_estimate object
f_merge_pred <- function(data, refit){
  y_pred <- array(NA, dim = c(dim(refit)[1],dim(refit)[3]))
  # loop over number of variables in refit object
  for(p in 1:dim(refit)[3]){
    y_pred[,p] <- refit[,,p][,"Post.mean"]
    dimnames(y_pred) <- list(NULL, paste0("V",1:6,"_pred"))
  }
  res <- cbind(data$Y, y_pred)
  return(res)
}

test <- f_merge_pred(data = l_dat[[1]], refit = l_refit[[1]])


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



