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


