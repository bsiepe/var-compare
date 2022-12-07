# Find ToDos --------------------------------------------------------------
# todor::todor(c("TODO"))

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


# Sim from posterior ------------------------------------------------------
# TODO save function arguments
# TODO make convert_bggm work even if data generation did not work
# TODO especially here, graphicalVARsim might be an issue? it has a slightly different model
sim_from_post_parallel <- function(fitobj, 
                                   n_datasets, 
                                   n,
                                   tp,
                                   iterations,
                                   means = 0,
                                   convert_bggm = TRUE){
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
  if(convert_bggm == TRUE){
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


# Cross compare two models with posterior ---------------------------------
#' Cross-compare two models with their posterior
#'
#' @param postdata_a List of data sampled from posterior of model/person A.
#' @param postdata_b List of data sampled from posterior of model/person B.
#' @param fitres_a List of model results for each participant under specific DGP.
#' @param fitres_b List of model results for each participant under specific DGP.
#' @param mod_a Numerical indicator of model/person A. 
#' @param mod_b Numerical indicator of model/person B. 
#' @param n_datasets Number of datasets sampled from posterior. 
#' @param comparison Which comparison to use. "Frob" for Frobenius-Norm, "maxdiff" for maximum edge difference. 
#' @param ... Currently ignored. 
#' 
#'
#' @return Dataframe with null distributions for Frobenius norm of both models under the Null and empirical Frobenius norm between both models 
#' @export
#'
#' 
cross_compare_emp <- function(postdata_a = l_data,
                          postdata_b = l_data,
                          fitres_a = l_res,
                          fitres_b = l_res,
                          mod_a = 1, 
                          mod_b = 2,
                          rho_prior = rho_sd, 
                          beta_prior = beta_sd,
                          n_datasets = 100,
                          iterations = n_iter, 
                          comparison = "frob",
                          ...){
  if(!is.numeric(mod_a) | !is.numeric(mod_b)){
    stop("Error: Model needs to have numerical index")
  }
  
  # Refit to posterior datasets of each model

  l_postpred_a <- foreach(i = seq(n_datasets), .packages = "BGGM") %dopar% {
    if(is.list(postdata_a[[mod_a]][[i]])){
      postpreda <- try(BGGM::var_estimate(postdata_a[[mod_a]][[i]]$Y,
                                          rho_sd = rho_prior,
                                          beta_sd = beta_prior,
                                          iter = iterations,
                                          progress = FALSE,
                                          seed = 2022))
    }
    else 
      postpreda <- NA
      

    postpreda
    }
    
    
  l_postpred_b <- foreach(i = seq(n_datasets), .packages = "BGGM") %dopar% {
    if(is.list(postdata_b[[mod_b]][[i]])){
      postpredb <- try(BGGM::var_estimate(postdata_b[[mod_b]][[i]]$Y,
                                          rho_sd = rho_prior,
                                          beta_sd = beta_prior,
                                          iter = iterations,
                                          progress = FALSE,
                                          seed = 2022))
      
    }
    else
      postpredb <- NA
      
    postpredb
  }
    

    

  # Frobenius Norm Comparison 
  if(comparison == "frob"){
    normtype = "F"
    # Compute Distance of empirical beta to posterior samples beta
    frob_null_a <- unlist(lapply(l_postpred_a, 
                                 function(x){
                                   if(!is.list(x))
                                     {NA}
                                   else 
                                     {norm(fitres_a[[mod_a]]$beta_mu-x$beta_mu, type = normtype)}}))
    frob_null_b <- unlist(lapply(l_postpred_b, 
                                 function(x){
                                   if(!is.list(x))
                                    {NA}
                                   else
                                   {norm(fitres_b[[mod_b]]$beta_mu-x$beta_mu, type = normtype)}}))

    
    # Compute Distance of empirical betas between a and b
    frob_emp <- norm(fitres_a[[mod_a]]$beta_mu - fitres_b[[mod_b]]$beta_mu)
    
    cc_res <- data.frame(model_ind = c(rep(mod_a, n_postds), rep(mod_b, n_postds)),
                         frob_null = c(frob_null_a, frob_null_b),
                         frob_emp = rep(frob_emp, n_postds*2))
    
    
  }
  # Max Difference Comparison
  if(comparison == "maxdiff"){
    # Compute maximum distance of empirical beta to posterior samples beta
    maxdiff_a <- unlist(lapply(l_postpred_a, 
                                 function(x){
                                   if(!is.list(x))
                                     {NA}
                                   else 
                                     {max(abs(fitres_a[[mod_a]]$beta_mu-x$beta_mu))}}))
    maxdiff_b <- unlist(lapply(l_postpred_b, 
                                 function(x){
                                   if(!is.list(x))
                                     {NA}
                                 else
                                   {max(abs(fitres_b[[mod_b]]$beta_mu-x$beta_mu))}})) 
    # Compute maxdiff of empirical betas between a and b
    maxdiff_emp <- max(abs(fitres_a[[mod_a]]$beta_mu - fitres_b[[mod_b]]$beta_mu))
    
    cc_res <- data.frame(model_ind = c(rep(mod_a, n_postds), rep(mod_b, n_postds)),
                         maxdiff_null = c(maxdiff_a, maxdiff_b),
                         maxdiff_emp = rep(maxdiff_emp, n_postds*2))
    

    } # end maxdiff
  
  return(cc_res)
} # end function

  










# Cross-compare all posterior samples -------------------------------------
# Don't look at empirical distance, but distances between 100 resamples
cross_compare_post <- function(postdata_a = l_dat_bggm,
                          postdata_b = l_dat_bggm,
                          fitres_a = l_res,
                          fitres_b = l_res,
                          mod_a = 1, 
                          mod_b = 2,
                          n_postds = 100,
                          normtype = "F",
                          ...){
  if(!is.numeric(mod_a) | !is.numeric(mod_b)){
    stop("Error: Model needs to have numerical index")
  }
  # Refit to posterior samples of each model
  l_postpred_a <- list()
  for(i in seq(n_postds)){
    l_postpred_a[[i]] <- try(BGGM::var_estimate(postdata_a[[mod_a]][[i]]$Y,
                                                    rho_sd = rho_sd,
                                                    beta_sd = beta_sd,
                                                    iter = n_iter,
                                                    progress = FALSE,
                                                    seed = 2022))
    
  }
  l_postpred_b <- list()
  for(i in seq(n_postds)){
    l_postpred_b[[i]] <- try(BGGM::var_estimate(postdata_b[[mod_b]][[i]]$Y,
                                                    rho_sd = rho_sd,
                                                    beta_sd = beta_sd,
                                                    iter = n_iter,
                                                    progress = FALSE,
                                                    seed = 2022))
    
  }
  

  
  # Compute Distance of empirical beta to posterior samples beta
  frob_null_a <- unlist(lapply(l_postpred_a, 
                               function(x){
                                 if(!is.list(x)){NA}
                                 else 
                                   norm(fitres_a[[mod_a]]$beta_mu-x$beta_mu, type = normtype)}))
  frob_null_b <- unlist(lapply(l_postpred_b, 
                               function(x)
                                 if(!is.list(x)){NA}
                               else{norm(fitres_b[[mod_b]]$beta_mu-x$beta_mu, type = normtype)}))
  
  # Compute distance between posterior samples
  frob_emp_post <- mapply(function(x,y){if(!is.list(x) | !is.list(y)){NA} else norm(x$beta_mu-y$beta_mu, type = normtype)} , l_postpred_a, l_postpred_b)
  
  
  # Compute Distance of empirical betas between a and b
  frob_emp <- norm(fitres_a[[mod_a]]$beta_mu - fitres_b[[mod_b]]$beta_mu)
  
  cc_res <- data.frame(model_ind = c(rep(paste0("null_", mod_a), n_postds), 
                                     rep(paste0("null", mod_b), n_postds),
                                     rep(paste0("post_", mod_a,"-", mod_b), n_postds)),
                       frob = c(frob_null_a, frob_null_b, frob_emp_post))
  cc_res
  
}




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

cross_compare_eval <- function(df_res,
                               plot = FALSE){
  
  # Obtain model indexes
  model_ind_a <- unique(df_res$model_ind)[1]
  model_ind_b <- unique(df_res$model_ind)[2]
  
  # Number of posterior difference > empirical difference
  teststat_a <- sum(df_res$frob_null[df_res$model_ind == model_ind_a] > df_res$frob_emp[df_res$model_ind == model_ind_a], na.rm = TRUE)
  teststat_b <- sum(df_res$frob_null[df_res$model_ind == model_ind_b] > df_res$frob_emp[df_res$model_ind == model_ind_b], na.rm = TRUE)
  
  
  
  if(plot == TRUE){
    print(
      df_res %>% 
        mutate(model_ind = as.factor(model_ind)) %>% 
        ggplot(aes(x = frob_null, fill = model_ind))+
        geom_density() +
        theme_minimal()+
        geom_vline(aes(xintercept = max(frob_emp)), col = "black",
                   linetype = "dashed")+
        geom_label(aes(x = max(frob_emp), y = 1),
                   label = "Emp. Norm", show.legend = FALSE)+
        labs(fill = "Model",
             x = "Frobenius Norm",
             y = "Density")
    )
    
    
  }
  testres <- list(res_a = teststat_a,
                  res_b = teststat_b)
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


















