# -------------------------------------------------------------------------
# POSTERIOR PREDICTIVE VAR ------------------------------------------------
# -------------------------------------------------------------------------
library(tidyverse)
library(BGGM)    # modelling
library(qgraph)  # plotting
library(cowplot) # combining graphs

# First part based on vignette here: 
# https://cran.r-project.org/web/packages/BGGM/vignettes/var_model.html

# Data --------------------------------------------------------------------
# Use data from ifit study, 98 obs. 
Y <- subset(ifit, id == 1)[,-1]


# VAR Estimation ----------------------------------------------------------
fit <- var_estimate(Y, beta_sd = 1)

# Summarize model coefficients
print(
  summary(fit,  cred = 0.95), 
  param = "pcor"
)

plts <- plot(summary(fit,  cred = 0.95))

cowplot::plot_grid(
  cowplot::plot_grid(
    plts$beta_plt$interested, 
    plts$beta_plt$disinterested, 
    plts$beta_plt$excited, 
    nrow = 1),
  cowplot::plot_grid(  
    plts$beta_plt$upset,
    plts$beta_plt$strong,
    plts$beta_plt$stressed, 
    nrow = 1
  ), 
  nrow = 2)



# select a Graph with Credibility intervals
sel <- select(fit, cred = 0.95)
par(mfrow=c(1,2))

# Plot the graphs 
# INTERESTING: I cannot replicate Donald Williams results exactly

# Network of undirected relations (so contemporaneous network?)
qgraph::qgraph(sel$pcor_weighted_adj, title = "Partials")

# Network of directed relations (lagged coefficients)
qgraph::qgraph(sel$beta_weighted_adj, title = "Coefficients")

# Compute predictability
# note that this uses all variables, also those not selected in the select
# function above
predictability(fit)
plot(predictability(fit), type = "ridgeline")



# Comparison with GraphicalVAR Results ------------------------------------



# Posterior Predictive Distribution ---------------------------------------
# For the GGMs, Williams calls the "BGGM_Theta_continuous" function
# or the "BGGM_mv_continuous" function 
# to sample from the posterior for continuous data

# The fit object contains the posterior samples
# but this should not be equivalent to the posterior predictive samples
data.frame(matrix(fit$fit$beta, nrow = 5050, byrow = TRUE))


arr_beta <- fit$fit$beta
dimnames(arr_beta) <- list(var1 = paste0("x", seq(1:7)), var2 = paste0("x", seq(1:7)),
                                         iter = seq(1:5050))

# dataframe of posterior distribution (long format)
df_pdist<- reshape2::melt(arr_beta)

df_pdist %>% 
  filter(var1 == "x1" & var2 == "x1") %>% 
  ggplot(aes(x = value)) + 
  ggdist::stat_halfeye() + 
  theme_light() 


# # just loop this out, this is very inefficient but just for now
# l_beta <- list()
# for(i in 1:dim(arr_beta)[3]){
#   l_beta[[i]] <- as_tibble(arr_beta[,,i])
# }


# to just approximate it, do I even need the posterior?
# or can I just compute the mean and the variance of the samples from the posterior distribution? 
#  as per Gelman & Hill p.143, I also need the covariance matrix 
# of regression coefficients
# can I just obtain im from the posterior samples???
#  maybe this helps? https://github.com/danheck/multinomineq/blob/master/R/posterior_predictive.R



# Posterior predictive distribution GGM -----------------------------------
# To understand the internal BGGM code, I use a cross-sectional 
# GGM to obtain some posterior predictive distributions
Y <- gss

fit <- estimate(as.matrix(Y),
                impute = TRUE,
                iter = 150, type = "mixed")
yrep <- posterior_predict(fit, iter = 100)



