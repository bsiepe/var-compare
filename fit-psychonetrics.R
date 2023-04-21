library(psychonetrics)
library(tidyverse)
graph <- graphicalVAR::randomGVARmodel(6, probKappaEdge = .3, probBetaPositive = )
psn <- graphicalVARsim(200, beta = graph$beta, kappa = graph$kappa)
psn2 <- graphicalVARsim(200, beta = graph$beta, kappa = graph$kappa)
lam <- diag(6)


tsres <- psychonetrics::tsdlvm1(data = psn,
                                lambda = lam)
tsres2 <- psychonetrics::tsdlvm1(data = psn2,
                                lambda = lam)

res <- tsres %>% 
  runmodel()
res2 <- tsres2 %>% 
  runmodel()

res %>% 
  fit()
res2 %>% 
  fit()


res_pruned <- res %>% 
  runmodel %>%
  prune(alpha = 0.01, adjust = "none", recursive = FALSE, matrices = c("beta", "omega_zeta"))
res_pruned2 <- res2 %>% 
  runmodel %>%
  prune(alpha = 0.01, adjust = "none", recursive = FALSE, matrices = c("beta", "omega_zeta"))

# modelsearch
res_pruned <- res_pruned %>% 
  modelsearch(prunealpha = 0.01, addalpha = 0.01)
res_pruned2 <- res_pruned2 %>% 
  modelsearch(prunealpha = 0.01, addalpha = 0.01)

# compare res against pruned
compare(
  original = res, 
  pruned = res_pruned
)


# compare two models
compare(
  mod1 = res_pruned, 
  mod2 = res_pruned2
)



?compare



# Equality constraints ----------------------------------------------------
# see https://osf.io/y58mey58me
# Combine data
psn <- as.data.frame(psn)
psn2 <- as.data.frame(psn2)
psn$id = 1
psn2$id = 2

comb_data <- rbind(psn, psn2)

mod_eq <- psychonetrics::tsdlvm1(data = comb_data,
                                 lambda = lam,
                                 contemporaneous = "ggm",
                                 residual = "ggm",
                                 estimator = "FIML",
                                 identification = c("loadings","variance"),
                                 beta = "full",
                                 groups = "id")

mod_eq_run <- mod_eq %>% 
  runmodel()

mod_eq_groupequal <- groupequal(mod_eq, matrix = "beta", runmodel = "TRUE")

mod_eq_groupequal %>% fit
