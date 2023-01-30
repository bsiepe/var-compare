library(psychonetrics)
psn <- graphicalVARsim(200, beta = graph$beta, kappa = graph$kappa)

lam <- diag(6)


tsres <- psychonetrics::tsdlvm1(data = psn,
                                lambda = lam)

res <- tsres %>% 
  runmodel()

res %>% 
  fit()

res_pruned <- res %>% 
  runmodel %>%
  prune(alpha = 0.01, adjust = "none", recursive = FALSE, matrices = c("beta", "omega_zeta"))


# modelsearch
res_pruned <- res_pruned %>% 
  modelsearch(prunealpha = 0.01, addalpha = 0.01)


compare(
  original = res, 
  pruned = res_pruned
)




