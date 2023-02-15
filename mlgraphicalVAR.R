# Try to fit multilevel VAR to data from two individuals 
# working with github version by  https://github.com/JPark93/graphicalVAR
# install with devtools::install_github()



ml_dat <- graphicalVAR::simMLgvar(100, 6, 2)
ml_dat


df_ml_dat <- as.data.frame(ml_dat$data)
df_ml_dat <- df_ml_dat %>% 
  mutate(beep = c(rep(seq(1,4),25),rep(seq(1,4),25)),
         day = c(rep(seq(1,25), each = 4), rep(seq(1,25), each = 4)))

var_names <- c("V1", "V2", "V3", "V4", "V5", "V6")

res_ml <- mlGraphicalVAR(df_ml_dat,
               vars = var_names,
               beepvar = "beep", 
               dayvar = "day",
               idvar = "ID",
               btwnNetworks = FALSE)




