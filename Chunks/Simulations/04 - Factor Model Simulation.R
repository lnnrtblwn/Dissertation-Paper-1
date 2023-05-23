library(tidyverse)
library(gridExtra)

rm(list = ls())

if (Sys.info()[6] == "jctoe"){
  setwd("C:/Promotion/SC_Paper/Documents/sc_sim Ferman/Ferman") 
  } else {
  setwd("~/Diss/Topics/Synthetic Control/Documents/sc_sim Ferman/Ferman") 
  }
source("my_functions.R")
set.seed(052023)

# 1. DATA GENERATING PROCESS: FACTOR MODEL WITHOUT COVARIATES ---- 

# Number of pre-and post-treatment periods
T1 = 20
T0 = 50

# AR-Term in Factor model. y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
# rho = 0.5
rho = 0

# Intercept. Set it equal to mean*(1-rho) to define mean of process
alpha = 0*(1-rho)

# Specify variance of u_t. Set it to (1 - rho^2) will lead to var(\lambda^k_t) = 1. Variance of the factors
var_u = (1-rho^2)

# Specify variance of transitory shocks in Factor model equation. Variance of the error terms 
var_epsilon = 1

# Post-treatment effects. Could be specified differently
post_effect = 10

# Number of factors
K = 2

# Adding a trend
#c = 0.02
c = 0

# Group distribution of each factor 
group_distribution = list(
  "lambda1" = c(1,0),
  "lambda2" = c(0,1))

# Specify intercept of treatment-unit. c(rnorm(1, mean = treat_inter, sd = 1), rnorm(J, mean = 0, sd = 1))
treat_inter = 0

iter = 1000
# J_max = min(round(T1 / 2.5,0), 70)
J_max = 30
CV_share = .5
my_by = 5
# J_seq = seq(5, J_max, by = my_by)
J_seq = c(5,10,15,20,25,30)

results = data.frame(matrix(NA, nrow = iter*length(J_seq), ncol = 1)) %>% 
    rename(Donors = c(1))

# 2. SIMULATION ---- 

for (J in J_seq) {
  
  for (i in 1:iter) {
    
    ID = (((J - J_seq[1]) / my_by) * iter) + i
    
    # ensure this is not overwritten for each J
    result_prelim = simulation_factor(J)
    
    results$Donors[ID] = J
    results$bound_check[ID] = result_prelim$bound_check
    
    results$PRE_SC_RMSPE[ID] = result_prelim$SC[1]
    results$PRE_SC_BIAS[ID] = result_prelim$SC[2]  
    results$PRE_SC_VAR[ID] = result_prelim$SC[3]      
    results$POST_SC_RMSFE[ID] = result_prelim$SC[4] 
    results$POST_SC_BIAS[ID] = result_prelim$SC[5]
    results$POST_SC_VAR[ID] = result_prelim$SC[6]
    
    results$PRE_OLS_RMSPE[ID] = result_prelim$OLS[1]
    results$PRE_OLS_BIAS[ID] = result_prelim$OLS[2]  
    results$PRE_OLS_VAR[ID] = result_prelim$OLS[3]      
    results$POST_OLS_RMSFE[ID] = result_prelim$OLS[4] 
    results$POST_OLS_BIAS[ID] = result_prelim$OLS[5]
    results$POST_OLS_VAR[ID] = result_prelim$OLS[6]
    
    results$PRE_REGOLS_RMSPE[ID] = result_prelim$REGOLS[1]
    results$PRE_REGOLS_BIAS[ID] = result_prelim$REGOLS[2]  
    results$PRE_REGOLS_VAR[ID] = result_prelim$REGOLS[3]      
    results$POST_REGOLS_RMSFE[ID] = result_prelim$REGOLS[4] 
    results$POST_REGOLS_BIAS[ID] = result_prelim$REGOLS[5]
    results$POST_REGOLS_VAR[ID] = result_prelim$REGOLS[6]
    
    results$PRE_NET_RMSPE[ID] = result_prelim$NET[1]
    results$PRE_NET_BIAS[ID] = result_prelim$NET[2]  
    results$PRE_NET_VAR[ID] = result_prelim$NET[3]      
    results$POST_NET_RMSFE[ID] = result_prelim$NET[4] 
    results$POST_NET_BIAS[ID] = result_prelim$NET[5]
    results$POST_NET_VAR[ID] = result_prelim$NET[6]

    results$PRE_FACTOR_RMSPE[ID] = result_prelim$FACTOR[1]
    results$PRE_FACTOR_BIAS[ID] = result_prelim$FACTOR[2]  
    results$PRE_FACTOR_VAR[ID] = result_prelim$FACTOR[3]      
    results$POST_FACTOR_RMSFE[ID] = result_prelim$FACTOR[4] 
    results$POST_FACTOR_BIAS[ID] = result_prelim$FACTOR[5]
    results$POST_FACTOR_VAR[ID] = result_prelim$FACTOR[6] 
    
    rm(result_prelim)
    svMisc::progress(ID, nrow(results))
  }
}

writexl::write_xlsx(results, 
                    "~/Diss/Topics/Synthetic Control/Chunks/Simulations/Results/Factor/Factor_results_50_20.xlsx")

results = readxl::read_excel("~/Diss/Topics/Synthetic Control/Chunks/Simulations/Results/Factor/Factor_results_50_20.xlsx")

results_mean = results %>% 
  group_by(Donors) %>% 
  summarise_at(.vars = dplyr::vars(PRE_SC_RMSPE:POST_FACTOR_VAR),
               .funs = mean) 
  
df_meta = results_mean %>%
  select(Donors, # bound_check,
         POST_SC_RMSFE, POST_SC_BIAS, 
         POST_OLS_RMSFE, POST_OLS_BIAS, 
         POST_REGOLS_RMSFE, POST_REGOLS_BIAS, 
         POST_NET_RMSFE, POST_NET_BIAS, 
         POST_FACTOR_RMSFE, POST_FACTOR_BIAS) %>% 
  gather(type, value, c(POST_SC_RMSFE, POST_SC_BIAS, 
                        POST_OLS_RMSFE, POST_OLS_BIAS, 
                        POST_REGOLS_RMSFE, POST_REGOLS_BIAS, 
                        POST_NET_RMSFE, POST_NET_BIAS, 
                        POST_FACTOR_RMSFE, POST_FACTOR_BIAS)) %>% 
  mutate(type = case_when(
    type == "POST_SC_RMSFE" ~ 1,
    type == "POST_OLS_RMSFE" ~ 2,
    type == "POST_REGOLS_RMSFE" ~ 3,
    type == "POST_NET_RMSFE" ~ 4,
    type == "POST_FACTOR_RMSFE" ~ 5,
    type == "POST_SC_BIAS" ~ 6,
    type == "POST_OLS_BIAS" ~ 7,
    type == "POST_REGOLS_BIAS" ~ 8,
    type == "POST_NET_BIAS" ~ 9,
    type == "POST_FACTOR_BIAS" ~ 10)) %>% 
  mutate(type = as.factor(type)) %>% 
  mutate(type = recode_factor(type,
                              `1` = "POST_SC_RMSFE",
                              `2` = "POST_OLS_RMSFE",
                              `3` = "POST_REGOLS_RMSFE",
                              `4` = "POST_NET_RMSFE",
                              `5` = "POST_FACTOR_RMSFE",
                              `6` = "POST_SC_BIAS",
                              `7` = "POST_OLS_BIAS",
                              `8` = "POST_REGOLS_BIAS",
                              `9` = "POST_NET_BIAS",
                              `10` = "POST_FACTOR_BIAS")) 

df_RMSFE = df_meta %>% 
  filter(type %in% c("POST_SC_RMSFE", "POST_OLS_RMSFE", "POST_REGOLS_RMSFE", "POST_NET_RMSFE", "POST_FACTOR_RMSFE"))

df_BIAS = df_meta %>% 
  filter(type %in% c("POST_SC_BIAS", "POST_OLS_BIAS", "POST_REGOLS_BIAS", "POST_NET_BIAS", "POST_FACTOR_BIAS"))

p_RMSFE = df_RMSFE %>%
  ggplot() +
  aes(x = Donors, y = value, colour = type) +
  geom_point(shape = "circle", size = 2.8, alpha = .7) +
  scale_color_hue(direction = 1) +
  labs(
    x = "Number of Donors",
    y = "MSFE-Mean",
    title = paste("Factor Data // ", iter, " Iterations"),
    subtitle = paste(T0, "Pre-Treatment, ", T1, "Post-Treatment Periods.")) +
  theme_minimal()

p_BIAS = df_BIAS %>%
  ggplot() +
  aes(x = Donors, y = value, colour = type) +
  geom_point(shape = "circle", size = 2.8, alpha = .7) +
  scale_color_hue(direction = 1) +
  labs(
    x = "Number of Donors",
    y = "BIAS-Mean",
    title = paste("Factor Data // ", iter, " Iterations"),
    subtitle = paste(T0, "Pre-Treatment, ", T1, "Post-Treatment Periods.")) +
  theme_minimal()

grid.arrange(p_RMSFE, p_BIAS, ncol=2)
