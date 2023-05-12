library(tidyverse)

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
T0 = T1

# AR-Term in Factor model. y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
rho = 0.5

# Intercept. Set it equal to mean*(1-rho) to define mean of process
alpha = 0*(1-rho)

# Specify variance of u_t. Set it to (1 - rho^2) will lead to var(\lambda^k_t) = 1. Variance of the factors
var_u = (1-rho^2)*20

# Specify variance of transitory shocks in Factor model equation. Variance of the error terms 
#var_epsilon = 0.1*5
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

id = 1
iter = 5
J_max = min(round(T1 / 2.5,0), 70)
CV_share = 3/4

results_mean = data.frame(matrix(NA, nrow = J_max-1, ncol = 7)) %>%
  rename(
    Donors = c(7),
    RMSPE_SC = c(1),
    RMSPE_OLS = c(2),
    RMSPE_REGOLS = c(3),
    RMSFE_SC = c(4),
    RMSFE_OLS = c(5),
    RMSFE_REGOLS = c(6))

# 2. SIMULATION ---- 

for (J in seq(2, J_max, by = 1)) {
  
  results = data.frame(matrix(NA, nrow = iter, ncol = 4)) %>%
    rename(
      RMSPE_SC = c(1),
      RMSPE_OLS = c(2),
      RMSFE_SC = c(3),
      RMSFE_OLS = c(4))
  
  for (i in 1:iter) {
    
    result_prelim = simulations(J)
    
    results$RMSPE_SC[i] = result_prelim[1]
    results$RMSFE_SC[i] = result_prelim[2]
    results$RMSPE_OLS[i] = result_prelim[3]
    results$RMSFE_OLS[i] = result_prelim[4]
    results$RMSPE_REGOLS[i] = result_prelim[5]
    results$RMSFE_REGOLS[i] = result_prelim[6]
    
    rm(result_prelim)
  }


  results_mean$Donors[id] = J
  # Predictions
  results_mean$RMSPE_SC[id] = mean(results$RMSPE_SC)
  results_mean$RMSPE_OLS[id] = mean(results$RMSPE_OLS)
  results_mean$RMSPE_REGOLS[id] = mean(results$RMSPE_REGOLS)
  
  # Forecasts
  results_mean$RMSFE_SC[id] = mean(results$RMSFE_SC)
  results_mean$RMSFE_OLS[id] = mean(results$RMSFE_OLS)
  results_mean$RMSFE_REGOLS[id] = mean(results$RMSFE_REGOLS)
  
  id = id + 1
  svMisc::progress(J, J_max)
}


# Visualization

df_meta = results_mean %>%
  gather(type, value, RMSPE_SC:RMSFE_REGOLS) %>% 
  mutate(type = case_when(
    type == "RMSPE_OLS" ~ 1,
    type == "RMSPE_SC" ~ 2,
    type == "RMSPE_REGOLS" ~ 3,
    type == "RMSFE_OLS" ~ 4,
    type == "RMSFE_SC" ~ 5,
    type == "RMSFE_REGOLS" ~ 6)) %>% 
  mutate(type = as.factor(type)) %>% 
  mutate(type = recode_factor(type,
                              `1` = "RMSPE_OLS",
                              `2` = "RMSPE_SC",
                              `3` = "RMSPE_REGOLS",
                              `4` = "RMSFE_OLS",
                              `5` = "RMSFE_SC",
                              `6` = "RMSFE_REGOLS"))

df_meta %>%
  filter(type %in% c("RMSFE_SC", "RMSFE_OLS", "RMSFE_REGOLS")) %>%
  ggplot() +
  aes(x = Donors, y = value, colour = type) +
  geom_line(size = 0.7) +
  #scale_color_viridis_b() +
  labs(
    x = "Number of Donors",
    y = "RMSFE",
    title = paste("Factor Data // ", T0, " Obs with Means based on ", iter, " Iterations")) +
  theme_minimal()

summary(results_mean)
