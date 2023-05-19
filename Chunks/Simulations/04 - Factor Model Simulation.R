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
T1 = 200
T0 = T1

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

treat_inter = 1

iter = 5
# J_max = min(round(T1 / 2.5,0), 70)
J_max = min(round(T1 / 2.5,0), 70)
CV_share = 3/4
my_by = 5
J_seq = seq(5, J_max, by = my_by)

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
    
    results$MSPE_SC[ID] = result_prelim$SC[1]
    results$MSPE_SC_BIAS[ID] = result_prelim$SC[2]  
    results$MSPE_SC_VAR[ID] = result_prelim$SC[3]      
    results$MSFE_SC[ID] = result_prelim$SC[4] 
    results$MSFE_SC_BIAS[ID] = result_prelim$SC[5]
    results$MSFE_SC_VAR[ID] = result_prelim$SC[6]
    
    results$MSPE_OLS[ID] = result_prelim$OLS[1]
    results$MSPE_OLS_BIAS[ID] = result_prelim$OLS[2]  
    results$MSPE_OLS_VAR[ID] = result_prelim$OLS[3]      
    results$MSFE_OLS[ID] = result_prelim$OLS[4] 
    results$MSFE_OLS_BIAS[ID] = result_prelim$OLS[5]
    results$MSFE_OLS_VAR[ID] = result_prelim$OLS[6] 
    
    results$MSPE_REGOLS[ID] = result_prelim$REGOLS[1]
    results$MSPE_REGOLS_BIAS[ID] = result_prelim$REGOLS[2]  
    results$MSPE_REGOLS_VAR[ID] = result_prelim$REGOLS[3]      
    results$MSFE_REGOLS[ID] = result_prelim$REGOLS[4] 
    results$MSFE_REGOLS_BIAS[ID] = result_prelim$REGOLS[5]
    results$MSFE_REGOLS_VAR[ID] = result_prelim$REGOLS[6]
    
    results$MSPE_NET[ID] = result_prelim$NET[1]
    results$MSPE_NET_BIAS[ID] = result_prelim$NET[2]  
    results$MSPE_NET_VAR[ID] = result_prelim$NET[3]      
    results$MSFE_NET[ID] = result_prelim$NET[4] 
    results$MSFE_NET_BIAS[ID] = result_prelim$NET[5]
    results$MSFE_NET_VAR[ID] = result_prelim$NET[6] 

    results$MSPE_FACTOR[ID] = result_prelim$FACTOR[1]
    results$MSPE_FACTOR_BIAS[ID] = result_prelim$FACTOR[2]  
    results$MSPE_FACTOR_VAR[ID] = result_prelim$FACTOR[3]      
    results$MSFE_FACTOR[ID] = result_prelim$FACTOR[4] 
    results$MSFE_FACTOR_BIAS[ID] = result_prelim$FACTOR[5]
    results$MSFE_FACTOR_VAR[ID] = result_prelim$FACTOR[6] 
    
    rm(result_prelim)
    svMisc::progress(ID, nrow(results))
  }
}

results = results %>% 
  select(Donors, bound_check,
         starts_with(c("MSPE_SC", "MSFE_SC",
                       "MSPE_OLS", "MSFE_OLS",
                       "MSPE_REGOLS", "MSFE_REGOLS",
                       "MSPE_NET", "MSFE_NET",
                       "MSPE_FACTOR", "MSFE_FACTOR")))

results_mean = results %>% 
  group_by(Donors, bound_check) %>% 
  summarise_at(.vars = vars(MSPE_SC:MSFE_FACTOR_VAR),
               .funs = mean)
  
df_meta = results_mean %>%
  gather(type, value, c(MSFE_SC, MSFE_OLS, MSFE_REGOLS, MSFE_NET, MSFE_FACTOR)) %>% 
  mutate(type = case_when(
    type == "MSFE_SC" ~ 1,
    type == "MSFE_OLS" ~ 2,
    type == "MSFE_REGOLS" ~ 3,
    type == "MSFE_NET" ~ 4,
    type == "MSFE_FACTOR" ~ 5)) %>% 
  mutate(type = as.factor(type)) %>% 
  mutate(type = recode_factor(type,
                              `1` = "MSFE_SC",
                              `2` = "MSFE_OLS",
                              `3` = "MSFE_REGOLS",
                              `4` = "MSFE_NET",
                              `5` = "MSFE_FACTOR"))

bound_1 = df_meta %>%
  filter(bound_check == 1) %>%
  ggplot() +
  aes(x = Donors, y = value, colour = type) +
  geom_point(shape = "circle", size = 2.8, alpha = .7) +
  scale_color_hue(direction = 1) +
  labs(
    x = "Number of Donors",
    y = "MSFE",
    title = paste("Factor Data // ", T0, " Obs with Means based on ", iter, " Iterations"),
    subtitle = "Treatment-Intercept within Donor-Intercept") +
  theme_minimal()

bound_0 = df_meta %>%
  filter(bound_check == 0) %>%
  ggplot() +
  aes(x = Donors, y = value, colour = type) +
  geom_point(shape = "circle", size = 2.8, alpha = .7) +
  scale_color_hue(direction = 1) +
  labs(
    x = "Number of Donors",
    y = "MSFE",
    title = paste("Factor Data // ", T0, " Obs with Means based on ", iter, " Iterations"),
    subtitle = "Treatment-Intercept not within Donor-Intercept") +
  theme_minimal()

grid.arrange(bound_1, bound_0, ncol=2)

writexl::write_xlsx(results, 
                    "~/Diss/Topics/Synthetic Control/Chunks/Simulations/Plots/Factor_results.xlsx")
