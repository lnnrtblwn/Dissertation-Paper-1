library(tidyverse)
library(gridExtra)

rm(list = ls())

if (Sys.info()[6] == "jctoe"){
  setwd("C:/Promotion/SC_Paper/") 
} else {
  setwd("~/Diss/Topics/Synthetic Control/") 
}
source("Documents/sc_sim Ferman/Ferman/my_functions.R")
#source("Chunks/Simulations/07 - VAR_simu_GDP.R")
#set.seed(052023)

# 1. DATA GENERATING PROCESS: FACTOR MODEL WITHOUT COVARIATES ---- 

# Number of pre-and post-treatment periods
T1 = 20
T0 = 50

# AR-Term in Factor model. y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
# rho = 0.8 // Non.Stationary => rho = 1.0
rho = 0.8

# Error AR-Term --> defined in my_functions
# rho_u = runif(1, .5, 0.95)

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

# Lag for univariate dynamic case 
p_uni = 3


# Group distribution of each factor 
group_distribution = list(
  "lambda1" = c(1,0),
  "lambda2" = c(0,1))

# Specify intercept of treatment-unit. c(rnorm(1, mean = treat_inter, sd = 1), rnorm(J, mean = 0, sd = 1))
treat_inter = 0

iter = 20
# J_max = min(round(T1 / 2.5,0), 70)
J_max = 30
CV_share = .5
my_by = 5
# J_seq = seq(5, J_max, by = my_by)
J_seq = c(5,10,15,20,25,30)
#J_seq = 3

# J = 5
# simu_type = 'Factor'


results = data.frame(matrix(NA, nrow = iter*length(J_seq), ncol = 1)) %>% 
  rename(Donors = c(1))

plots_UNIDYN1 = list()
plots_UNIDYN2 = list()
plots_REGOLS = list()

# 2. SIMULATION ---- 

#simu_type = "VAR"
simu_type = "Factor"
# J = 5

for (J in J_seq) {
  
  for (i in 1:iter) {
    
    ID = (((J - J_seq[1]) / my_by) * iter) + i
    
    # ensure this is not overwritten for each J
    result_prelim = simulation_factor(J, simu_type = simu_type)
    
    results$Donors[ID] = J
    results$bound_check[ID] = result_prelim$bound_check
    
    results$rho_factor[ID] = result_prelim$rho_factor
    results$rho_error[ID] = result_prelim$rho_error

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
    
    results$PRE_UNIDYN1_RMSPE[ID] = result_prelim$UNIDYN1[1]
    results$PRE_UNIDYN1_BIAS[ID] = result_prelim$UNIDYN1[2]  
    results$PRE_UNIDYN1_VAR[ID] = result_prelim$UNIDYN1[3]      
    results$POST_UNIDYN1_RMSFE[ID] = result_prelim$UNIDYN1[4] 
    results$POST_UNIDYN1_BIAS[ID] = result_prelim$UNIDYN1[5]
    results$POST_UNIDYN1_VAR[ID] = result_prelim$UNIDYN1[6] 
    
    results$PRE_UNIDYN2_RMSPE[ID] = result_prelim$UNIDYN2[1]
    results$PRE_UNIDYN2_BIAS[ID] = result_prelim$UNIDYN2[2]  
    results$PRE_UNIDYN2_VAR[ID] = result_prelim$UNIDYN2[3]      
    results$POST_UNIDYN2_RMSFE[ID] = result_prelim$UNIDYN2[4] 
    results$POST_UNIDYN2_BIAS[ID] = result_prelim$UNIDYN2[5]
    results$POST_UNIDYN2_VAR[ID] = result_prelim$UNIDYN2[6] 
    
    plots_REGOLS[[ID]] = result_prelim$Plots_REGOLS
    plots_UNIDYN1[[ID]] = result_prelim$Plots_UNIDYN1
    plots_UNIDYN2[[ID]] = result_prelim$Plots_UNIDYN2
    
    rm(result_prelim)
    svMisc::progress(ID, nrow(results))
  }
}

ggsave(
  filename = "Chunks/Simulations/Results/VAR/hybrid/plots/TS_plots_REGOLS.pdf",
  plot = marrangeGrob(plots_REGOLS, nrow =1, ncol = 1),
  width = 15, height = 10)

ggsave(
  filename = "Chunks/Simulations/Results/VAR/hybrid/plots/TS_plots_UNIDYN1.pdf",
  plot = marrangeGrob(plots_UNIDYN1, nrow =1, ncol = 1),
  width = 15, height = 10)

ggsave(
  filename = "Chunks/Simulations/Results/VAR/hybrid/plots/TS_plots_UNIDYN2.pdf",
  plot = marrangeGrob(plots_UNIDYN2, nrow =1, ncol = 1),
  width = 15, height = 10)

writexl::write_xlsx(results, "Chunks/Simulations/Results/VAR/hybrid/Results_50_20_Factor.xlsx")

results_mean = results %>% 
  group_by(Donors) %>% 
  summarise_at(.vars = dplyr::vars(PRE_SC_RMSPE:POST_UNIDYN2_VAR),
               .funs = mean) %>% 
  select(Donors,
         ends_with("RMSFE")) %>% 
  select(Donors, 
         POST_FACTOR_RMSFE, POST_REGOLS_RMSFE, POST_NET_RMSFE, POST_UNIDYN1_RMSFE, POST_UNIDYN2_RMSFE)


t_0 = results %>% dplyr::select(POST_SC_BIAS, POST_REGOLS_BIAS)

t_0 = results %>% as.data.frame %>% select(POST_SC_BIAS)
colnames(results)
# writexl::write_xlsx(results, 
#                     "~/Diss/Topics/Synthetic Control/Chunks/Simulations/Results/Factor/Factor_results_50_30.xlsx")


results %>%  rowwise() %>%
  mutate(colnames(min(POST_SC_RMSFE, POST_OLS_RMSFE,POST_REGOLS_RMSFE,POST_NET_RMSFE,POST_FACTOR_RMSFE,POST_UNIDYN_RMSFE)))

df_check = results %>% 
  dplyr::select(c(POST_SC_RMSFE, POST_SC_BIAS, 
                  POST_OLS_RMSFE, POST_OLS_BIAS, 
                  POST_REGOLS_RMSFE, POST_REGOLS_BIAS, 
                  POST_NET_RMSFE, POST_NET_BIAS, 
                  POST_FACTOR_RMSFE, POST_FACTOR_BIAS,
                  POST_UNIDYN_RMSFE, POST_UNIDYN_BIAS)) %>% 
  gather(type, value, c(POST_SC_RMSFE, POST_SC_BIAS, 
                        POST_OLS_RMSFE, POST_OLS_BIAS, 
                        POST_REGOLS_RMSFE, POST_REGOLS_BIAS, 
                        POST_NET_RMSFE, POST_NET_BIAS, 
                        POST_FACTOR_RMSFE, POST_FACTOR_BIAS,
                        POST_UNIDYN_RMSFE, POST_UNIDYN_BIAS))

# writexl::write_xlsx(results, 
#                     "~/Diss/Topics/Synthetic Control/Chunks/Simulations/Results/VAR/hybrid/VAR_results_100_30_neg.xlsx")



results = readxl::read_excel("C:/Promotion/SC_Paper/Chunks/Simulations/Results/VAR/hybrid/VAR_results_100_20_neg.xlsx")

results_mean = results %>% 
  group_by(Donors)  %>% 
  summarise_at(.vars = dplyr::vars(PRE_SC_RMSPE:POST_UNIDYN_VAR),
               .funs = mean) %>% 
  dplyr::select(c(Donors, POST_SC_RMSFE, PRE_SC_RMSPE, POST_SC_BIAS,
                  POST_OLS_RMSFE, PRE_OLS_RMSPE, POST_OLS_BIAS,
                  POST_REGOLS_RMSFE, PRE_REGOLS_RMSPE,POST_REGOLS_BIAS,
                  POST_NET_RMSFE, PRE_NET_RMSPE, POST_NET_BIAS,
                  POST_FACTOR_RMSFE, PRE_FACTOR_RMSPE,POST_FACTOR_BIAS,
                  POST_UNIDYN_RMSFE, PRE_UNIDYN_RMSPE, POST_UNIDYN_BIAS))

df_meta = results_mean %>%
  dplyr::select(Donors, 
         POST_SC_RMSFE, POST_SC_BIAS, 
         POST_OLS_RMSFE, POST_OLS_BIAS, 
         POST_REGOLS_RMSFE, POST_REGOLS_BIAS, 
         POST_NET_RMSFE, POST_NET_BIAS, 
         POST_FACTOR_RMSFE, POST_FACTOR_BIAS,
         POST_UNIDYN_RMSFE, POST_UNIDYN_BIAS) %>% 
  gather(type, value, c(POST_SC_RMSFE, POST_SC_BIAS, 
                        POST_OLS_RMSFE, POST_OLS_BIAS, 
                        POST_REGOLS_RMSFE, POST_REGOLS_BIAS, 
                        POST_NET_RMSFE, POST_NET_BIAS, 
                        POST_FACTOR_RMSFE, POST_FACTOR_BIAS,
                        POST_UNIDYN_RMSFE, POST_UNIDYN_BIAS)) %>% 
  mutate(type = case_when(
    type == "POST_SC_RMSFE" ~ 1,
    type == "POST_OLS_RMSFE" ~ 2,
    type == "POST_REGOLS_RMSFE" ~ 3,
    type == "POST_NET_RMSFE" ~ 4,
    type == "POST_FACTOR_RMSFE" ~ 5,
    type == "POST_UNIDYN_RMSFE" ~ 6,
    type == "POST_SC_BIAS" ~ 7,
    type == "POST_OLS_BIAS" ~ 8,
    type == "POST_REGOLS_BIAS" ~ 9,
    type == "POST_NET_BIAS" ~ 10,
    type == "POST_FACTOR_BIAS" ~ 11,
    type == "POST_UNIDYN_BIAS" ~ 12)) %>% 
  mutate(type = as.factor(type)) %>% 
  mutate(type = recode_factor(type,
                              `1` = "POST_SC_RMSFE",
                              `2` = "POST_OLS_RMSFE",
                              `3` = "POST_REGOLS_RMSFE",
                              `4` = "POST_NET_RMSFE",
                              `5` = "POST_FACTOR_RMSFE",
                              `6` = "POST_UNIDYN_RMSFE",
                              `7` = "POST_SC_BIAS",
                              `8` = "POST_OLS_BIAS",
                              `9` = "POST_REGOLS_BIAS",
                              `10` = "POST_NET_BIAS",
                              `11` = "POST_FACTOR_BIAS",
                              `12` = "POST_UNIDYN_BIAS")) 

df_RMSFE = df_meta %>% 
  filter(type %in% c("POST_SC_RMSFE", "POST_OLS_RMSFE", "POST_REGOLS_RMSFE", "POST_NET_RMSFE", "POST_FACTOR_RMSFE", "POST_UNIDYN_RMSFE"))

df_BIAS = df_meta %>% 
  filter(type %in% c("POST_SC_BIAS", "POST_OLS_BIAS", "POST_REGOLS_BIAS", "POST_NET_BIAS", "POST_FACTOR_BIAS", "POST_UNIDYN_BIAS"))

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






