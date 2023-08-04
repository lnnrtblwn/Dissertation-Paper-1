# 00 PACKAGES AND WD ----

library(tidyverse)
library(gridExtra)

rm(list = ls())

if (Sys.info()[6] == "jctoe"){
  setwd("C:/Promotion/SC_Paper/") 
} else {
  setwd("~/Diss/Topics/Synthetic Control/") 
}

source("R-Scripts/Simulations/Functions/my_functions.R")
source("R-Scripts/Simulations/07 - VAR_simu_GDP.R")
#set.seed(052023)

# 01 DGP: FACTOR/VAR ---- 

## 01.1 Joint Settings ----

# Number of pre-and post-treatment periods
T0 = 50
T1 = 20

# Treatment Effect
post_effect = 10

# Lag for univariate and multivariate dynamic case 
p = 2 
p_uni = p
p_multi = p

## 01.2 FACTOR Settings ----

# AR-Term of Factors
rho_factor = 0

# AR-Term of Errors: rho_u =  runif(1, rho_u_left, rho_u_right)
rho_u = 0
# rho_u_left = 0.5
# rho_u_right = 0.95

# Factor-Intercept
alpha = 0*(1-rho_factor)

# Factor-Variance
var_factor = 1

# Error-Variance
var_u = 1

# Number of Factors
K = 2

# Group Distribution of Factors 
group_distribution = list(
  "lambda1" = c(1,0),
  "lambda2" = c(0,1))

# Adding a Trend
c = 0

# Specify Treatment-unit. c(rnorm(1, mean = treat_inter, sd = 1), rnorm(J, mean = 0, sd = 1))
treat_inter = 0

## 01.3 VAR Settings ----

# Error-Variance
var_error_VAR = 1

# 02 SIMULATION ---- 

## 02.1 Settings ----

iter = 50
CV_share = .5

# J = 4

# Factor
my_by = 1
# J_seq = c(5,10,15,20,25,30)
J_seq = c(5:20)
simu_type = "Factor"
dynamic = "no"

# VAR
# my_by = 2
# J_seq = c(2,4,6,8)
# simu_type = "VAR"
# dynamic = "yes"

results = data.frame(matrix(NA, nrow = iter*length(J_seq), ncol = 1)) %>% 
  rename(Donors = c(1))

plots_UNIDYN = list()
plots_REGOLS = list()
plots_MULTIDYN1 = list()
plots_MULTIDYN2 = list()
plots_MULTIDYN3 = list()
plots_VAR = list()

## 02.2 Simulation ----

for (J in J_seq) {
  
  for (i in 1:iter) {
    
    ID = (((J - J_seq[1]) / my_by) * iter) + i
    
    # ensure this is not overwritten for each J
    result_prelim = simulation_factor(J, simu_type = simu_type)
    
    results$Donors[ID] = J
    results$bound_check[ID] = result_prelim$bound_check
    
    results$rho_factor[ID] = result_prelim$rho_factor
    results$rho_error[ID] = result_prelim$rho_error
    results$l1[ID] = result_prelim$REGOLS[7] # Ridge
    results$l2[ID] = result_prelim$REGOLS[8] # Inverse Ridge   

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
    
    if (dynamic == "yes"){
      
    results$PRE_UNIDYN_RMSPE[ID] = result_prelim$UNIDYN[1]
    results$PRE_UNIDYN_BIAS[ID] = result_prelim$UNIDYN[2]  
    results$PRE_UNIDYN_VAR[ID] = result_prelim$UNIDYN[3]      
    results$POST_UNIDYN_RMSFE[ID] = result_prelim$UNIDYN[4] 
    results$POST_UNIDYN_BIAS[ID] = result_prelim$UNIDYN[5]
    results$POST_UNIDYN_VAR[ID] = result_prelim$UNIDYN[6] 
    
    results$PRE_MULTIDYN1_RMSPE[ID] = result_prelim$MULTIDYN1[1]
    results$PRE_MULTIDYN1_BIAS[ID] = result_prelim$MULTIDYN1[2]  
    results$PRE_MULTIDYN1_VAR[ID] = result_prelim$MULTIDYN1[3]      
    results$POST_MULTIDYN1_RMSFE[ID] = result_prelim$MULTIDYN1[4] 
    results$POST_MULTIDYN1_BIAS[ID] = result_prelim$MULTIDYN1[5]
    results$POST_MULTIDYN1_VAR[ID] = result_prelim$MULTIDYN1[6] 
    
    results$PRE_MULTIDYN2_RMSPE[ID] = result_prelim$MULTIDYN2[1]
    results$PRE_MULTIDYN2_BIAS[ID] = result_prelim$MULTIDYN2[2]  
    results$PRE_MULTIDYN2_VAR[ID] = result_prelim$MULTIDYN2[3]      
    results$POST_MULTIDYN2_RMSFE[ID] = result_prelim$MULTIDYN2[4] 
    results$POST_MULTIDYN2_BIAS[ID] = result_prelim$MULTIDYN2[5]
    results$POST_MULTIDYN2_VAR[ID] = result_prelim$MULTIDYN2[6]
    
    results$PRE_MULTIDYN3_RMSPE[ID] = result_prelim$MULTIDYN3[1]
    results$PRE_MULTIDYN3_BIAS[ID] = result_prelim$MULTIDYN3[2]  
    results$PRE_MULTIDYN3_VAR[ID] = result_prelim$MULTIDYN3[3]      
    results$POST_MULTIDYN3_RMSFE[ID] = result_prelim$MULTIDYN3[4] 
    results$POST_MULTIDYN3_BIAS[ID] = result_prelim$MULTIDYN3[5]
    results$POST_MULTIDYN3_VAR[ID] = result_prelim$MULTIDYN3[6] 
    
    results$PRE_VAR_RMSPE[ID] = result_prelim$VAR[1]
    results$PRE_VAR_BIAS[ID] = result_prelim$VAR[2]  
    results$PRE_VAR_VAR[ID] = result_prelim$VAR[3]      
    results$POST_VAR_RMSFE[ID] = result_prelim$VAR[4] 
    results$POST_VAR_BIAS[ID] = result_prelim$VAR[5]
    results$POST_VAR_VAR[ID] = result_prelim$VAR[6] 
    
    plots_REGOLS[[ID]] = result_prelim$Plots_REGOLS
    plots_UNIDYN[[ID]] = result_prelim$Plots_UNIDYN
    plots_MULTIDYN1[[ID]] = result_prelim$Plots_MULTIDYN1
    plots_MULTIDYN2[[ID]] = result_prelim$Plots_MULTIDYN2
    plots_MULTIDYN3[[ID]] = result_prelim$Plots_MULTIDYN3
    plots_VAR[[ID]] = result_prelim$Plots_VAR
    
    }
    
    rm(result_prelim)
    svMisc::progress(ID, nrow(results))
  }
  
}

## 02.3 Results ----

# results = results %>% 
#   filter(POST_VAR_RMSFE < 50)

results_mean = results %>% 
  group_by(Donors) %>% 
  summarise_at(.vars = dplyr::vars(PRE_SC_RMSPE:POST_FACTOR_VAR),
               .funs = mean) %>% 
  dplyr::select(Donors,
                ends_with("RMSFE")) 

t(results_mean)

boxplot(results %>%  
          dplyr::select(c(POST_UNIDYN_RMSFE, POST_MULTIDYN1_RMSFE, POST_MULTIDYN2_RMSFE, POST_MULTIDYN3_RMSFE, POST_VAR_RMSFE)) %>% 
          rename_with(~str_remove(.x, "^.{1,5}"), everything()))
          # dplyr::select(ends_with("RMSFE")) %>% 
          

# ggsave(
#   filename = "Chunks/Simulations/Results/Factor/20230730/TS_plots_REGOLS.pdf",
#   plot = marrangeGrob(plots_REGOLS, nrow =1, ncol = 1),
#   width = 15, height = 10)
# 
# ggsave(
#   filename = "Chunks/Simulations/Results/Factor/20230730/TS_plots_UNIDYN1.pdf",
#   plot = marrangeGrob(plots_UNIDYN1, nrow =1, ncol = 1),
#   width = 15, height = 10)
# 
# ggsave(
#   filename = "Chunks/Simulations/Results/Factor/20230730/TS_plots_UNIDYN2.pdf",
#   plot = marrangeGrob(plots_UNIDYN2, nrow =1, ncol = 1),
#   width = 15, height = 10)
# 
# ggsave(
#   filename = "Chunks/Simulations/Results/Factor/20230730/TS_plots_OLSDIST.pdf",
#   plot = marrangeGrob(plots_OLSDIST, nrow =1, ncol = 1),
#   width = 15, height = 10)

# ggsave(
#   filename = "Chunks/Simulations/Results/VAR/fulll dynamic/20230803/TS_plots_VAR.pdf",
#   plot = marrangeGrob(plots_VAR, nrow =1, ncol = 1),
#   width = 15, height = 10)
# 
# ggsave(
#   filename = "Chunks/Simulations/Results/VAR/fulll dynamic/20230803/TS_plots_MULTIDYN.pdf",
#   plot = marrangeGrob(plots_MULTIDYN, nrow =1, ncol = 1),
#   width = 15, height = 10)

# writexl::write_xlsx(results, "Chunks/Simulations/Results/Factor/20230730/Results_50_20_Factor.xlsx")




# 03 OLD STUFF ----
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






