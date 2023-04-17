library(tidyverse)

if (Sys.info()[6] == "jctoe"){
  setwd("C:/Promotion/SC_Paper/Documents/sc_sim Ferman/Ferman") 
  } else {
    setwd("~/Diss/Topics/Synthetic Control/Documents/sc_sim Ferman/Ferman") 
}


source("my_functions.R")

# 1. DATA GENERATING PROCESS: FACTOR MODEL WITHOUT COVARIATES ---- 

iter = 100

# Number of donors
J = 35

# Number of pre-and post-treatment periods
T1 = 50
T0 = T1

# AR-Term in Factor model. y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
rho = 0.5

# Intercept. Set it equal to mean*(1-rho) to define mean of process
# alpha = 3
alpha = 0*(1-rho)

# Specify variance of u_t. Set it to (1 - rho^2) will lead to var(\lambda^k_t) = 1. Variance of the factors
var_u = (1-rho^2)

# Specify variance of transitory shocks in Factor model equation. Variance of the error terms 
var_epsilon = 1

# Post-treatment effects. Could be specified differently. Easier to assesss counterfactual with flat effect
post_effect = 10

# Number of factors
K = 2

# Adding a trend
c = 0

# Group distribution of each factor 
group_distribution = list(
  "lambda1" = c(1,0),
  "lambda2" = c(0,1))

Mu = matrix(0, nrow = J+1, ncol = K)

for(k in 1:K) {
  fac_distr = group_distribution[[k]]
  for(pos in 1:length(fac_distr))
    if(pos==1)
      Mu[1:(1 + J %/% length(fac_distr)),k] = fac_distr[pos] else 
        if (pos < length(fac_distr))
          Mu[(2 + (pos - 1)*J %/% length(fac_distr)):(1 + pos * J %/% length(fac_distr)),k] = fac_distr[pos] else
            Mu[(2+(pos-1)*J%/%length(fac_distr)):nrow(Mu),k] = fac_distr[pos]
          
}

results = data.frame(matrix(NA, nrow = iter, ncol = 4)) %>%
  rename(
    RMSPE_SC = c(1),
    RMSPE_OLS = c(2),
    RMSFE_SC = c(3),
    RMSFE_OLS = c(4))

my_func = function(J){
  
  Factors = sapply(1:K, function(k){simulate_ar1(rho = rho, var_shock = var_u, T0 = T0+T1, intercept = alpha)}) 
  
  transitory_shocks = matrix(rnorm((J+1)*(T0+T1), sd = sqrt(var_epsilon)), nrow = (T0+T1), ncol = J+1)
  
  y = Factors%*%t(Mu) + transitory_shocks
  y = y + (1:nrow(y))*c
  
  # Adding effect
  y[(T0+1):(T0+T1),1] = post_effect + y[(T0+1):(T0+T1),1] 
  
  y_pre = y[1:T0,1]
  y_post = y[(T0+1):(T0+T1),1]
  
  x_pre = y[1:T0, -1]
  x_post = y[(T0+1):(T0+T1), -1]
  
  matplot(ts(y),
          type = "l",
          lty = 1,
          lwd = 2,
          main = "Different Simulations",
          xlab = "Time",
          ylab = "Value")
  
  #rm(list = setdiff(ls(), c("y", "y_pre", "y_post", "x_pre", "x_post", "Mu", "T0", "T1")))
  
  # 2. ESTIMATION ----
  
  # SC
  Dmat = t(x_pre) %*% x_pre
  dvec = t(x_pre) %*% y_pre
  Amat = t(rbind(rep(1, ncol(x_pre)), diag(ncol(x_pre)), -1*diag(ncol(x_pre))))
  bvec = c(1, rep(0, ncol(x_pre)), rep(-1,ncol(x_pre)))
  synth_model = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  w_sc = synth_model$solution
  y_sc_pre = x_pre %*% w_sc
  y_sc_post = x_post %*% w_sc
  
  # y_treat_sc = as.data.frame(c(y_pre, y_post)) %>% 
  #   rename(y = c(1))
  # y_treat_sc$y_hat = c(y_ols_pre, y_ols_post)
  # 
  # matplot(ts(y_treat_sc),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "Different Simulations",
  #         xlab = "Time",
  #         ylab = "Value")
  
  results = c()
  results[1] = sqrt(mean((y_pre - y_sc_pre)^2))
  results[2] = sqrt(mean(((y_post-post_effect) - y_sc_post)^2))
  
  # OLS
  
  w_ols = summary(lm(y_pre ~ x_pre))$coefficients[,1]
  y_ols_pre = cbind(rep(1, (T0)), x_pre) %*% w_ols
  y_ols_post = cbind(rep(1, (T1)), x_post) %*% w_ols
  
  # y_treat_ols = as.data.frame(c(y_pre, y_post)) %>%
  #   rename(y = c(1))
  # y_treat_ols$y_hat = c(y_ols_pre, y_ols_post)
  # 
  # matplot(ts(y_treat_ols),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "Different Simulations",
  #         xlab = "Time",
  #         ylab = "Value")
  
  results[3] = sqrt(mean((y_pre - y_ols_pre)^2))
  results[4] = sqrt(mean(((y_post-post_effect) - y_ols_post)^2))
  
  
  return(results)
}

for (i in 1:iter) {
  
  result_prelim = my_func(J)
  
  results$RMSPE_SC[i] = result_prelim[1]
  results$RMSFE_SC[i] = result_prelim[2]
  results$RMSPE_OLS[i] = result_prelim[3]
  results$RMSFE_OLS[i] = result_prelim[4]
  
  rm(result_prelim)
  svMisc::progress(i,iter)
}

# Visualization

df_meta = results %>%
  gather(type, value, RMSPE_SC:RMSFE_OLS) %>% 
  mutate(type = case_when(
    type == "RMSPE_OLS" ~ 1,
    type == "RMSPE_SC" ~ 2,
    type == "RMSFE_OLS" ~ 3,
    type == "RMSFE_SC" ~ 4)) %>% 
  mutate(type = as.factor(type)) %>% 
  mutate(type = recode_factor(type,
                              `1` = "RMSPE_OLS",
                              `2` = "RMSPE_SC",
                              `3` = "RMSFE_OLS",
                              `4` = "RMSFE_SC"))

means <- aggregate(value ~  type, df_meta, mean)
means$value = round(means$value, 3)

ggplot(df_meta) +
  aes(x = type, y = value, fill = type, color = type) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "viridis", direction = 1) +
  stat_summary(fun=mean, colour="darkred", geom="point", 
               shape=18, size=3, show.legend=FALSE) + 
  geom_text(data = means, aes(label = value, y = value + 0.035), color = "black")+
  scale_color_viridis_d(option = "viridis", direction = 1) +
  labs(
    x = "Model: OLS / SC",
    y = "RMSPE / RMSFE") +
  theme_minimal()



