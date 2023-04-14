library(quadprog)
library(MASS)
library(tidyverse)

# mu = c(1,1,1)
# A =  matrix(c(1, 0.1, 0.4, 0.1, 1, 0.5, 0.4, 0.5, 1), byrow = T, nrow = 3)
#Simulate stationary coefficient matrix

k = 20
mu = c(rep(1, times = k))
a = matrix(runif(k^2)*2-1, ncol = k)
A = crossprod(a, a)
diag(A) = diag(A)*20
  
iter = 50
obs = 500

sc_fun <- function(obs, mu, sig) {
  
  # Preparing output data frame
  result_prelim = data.frame(matrix(NA, nrow = 1, ncol = 4)) %>%
    rename(
      RMSPE_SC = c(1),
      RMSPE_OLS = c(2),
      RMSFE_SC = c(3),
      RMSFE_OLS = c(4))
  
  # Generating, de-mean data and splitting it
  df = MASS::mvrnorm(n = obs, mu = mu, Sigma = A, ) 
  # df[,1] = df[,1] - mean(df[,1])
  # df[,2] = df[,2] - mean(df[,2])
  # df[,3] = df[,3] - mean(df[,3])
  df_pre = df[1:(obs/2),]
  df_post = df[((obs/2)+1):obs,]
  X_pre = df_pre[,-1]
  X_post = df_post[,-1]
  y_pre = df_pre[,1]
  y_post = df_post[,1]
  
  # OLS estimation
  w_ols = summary(lm(y_pre ~ X_pre))$coefficients[,1]
  y_ols_pre = cbind(rep(1, (obs/2)), X_pre) %*% w_ols
  y_ols_post = cbind(rep(1, (obs/2)), X_post) %*% w_ols
  # summary(lm(y_pre ~ X_pre))
  # w_ols = solve(A[2:ncol(A), 2:ncol(A)]) %*% A[2:ncol(A),1]
  # w_ols = as.matrix(c(w_ols, 1-w_ols[1,1] - w_ols[2,1]))
  # y_ols_pre = cbind(X_pre, rep(1, (obs/2))) %*% w_ols
  # y_ols_post = cbind(X_post, rep(1, (obs/2))) %*% w_ols
  result_prelim$RMSPE_OLS = sqrt(mean((y_pre - y_ols_pre)^2))
  result_prelim$RMSFE_OLS = sqrt(mean((y_post - y_ols_post)^2))
  
  # SC estimation
  Dmat = t(X_pre) %*% X_pre
  dvec = t(X_pre) %*% y_pre
  Amat = t(rbind(rep(1, ncol(X_pre)), diag(ncol(X_pre)), -1*diag(ncol(X_pre))))
  bvec = c(1, rep(0, ncol(X_pre)), rep(-1,ncol(X_pre)))
  synth_model = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  w_sc = synth_model$solution
  y_sc_pre = X_pre %*% w_sc
  y_sc_post = X_post %*% w_sc
  result_prelim$RMSPE_SC = sqrt(mean((y_pre - y_sc_pre)^2))
  result_prelim$RMSFE_SC = sqrt(mean((y_post - y_sc_post)^2))
  
  return(result_prelim)
}

results = data.frame(matrix(NA, nrow = iter, ncol = 4)) %>%
  rename(
    RMSPE_SC = c(1),
    RMSPE_OLS = c(2),
    RMSFE_SC = c(3),
    RMSFE_OLS = c(4))

for (i in 1:iter) {
  
  result_prelim = sc_fun(obs = obs, mu = mu, sig = A)
  
  results$RMSPE_SC[i] = result_prelim$RMSPE_SC[1]
  results$RMSFE_SC[i] = result_prelim$RMSFE_SC[1]
  results$RMSPE_OLS[i] = result_prelim$RMSPE_OLS[1]
  results$RMSFE_OLS[i] = result_prelim$RMSFE_OLS[1]
  
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

  