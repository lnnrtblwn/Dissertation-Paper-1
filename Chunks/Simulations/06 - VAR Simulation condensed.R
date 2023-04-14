library(quadprog)
library(MASS)
library(tidyverse)

# mu = c(1,1,1)
# A =  matrix(c(1, 0.1, 0.4, 0.1, 1, 0.5, 0.4, 0.5, 1), byrow = T, nrow = 3)
#Simulate stationary coefficient matrix

iter = 100
obs = 1000
my_var = 1.5

results_mean = data.frame(matrix(NA, nrow = 46, ncol = 5)) %>%
  rename(
    Donors = c(5),
    RMSPE_SC = c(1),
    RMSPE_OLS = c(2),
    RMSFE_SC = c(3),
    RMSFE_OLS = c(4))

id = 1

for (k in seq(5,50, by = 1)) {
  
  results = data.frame(matrix(NA, nrow = iter, ncol = 4)) %>%
    rename(
      RMSPE_SC = c(1),
      RMSPE_OLS = c(2),
      RMSFE_SC = c(3),
      RMSFE_OLS = c(4))
  
  sc_fun <- function(obs, mu, sig) {
    
    # Preparing output data frame
    result_prelim = data.frame(matrix(NA, nrow = 1, ncol = 4)) %>%
      rename(
        RMSPE_SC = c(1),
        RMSPE_OLS = c(2),
        RMSFE_SC = c(3),
        RMSFE_OLS = c(4))
    
    n = k
    mu = c(rep(1, times = n))
    A = matrix(runif(n*n), nrow = n, ncol = n)
    A = t(A) %*% A
    diag(A) = n/my_var
    
    # Generating, de-mean data and splitting it
    df = MASS::mvrnorm(n = obs, mu = mu, Sigma = A, tol = T) 
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
  
  for (i in 1:iter) {
    
    result_prelim = sc_fun(obs = obs, mu = mu, sig = A)
    
    results$RMSPE_SC[i] = result_prelim$RMSPE_SC[1]
    results$RMSFE_SC[i] = result_prelim$RMSFE_SC[1]
    results$RMSPE_OLS[i] = result_prelim$RMSPE_OLS[1]
    results$RMSFE_OLS[i] = result_prelim$RMSFE_OLS[1]
    
    rm(result_prelim)
  }
  
  results_mean$Donors[id] = k
  results_mean$RMSPE_SC[id] = mean(results$RMSPE_SC)
  results_mean$RMSPE_OLS[id] = mean(results$RMSPE_OLS)
  results_mean$RMSFE_SC[id] = mean(results$RMSFE_SC)
  results_mean$RMSFE_OLS[id] = mean(results$RMSFE_OLS)

  id = id + 1
  svMisc::progress(k,50)
}

rm(id)

# Visualization

df_meta = results_mean %>%
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

p = df_meta %>%
  filter(type %in% c("RMSFE_SC", "RMSFE_OLS")) %>%
  ggplot() +
  aes(x = Donors, y = value, colour = type) +
  geom_line(size = 0.7) +
  scale_color_manual(
    values = c(RMSFE_OLS = "#E41A1C",
               RMSFE_SC = "#2B7B50")
  ) +
  labs(
    x = "Number of Donors",
    y = "RMSFE",
    title = "VAR-Data // 1,000 Obs with 100 Simulations",
    subtitle = paste0("Variance of Donors = frac{n_donors}{",my_var,"}")
  ) +
  theme_minimal()

setwd("C:/Users/lbolwin/Documents/Diss/Topics/Synthetic Control/Chunks/Simulations/Plots")
png(paste0("plot",my_var,".png"), width = 7.83, height = 4.96, units = 'in', res = 600)
p
dev.off()
