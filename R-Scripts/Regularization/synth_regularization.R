# PACKAGE INSTALLATION ####

suppressPackageStartupMessages({
  list.of.packages = c("tidyverse")
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages, require, character.only = T)
})

gen_data = function(obs, k){
  n = k
  mu = c(rep(1, times = n))
  A = matrix(runif(n * n), nrow = n, ncol = n)
  A = t(A) %*% A
  diag(A) = my_var
  
  df = MASS::mvrnorm(n = obs,
                     mu = mu,
                     Sigma = A,
                     tol = T) %>% 
    as.data.frame()
  
  colnames(df) = c("y", paste0("x", 1:(k-1)))
  
  return(df)
}

obs = 200
k = 100
my_var = k * (200/5)

df = gen_data(obs, k)

ggplot(df) +
  aes(x = x6, y = y) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  theme_minimal()

df_test = df %>% 
  slice(1:(obs/2))

df_train = df %>% 
  slice((obs/2)+1:obs)

# Function to formalize regularized SC-approach
# x, y as data frame
# l1, l2 as numerics
# intercept as boolean

regularized_synth <- function(x, y, l1, l2, intercept) {
  
  y = y %>% 
    as.matrix()
  
  if (intercept == TRUE){
    
    I = diag(k)
    U = matrix(data = 1, nrow = k, ncol = 1)
    
    x = x %>% 
      mutate(x0 = 1) %>% 
      select(x0, everything()) %>% 
      as.matrix()
    
    I[1,1] = 0
    U[1,1] = 0
    
  } else {
    
    I = diag(k-1)
    U = matrix(data = 1, nrow = k-1, ncol = 1)
    
    x = x %>% 
      as.matrix()
    
  }
  
  
  beta = solve(t(x) %*% x + l1 * I + l2 * (U %*% t(U))) %*% ((t(x) %*% y) + l2 * U)
  y_hat = x %*% beta
  
  synth_out = list()
  
  synth_out$y_hat = y_hat
  synth_out$beta = beta
  
  precision = c(sqrt(mean((y - y_hat)^2)),
                mean(abs(y - y_hat)))
  names(precision) = c("RMSPE", "MAPE")
  
  synth_out$precision = precision
  
  return(synth_out)
}

# hyperparameter grid search
# can experiment with: my_var (determines var-cov of data frame, n/k (number of obs/covariates)
# above settings impact paramter grid

param_grid <- expand.grid(
  l1 = seq(0,1/4*10^4, by = 10),
  #l2 = seq(0,10^5, by = 5000), # commented out for sake of simplicity
  intercept = c(TRUE, FALSE))

param_grid$coeff_sum = NA
param_grid$RMSPE = NA
param_grid$MAPE = NA
param_grid$RMSFE = NA
param_grid$MAFE = NA

for (i in 1:nrow(param_grid)) {
  
  current_params = param_grid[i, ]
  
  model = regularized_synth(
    x = df_train[,2:k] %>% 
      scale(center = TRUE, scale = TRUE) %>% 
      as.data.frame(),
    y = df_train[,1],
    l1 = current_params[["l1"]],
    # l2 = current_params[["l2"]] # same as above. If commented in, have to compute n^2 instead of n loops
    l2 = current_params[["l1"]],
    intercept = current_params[["intercept"]])
  
  param_grid$coeff_sum[i] = sum(model$beta[rownames(model$beta) != "x0"])
  
  # test performance
  
  param_grid$RMSPE[i] = model$precision["RMSPE"]
  param_grid$MAPE[i] = model$precision["MAPE"]
  
  # train performance
  
  y_test = df_test[,1] %>% 
    as.matrix()
  x_test = df_test[,2:k] %>% 
    scale(center = TRUE, scale = TRUE) %>% 
    as.matrix()
  
  if ("x0" %in% rownames(model$beta)) {
    y_hat = as.matrix(cbind(1, x_test)) %*% model$beta
    } else {
      y_hat = x_test %*% model$beta
    }
  
  param_grid$RMSFE[i] = sqrt(mean((y_test - y_hat)^2))
  param_grid$MAFE[i] =  mean(abs(y_test - y_hat))
  
  svMisc::progress(i, nrow(param_grid))
  
  }

# Visualizing

df_meta = param_grid %>%
  filter(intercept == TRUE) %>% 
  gather(type, value, RMSPE:MAFE) %>% 
  select(l1, type, value) %>% 
  mutate(type = case_when(
    type == "RMSPE" ~ 1,
    type == "MAPE" ~ 2,
    type == "RMSFE" ~ 3,
    type == "MAFE" ~ 4)) %>% 
  mutate(type = as.factor(type)) %>% 
  mutate(type = recode_factor(type,
                              `1` = "RMSPE",
                              `2` = "MAPE",
                              `3` = "RMSFE",
                              `4` = "MAFE")) %>% 
  filter(type %in% c("RMSPE", "RMSFE"),
         l1 != 0) 


df_meta %>%
  ggplot() +
  aes(x = l1, y = value, colour = type) +
  geom_line() +
  geom_vline(xintercept  = df_meta$l1[df_meta$value == min(df_meta$value[df_meta$type == "RMSFE"])],
             colour="cyan", linetype = "longdash", size = 1.5)+
  scale_color_hue(direction = 1) + 
  theme_minimal()


