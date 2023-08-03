library(tidyverse)
library(dynamac)
library(forecast)
library(tseries)
library(urca)
library(TSstudio)
library(dLagM)
library(vars)

# Simulations: Random walk plus drift is simulated by the cumulative sum of a random variable

N = 250
delta = 0.99 # drift parameter. The higher, the steeper and more linear the process
variance = 0.5 # noise
treat = 100 # size of treatment

# RW1: Simulate for y-variable observed and counterfactual value
RW1 <- function(N, delta, variance, treat) {
  z = cumsum(rnorm(n=N/2, mean=0, sd=sqrt(variance)))
  t = 1:(N/2)
  y1 = t*delta + z
  
  t = ((N/2)+1) :N
  y2 = t*delta + z + treat
  y3 = t*delta + z
  y_treat = c(y1, y2)
  y_counter = c(y1, y3)
  
  df_treat = 1:N %>% 
    as.data.frame() %>% 
    rename(time = c(1)) %>% 
    mutate(y_treat = y_treat,
           y_counter = y_counter)
  
  return(df_treat)}

# RW2: Simulate for x-variables (donors) observed value
RW2 <- function(N, delta, variance) {
  z = cumsum(rnorm(n=N, mean=0, sd=sqrt(variance)))
  t = 1:N
  x = t*delta + z
  return(x)}

df_treat = RW1(N = N, delta = delta, variance = variance, treat = treat) %>% 
  gather(type, value, y_treat:y_counter)

ggplot(df_treat) +
  aes(x = time, y = value, linetype = type) +
  geom_line(size = 1) +
  scale_color_hue(direction = 1) +
  scale_linetype_manual(values=c("dotted", "solid")) +
  geom_vline(xintercept = 125, linetype="dotted", colour = "red", size = 1) +
  theme_minimal()

# adding covariates

df_long = df_treat %>% 
  spread(key = type, value = value) %>% 
  bind_cols(replicate(5, RW2(N = N, delta = delta, variance = variance))) %>% 
  rename_with(.fn = ~paste0("x", 1:5), 
              .cols = c(4:8)) %>% 
  gather(type, value, y_counter:x5)

df_wide = df_long %>% 
  spread(key = type, value = value)
  

# ggplot(df_long) +
#   aes(x = time, y = value, colour = type, linetype = type) +
#   geom_line(size = 1) +
#   scale_color_hue(direction = 1) +
#   scale_linetype_manual(values=c("dotted","dotted","dotted","dotted","dotted",
#                                  "solid", "dashed")) +
#   geom_vline(xintercept = 125, linetype="dotted", colour = "red", size = 1) +
#   theme_minimal()

# Estimating the counterfactual (= model forecast after treatment)

# OLS

model_OLS = lm(y_counter ~ x1 + x2 + x3 + x4 + x5, data = df_wide %>% 
                    filter(time %in% c(1:(N/2))))

sqrt(c(crossprod(model_OLS$residuals)) / length(model_OLS$residuals)) # RMSPE
sqrt(c(crossprod(predict(model_OLS, df_wide %>% 
                      filter(time %in% c((1 + N/2):N))) - df_wide$y_counter[(1 + N/2):N])) / length(c(1:(N/2)))) # RMSFE

test = df_wide %>%
  mutate(y_forecast_OLS = c(rep(NA, times = N/2), predict(model_OLS, df_wide %>%
                                    filter(time %in% c((1 + N/2):N))))) %>%
  gather(type, value, x1:y_forecast_OLS)

ggplot(test) +
  aes(x = time, y = value, colour = type, linetype = type) +
  geom_line(size = 1) +
  scale_color_hue(direction = 1) +
  scale_linetype_manual(values=c("dotted","dotted","dotted","dotted","dotted",
                                 "solid", "solid", "dashed")) +
  geom_vline(xintercept = 125, linetype="dotted", colour = "red", size = 1) +
  theme_minimal()

# ARDL Model. Here only quick and dirty. Actual fitting process is more complicated (see ARDL script)

model_ARDL = dLagM::ardlDlm(formula = y_counter ~ x1 + x2 + x3 + x4 + x5,
                             data = df_wide %>%
                               filter(time %in% c(1:(N/2))),
                             p = 9, q = 9)
summary(model_ARDL)

# come up with automatic way to get rid of insignificant covariates.

vars_remove = summary(model_ARDL)$coefficients[,4] %>% 
  as.data.frame() %>% 
  rename(p_value = c(1)) %>% 
  filter(p_value > 0.1)

vars_remove$vars = rownames(vars_remove)
rownames(vars_remove) = NULL

vars_remove = vars_remove %>% 
  mutate(count = str_count(vars)) %>% 
  mutate(v1 = ifelse(count == 4, substr(vars, 1, 2), substr(vars, 1, 1))) %>% 
  mutate(v2 = ifelse(count == 4, substr(vars, 4, 4), substr(vars, 11, 11))) %>% 
  slice(-1) %>% 
  mutate(v2 = as.numeric(ifelse(v2 == "t", 0, v2))) %>%
  filter(v2 != 0) %>% 
  dplyr::select(v1, v2)

for (i in c("y", "x1", "x2", "x3", "x4", "x5")) {
  name = paste("rm", i, sep = "_")
  assign(name, vars_remove %>% 
           filter(v1 == i) %>% 
           dplyr::select(v2) %>% 
           pull())
}

my.remove = list(p = list(x1 = rm_x1,
                          x2 = rm_x2,
                          x3 = rm_x3,
                          x4 = rm_x4,
                          x5 = rm_x5),
                 q = c(y_counter = rm_y))

model_ARDL = dLagM::ardlDlm(formula = y_counter ~ x1 + x2 + x3 + x4 + x5,
                             data = df_wide %>%
                               filter(time %in% c(1:(N/2))),
                             p = 9, q = 9,
                             remove = my.remove)
summary(model_ARDL)

sqrt(c(crossprod(model_ARDL$model$residuals)) / length(model_ARDL$model$residuals))

# forecast recursively:
x = df_wide %>%
   filter(time > 125) %>%
   dplyr::select(x1:x5) %>%
   t() %>%
   as.matrix()

fc = forecast(model_ARDL, x, h = 125)

sqrt(c(crossprod(fc$forecasts - df_wide$y_counter[(1 + N/2):N])) / length(c(1:(N/2)))) # RMSFE

test = df_wide %>%
   mutate(y_forecast_ARDL = c(rep(NA, times = (N/2)), fc$forecasts)) %>%
   gather(type, value, x1:y_forecast_ARDL)

ggplot(test) +
  aes(x = time, y = value, colour = type, linetype = type) +
  geom_line(size = 1) +
  scale_color_hue(direction = 1) +
  scale_linetype_manual(values=c("dotted","dotted","dotted","dotted","dotted",
                                 "solid", "solid", "dashed")) +
  geom_vline(xintercept = 125, linetype="dotted", colour = "red", size = 1) +
  theme_minimal()
                            
# MCMC simulation. analog zu bisherigen Berechnungen aber jetzt im Loop für mehrere Datenziehungen 

MCMC = data.frame(matrix(NA, nrow = 50, ncol = 5)) %>% 
  rename(Iteration = c(1),
         RMSPE_OLS = c(2),
         RMSFE_OLS = c(3),
         RMSPE_ARDL = c(4),
         RMSFE_ARDL = c(5))


for (iteration in c(1:50)) {
  
  MCMC$Iteration[iteration] = iteration
  
  df_treat = RW1(N = N, delta = delta, variance = variance, treat = treat) %>% 
    gather(type, value, y_treat:y_counter)
  
  df_long = df_treat %>% 
    spread(key = type, value = value) %>% 
    bind_cols(replicate(5, RW2(N = N, delta = delta, variance = variance))) %>% 
    rename_with(.fn = ~paste0("x", 1:5), 
                .cols = c(4:8)) %>% 
    gather(type, value, y_counter:x5)
  
  df_wide = df_long %>% 
    spread(key = type, value = value)
  
  model_OLS = lm(y_counter ~ x1 + x2 + x3 + x4 + x5, data = df_wide %>% 
                   filter(time %in% c(1:(N/2))))
  
  model_ARDL = dLagM::ardlDlm(formula = y_counter ~ x1 + x2 + x3 + x4 + x5,
                              data = df_wide %>%
                                filter(time %in% c(1:(N/2))),
                              p = 9, q = 9)
  
  # Automatic vars remove
  vars_remove = summary(model_ARDL)$coefficients[,4] %>% 
    as.data.frame() %>% 
    rename(p_value = c(1)) %>% 
    filter(p_value > 0.1)
  
  vars_remove$vars = rownames(vars_remove)
  rownames(vars_remove) = NULL
  
  vars_remove = vars_remove %>% 
    mutate(count = str_count(vars)) %>% 
    mutate(v1 = ifelse(count == 4, substr(vars, 1, 2), substr(vars, 1, 1))) %>% 
    mutate(v2 = ifelse(count == 4, substr(vars, 4, 4), substr(vars, 11, 11))) %>% 
    slice(-1) %>% 
    mutate(v2 = as.numeric(ifelse(v2 == "t", 0, v2))) %>%
    filter(v2 != 0) %>% 
    dplyr::select(v1, v2)
  
  for (i in c("y", "x1", "x2", "x3", "x4", "x5")) {
    name = paste("rm", i, sep = "_")
    assign(name, vars_remove %>% 
             filter(v1 == i) %>% 
             dplyr::select(v2) %>% 
             pull())
  }
  
  my.remove = list(p = list(x1 = rm_x1,
                            x2 = rm_x2,
                            x3 = rm_x3,
                            x4 = rm_x4,
                            x5 = rm_x5),
                   q = c(y_counter = rm_y))
  
  model_ARDL = dLagM::ardlDlm(formula = y_counter ~ x1 + x2 + x3 + x4 + x5,
                              data = df_wide %>%
                                filter(time %in% c(1:(N/2))),
                              p = 9, q = 9,
                              remove = my.remove)
  
  MCMC$RMSPE_OLS[iteration] = sqrt(c(crossprod(model_OLS$residuals)) / length(model_OLS$residuals))
  MCMC$RMSFE_OLS[iteration] = sqrt(c(crossprod(predict(model_OLS, df_wide %>% 
                                                         filter(time %in% c((1 + N/2):N))) - df_wide$y_counter[(1 + N/2):N])) / length(c(1:(N/2))))
  
  MCMC$RMSPE_ARDL[iteration] = sqrt(c(crossprod(model_ARDL$model$residuals)) / length(model_ARDL$model$residuals))
  
  x = df_wide %>%
    filter(time > 125) %>%
    dplyr::select(x1:x5) %>%
    t() %>%
    as.matrix()
  
  fc = forecast(model_ARDL, x, h = 125)
  
  MCMC$RMSFE_ARDL[iteration] = sqrt(c(crossprod(fc$forecasts - df_wide$y_counter[(1 + N/2):N])) / length(c(1:(N/2))))
}

mean(MCMC$RMSPE_OLS)
mean(MCMC$RMSPE_ARDL)
mean(MCMC$RMSFE_OLS)
mean(MCMC$RMSFE_ARDL)

rm(MCMC, model_ARDL, model_OLS)

