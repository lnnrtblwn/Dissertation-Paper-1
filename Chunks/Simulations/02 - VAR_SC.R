library(tidyverse)
library(CVXR) # for convex optimization in SC

#set.seed(123)

rm(list = ls())
graphics.off()

setwd("~/Diss/Topics/Synthetic Control/Chunks/Simulations")

# Functions for coefficient matrix and SC-Process
source("functions/gmvarkit.R")
source("functions/SC_simulation.R")
rm(list = setdiff(ls(), c("random_coefmats2", "SC_simulation")))

# Generate sample
t = 1000 # half sample. At least 3 times k 
t_full = 2*t # full sample
k = 5 # number of time series
p = 1 # AR specification of VAR-Process
iter =  1000 # iterations per MC-simulation
# series_share = 0 # share of VAR, SC is 1-(share of VAR)

# Start of simulation

MCMC_meta = list()

for (series_share in seq(0,1, by = 1/20)) {
  
  MCMC = data.frame(matrix(NA, nrow = iter, ncol = 10)) %>%
    rename(Iteration = c(1),
           RMSPE_SC = c(2),
           RMSPE_OLS = c(3),
           RMSPE_VAR = c(4),
           RMSPE_VAR_SC = c(5),
           RMSFE_SC = c(6),
           RMSFE_OLS = c(7),
           RMSFE_VAR = c(8),
           RMSFE_VAR_SC = c(9),
           VAR_Share = c(10))
  
  MCMC$VAR_Share = series_share
  
for (iteration in c(1:iter)) {
  # iteration = 1
  share_help = seq(0,1, by = 1/20)
  
  MCMC$Iteration[iteration] = iteration

# Simulate stationary coefficient matrix
A = matrix(c(random_coefmats2(p = 1, d = k, ar_scale = 1)), k) 

# Generate VAR Series
series_VAR = matrix(0, k, t_full + 2*p) # Raw series with zeros
for (i in (p + 1):(t_full + 2 * p)) {
  series_VAR[, i] <- A %*% series_VAR[, i - 1] + rnorm(k, 0, 0.5) # Generate series with e ~ N(0,0.5)
}

# Generate SC_Series
series_SC = t(SC_simulation(t_full+2,k-1))

# Series full to assess forecasts
series_full = t(series_share * series_VAR[, -(1:p)]) + t((1-series_share) * series_SC[, -(1:p)])
series = ts(series_full[1:(t+1),])

rm(series_SC, series_VAR)

#plot.ts(series) 

predicts = tibble(
  y_true = series[2:(t+1),1],
  pred_SC = rep(NA, t),
  pred_OLS = rep(NA, t),
  pred_VAR = rep(NA, t),
  pred_VAR_SC = rep(NA, t))

forecasts = tibble(
  y_true = series_full[(t+2):t_full,1],
  fore_SC = rep(NA, t-1),
  fore_OLS = rep(NA, t-1),
  fore_VAR = rep(NA, t-1),
  fore_VAR_SC = rep(NA, t-1))


####
# ESTIMATION 1: SC
####

  # prediction

df_SC_pred = data.frame(series[2:(t+1),1:k]) 
colnames(df_SC_pred) = c(paste("y"), paste0("x", 1:(k-1)))

# comvex optimization via CVXR package.
b = Variable(ncol(df_SC_pred)-1)
X = cbind(as.matrix(df_SC_pred[,-1]))
obj = Minimize(sum((df_SC_pred$y - X %*% b)^2))
constraints = list(sum(b[1:(k-1)]) == 1,
                   b[1:(k-1)] >= 0)
problem = Problem(obj, constraints)
soln = solve(problem)
bval = soln$getValue(b)

# grid search to confirm minization results. all good.

# check <- do.call(rbind, lapply(0:100, function(i) data.frame(x=i, y=0:(100-i))))
# check$x = check$x/100
# check$y = check$y/100
# check$z <- 1-check$x-check$y
# 
# check$RMSE = NA
# 
# for (i in 1:nrow(check)) {
#   y_hat = t(c(check$x[i], check$y[i], check$z[i])  %*% t(df_SC_pred[,2:4]))
#   check$RMSE[i] = sqrt(c(crossprod(y_hat - df_SC_pred$y)) / t)
#   rm(y_hat)
# }
# 
# min(check$RMSE) < sqrt(c(crossprod(t(c(bval) %*% t(df_SC_pred[,2:4])) - df_SC_pred$y)) / t)

rm(b, constraints, obj, problem, soln, X, i)
predicts$pred_SC = c(c(bval)  %*% t(df_SC_pred[,2:k])) 

  # forecast

df_SC_fore = data.frame(series_full[(t+2):t_full,1:k]) 
colnames(df_SC_fore) = c(paste("y"), paste0("x", 1:(k-1)))


forecasts$fore_SC = c(c(bval)  %*% t(df_SC_fore[,2:k])) 

rm(bval)

####
# ESTIMATION 2: OLS
####

  # prediction

df_OLS_pred = df_SC_pred
df_OLS_fore = df_SC_fore
rm(df_SC_pred, df_SC_fore)

model_OLS = lm(y ~ ., data = df_OLS_pred)
# summary(model_OLS)
predicts$pred_OLS = model_OLS$fitted.values

  # forecast

forecasts$fore_OLS = c(model_OLS$coefficients %*% t(cbind(1, df_OLS_fore[,2:k])))

rm(df_OLS_fore, df_OLS_pred, model_OLS)

####
# ESTIMATION 3: VAR
####

  # prediction

  # VAR Estimation
# var <- vars::VAR(series, 1, type = "none")
# var$varresult
# summary(var)
# rm(var)

  # manually using OLS
# X = dplyr::lag(as.data.frame(as.matrix(series))) %>%
#   slice(2:(t+1)) %>%
#   as.matrix()
# 
# y = as.matrix(series[2:(t+1),])
# 
# solve(t(X) %*% X) %*%  t(X) %*% y

  # using lm()

X = dplyr::lag(as.data.frame(as.matrix(series))) %>%
  slice(2:(t+1)) %>%
  as.matrix()

y = as.matrix(series[2:(t+1),])

df_VAR_pred = data.frame(y[,1], X) %>% 
  rename(y = c(1))

model_VAR = lm(y ~ ., data = df_VAR_pred)
# summary(model_VAR)
predicts$pred_VAR = model_VAR$fitted.values

  # forecast

X = dplyr::lag(as.data.frame(as.matrix(series_full))) %>%
  slice((t+2):t_full) %>%
  as.matrix()

df_VAR_fore = data.frame(X)
df_VAR_fore[(2:t),1] = NA

for (i in 2:nrow(df_VAR_fore)) {
  df_VAR_fore[i,1] = c(model_VAR$coefficients %*% t(cbind(1, df_VAR_fore)[(i-1),]))
}

forecasts$fore_VAR = df_VAR_fore[2:t,1]

rm(X, df_VAR_pred, df_VAR_fore, y, model_VAR, i)

####
# ESTIMATION 1: VAR_SC
####

  # Estimation. Similar to OLS 
  # y = constant + current donors values (d) + all lagged values (a).

  # prediction
d = as.data.frame(series[,-1]) %>% 
  slice(2:(t+1)) %>%
  as.matrix()

a = dplyr::lag(as.data.frame(series)) %>% 
  slice(2:(t+1)) %>%
  as.matrix()

X = as.matrix(cbind(d, a))
y = as.matrix(series[2:(t+1),])

df_VAR_SC_pred = data.frame(y[,1], X) %>% 
  rename(y = c(1))

model_VAR_SC = lm(y ~ ., data = df_VAR_SC_pred)
# summary(model_VAR_SC)
predicts$pred_VAR_SC = model_VAR_SC$fitted.values


  # forecast
d = as.data.frame(series_full[,-1]) %>% 
  slice((t+2):(t_full+1)) %>%
  as.matrix()

a = dplyr::lag(as.data.frame(series_full)) %>% 
  slice((t+2):(t_full+1)) %>%
  as.matrix()

X = as.matrix(cbind(1, d, a))

df_VAR_SC_fore = data.frame(X) 
df_VAR_SC_fore[2:t,(k+1)] = NA

for (i in 2:nrow(df_VAR_SC_fore)) {
  df_VAR_SC_fore[i,(k+1)] = c(model_VAR_SC$coefficients %*% t(df_VAR_SC_fore[(i-1),]))
}


forecasts$fore_VAR_SC = df_VAR_SC_fore[2:t,(k+1)]
rm(df_VAR_SC_fore, df_VAR_SC_pred, a, d, model_VAR_SC, X, y)


MCMC$RMSPE_SC[iteration] = sqrt(c(crossprod(predicts$y_true - predicts$pred_SC)) / nrow(predicts))
MCMC$RMSPE_OLS[iteration] = sqrt(c(crossprod(predicts$y_true - predicts$pred_OLS)) / nrow(predicts))
MCMC$RMSPE_VAR[iteration] = sqrt(c(crossprod(predicts$y_true - predicts$pred_VAR)) / nrow(predicts))
MCMC$RMSPE_VAR_SC[iteration] = sqrt(c(crossprod(predicts$y_true - predicts$pred_VAR_SC)) / nrow(predicts))

MCMC$RMSFE_SC[iteration] = sqrt(c(crossprod(forecasts$y_true - forecasts$fore_SC)) / nrow(forecasts))
MCMC$RMSFE_OLS[iteration] = sqrt(c(crossprod(forecasts$y_true - forecasts$fore_OLS)) / nrow(forecasts))
MCMC$RMSFE_VAR[iteration] = sqrt(c(crossprod(forecasts$y_true - forecasts$fore_VAR)) / nrow(forecasts))
MCMC$RMSFE_VAR_SC[iteration] = sqrt(c(crossprod(forecasts$y_true - forecasts$fore_VAR_SC)) / nrow(forecasts))

MCMC_meta[[which(share_help == series_share)]] = MCMC
}
}

# Results of Simulation
df_meta = data.frame(matrix(NA, nrow = 0, ncol = 5)) %>% 
  rename(VAR_Share = c(1),
         RMSFE_SC = c(2),
         RMSFE_OLS = c(3),
         RMSFE_VAR = c(4),
         RMSFE_VAR_SC = c(5))

for (i in 1:length(share_help)) {
  name = paste("data")
  assign(
    name,
    data.frame(matrix(NA, nrow = iter, ncol = 5)) %>%
      rename(
        VAR_Share = c(1),
        RMSFE_SC = c(2),
        RMSFE_OLS = c(3),
        RMSFE_VAR = c(4),
        RMSFE_VAR_SC = c(5)))
  
  data$VAR_Share = MCMC_meta[[i]]$VAR_Share
  data$RMSFE_SC = MCMC_meta[[i]]$RMSFE_SC
  data$RMSFE_OLS = MCMC_meta[[i]]$RMSFE_OLS
  data$RMSFE_VAR = MCMC_meta[[i]]$RMSFE_VAR
  data$RMSFE_VAR_SC = MCMC_meta[[i]]$RMSFE_VAR_SC
  
  df_meta = df_meta %>%
    bind_rows(data)
  rm(data)
}

df_meta = df_meta %>% 
  gather(type, value, RMSFE_SC:RMSFE_VAR_SC, -VAR_Share) %>% 
  mutate(VAR_Share = as.factor(round(VAR_Share,2)))

ggplot(df_meta) +
  aes(x = VAR_Share, y = value, fill = type, colour = type) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "viridis", direction = 1) +
  scale_color_viridis_d(option = "viridis", direction = 1) +
  theme_minimal() +
  facet_wrap(vars(type))

save(MCMC_meta ,file="~/Diss/Topics/Synthetic Control/Chunks/Simulations/MC_meta.Rda")



# MCMC_pred = MCMC %>% 
#   select(-c(Iteration,RMSFE_SC:RMSFE_VAR_SC)) %>% 
#   gather(type, value, RMSPE_SC:RMSPE_VAR_SC)
# 
# p_pred = ggplot(MCMC_pred) +
#   aes(x = value, fill = type) +
#   geom_density(adjust = 1L, alpha = 0.5) +
#   scale_fill_viridis_d(option = "viridis", direction = 1)+
#   labs(title = "RMSPE")+
#   theme_minimal()
# 
# MCMC_fore = MCMC %>% 
#   select(-c(Iteration:RMSPE_VAR_SC)) %>% 
#   gather(type, value, RMSFE_SC:RMSFE_VAR_SC)
# 
# p_fore = ggplot(MCMC_fore) +
#   aes(x = value, fill = type) +
#   geom_density(adjust = 1L, alpha = 0.5) +
#   scale_fill_viridis_d(option = "viridis", direction = 1)+
#   labs(title = "RMSFE")+
#   theme_minimal()
# 
# ggpubr::ggarrange(p_pred, p_fore, ncol=2)
# 
# summary(MCMC)

# Plots for RMSPE and RMSFE
# test = predicts %>% 
#   mutate(time = 1:n()) %>% 
#   gather(type, value, y_true:pred_VAR_SC)
# 
# ggplot(test) +
#   aes(x = time, y = value, colour = type) +
#   geom_line(size = 1) +
#   scale_color_hue(direction = 1) +
#   theme_minimal()
# 
# test = forecasts %>% 
#   mutate(time = 1:n()) %>% 
#   gather(type, value, y_true:fore_VAR_SC)
# 
# ggplot(test) +
#   aes(x = time, y = value, colour = type) +
#   geom_line(size = 1) +
#   scale_color_hue(direction = 1) +
#   theme_minimal()



