library(tidyverse)
library(tsm)
library(dynlm)

# Estimate TS-Model from simulated data

rm(list=ls())
graphics.off()

Obs = 300
Win = 50

df_total <- tibble(
  time = ts(1:Obs),
  y = ts(arima.sim(n=Obs, list(order=c(1,0,1), ar = c(.4), ma = -.2), mean = 1, sd = 5) + .6 * seq(1,Obs)),
  x1 = ts(arima.sim(n=Obs, list(order=c(1,0,1), ar = c(.4), ma = -.2), mean = 3, sd = 3.5) + .6 * seq(1,Obs)),
  x2 = ts(arima.sim(n=Obs, list(order=c(1,0,1), ar = c(.2), ma = -0.2), mean = 1, sd = 2) + -0.8 * seq(1,Obs)),
  x3 = ts(arima.sim(n=Obs, list(order=c(2,0,1), ar = c(1, -0.5), ma = 0.5), mean = 1, sd = 5) + -0.5 * seq(1,Obs)))

# df_total <- tibble(
#   time = 1:Obs,
#   y = arima.sim(n=Obs, list(order=c(1,0,1), ar = c(.4), ma = -.2), mean = 1, sd = 5) + .6 * seq(1,Obs),
#   x1 = arima.sim(n=Obs, list(order=c(1,0,0), ar = c(.2)), mean = 1, sd = 2) + 0.8 * seq(1,Obs),
#   x2 = arima.sim(n=Obs, list(order=c(1,0,0), ar = c(.2)), mean = 1, sd = 2) + -0.8 * seq(1,Obs),
#   x3 = arima.sim(n=Obs, list(order=c(2,0,1), ar = c(1, -0.5), ma = 0.5), mean = 1, sd = 5) + -0.5 * seq(1,Obs))

df = df_total[1:250,] %>% 
  mutate(time = ts(time),
         y = ts(y),
         x1 = ts(x1),
         x2 = ts(x2),
         x3 = ts(x3))

matplot(ts(df), 
        type = "l", 
        col = c("steelblue", "darkgreen", "darkred", "orange"), 
        lty = 1, 
        lwd = 2,
        main = "Different Simulations",
        xlab = "Time",
        ylab = "Value")

# 1. ARIMA-Model
plot(df$y)

model_ARIMA = forecast::auto.arima(df$y)
print(summary(model_ARIMA))
# forecast::checkresiduals(model_ARIMA)

plot(forecast::forecast(model_ARIMA,h=50))

fc = tibble(
  ARIMA = forecast::forecast(model_ARIMA,h=50)$mean)

# 2. ADL-Model

  # Order of ADL-Model: Check PACF to get AR-Specification for each Covariate
  
significant_pacf = function(x){
  
  options(warn=-1)
  
  df_pacf = ((pacf(df[[x]], type='correlation', plot = FALSE)$acf / 
                sqrt(1/pacf(df[[x]], type='correlation', plot = FALSE)$n.used)) > 1.96) %>% 
    as.data.frame() %>% 
    pull(c(1)) 
  
  lags = max(which(df_pacf == TRUE))
  options(warn=0)
  return(lags)
  
  }

df_lags = tibble(
  Variable = c("y", paste0("x", 1:3)),
  n_lags = NA)

for (var in c("y", paste0("x", 1:3))) {
  df_lags$n_lags[df_lags$Variable == var] = significant_pacf(var)
  rm(var)
}

  # Run model with all lags in first place, keep only significant lags in second step

model_ADL_prelim = summary(dynlm(y ~ L(y, 1:df_lags$n_lags[df_lags$Variable == "y"]) + 
                L(x1, 0:df_lags$n_lags[df_lags$Variable == "x1"]) +
                L(x2, 0:df_lags$n_lags[df_lags$Variable == "x2"]) + 
                L(x3, 0:df_lags$n_lags[df_lags$Variable == "x3"]),
                data = df))

df_sig_lags = tibble(
  Variable = c(rep(df_lags$Variable[1], df_lags$n_lags[1]),
               rep(df_lags$Variable[2], df_lags$n_lags[2]+1),
               rep(df_lags$Variable[3], df_lags$n_lags[3]+1),
               rep(df_lags$Variable[4], df_lags$n_lags[4]+1)),
  n_lags = c(1:df_lags$n_lags[1],
            0:df_lags$n_lags[2],
            0:df_lags$n_lags[3],
            0:df_lags$n_lags[4]),
  p_value = as.numeric(model_ADL_prelim$coefficients[,4])[2:length(model_ADL_prelim$coefficients[,4])]) %>% 
  filter(p_value <= 0.1)

rm(model_ADL_prelim, df_lags)

y_lag  = df_sig_lags$n_lags[df_sig_lags$Variable == "y"]
x1_lag = df_sig_lags$n_lags[df_sig_lags$Variable == "x1"]
x2_lag = df_sig_lags$n_lags[df_sig_lags$Variable == "x2"]
x3_lag = df_sig_lags$n_lags[df_sig_lags$Variable == "x3"]

reg_formula = as.formula(
  paste("y ~", 
        ifelse(length(y_lag) >= 1, "L(y, c(y_lag))", ""),
        ifelse(length(x1_lag) >= 1, " + ", ""),
        ifelse(length(x1_lag) >= 1, "L(x1, c(x1_lag))", ""),
        ifelse(length(x2_lag) >= 1, " + ", ""),
        ifelse(length(x2_lag) >= 1, "L(x2, c(x2_lag))", ""),
        ifelse(length(x3_lag) >= 1, " + ", ""),
        ifelse(length(x3_lag) >= 1, "L(x3, c(x3_lag))", ""),
        sep = ""))

model_ADL = dynlm(print(reg_formula),
                                 data=df)

rm(reg_formula)

  # Come up with logic to get Variables of df_sig_lags into recursive forecasts procedure.
  # LÖSUNG FÜR MEHRERE Y NOCH FINDEN.
new_df = data.frame(matrix(NA,
                           nrow = Obs,
                           ncol = nrow(df_sig_lags))) %>% 
  mutate(constant = 1) %>% 
  select(constant, everything())

for(i in 1:nrow(df_sig_lags)) {  
  if (df_sig_lags$Variable[i] == "y") {
    new_df[i+1] <- c(df[[df_sig_lags$Variable[i]]], rep(NA, Win))  
    colnames(new_df)[i+1] <- paste0(df_sig_lags$Variable[i]) 
  } else {
    new_df[i+1] <- df_total[[df_sig_lags$Variable[i]]]                
    colnames(new_df)[i+1] <- paste0(df_sig_lags$Variable[i])   
  }
}

# Lagging
for (i in 1:nrow(df_sig_lags)) {
  new_df[[i+1]] = Hmisc::Lag(new_df[[i+1]], df_sig_lags$n_lags[i])
}

for (i in 1:ncol(new_df)) {
  if (colnames(new_df)[i] == "y"){
    new_df[251:300, i] = NA
  }
}

fc = fc %>% 
  mutate(ADL = NA)

for (i in 1:Win) {  
  # new_df$y[250 + i] = coef(model_ADL)  %*% as.numeric(new_df[249 + i,])
  fc$ADL[i] = coef(model_ADL)  %*% as.numeric(new_df[249 + i,])
  if (colnames(new_df)[2] == "y") {
    new_df$y[250+i] = coef(model_ADL)  %*% as.numeric(new_df[249 + i,])
  }
}

# rm(new_df, df_sig_lags, i, Obs, Win, x1_lag, x2_lag, x3_lag, y_lag)

output = tibble(
  time = 1:50,
  y_counter = df_total$y[251:300],
  y_ARIMA = fc$ARIMA,
  y_ADL = fc$ADL)

sqrt(c(crossprod(output$y_ARIMA - output$y_counter)) / 50)
sqrt(c(crossprod(output$y_ADL - output$y_counter)) / 50)

plot = output %>% 
  gather(type, value, y_counter:y_ADL) 

ggplot(plot) +
  aes(x = time, y = value, colour = type) +
  geom_line(size = 1) +
  scale_color_hue(direction = 1) +
  theme_minimal()







# 3. SC-Model

# 4. OLS-Model


model_OLS = lm(y~ x1 + x2 + x3, data = df)
summary(model_OLS)

# 
# # Further Simulations
# # White Noise
# 
# WN <- arima.sim(model = list(order = c(0, 0, 0)), n = 200, mean = 2, sd = 1)
# plot(WN)
# 
# # Moving Average
# 
# MA <- arima.sim(model = list(order = c(0, 0, 1), ma = .9), n = 200, mean = 2, sd = 1)  
# plot(MA)
# 
# # Autoregressive
# 
# AR <- arima.sim(model = list(order = c(2, 0, 0), ar = c(1.74, -.75)), n = 200) 
# plot(AR)
# 
# # ARMA
# 
# ARMA <- arima.sim(model = list(order = c(2, 0, 1), ar = c(1, -.9), ma = .8), n = 250)
# plot(ARMA)
# 
# # ARIMA
# 
# ARIMA <- arima.sim(model = list(order = c(2, 1, 1), ar = c(1, -.9), ma = .8), n = 250)
# plot(ARIMA)
# 
# # Seasonal ARMA's
# 
# library(astsa)
# plot(chicken)
# 
# sarima(chicken, p = 2, d = 1, q = 0, P = 1, D = 0, Q = 0, S = 12)
# sarima.for(chicken, n.ahead = 60, p = 2, d = 1, q = 0, P = 1, D = 0, Q = 0, S = 12)
# 
# library(sarima)
# tsplot(sim_sarima(n=144, model = list(ma=0.8)))
# tsplot(sim_sarima(n=144, model = list(ar=c(rep(0,11),0.95))))  # SAR(1), 12  seasons
# tsplot(sim_sarima(n=144, model = list(ma=c(rep(0,11),0.8))))  # SMA(1)
# tsplot(sim_sarima(n=144,model=list(sar=0.8, nseasons=12, sigma2 = 1)))  # SAR(1), 12 seasons
# 
# # Deterministic Trend
# 
# x = arima.sim(n=50, list(order=c(1,0,1), ar = c(.9), ma = -.2)) + -.8 * seq(1,50)
# plot(x)
# 
# # Stochastic Trend
# 
# x = arima.sim(n=50, list(order=c(1,0,1),ar=c(.9), ma=-.2)) + 2
# plot(cumsum(x))

