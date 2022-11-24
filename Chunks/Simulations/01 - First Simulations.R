library(tidyverse)
library(dLagM)

# Estimate TS-Model from simulated data

rm(list = setdiff(ls(), "MCMC"))
graphics.off()

MCMC = data.frame(matrix(NA, nrow = 100, ncol = 7)) %>%
  rename(Iteration = c(1),
         RMSPE_OLS = c(2),
         RMSPE_ADL = c(3),
         RMSPE_ARIMA = c(4),
         RMSFE_OLS = c(5),
         RMSFE_ADL = c(6),
         RMSFE_ARIMA = c(7))

for (iteration in c(1:100)) {
  MCMC$Iteration[iteration] = iteration
  
  Obs = 300
  Win = 50
  
  
  df_total <- tibble(
    time = ts(1:Obs),
    y = ts(arima.sim(n = Obs, list(order = c(1, 0, 1), ar = c(.8), ma = 0.2), mean = 1, sd = 5) + 1.015^seq(1, Obs)),
    x1 = ts(arima.sim(n = Obs, list(order = c(1, 0, 0), ar = c(.8)), mean = -5, sd = 5) + 0.6 * y + rnorm(300, sd = 5)),
    x2 = 0.1 * y + rnorm(300, sd = 1),
    x3 = ts(arima.sim(n = Obs, list(order = c(1, 0, 1), ar = c(-0.5), ma = 0.5), mean = 1, sd = 10) +-0.5 * seq(1, Obs)))
  
  # linear time trend: exponent in seq  = 1
  # concave time trend: exponent in seq < 1
  # convex time trend: exponent in seq > 1
  # alternating time trend: -x^seq(), e.g. (0.6 * (-1.02)^seq(1,Obs))
  
  df = df_total[1:250, ] %>%
    mutate(
      time = ts(time),
      y = ts(y),
      x1 = ts(x1),
      x2 = ts(x2),
      x3 = ts(x3)
    )
  
  matplot(ts(df[-1]),
          type = "l",
          col = c("steelblue", "darkgreen", "darkred", "orange"),
          lty = 1,
          lwd = 2,
          main = "Different Simulations",
          xlab = "Time",
          ylab = "Value")
  
  # 1. ARIMA-Model
  # plot(df$y)
  
  model_ARIMA = forecast::auto.arima(df$y)
  # print(summary(model_ARIMA))
  # forecast::checkresiduals(model_ARIMA)
  
  plot(forecast::forecast(model_ARIMA,h=50))
  
  fc = tibble(ARIMA = forecast::forecast(model_ARIMA, h = 50)$mean)
  
  # 2. ADL-Model
  
  # Order of ADL-Model: Check PACF to get AR-Specification for each Covariate
  
  significant_pacf = function(x) {
    options(warn = -1)
    
    df_pacf = ((
      pacf(df[[x]], type = 'correlation', plot = FALSE)$acf /
        sqrt(1 / pacf(
          df[[x]], type = 'correlation', plot = FALSE
        )$n.used)
    ) > 1.96) %>%
      as.data.frame() %>%
      pull(c(1))
    
    lags = max(which(df_pacf == TRUE))
    options(warn = 0)
    return(lags)
    
  }
  
  df_lags = tibble(Variable = c("y", paste0("x", 1:3)),
                   n_lags = NA)
  
  for (var in c("y", paste0("x", 1:3))) {
    df_lags$n_lags[df_lags$Variable == var] = significant_pacf(var)
    rm(var)
  }
  
  df_lags = df_lags %>%
    mutate(n_lags = ifelse(n_lags > 5, 5, n_lags)) # truncate at 10 lags.
  
  # Run model with all lags in first place, keep only significant lags in second step
  
  model_ADL_prelim = dLagM::ardlDlm(
    formula = y ~ x1 + x2 + x3,
    data = df,
    p = df_lags %>%
      filter(Variable != "y") %>%
      filter(n_lags == max(n_lags)) %>%
      slice(1) %>%
      pull(),
    q = df_lags %>%
      filter(Variable == "y") %>%
      select(n_lags) %>%
      pull()
  )
  # Remove insignificant lags
  vars_remove = summary(model_ADL_prelim)$coefficients[, 4] %>%
    as.data.frame() %>%
    rename(p_value = c(1)) %>%
    filter(p_value > 0.05)
  
  vars_remove$vars = rownames(vars_remove)
  rownames(vars_remove) = NULL
  
  vars_remove = vars_remove %>%
    mutate(count = str_count(vars)) %>%
    mutate(v1 = ifelse(startsWith(vars, "x"), substr(vars, 1, 2), substr(vars, 1, 1))) %>%
    mutate(v2 = gsub(".*\\.", "", vars)) %>%
    mutate(v2 = as.numeric(ifelse(v2 == "t", 0, v2))) %>%
    dplyr::select(v1, v2)
  
  for (i in c("y", "x1", "x2", "x3")) {
    name = paste("rm", i, sep = "_")
    assign(name,
           vars_remove %>%
             filter(v1 == i) %>%
             dplyr::select(v2) %>%
             pull())
  }
  
  my.remove = list(p = list(x1 = rm_x1,
                            x2 = rm_x2,
                            x3 = rm_x3),
                   q = c(y = rm_y))
  
  my.remove_v = c(rep(NA, times = 3)) # hardcoden for loop
  
  # remove series with insignificant only lags.
  for (i in 1:length(names(my.remove$p))) {
    if (all(as.numeric(unlist(my.remove$p[i])) == c(0:as.numeric(
      df_lags %>%
      filter(Variable != "y") %>%
      filter(n_lags == max(n_lags)) %>%
      pull()
    ))))
    {
      my.remove_v[i] = T
    } else {
      my.remove_v[i] = F
    }
  }
  
  my.remove$p[my.remove_v == T] = NULL
  rm(my.remove_v)
  
  # automatic regression formula
  n_vars = length(names(my.remove$p))
  
  if (n_vars > 0) {
    
    reg_formula = as.formula(paste("y ~",
                                   if (n_vars == 1) {
                                     names(my.remove$p)
                                   } else if (n_vars == 2) {
                                     paste(names(my.remove$p)[1], names(my.remove$p)[2], sep = " + ")
                                   } else if (n_vars == 3) {
                                     paste(names(my.remove$p)[1],
                                           names(my.remove$p)[2],
                                           names(my.remove$p)[3],
                                           sep = " + ")
                                   }, sep = ""))
    
    model_ADL = dLagM::ardlDlm(
      formula = reg_formula,
      data = df,
      p = df_lags %>%
        filter(Variable != "y") %>%
        filter(n_lags == max(n_lags)) %>%
        slice(1) %>%
        pull(),
      q = df_lags %>%
        filter(Variable == "y") %>%
        select(n_lags) %>%
        pull(),
      remove = my.remove
    )
    
    
    # summary(model_ADL)
    
    # forecast recursively: select only vars that are used in forecast
    x = df_total %>%
      filter(time > 250) %>%
      dplyr::select(x1:x3) %>%
      t() %>%
      as.matrix()
    
    if (reg_formula == "y ~ x1") {
      x = t(as.matrix(x[1, ]))
    } else if (reg_formula == "y ~ x2") {
      x = x = t(as.matrix(x[2, ]))
    } else if (reg_formula == "y ~ x3") {
      x = x = t(as.matrix(x[3, ]))
    }
    
    #forecast(model_ADL, x, h = 50)$forecasts
    
    fc = fc %>%
      mutate(ADL = forecast(model_ADL, x, h = 50)$forecasts)
      
  } else {
    fc = fc %>%
      mutate(ADL = NA)
  }
  
    
  
    # 3. OLS-Model
  
  model_OLS = lm(y ~ x1 + x2 + x3, data = df)
  summary(model_OLS)
  
  fc = fc %>%
    mutate(OLS = predict(model_OLS, newdata = df_total)[251:300])
  
  # rm(new_df, df_sig_lags, i, Obs, Win, x1_lag, x2_lag, x3_lag, y_lag)
  
  output = tibble(
    time = 1:50,
    y_counter = df_total$y[251:300],
    y_ARIMA = fc$ARIMA,
    y_ADL = fc$ADL,
    y_OLS = fc$OLS
  )
  
  plot = output %>%
    gather(type, value, y_counter:y_OLS)

  ggplot(plot) +
    aes(x = time, y = value, colour = type) +
    geom_line(size = 1) +
    scale_color_hue(direction = 1) +
    theme_minimal()

  MCMC$RMSFE_OLS[iteration] = sqrt(c(crossprod(output$y_OLS - output$y_counter)) / 50)
  sqrt(c(crossprod(output$y_OLS - output$y_counter)) / 50)
  MCMC$RMSFE_ADL[iteration] = sqrt(c(crossprod(output$y_ADL - output$y_counter)) / 50)
  sqrt(c(crossprod(output$y_ADL - output$y_counter)) / 50)
  MCMC$RMSFE_ARIMA[iteration] = sqrt(c(crossprod(output$y_ARIMA - output$y_counter)) / 50)
  sqrt(c(crossprod(output$y_ARIMA - output$y_counter)) / 50)
}


summary(MCMC$RMSFE_OLS)
summary(MCMC$RMSFE_ADL)
summary(MCMC$RMSFE_ARIMA)



MCMC_long = MCMC %>% 
  select(-c(Iteration:RMSPE_ARIMA)) %>% 
  gather(type, value, RMSFE_OLS:RMSFE_ARIMA)

ggplot(MCMC_long) +
  aes(x = value, fill = type) +
  geom_density(adjust = 1L, alpha = 0.7) +
  scale_fill_hue(direction = 1) +
  theme_minimal()


# 3. SC-Model

# als nächstes Implementieren

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

