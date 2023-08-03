rm(list=ls()) 
graphics.off() 

# 1 THE BASIC SETUP

devtools::install_gitlab("KevinKotze/sarb2020q1")
devtools::install_github("KevinKotze/tsm")
install.packages('strucchange', repos='https://cran.rstudio.com/', dependencies=TRUE)
install.packages('forecast', repos='https://cran.rstudio.com/', dependencies=TRUE)

library(tsm)
library(forecast)
library(strucchange)
library(tidyverse)
library(lubridate)
library(sarb2020q1)

# 2 TRANSFORMING THE DATA

gdp <- sarb_quarter %>%
  select(date, KBP6006D) %>%
  mutate(growth = 100 * ((KBP6006D / lag(KBP6006D)) - 1),
         grow_lag = lag(growth)) %>%
  drop_na()

gdp %>%
  pull(growth) %>%
  plot.ts()

sa_bp <- breakpoints(growth ~ grow_lag, data = gdp, breaks = 5)
summary(sa_bp)

plot(sa_bp, breaks = 5)

gdp %>%
  slice(sa_bp$breakpoint)

# 3 ESTIMATING THE MODEL PARAMETERS

y <- gdp %>%
  filter(date > ymd('1981-01-01')) %>%
  select(date, growth)

y %>%
  pull(growth) %>%
  ac()

# 4 RECURSIVE FORECAST GENERATION AND EVALUATION

H <- 8
P <- 40
R <- dim(y)[1] - (P + H) + 1

train <- y %>%
  slice_head(n = R)

test <- y %>%
  slice_tail(n = (P + H - 1))

train %>%
  tail()

test %>%
  head()

actual <- matrix(data = rep(0, P * H), ncol = H)

illust <- y %>%
  mutate(num = seq(1, n()))

for (i in 1:P) {
  first <- R + i
  last <- first + H - 1
  actual[i, 1:H] <- illust$num[first:last]
}

actual <- matrix(data = rep(0, P * H), ncol = H)

for (i in 1:P) {
  first <- R + i
  last <- first + H - 1
  actual[i, 1:H] <- y$growth[first:last]
}

arma_fore <- matrix(data = rep(0, P * H), ncol = H)

for (i in 1:P) {
  arma_fore[i, 1:H] <- y %>%
    slice_head(n = (R + i - 1)) %>%
    pull(growth) %>%
    arima(., order = c(1, 0, 1)) %>%
    predict(., n.ahead = H) %>%
    pluck(1)
}

ar_fore <- matrix(data = rep(0, P * H), ncol = H)

for (i in 1:P) {
  ar_fore[i, 1:H] <- y %>%
    slice_head(n = (R + i - 1)) %>%
    pull(growth) %>%
    arima(., order = c(1, 0, 0)) %>%
    predict(., n.ahead = H) %>%
    pluck(1)
}

rowf <- 30

par(mfrow = c(1, 1))
plot.ts(
  actual[rowf,],
  col = "black",
  ylab = "Example",
  xlab = "Steps",
  ylim = c(min(actual[rowf,] * (-2)), max(actual[rowf,] * 2))
)
lines(ar_fore[rowf,],
      col = "red",
      lty = 2,
      lwd = 2)
legend(
  "topright",
  legend = c("actual", "forecast"),
  lty = c(1, 2),
  col = c("black", "red"),
  bty = "n"
)

error_arma <- actual - arma_fore
error_ar <- actual - ar_fore

bias_arma_step <-  colMeans(error_arma)
bias_ar_step <- colMeans(error_ar)

bias_arma_time <- rowMeans(error_arma)
bias_ar_time <- rowMeans(error_ar)

RMSE_arma_step <- sqrt(colMeans(error_arma ^ 2))
RMSE_ar_step <- sqrt(colMeans(error_ar ^ 2))

RMSE_arma_time <- sqrt(rowMeans(error_arma ^ 2))
RMSE_ar_time <- sqrt(rowMeans(error_ar ^ 2))

RMSE_arma_step %>%
  plot()

RMSE_ar_step %>%
  plot()

RMSE_arma_time %>%
  plot()

RMSE_ar_time %>%
  plot()

# 5 MODEL COMPARISON

dm_steps <- rep(0, H)

for (i in 1:H) {
  dm_steps[i] <-
    dm.test(error_arma[, i], error_ar[, i], h = i)$statistic
}

plot(dm_steps,
     ylim = c(-3, 3),
     xlab = "Steps",
     ylab = "DM statistic")
lines(rep(0, H), col = "grey", lty = 1)
lines(rep(1.96, H), col = "red", lty = 2)
lines(rep(-1.96, H), col = "red", lty = 2)

RMSE_arma_step[2]
RMSE_ar_step[2]

dm_time <- rep(0, P)

for (i in 1:P) {
  dm_time[i] <-
    dm.test(error_arma[i,], error_ar[i,], h = 1)$statistic
}

plot(dm_time,
     ylim = c(-3, 3),
     xlab = "Steps",
     ylab = "DM statistic")
lines(rep(0, P), col = "grey", lty = 1)
lines(rep(1.96, P), col = "red", lty = 2)
lines(rep(-1.96, P), col = "red", lty = 2)

RMSE_arma_time[5]
RMSE_ar_time[5]

# 5.1 Comparing nested models

cwstat_steps <- rep(0, H)
cwpval_steps <- rep(0, H)

for (i in 1:H) {
  cwstat_steps[i] <-
    cw(error_arma[, i], error_ar[, i], arma_fore[, i], ar_fore[, i])$test
  cwpval_steps[i] <-
    cw(error_arma[, i], error_ar[, i], arma_fore[, i], ar_fore[, i])$pvalue
}


plot(cwstat_steps,
     ylim = c(-3, 3),
     xlab = "Steps",
     ylab = "CW statistic")
lines(rep(0, H), col = "grey", lty = 1)
lines(rep(1.96, H), col = "red", lty = 2)
lines(rep(-1.96, H), col = "red", lty = 2)

cwstat_time <- rep(0, P)
cwpval_time <- rep(0, P)

for (i in 1:P) {
  cwstat_time[i] <-
    cw(error_arma[i,], error_ar[i,], arma_fore[i,], ar_fore[i,])$test
  cwpval_time[i] <-
    cw(error_arma[i,], error_ar[i,], arma_fore[i,], ar_fore[i,])$test
}

plot(cwstat_time,
     ylim = c(-5, 5),
     xlab = "Time",
     ylab = "CW statistic")
lines(rep(0, P), col = "grey", lty = 1)
lines(rep(1.96, P), col = "red", lty = 2)
lines(rep(-1.96, P), col = "red", lty = 2)

# 5.2 Other forecast plots

h_start <- R - 12

y_plot <-
  ts(y$growth[h_start:length(y$growth)], start = c(2005, 1), frequency =
       4)

n_obs <- length(y_plot)
no_fore <- n_obs - P - H
nons <- n_obs - H

fore_plot <- matrix(data = rep(0, n_obs * P), nrow = P)

for (i in 1:P) {
  nan_beg <- rep(NaN, no_fore + i)
  nan_end <- rep(NaN, P - i)
  fore_plot[i, 1:n_obs] <- c(nan_beg, ar_fore[i,], nan_end)
}

length(fore_plot[i, 1:n_obs])

length(c(nan_beg, ar_fore[i,], nan_end))

plot.ts(y_plot)

for (i in 1:P) {
  lines(ts(fore_plot[i,], start = c(2005, 1), frequency = 4),
        col = "red",
        lty = 2)
}

# 6 CONSTRUCTING A FORECAST DENSITY

H <- 16
percentiles <- c(0.10, 0.15, 0.25, 0.35)
num <- dim(y)[1]
inter <- length(percentiles) * 2

ar_mod <- tibble(forecast = y %>%
                   pull(growth) %>%
                   arima(., order = c(1, 0, 0)) %>%
                   predict(., n.ahead = H) %>%
                   pluck(1))

ar_mod <- ar_mod %>%
  mutate(mse = rep(0, H))

for (i in 1:H) {
  ar_mod$mse[i] <- var(c(y$growth, ar_mod$forecast[i]))
}

fanChart <- matrix(rep(0, H * inter), nrow = H)

for (i in 1:length(percentiles)) {
  z2 <- abs(qnorm(percentiles[i] / 2, 0, 1))
  forcInt <- matrix(rep(0, H * 2), ncol = 2)
  
  for (h in 1:H) {
    forcInt[h, 1] <- ar_mod$forecast[h] - z2 * sqrt(ar_mod$mse[h])
    forcInt[h, 2] <- ar_mod$forecast[h] + z2 * sqrt(ar_mod$mse[h])
    
    fanChart[h, i] <- forcInt[h, 1]
    fanChart[h, inter - i + 1] <- forcInt[h, 2]
  }
}

plot_len <- num + H

plot_yt <- c(y$growth, rep(NaN, H))
plot_ytF <- c(rep(NaN, num), ar_mod$forecast)

plot.ts(plot_yt, type = "o")
lines(plot_ytF, col = "red", lty = 1)

plot_ytD <- rbind(matrix(rep(NaN, num * inter), ncol = inter), fanChart)

no_int <- dim(plot_ytD)

for (i in 1:no_int[2]) {
  lines(plot_ytD[, i], col = "darkgrey", lty = 2)
}

start_obs <- num + 1

plot(
  plot_yt,
  type = "o",
  ylim = c(-6, 6),
  ylab = "",
  xlab = ""
)

polygon(c(c(start_obs:plot_len), rev(c(start_obs:plot_len))),
        c(plot_ytD[start_obs:plot_len, 1], rev(plot_ytD[start_obs:plot_len, 8])),
        col = "grey90",
        border = NA)

polygon(c(c(start_obs:plot_len), rev(c(start_obs:plot_len))),
        c(plot_ytD[start_obs:plot_len, 2], rev(plot_ytD[start_obs:plot_len, 7])),
        col = "grey80",
        border = NA)

polygon(c(c(start_obs:plot_len), rev(c(start_obs:plot_len))),
        c(plot_ytD[start_obs:plot_len, 3], rev(plot_ytD[start_obs:plot_len, 6])),
        col = "grey70",
        border = NA)

polygon(c(c(start_obs:plot_len), rev(c(start_obs:plot_len))),
        c(plot_ytD[start_obs:plot_len, 4], rev(plot_ytD[start_obs:plot_len, 5])),
        col = "grey60",
        border = NA)

lines(plot_ytF, col = "red", lty = 1)

