
 # 1 BASIC SETUP FOR MOST EMPIRICAL WOR

rm(list=ls())
graphics.off()

devtools::install_github("KevinKotze/tsm")
devtools::install_gitlab("KevinKotze/sarb2020q1")
install.packages("fArma")
install.packages("forecast")

library(tidyverse)
library(tsm)

# 2 SIMULATED PROCESSES

set.seed(123) 

dat <- tibble(
  ar1 = arima.sim(model = list(ar = 0.8), n = 200),
  ma1 = arima.sim(model = list(ma = 0.7), n = 200),
  arma10 = arima.sim(model = list(ar  = 0.4), n = 200),
  arma01 = arima.sim(model = list(ma = 0.5), n = 200),
  arma11 = arima.sim(model = list(ar = 0.4, ma = 0.5), n = 200),
  arma21 = arima.sim(model = list(ar = c(0.6,-0.2), ma = c(0.4)), n = 200))

# 3 SIMULATED AUTOREGRESSIVE PROCESS

dat %>%
  pull(ar1) %>%
  plot.ts()

dat %>%
  pull(ar1) %>%
  ac(., max.lag = 20)

Box_res <- dat %>%
  pull(ar1) %>%
  Box.test(., lag = 1, type = "Ljung-Box") 

Box_res  

names(Box_res)

Q_stat <- rep(0, 10)
Q_prob <- rep(0, 10)

Q_stat[1] <- Box.test(dat$ar1, lag=1, type="Ljung-Box")$statistic
Q_prob[1] <- Box.test(dat$ar1, lag=1, type="Ljung-Box")$p.value

Q_stat
Q_prob

for (i in 1:10) {
  Q_stat[i] <- Box.test(dat$ar1, lag=i, type="Ljung-Box")$statistic
  Q_prob[i] <- Box.test(dat$ar1, lag=i, type="Ljung-Box")$p.value
}

op <- par(mfrow = c(1, 2)) # create plot area of (1 X 2)
plot(Q_stat, ylab = "", main = 'Q-Statistic')
plot(Q_prob, ylab = "", ylim = c(0,1), main = 'Probability values')

par(op) # close plot area

arma10 <- dat %>%
  pull(ar1) %>%
  arima(., order = c(1,0,0), include.mean = FALSE) # uses ARIMA(p,d,q) specification
arma10

par(mfrow=c(1, 1))
arma10$residuals %>%
  plot()

arma10$residuals %>%
  ac(., max.lag=20)

# 4 SIMULATED MOVING AVERAGE PROCESS

par(mfrow=c(1, 1))
dat %>%
  pull(ma1) %>%
  plot.ts()

dat %>%
  pull(ma1) %>%
  ac(., max.lag=20)

arma01 <- dat %>%
  pull(ma1) %>%
  arima(., order=c(0,0,1)) # uses ARIMA(p,d,q) with constant
arma01

par(mfrow=c(1, 1))
arma01$residuals %>%
  plot()

par(mfrow=c(1, 1))
arma01$residuals %>%
  ac(., max.lag=20)

# 5 SIMULATED ARMA PROCESS

# autocorrelation function for ARMA(1,0)
dat %>%
  pull(arma10) %>%
  ac(., max.lag=20)

# autocorrelation function for ARMA(0,1)
dat %>%
  pull(arma01) %>%
  ac(., max.lag=20)

# autocorrelation function for ARMA(1,1)
dat %>%
  pull(arma11) %>%
  ac(., max.lag=20)

# autocorrelation function for ARMA(2,1)
dat %>%
  pull(arma21) %>%
  ac(., max.lag=20)

arma_res <- rep(0,16)
arma_res[1] <- arima(dat$arma21, order=c(3,0,2))$aic # fit arma(3,2) and save aic value
arma_res[2] <- arima(dat$arma21, order=c(2,0,2))$aic
arma_res[3] <- arima(dat$arma21, order=c(2,0,1))$aic
arma_res[4] <- arima(dat$arma21, order=c(1,0,2))$aic
arma_res[5] <- arima(dat$arma21, order=c(1,0,1))$aic
arma_res[6] <- arima(dat$arma21, order=c(3,0,0))$aic
arma_res[7] <- arima(dat$arma21, order=c(2,0,0))$aic
arma_res[8] <- arima(dat$arma21, order=c(0,0,2))$aic
arma_res[9]  <- arima(dat$arma21, order=c(3,0,2), include.mean=FALSE)$aic
arma_res[10] <- arima(dat$arma21, order=c(2,0,2), include.mean=FALSE)$aic
arma_res[11] <- arima(dat$arma21, order=c(2,0,1), include.mean=FALSE)$aic
arma_res[12] <- arima(dat$arma21, order=c(1,0,2), include.mean=FALSE)$aic
arma_res[13] <- arima(dat$arma21, order=c(1,0,1), include.mean=FALSE)$aic
arma_res[14] <- arima(dat$arma21, order=c(3,0,0), include.mean=FALSE)$aic
arma_res[15] <- arima(dat$arma21, order=c(2,0,0), include.mean=FALSE)$aic
arma_res[16] <- arima(dat$arma21, order=c(0,0,2), include.mean=FALSE)$aic

which(arma_res == min(arma_res))

arma11 <- arima(dat$arma21, order=c(1, 0, 1), include.mean=FALSE)

arma11$residuals %>%
  ac(., max.lag=20)

Q_stat <- rep(NA,10) # vector of ten observations
Q_prob <- rep(NA,10)

for (i in 4:10) {
  Q_stat[i] <- Box.test(arma11$residuals, lag=i, type="Ljung-Box", fitdf=3)$statistic
  Q_prob[i] <- Box.test(arma11$residuals, lag=i, type="Ljung-Box", fitdf=3)$p.value
}

op <- par(mfrow = c(1, 2))
plot(Q_stat, ylab = "", main='Q-Statistic')
plot(Q_prob, ylab = "", ylim=c(0,1), main='Probability values')

par(op)

# 6 UNIVARIATE MODELS FOR REAL DATA

  # 6.1 South African gross domestic product

rm(list=ls()) 
graphics.off() 

library(fArma)
library(forecast)
library(strucchange)
library(tidyverse)
library(lubridate)
library(tsm)
library(sarb2020q1)

gdp <- sarb_quarter %>%
  select(date, KBP6006D) %>%
  mutate(growth = 100 * ((KBP6006D / lag(KBP6006D)) - 1),
         grow_lag = lag(growth)) %>%
  drop_na()

gdp %>%
  pull(growth) %>%
  plot.ts()

sa_bp <- breakpoints(growth ~ grow_lag, data = gdp, breaks = 5)
summary(sa_bp)  # where the breakpoints are

plot(sa_bp, breaks = 5)

gdp %>%
  slice(sa_bp$breakpoint)

y <- gdp %>%
  filter(date > ymd('1981-01-01')) %>%
  pull(growth)

y %>% ac() 

arma_res <- rep(0,5)

arma_res[1] <- arima(y, order=c(1, 0, 2))$aic
arma_res[2] <- arima(y, order=c(1, 0, 1))$aic
arma_res[3] <- arima(y, order=c(1, 0, 0))$aic
arma_res[4] <- arima(y, order=c(0, 0, 2))$aic
arma_res[5] <- arima(y, order=c(0, 0, 1))$aic

which(arma_res == min(arma_res))

arma <- y %>%
  arima(., order=c(1,0,0))

par(mfrow=c(1, 1))
arma$residuals %>%
  plot()

arma$residuals %>%
  ac()

Q_stat <- rep(NA,12) # vector of twelve observations
Q_prob <- rep(NA,12)

for (i in 6:12) {
  Q_stat[i] <- Box.test(arma$residuals, lag = i, type =  "Ljung-Box",fitdf = 2)$statistic
  Q_prob[i] <- Box.test(arma$residuals, lag = i, type =  "Ljung-Box",fitdf = 2)$p.value
}

op <- par(mfrow = c(1, 2))
plot(Q_stat, ylab = "", main='Q-Statistic')
plot(Q_prob, ylab = "", ylim=c(0,1), main='Probability values')

par(op)

# arma_fit1 <- fArma::armaFit(~ arma(1,2), data=y, include.mean=FALSE)
# summary(arma_fit1)
auto_mod <- forecast::auto.arima(y, max.p = 2, max.q = 2, start.p = 0, start.q = 0,
                                 stationary = TRUE, seasonal = FALSE, allowdrift = FALSE,
                                 ic = "aic")

auto_mod # AR(1) with intercept

y_sub <- gdp %>%
  filter(date > ymd('1981-01-01')) %>%
  slice(1:(n() - 12)) %>%
  pull(growth)


arma10 <- y_sub %>%
  arima(., order=c(1, 0, 0))
arma12 <- y_sub %>%
  arima(., order=c(1, 0, 2))

arma10_pred <- predict(arma10, n.ahead=12)
arma12_pred <- predict(arma12, n.ahead=12)

par(mfrow=c(1, 1))
plot.ts(y)
lines(arma10_pred$pred, col="blue", lty=2)
lines(arma12_pred$pred, col="green", lty=2)

# 6.2 South African consumer price index

rm(list=ls()) 
graphics.off() 

library(forecast)
#library(astsa)
library(strucchange)
library(tidyverse)
library(lubridate)
library(tsm)
library(sarb2020q1)

dat <- sarb_month %>%
  select(date, KBP7170N) %>%
  mutate(inf_yoy  = 100 * ((KBP7170N / lag(KBP7170N, n = 12)) - 1),
         infl_lag = lag(inf_yoy)) %>%
  drop_na()

dat %>%
  pull(inf_yoy) %>%
  plot.ts()

dat <- dat %>%
  filter(date >= ymd("2000-01-01"))

dat %>%
  pull(inf_yoy) %>%
  plot.ts()

sa_bp <- breakpoints(inf_yoy ~ infl_lag, data = dat, breaks = 5)
summary(sa_bp)  # where the breakpoints are

plot(sa_bp, breaks = 5)

dat %>%
  slice(sa_bp$breakpoint)

y <- dat %>%
  filter(date > ymd('2006-12-01')) %>%
  pull(inf_yoy)

y %>% ac()

arma_res <- rep(0,8)
arma_res[1]  <- arima(y, order=c(2, 0, 4))$aic
arma_res[2]  <- arima(y, order=c(2, 0, 3))$aic
arma_res[3]  <- arima(y, order=c(2, 0, 2))$aic
arma_res[4]  <- arima(y, order=c(2, 0, 1))$aic
arma_res[5]  <- arima(y, order=c(2, 0, 0))$aic
arma_res[6]  <- arima(y, order=c(1, 0, 2))$aic
arma_res[7]  <- arima(y, order=c(1, 0, 1))$aic
arma_res[8] <- arima(y, order=c(1, 0, 0))$aic

which(arma_res == min(arma_res))

arma <- y %>%
  arima(., order=c(2,0,3))

summary(arma)

par(mfrow=c(1, 1))
arma$residuals %>%
  plot()

arma$residuals %>%
  ac()

sarima_fit <- astsa::sarima(y,2,0,3,1,0,1,12) # sarima(x,p,d,q,P,D,Q,S) will fit SARIMA(2,0,3)*(1,0,1)_{12}

ac(as.numeric(sarima_fit$fit$residuals))
