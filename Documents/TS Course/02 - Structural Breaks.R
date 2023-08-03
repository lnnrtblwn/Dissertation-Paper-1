install.packages("devtools")
install.packages("changepoint")
install.packages("strucchange")
install.packages("tidyverse")
install.packages("lubridate")

devtools::install_gitlab("KevinKotze/sarbcurrent")

rm(list=ls())  # remove variables
graphics.off() # close figures

library(changepoint)
library(sarbcurrent)
library(tidyverse)
library(lubridate)

# 2 CHANGE POINT TESTS

  # 2.1 Change point in mean - simulated data

set.seed(123) 

sim_mean <- c(rnorm(100, 0, 1),
              rnorm(50, 1.5, 1),
              rnorm(90, 0, 1),
              rnorm(120, -0.8, 1))

plot.ts(sim_mean)

  # binary segmentation method

m_binseg <- cpt.mean(sim_mean, penalty = "BIC", method = "BinSeg", Q = 5)

plot(m_binseg, type = "l", xlab = "Index", cpt.width = 4)

cpts(m_binseg) 

  # segmented neighbour method

m_segneigh <- cpt.mean(sim_mean, penalty = "BIC", method = "SegNeigh", Q = 5)

plot(m_segneigh, type = "l", xlab = "Index", cpt.width = 4)

cpts(m_segneigh)

  # PELT algorithm

m_pelt <- cpt.mean(sim_mean, penalty = "BIC", method = "PELT")

plot(m_pelt, type = "l", cpt.col = "blue", xlab = "Index", cpt.width = 4)

cpts(m_pelt)

  # adding penalty

m_pm <- cpt.mean(sim_mean, penalty = "Manual", pen.value = "1.5 * log(n)",
                 method = "PELT")

plot(m_pm, type = "l", cpt.col = "blue", xlab = "Index", cpt.width = 4)

cpts(m_pm)

  # 2.2 Change point in mean - Real Gross Domestic Product

data(package = 'sarbcurrent')

View(sarb_quarter, title = "sarb_data")

sarb_quarter %>%
  select(date, KBP6006D)

m_gdp <- sarb_quarter %>%
  select(date, KBP6006D) %>%
  mutate(growth = 100 * ((KBP6006D / lag(KBP6006D)) - 1)) %>%
  filter(date > "1960-01-01") %>%
  pull(growth) %>%
  cpt.mean(., penalty = "SIC", method = "PELT")

plot(m_gdp, type = "l", cpt.col = "blue", xlab = "Index", cpt.width = 4)

sarb_quarter %>%
  select(date, KBP6006D) %>%
  mutate(growth = 100 * ((KBP6006D / lag(KBP6006D)) - 1)) %>%
  filter(date > "1960-01-01") %>%
  slice(cpts(m_gdp))

  # 2.3 Change point in variance - simulated data

sim_var <- c(rnorm(100, 0, 1),
             rnorm(50, 0, 2),
             rnorm(90, 0, 1),
             rnorm(120, 0, 0.5))
plot.ts(sim_var)

v_pelt <- cpt.var(sim_var, method = "PELT")
plot(v_pelt, type = "l", cpt.col = "blue", xlab = "Index", cpt.width = 4)
cpts(v_pelt)

  # 2.4 Change point in variance - Real Gross Domestic Product

v_gdp <- sarb_quarter %>%
  select(date, KBP6006D) %>%
  mutate(growth = 100 * ((KBP6006D / lag(KBP6006D)) - 1)) %>%
  filter(date > "1960-01-01") %>%
  pull(growth) %>%
  cpt.var(., method = "PELT")

plot(v_gdp, type = "l", cpt.col = "blue", xlab = "Index", cpt.width = 4)

sarb_quarter %>%
  select(date, KBP6006D) %>%
  mutate(growth = 100 * ((KBP6006D / lag(KBP6006D)) - 1)) %>%
  filter(date > "1960-01-01") %>%
  slice(cpts(v_gdp))

   # 2.5 Change points in the mean and variance - simulated data

sim_mv <- c(rnorm(100, 0, 1),
            rnorm(50, 1, 2),
            rnorm(90, 0, 1),
            rnorm(120, -0.8, 0.5))
plot.ts(sim_mv)

mv_pelt <- cpt.meanvar(sim_mv, method = "PELT")
plot(mv_pelt)

  # 2.6 Change point in mean and variance - Real Gross Domestic Product

mv_gdp <- sarb_quarter %>%
  select(date, KBP6006D) %>%
  mutate(growth = 100 * ((KBP6006D / lag(KBP6006D)) - 1)) %>%
  filter(date > "1960-01-01") %>%
  pull(growth) %>%
  cpt.meanvar(., method = "PELT")

plot(mv_gdp, type = "l", cpt.col = "blue", xlab = "Index", cpt.width = 4)

sarb_quarter %>%
  select(date, KBP6006D) %>%
  mutate(growth = 100 * ((KBP6006D / lag(KBP6006D)) - 1)) %>%
  filter(date > "1960-01-01") %>%
  slice(cpts(mv_gdp))

# 3 STRUCTURAL BREAK TESTS

rm(list=ls())
graphics.off

library(strucchange)
library(sarbcurrent)
library(tidyverse)
library(lubridate)

  # 3.1 Structural break tests - simulated data

set.seed(123) 

x1 <- arima.sim(model = list(ar = 0.9), n = 100)
x2 <- arima.sim(model = list(ma = 0.1), n = 100)
x3 <- arima.sim(model = list(ar = 0.5, ma = 0.3), n = 100)

y <- c((1 + x1),
       x2,
       (0.5 - x3))
plot.ts(y)

dat <- tibble(ylag0 = y,
              ylag1 = lag(y)) %>%
  drop_na()

qlr <- Fstats(ylag0 ~ ylag1, data = dat)

breakpoints(qlr)

sctest(qlr, type = "supF")

plot(qlr)

cusum <- efp(ylag0 ~ ylag1, type = "OLS-CUSUM", data = dat)
plot(cusum)

  # 3.2 Structural break tests - Real Gross Domestic Product

sa_dat <- sarb_quarter %>%
  select(date, KBP6006D) %>%
  mutate(growth = 100 * ((KBP6006D / lag(KBP6006D)) - 1),
         grow_lag = lag(growth)) %>%
  drop_na()

sa_qlr <- Fstats(growth ~ grow_lag, data = sa_dat)
breakpoints(sa_qlr)

sctest(sa_qlr, type = "supF")

plot(sa_qlr)

sarb_quarter %>%
  select(date, KBP6006D) %>%
  mutate(growth = 100 * ((KBP6006D / lag(KBP6006D)) - 1),
         grow_lag = lag(growth)) %>%
  drop_na() %>%
  slice(sa_qlr$breakpoint)

sa_bp <- breakpoints(growth ~ grow_lag, data = sa_dat, breaks = 5)
summary(sa_bp)  # where the breakpoints are

plot(sa_bp, breaks = 5)

sarb_quarter %>%
  select(date, KBP6006D) %>%
  mutate(growth = 100 * ((KBP6006D / lag(KBP6006D)) - 1),
         grow_lag = lag(growth)) %>%
  drop_na() %>%
  slice(sa_bp$breakpoint)

sa_cusum <- efp(growth ~ grow_lag, data = sa_dat, type = "OLS-CUSUM")
plot(sa_cusum)
