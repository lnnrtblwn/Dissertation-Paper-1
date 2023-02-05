options(scipen = 50)

suppressPackageStartupMessages({
  list.of.packages = c("tidyverse", "magrittr", "expss", "writexl", "readxl", "sjlabelled", "haven", "urca", "vars", "mFilter", "tseries", "forecast", "faux")
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages, require, character.only = T,)
})

rm(list = ls(all = T))

# df = read_excel("Diss/Topics/Synthetic Control/Chunks/Simulations/Data/P_Data_Extract_From_World_Development_Indicators.xlsx")
df = read_csv("Diss/Topics/Synthetic Control/Chunks/Simulations/Data/DP_LIVE_03022023114159743.csv")

oecd = df %>% 
  filter(LOCATION %in% c("AUT", "BEL", "CAN", "DNK", "FRA", "DEU", "GRC", "ISL", "IRL", "ITA", "LUX", "NLD", "NOR", "PRT", "ESP", "SWE", "CHE", "TUR", "GBR", "USA"),
         MEASURE == "MLN_USD") %>% 
  mutate(LOCATION = as.factor(LOCATION)) %>% 
  mutate(TIME = as.numeric(TIME)) %>% 
  filter(TIME >= 1980) 

ggplot(oecd) +
  aes(x = TIME, y = Value, colour = LOCATION) +
  geom_line() +
  scale_color_hue(direction = 1) +
  theme_minimal()

oecd %<>% 
  dplyr::select(-c(8)) %>% 
  spread(LOCATION, Value) %>% 
  dplyr::select(-c(INDICATOR:FREQUENCY)) %>% 
  dplyr::select(-TIME)

# Time series interpolation to monthly data 

oecd_inter = oecd 

oecd_inter[43:(42*12),] = NA

oecd_inter = oecd_inter %>% 
  mutate(Group = rep(1:12, each = 42)) %>% 
  group_by(Group) %>% 
  mutate(Counter = row_number()) %>% 
  ungroup() %>% 
  arrange(Counter) %>% 
  dplyr::select(-Counter) %>% 
  head(41*12+1) %>% 
  mutate_at(dplyr::vars(c(AUT:USA)), 
            .funs = list(~na.approx(.))) %>% 
  dplyr::select(-Group)

# Assign TS-Structure

oecd_ts = ts(oecd_inter$AUT, start = 1980, end = 2021, frequency = 12)

for (i in 2:ncol(oecd_inter)) {
  oecd_ts = cbind(oecd_ts, ts(oecd_inter[,i], start = 1980, end = 2021, frequency = 12))
}

colnames(oecd_ts) = colnames(oecd_inter)

oecd5 = oecd_inter %>% 
  dplyr::select(USA, DEU, FRA, GBR, ITA)

oecd5_ts = ts(oecd5$USA, start = 1980, end = 2021, frequency = 12)

for (i in 2:ncol(oecd5)) {
  oecd5_ts = cbind(oecd5_ts, ts(oecd5[,i], start = 1980, end = 2021, frequency = 12))
}

colnames(oecd5_ts) = colnames(oecd5)

autoplot(oecd_ts)
autoplot(oecd5_ts)

#========================================================
# VAR model in level
#========================================================

# forecasting horizon
nhor = 48 

# lag length
VARselect(oecd5_ts, lag.max = 12,
          type = "const")

# estimatio
var.model_lev = VAR(oecd5_ts, p = 2,
                    type = "const")

# forecast of lev data
var.pred = predict(var.model_lev, n.ahead = nhor)

test = 
















# Augmented Dickey Fuller to test Stationarity

adf = adf.test(oecd_inter$AUT)
adf

# cant take first differences if data is interpolated ... 
# oecd_diff = ts(diff(oecd_inter$AUT), start = 1980, end = 2021, frequency = 12)
# 
# for (i in 2:ncol(oecd_inter)) {
#   oecd_diff = cbind(oecd_diff, ts(diff(oecd_inter[[i]]), start = 1980, end = 2021, frequency = 12))
# }
# 
# colnames(oecd_diff) = colnames(oecd_inter)

# Tutorial: https://www.youtube.com/watch?v=xOK0UlrR3Ks&t=642s&ab_channel=Hanomics

autoplot(oecd_ts)

lag = VARselect(oecd_ts, type = "both", lag.max = 10)  

m1 = VAR(oecd_ts, p = 2, type = "both")
summary(m1)

stargazer::stargazer(m1[["varresult"]], type = 'html', 
                     out = "Diss/Topics/Synthetic Control/Chunks/Simulations/Data/table.html")

# Stable model
roots(m1, modulus = T)

# Granger Causality
granger = causality(m1, cause = "AUT")
granger$Granger
# small p-valu, reject H0 --> AUT does granger cause the remaining variables.

# Impulse response functions
irf1 = irf(m1, impulse = "AUT", response = "BEL", n.ahead = 12, boot = T, run = 50, ci = 0.95)
plot(irf1)

# Variance decomposition
vd = fevd(m1, n.aheda = 10)
plot(vd)

# Ideen:
# an Born orientieren und nur Länder verwenden, die sie auch verwendet haben
# für diese Länder VAR fitten und auch SC Modell fitten
# prüfen ob SC mit erklärenden Variablen besser performt als VAR, das nur GDP bekommt
# Wie verwende ich Varianz-Covarianz-Matrix um Simulation darauf aufzubauen
# Was passiert wenn man nicht-stationären Prozess mit VAR fittet? Inkonsistent?
