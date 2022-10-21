library(dynamac)
library(forecast)
library(tidyverse)
library(tseries)
library(urca)
library(TSstudio)
library(dLagM)

# dynamac vignette
# https://cran.r-project.org/web/packages/dynamac/vignettes/dynamac-vignette.html

rm(list = ls())
head(ineq)

df = ineq %>% 
  filter(year <= 2000)

model1 = dLagM::ardlDlm(formula = concern ~ urate + incshare10, data = df, p = 1, q = 1)
summary(model1)

df = df %>% 
  bind_cols(as.data.frame(c(NA, (model1$model$fitted.values)))) %>% 
  rename(fitted = c(9)) %>% 
  filter(!is.na(fitted))

# prediction
df_plot = df %>% 
  select(year, concern) %>% 
  rename(value = concern) %>% 
  mutate(type = "concern") %>% 
  bind_rows(select(df, year, fitted) %>% 
              rename(value = fitted) %>% 
              mutate(type = "fitted"))

ggplot(df_plot) +
  aes(x = year, y = value, colour = type) +
  geom_line(size = 0.5) +
  scale_color_hue(direction = 1) +
  theme_minimal()

# forecast
x = ineq %>% 
  filter(year > 2000) %>% 
  select(urate, incshare10) %>% 
  t() %>% 
  as.matrix()

fc = forecast(model1, x, h = 14)

df_help = ineq %>% 
  filter(year > 2000) %>% 
  select(year, concern) %>% 
  mutate(forecast = fc$forecasts) 

df_plot %<>% 
  bind_rows(select(df_help, year, concern) %>% 
              rename(value = concern) %>% 
              mutate(type = "concern")) %>% 
  bind_rows(select(df_help, year, forecast) %>% 
              rename(value = forecast) %>% 
              mutate(type = "fitted")) 
  
ggplot(df_plot) +
  aes(x = year, y = value, colour = type) +
  geom_line(size = 0.5) +
  scale_color_hue(direction = 1) +
  geom_vline(xintercept = 2000, linetype="dashed", size=.5)+
  theme_minimal()

# Jetzt Vergleich con SC und ARDL-Modell
rm(list = ls())

library(Synth)
data("basque")

dataprep.out <- dataprep(
  foo = basque,
  predictors = c("school.illit", "school.prim", "school.med",
                 "school.high", "school.post.high", "invest"),
  predictors.op = "mean",
  time.predictors.prior = 1964:1969,
  special.predictors = list(
    list("gdpcap", 1960:1969 ,"mean"),
    list("sec.agriculture", seq(1961, 1969, 2), "mean"),
    list("sec.energy", seq(1961, 1969, 2), "mean"),
    list("sec.industry", seq(1961, 1969, 2), "mean"),
    list("sec.construction", seq(1961, 1969, 2), "mean"),
    list("sec.services.venta", seq(1961, 1969, 2), "mean"),
    list("sec.services.nonventa", seq(1961, 1969, 2), "mean"),
    list("popdens",               1969,               "mean")),
  dependent = "gdpcap",
  unit.variable = "regionno",
  unit.names.variable = "regionname",
  time.variable = "year",
  treatment.identifier = 17,
  controls.identifier = c(2:16, 18),
  time.optimize.ssr = 1960:1969,
  time.plot = 1955:1997)

synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS")

df_out = data.frame(
  year = as.numeric(rownames(dataprep.out$Y1plot))) %>% 
  mutate(Basque = as.numeric(dataprep.out$Y1plot),
         SC_Basque = as.numeric(dataprep.out$Y0plot %*% synth.out$solution.w)) 

df_plot = df_out %>% 
  select(year, Basque) %>% 
  rename(value = Basque) %>% 
  mutate(type = "Basque") %>% 
  bind_rows(select(df_out, year, SC_Basque) %>% 
              rename(value = SC_Basque) %>% 
              mutate(type = "SC_Basque"))

ggplot(df_plot) +
  aes(x = year, y = value, colour = type) +
  geom_line(size = 0.5) +
  scale_color_hue(direction = 1) +
  geom_vline(xintercept = 1975, linetype="dashed", size=.5)+
  ylim(0,12)+
  theme_minimal()

rm(dataprep.out, synth.out)

# alternative estimation using ARDL. For ease of computation: include only Catalonia and Madrid
table(basque$regionname)
table(basque$year)

df_help = basque %>% 
  filter(regionname %in% c("Basque Country (Pais Vasco)", "Cataluna", "Madrid (Comunidad De)"),
         year <= 1975)  %>% 
  select(year, regionname, gdpcap)

df = df_help %>% 
  filter(regionname == "Basque Country (Pais Vasco)") %>% 
  select(year, gdpcap) %>% 
  rename(gdpcap_Basque = gdpcap) %>% 
  bind_cols(df_help %>% 
              filter(regionname == "Cataluna") %>% 
              select(gdpcap) %>% 
              rename(gdpcap_Cataluna = gdpcap)) %>% 
  bind_cols(df_help %>% 
              filter(regionname == "Madrid (Comunidad De)") %>% 
              select(gdpcap) %>% 
              rename(gdpcap_Madrid = gdpcap)) 


model1 = dLagM::ardlDlm(formula = gdpcap_Basque ~ gdpcap_Cataluna + gdpcap_Madrid, data = df, p = 1, q = 1)
summary(model1)

df = df %>% 
  bind_cols(as.data.frame(c(NA, (model1$model$fitted.values)))) %>% 
  rename(fitted = c(5)) %>% 
  filter(!is.na(fitted))

# add forecast to df.

df_help = basque %>% 
  filter(regionname %in% c("Basque Country (Pais Vasco)", "Cataluna", "Madrid (Comunidad De)"),
         year > 1975)  %>% 
  select(year, regionname, gdpcap)

x = df_help %>% 
  filter(regionname == "Cataluna") %>% 
  select(gdpcap) %>% 
  rename(gdpcap_Cataluna = gdpcap) %>% 
  bind_cols(df_help %>% 
              filter(regionname == "Madrid (Comunidad De)") %>% 
              select(gdpcap) %>% 
              rename(gdpcap_Madrid = gdpcap)) %>% 
  t() %>% 
  as.matrix()

colnames(x) = paste0("V", c(1:22))

fc = forecast(model1, x, h = 22)

df_help = data.frame(
  year = c(1956:1997)) %>% 
  left_join(select(df, year, fitted), by = "year") %>% 
  mutate(forecast = (c(rep(NA, times = 20), fc$forecasts))) %>% 
  mutate(combined = ifelse(is.na(forecast), fitted, forecast))
  
df_plot  = df_plot %>% 
  bind_rows(select(df_help, year, combined) %>% 
              rename(value = combined) %>% 
              mutate(type = "ARDL_Basque")) 

ggplot(df_plot) +
  aes(x = year, y = value, colour = type) +
  geom_line(size = 0.75) +
  scale_color_hue(direction = 1) +
  geom_vline(xintercept = 1975, linetype="dashed", size=.5)+
  theme_minimal()

df_out = df_out %>% 
  left_join(select(df_help, year, combined), by = "year") %>% 
  rename(ARDL_Basque = combined)

# calculating RMSE

sqrt(mean((df_out$Basque[df_out$year <= 1975] - df_out$SC_Basque[df_out$year <= 1975])^2))
sqrt(mean((df_out$Basque[df_out$year <= 1975] - df_out$ARDL_Basque[df_out$year <= 1975])^2, na.rm = T))

