





# Packages ----


library(dplyr)
library(tidyr)
library(ggplot2)
library(WDI)
library(lubridate)
library(corrplot)
library(vars)


# Data Input via World Bank Database ----


#Country_List = WDI(indicator='NY.GDP.MKTP.KD.ZG', country="all", start=2020, end=2020)

#Indicator Codes for GDP growth rate, GDP and GDP per Capita 
indicator_code = c('NY.GDP.MKTP.KD.ZG','NY.GDP.MKTP.CD','NY.GDP.PCAP.CD','NY.GDP.PCAP.KD.ZG')

# Choose Indicator Code: Here 3: GDP per Capita

#ic = 1 # GDP growth (annual %)
#ic = 2 # GDP Current in USD
#ic = 3 # GDP Per Capita in USD 
ic = 4 # GDP Per Capita Growth Rate


# Country Presets Like G20, EU, Abadie2015 Country List ----

g20=c('CN','IN','US','ID','BR','RU','JP','MX','DE','TR','FR','GB','IT','ZA','KR'
      ,'AR','CA','SA','AU')
#euro_area = 'EMU'

EU = c('BE','BG','CZ','DK','DE','EE','IE','GR','ES','FR','HR','IT','CY','LV','LT','LU','HU','MT','NL',
       'AT','PL','PT','RO','SI','SK','FI','SE')

#
Abadie15 = c('AT','NL','JP','US','CH', 'DE')
Abadie15 = c('AT','NL','CH', 'DE')

# Obtain Data unsing Country List and Indicator Code----
#df = WDI(indicator=indicator_code[ic], country=Abadie15, start=1960, end=NULL)
df_full = WDI(indicator=indicator_code[ic], country=c(EU,g20), start=1960, end=NULL)
#df_full = WDI(indicator=indicator_code[ic], start=1960, end=NULL)

#----

#Define Top Element such that the correlation Matrix has the desired structure later on

# top_element <- "Germany"
# 
# df_spread = df %>%
#   filter(year >= 1981) %>%
#   filter(year <= 2021) %>%
#   dplyr::select(year, country, indicator_code[ic]) %>%
#   pivot_wider(names_from = country, values_from = indicator_code[ic]) %>%
#   relocate(paste(top_element), .after = year)


#Choose mix from G20 and EU Countries as Base

df_spread_full = df_full %>%
  filter(year >= 1981) %>% 
  filter(year <= 2021) %>% 
  dplyr::select(year, country, indicator_code[ic]) %>% 
  pivot_wider(names_from = country, values_from = indicator_code[ic]) %>% 
  select_if(~all(complete.cases(.))) 

#Check Stationarity and remove all TS with a ADF p Value > 0.1
adf = c()
names =c()
for (i in 2:dim(df_spread_full)[2]){
  adf[i]=as.numeric(tseries::adf.test(ts(df_spread_full[i]))$p.value)
  names[i] = colnames(df_spread_full[i])
}

resid = as.data.frame(cbind(names[-1],adf[-1])) %>% 
  filter(V2 >=0.1) 

df_spread_resid <- df_spread_full[, colnames(df_spread_full) %in% resid[,1]]



VAR_est = function(J = 3, J_max =14, p = 2){

  #size = J+1
  size = (J_max+1)
  
#Choose random subgroup
df_spread = cbind(df_spread_full[1], sample((df_spread_resid[-1]), size = size))

#Randomly Relocate the Treatment Country
df_spread = df_spread %>% 
  relocate(paste(sample(colnames(df_spread), size = (1))), .after = year)

# perform unit root tests on output and interest rates ----
#tseries::adf.test(ts(df_spread[2]))
 
  

#Optimal Lag Number
# varselect = vars::VARselect(df_spread[-1])
# p= as.numeric(varselect$selection[1])

p = p
k = dim(df_spread)[2]-1

#VAR Model
VAR_mod= vars::VAR(df_spread[-1], p = p, type = "const")
#summary(VAR_mod)
#plot(vars::fevd(VAR_mod, n.ahead = 10))
#VAR_mod$varresult$Germany
#synth = predict(VAR_mod, n.ahead = 50)$fcst


#Coef extraction
est_coef <- coef(VAR_mod)

# Extract only the coefficients for both dependend variables
# and combine them to a single matrix

est_coefs <- est_coef[[1]][, 1]

for (i in 1:(k-1)){
  est_coefs <- rbind(est_coefs, est_coef[[i+1]][, 1])
}

est_coefs=(as.matrix(round(est_coefs, 4)))

return(est_coefs)
}

#rm(list= ls()[!(ls() %in% c('est_coefs'))])


###### ### ### ###
### VAR Simu ###
################

VAR_simu = function(est_coefs){

#set.seed(123)


t = 2000 # Number of time series observations
k <- dim(est_coefs)[1] # Number of variables
p <- (dim(est_coefs)[2]-1)/k # Number of lags (Note: here the case with const)


A = est_coefs

## Generate coefficient matrices if const
A = est_coefs[,-(dim(est_coefs)[2])]
const = est_coefs[,(dim(est_coefs)[2])]

submatrices <- list()
# Generate submatrices using a for loop
for (i in 1:p) {
  submatrix <- matrix(A[ ((i - 1) * k * k + 1):(i * k * k)],k)
  submatrices[[i]] <- (submatrix)
  assign(paste0("A_", i), submatrices[[i]])
}



#series <- matrix(0, k, t + 2 * p )  # Raw series with zeros
series <- matrix(const, k, t + 2 * p )  # Raw series with const

for (i in (p + 1):(t + 2 * p )) {
  # Generate series with e ~ N(0,1)
  for (j in 1:p) {
    series[, i] <- series[, i] + (submatrices[[j]]) %*% (series[, i - j])
  }
  series[, i] <- series[, i] + rnorm(k, mean=0, sd= var_error_VAR)
}




series <- ts(t(series[, -(1:k)])) # Convert to time series format



#plot.ts(series) # Plot the series

#companion = rbind(est_coefs[,-(dim(est_coefs)[2])],cbind(diag((p-1)*k),matrix(0, nrow = (p-1)*k, ncol = k)))

#eigen_val = eigen(companion)$values

return(series)

}



Stat_test = function(est_coefs) {
  
  k = dim(est_coefs)[1] # Number of variables
  p = (dim(est_coefs)[2]-1)/k # Number of lags (Note: here the case with const)
  
  companion = rbind(est_coefs[,-(dim(est_coefs)[2])],cbind(diag((p-1)*k),matrix(0, nrow = (p-1)*k, ncol = k)))
  
  eigen_val = eigen(companion)$values
  
  return(max(abs(eigen_val)))
  
}


rm(list=c("Abadie15", "EU", "g20", "i", "ic","indicator_code"))


# rm("df_full", "df_spread_full", "df_spread_resid", "resid", "adf", "names")


#est_coefs = VAR_est(J=4,p=4)

#Stat_test(est_coefs)

#VAR_simu(est_coefs)










# #Austesten
# y = VAR_simu(VAR_est(J=5, p =2))
# 
# 
# #VAR Model
# VAR_mod= vars::VAR(VAR_series, p =p, type = "none")
# summary(VAR_mod)
# 
# 
# #Coef extraction
# est_coef <- coef(VAR_mod)
# 
# # Extract only the coefficients for both dependend variables
# # and combine them to a single matrix
# 
# est_coefs <- est_coef[[1]][, 1]
# 
# for (i in 1:(k-1)){
#   est_coefs <- rbind(est_coefs, est_coef[[i+1]][, 1])
# }
# 
# est_coefs=(as.matrix(round(est_coefs, 4)))
# 
# est_coefs
# A
# A-est_coefs
# 
# 
# VAR_series = as.data.frame(VAR_series) %>% 
#   mutate(Series1_treat = VAR_series[,1] + sd(VAR_series[,1]))
# 
# 
# 
# VAR_series = as.data.frame(VAR_series)







































# #### Dynamic
# 
# w_init = as.vector(1/(dim(VAR_series)[2]-2))
# 
# t(w_init) %*% y %*% w_init
# 
# y_0n = VAR_series[,1]
# y =  as.data.frame(VAR_series[,2:4])
# 
# lm(y_0n ~ VAR_series[,2:4])
