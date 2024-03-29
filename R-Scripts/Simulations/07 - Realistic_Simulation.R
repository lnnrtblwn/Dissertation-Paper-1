

# Packages ----


library(dplyr)
library(tidyr)
library(ggplot2)
library(WDI)
library(lubridate)
library(corrplot)
library(vars)


# Data Input via World Bank Database ----


Country_List = WDI(indicator='NY.GDP.MKTP.KD.ZG', country="all", start=2020, end=2020)

#Indicator Codes for GDP growth rate, GDP and GDP per Capita 
indicator_code = c('NY.GDP.MKTP.KD.ZG','NY.GDP.MKTP.CD','NY.GDP.PCAP.CD')

# Choose Indicator Code: Here 3: GDP per Capita
ic = 3


# Country Presets Like G20, EU, Abadie2015 Country List ----

#g20=c('CN','IN','US','ID','BR','RU','JP','MX','DE','TR','FR','GB','IT','ZA','KR'
#      ,'AR','CA','SA','AU')
#euro_area = 'EMU'

#EU = c('BE','BG','CZ','DK','DE','EE','IE','GR','ES','FR','HR','IT','CY','LV','LT','LU','HU','MT','NL',
#       'AT','PL','PT','RO','SI','SK','FI','SE','XC')

#
Abadie15 = c('AT','NL','JP','US','CH', 'DE')


# Obtain Data unsing Country List and Indicator Code----
df = WDI(indicator=indicator_code[ic], country=Abadie15, start=1960, end=NULL)
#----

#Define Top Element such that the correlation Matrix has the desired structure later on
top_element <- "Germany"


df_spread = df %>%
  filter(year >= 1980) %>% 
  dplyr::select(year, country, indicator_code[ic]) %>% 
  pivot_wider(names_from = country, values_from = indicator_code[ic]) %>% 
  relocate(paste(top_element), .after = year)


cor=cor(as.matrix(df_spread[,-1]))
cov=cov(as.matrix(df_spread[,-1]))
corrplot(cor(as.matrix(df_spread[,-1])))



ggplot(df %>% 
         filter(year >= 1980), 
       aes(x = year, y = NY.GDP.PCAP.CD, group = country, color = country)) +
  geom_line() +
  theme_bw()

ggplot(df %>% 
         filter(year >= 1980), 
       aes(x = year, y = colnames(df)[5], group = country, color = country)) +
  geom_line() +
  theme_bw()





#Granger Causality
granger_matrix <- matrix(0, nrow = ncol(df_spread[,-1]), ncol = ncol(df_spread[,-1]))
for (i in 1:(ncol(df_spread[,-1])-1)) {
  for (j in (i+1):ncol(df_spread[,-1])) {
    # Subset the data to the two time series being tested
    subset_ts <- df_spread[,-1][,c(i,j)]
    # Run the VAR model with the maximum lag order
    var_fit <- VAR(subset_ts, p = 2, type = "const")
    # Perform the Granger causality test and extract the p-value
    gc_test <- causality(var_fit)
    p_value <- gc_test$Granger$p.value
    # Store the p-value in the appropriate position in the output matrix
    granger_matrix[i,j] <- p_value
    granger_matrix[j,i] <- p_value
  }
}

# View the resulting matrix of p-values
granger_matrix









#Optimal Lag Number
vars::VARselect(df_spread[-1])

#VAR Model
VAR_mod= vars::VAR(df_spread[-1], p =1, type = "const")
summary(VAR_mod)
plot(vars::fevd(VAR_mod, n.ahead = 10))
VAR_mod$varresult$Germany
synth = predict(VAR_mod, n.ahead = 50)$fcst



acf(resid(VAR_mod), 52)$acf
serial.test(VAR_mod, lags.pt=12, type="PT.adjusted")

(fit.pr = predict(VAR_mod, n.ahead = 24, ci = 0.95)) 
fanchart(fit.pr)



#Obtain Residuals
u=resid(VAR_mod)

#Obtaining Sigma as positive definite covariance matrix
Sigma = t(u)%*%u

#First row of Cholesky Factor of the inverse covariance matrix
q=chol(solve(Sigma))[1,]


#Obtaining Weights vector
w=c()

for (i in 2:length(q)) {w[i-1]=-q[i]/q[1]}





#Forecast Formula (1)


formula1 <- function(x, a1, b1) {
  return(mu +d )
}




#Forecast Formula (2) Distributed Lag of Synthetic control


install.packages("dLagM")

library(dLagM)

#Switching Algo Part1_1

data(M1Germany)
data = M1Germany[1:144,]
model.dlm  = dlm(formula = logprice ~ interest + logm1, 
                 data = data.frame(data) , q = 4)
summary(model.dlm)


formula = dlm(formula = y ~ x1 + x2 + x3, data = df, q = 4)



model.dlm  = dlm(formula = y ~ x1 + x2 + x3, 
                 data = data.frame(df) , q = 4)
summary(model.dlm)



model.dlm = df$y ~ mu + sum(delta_0 (1 + exp(-r * (x - t))))


nls(model.dlm, data = df, start = list(K = 10, r = 0.5, t = 5))




# Create the dummy variable
d <- c(rep(0, nrow(data) / 2), rep(1, nrow(data) / 2))

# Fit the regression model
model <- lm(y ~ x1 + x2 + x3 -1, data = df)

summary(model)





###### ### ### ###
### VAR Simu ###
################


## parameter and error terms used throughout
A_1 <- matrix(c(0.8,0.3,0.1,0.5),nrow=2)
A_2 <- matrix(c(0.4,0.2,0.5,0.9),nrow=2)
e <- matrix(rnorm(10000),ncol=2)

## Let's start with the R version
VAR_sim <- function(coeff, errors) {
  simdata <- matrix(0, nrow(errors), ncol(errors))
  for (row in 2:nrow(errors)) {
    simdata[row,] = coeff %*% simdata[(row-1),] + errors[row,]
  }
  return(simdata)
}

Data_VAR_sim <- VAR_sim(A_1, e)    

vars::VAR(Data_VAR_sim, p = 1)



compRsim <- compiler::cmpfun(VAR_sim)
compRData <- compRsim(A_1,e)              # generated by R 'compiled'
stopifnot(all.equal(Data_VAR_sim, compRData))




### Alternative using simulateVAR Package
Data = simulateVAR(pars = list(A_1,A_2), lags = 2, Nt = 100000)

vars::VAR(Data, p = 2)


Model <- mlVARsim(nPerson = 100, nNode = 3, nTime = 50, lag=1)
# Estimate using correlated random effects:
fit <- mlVAR(Model$Data, vars = Model$vars,
             idvar = Model$idvar, lags = 1,
             temporal = "correlated")
# Sample from fitted model:
samples <- mlVARsample(fit, nTime = 50, nSample = 50, pMissing = 0.1,
                       nReps = 5, nCores = 1)
# Summarize result







### Alternative usind VAR.sim
B1<-matrix(c(0.7, 0.2, 0.2, 0.7), 2)
var1 <- VAR.sim(B=B1, n=100, include="none")
ts.plot(var1, type="l", col=c(1,2))


B2<-rbind(c(0.5, 0.5, 0.5), c(0, 0.5, 0.5))
varcov<-matrix(c(1,0.2, 0.3, 1),2)
var2 <- VAR.sim(B=B2, n=100, include="const", varcov=varcov)
ts.plot(var2, type="l", col=c(1,2))

## VAR.boot: Bootstrap a VAR 
data(zeroyld)
mod <- lineVar(data=zeroyld,lag=1)
VAR.boot(mod)





# Ablage ----




colMeans(df_spread)[-1]
n = 50
T0 = 30                       
i = 1

mean = colMeans(as.matrix(df_spread[,-1]))
cor = cor(as.matrix(df_spread[,-1]))
cov = cov(as.matrix(df_spread[,-1]))


mean = c(1,1,1)
cov = matrix(c(1,0.1,0.4,0.1,1,0.5,0.4,0.5,1), nrow= 3, ncol = 3)

df_synth = as.data.frame(MASS::mvrnorm(n = n, mu = mean, Sigma = cov))

treat = c(rep(0, T0), rep(sqrt(cov[2,2]), n-T0))



eigen(cov)$values



q=chol(solve(cov))[1,]


w=c()

for (i in 2:length(q)) {w[i-1]=-q[i]/q[1]}



mu_s_j = 0.667

df_synth = df_synth %>% 
  mutate(T = as.numeric(seq(1:n))) %>% 
  mutate(treat = c(rep(0, T0), rep(sqrt(cov[i,i]), n-T0))) %>% 
  mutate(Germany_treat = Germany + treat) %>%
  dplyr::select(-treat)


test = df_synth %>% gather("Country", "GDP", -T)


ggplot(test) +
  aes(x = T, y = GDP, colour = Country) +
  geom_line() +
  scale_color_hue(direction = 1) +
  theme_bw()

df_synth = df_synth %>% 
  mutate(v1_cf = w1*V2 + w2*V3 + mu_s) %>% 
  mutate(v1_cf_j = w1_j*V2 + w2_j*V3 + mu_s_j)

test2 = df_synth %>% select(-c("V2","V3")) %>% gather("Country", "GDP", -T) 

plot_joergweights = {
  ggplot(test2) +
    aes(x = T, y = GDP, colour = Country) +
    geom_line() +
    scale_color_hue(direction = 1) +
    theme_bw()
}


mean(df_synth$v1_cf)
mean(df_synth$V1)

vector1[i] = var(df_synth$V1-df_synth$v1_cf)
vector2[i] = var(df_synth$V1-df_synth$v1_cf_j)


df_synth_after = df_synth %>% 
  filter(T>=T0) %>% 
  mutate(diff = v_treat - v1_cf) %>% 
  mutate(diff_j = v_treat - v1_cf_j )

vector3[i] = mean(df_synth_after$diff)
vector4[i] = mean(df_synth_after$diff_j)

