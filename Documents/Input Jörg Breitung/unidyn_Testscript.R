library(dplyr)
library(ggplot2)



#### Input Param ####



# Number of pre-and post-treatment periods
T1 = 20
T0 = 50

# AR-Term in Factor model. y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
# rho = 0.8 // Non.Stationary => rho = 1.0
rho = 0.8

# Error AR-Term
rho_u = runif(1, .8, 0.99)

# Intercept. Set it equal to mean*(1-rho) to define mean of process
alpha = 0*(1-rho)

# Specify variance of u_t. Set it to (1 - rho^2) will lead to var(\lambda^k_t) = 1. Variance of the factors
var_u = (1-rho^2)

# Specify variance of transitory shocks in Factor model equation. Variance of the error terms 
var_epsilon = 1

# Post-treatment effects. Could be specified differently
post_effect = 10

# Number of factors
K = 2

#trend
c = 0

# Number of Donors
J=4


# Group distribution of each factor 
group_distribution = list(
  "lambda1" = c(1,0),
  "lambda2" = c(0,1))


treat_inter = 1



##### AR Simulations Funktion - Wird in der Faktor Simulation verwendet #####


simulate_ar1 = function(rho, var_shock, T0, intercept = 0, start = NA) {
  y = ifelse(is.na(start), rnorm(1,mean=intercept/(1-rho),sd = sqrt(var_shock/(1-rho^2))), start)
  
  for(t in 1:T0)
    y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
  
  return(y[-1])
}




#### Factor Simulation ####



  Mu = matrix(0, nrow = J+1, ncol = K)
  
  for(k in 1:K) {
    fac_distr = group_distribution[[k]]
    for(pos in 1:length(fac_distr))
      if(pos==1)
        Mu[1:(1 + J %/% length(fac_distr)),k] = fac_distr[pos] else 
          if (pos < length(fac_distr))
            Mu[(2 + (pos - 1)*J %/% length(fac_distr)):(1 + pos * J %/% length(fac_distr)),k] = fac_distr[pos] else
              Mu[(2+(pos-1)*J%/%length(fac_distr)):nrow(Mu),k] = fac_distr[pos]
            
  }
  
  Factors = sapply(1:K, function(k){simulate_ar1(rho = rho, var_shock = var_u, T0 = T0+T1, intercept = alpha)}) 
  
  # AR(1) shocks
  transitory_shocks = sapply(1:dim(Mu)[1], function(k){simulate_ar1(rho = rho_u, var_shock = 1, T0 = T0+T1, intercept = 0)}) 
  
  #Stationary Shocks
  #transitory_shocks = matrix(rnorm((J+1)*(T0+T1), sd = sqrt(var_epsilon)), nrow = (T0+T1), ncol = J+1)
  
  y = Factors%*%t(Mu) + transitory_shocks
  y = y + (1:nrow(y))^1.5*c
  

results = list()



# Define Series specific intercepts for Factor Simu
  my_means = c(rnorm(1, mean = treat_inter, sd = 1), rnorm(J, mean = 0, sd = 1))
  
  
  
  # interval check: 1 if treated series within donor series, 0 if not.
  # results[["bound_check"]] = ifelse(my_means[1] > min(my_means[c(FALSE, Mu[-1,1] == 1)]) & 
  #                                my_means[1] < max(my_means[c(FALSE, Mu[-1,1] == 1)]),
  #                              1, 0)
  
  results[["bound_check"]] = as.numeric(rank(colMeans(y))[1])
  
  
  for (i in 1:J) {
    y[,i] = y[,i] + my_means[i]
  }


# Adding effect
y[(T0+1):(T0+T1),1] = post_effect + y[(T0+1):(T0+T1),1] 

y_pre = y[1:T0,1]
y_post = y[(T0+1):(T0+T1),1]

x_pre = y[1:T0, -1]
x_post = y[(T0+1):(T0+T1), -1]







#### Dyn-Uni Model ####



# if(p_unidyn == "auto"){
#   p_uni = as.numeric(min(VARselect(y_pre)$selection))
# } else { p_uni = p_unidyn }

p_uni = 4

w0=matrix(1/J,J,1)
S0=100000000000000000000
S1=90000000000000000000
y0=y_pre[(p_uni+1):T0]

iter=0
while ((S0-S1)>0.00001) {
  z = x_pre %*% w0
  lagz = z[(p_uni+1):T0]
  for (i in (1:p_uni)) {
    lagz = cbind(lagz,z[(p_uni + 1 - i):(T0-i)])
  }
  outreg = lm(y0 ~ lagz)
  alphas = outreg$coefficients[2:(p_uni + 2)]
  
  x1 = alphas[1] * x_pre[(p_uni+1):T0,]
  for (i in (1:p_uni)) {
    x1 = x1 + alphas[i+1] * x_pre[(p_uni+1-i):(T0-i),]
  } 
  outreg = lm(y0~x1)
  w0 = outreg$coefficients[2:(J+1)]
  S0 = S1
  S1 = sd(outreg$residuals)
  iter = iter+1
}

y_unidyn_pre = outreg$fitted.values
# w_unidyn = w0  #  ohne Konstante!  
w_unidyn = outreg$coefficients # mit Konstante

# building x1_post
x1_post = alphas[1] * x_post

for (i in (1:p_uni)) {
  
  x1_post = x1_post + alphas[i+1] * rbind(x_pre[(T0+1-i):T0,], 
                                          x_post[1:(T1-i),])
  
} 

#y_unidyn_post = x1_post %*% w_unidyn   # Ohne Konstante
y_unidyn_post = cbind(1,x1_post) %*% w_unidyn # Mit Konstante

y_treat_unidyn = as.data.frame(c(y_pre, y_post)) %>%
  rename(y = c(1))

y_treat_unidyn$y_hat = c(rep(NA, p_uni), y_unidyn_pre, y_unidyn_post)


matplot(ts(y_treat_unidyn),
        type = "l",
        lty = 1,
        lwd = 2,
        main = paste0("Dynamic (Univar) Path. \n","rho_u = ", round(rho_u,4)),
        xlab = "Time",
        ylab = "Value")


sqrt(mean(((y_post-post_effect) - y_unidyn_post)^2))

seed = seed +1
set.seed(seed)
seed

# 24 57