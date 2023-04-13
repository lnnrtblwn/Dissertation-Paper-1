setwd("~/Diss/Topics/Synthetic Control/Documents/sc_sim Ferman/Ferman")
source("my_functions.R")

# 1. DATA GENERATING PROCESS: FACTOR MODEL WITHOUT COVARIATES ---- 

# Number of donors
J = 3

# Number of pre-and post-treatment periods
T1 = 50
T0 = T1

# AR-Term in Factor model. y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
rho = 0.5

# Intercept. Set it equal to mean*(1-rho) to define mean of process
#alpha = 3
alpha = 0*(1-rho)

# Specify variance of u_t. Set it to (1 - rho^2) will lead to var(\lambda^k_t) = 1. Variance of the factors
var_u = (1-rho^2)

# Specify variance of transitory shocks in Factor model equation. Variance of the error terms 
var_epsilon = 1

# Post-treatment effects. Could be specified differently. Easier to assesss counterfactual with flat effect
post_effect = 10

# Number of factors
K = 2

# Adding a trend
c = 0

# Group distribution of each factor 
group_distribution = list(
  "lambda1" = c(1,0),
  "lambda2" = c(0,1))

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

transitory_shocks = matrix(rnorm((J+1)*(T0+T1), sd = sqrt(var_epsilon)), nrow = (T0+T1), ncol = J+1)

y = Factors%*%t(Mu) + transitory_shocks
y = y + (1:nrow(y))*c

# Adding effect
y[(T0+1):(T0+T1),1] = post_effect + y[(T0+1):(T0+T1),1] 

y_pre = y[1:T0,1]
y_post = y[(T0+1):(T0+T1),1]

x_pre = y[1:T0, -1]
x_post = y[(T0+1):(T0+T1), -1]

matplot(ts(y),
        type = "l",
        lty = 1,
        lwd = 2,
        main = "Different Simulations",
        xlab = "Time",
        ylab = "Value")

rm(list = setdiff(ls(), c("y", "y_pre", "y_post", "x_pre", "x_post", "Mu")))

# 2. DATA GENERATING PROCESS: VAR-MODEL ---- 

# to be done.

# 3. ESTIMATION ----

Dmat = t(x_pre) %*% x_pre
dvec = t(x_pre) %*% y_pre
Amat = t(rbind(rep(1, ncol(x_pre)), diag(ncol(x_pre)), -1*diag(ncol(x_pre))))
bvec = c(1, rep(0, ncol(x_pre)), rep(-1,ncol(x_pre)))
round(quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)$solution,3)

summary(lm(y_pre ~ x_pre))
