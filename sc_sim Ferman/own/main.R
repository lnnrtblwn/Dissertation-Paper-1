#Clears previously loaded files
rm(list=ls())

#Sets working directory
# setwd("/Users/luisalvarez/Dropbox/EESP/RA/synthetic_control/Simulations_Final/Tables1-3_new")


for(pack in c("Synth","nnls"))
  if(!require(pack, character.only = T))
  {
    install.packages(pack)
    require(pack, character.only = T)
  }

# source("first.R")


#Sets seed to allow for replication
set.seed(47)

#Number of replications
nreps = 5000

#Outcomes in pre-treatment are written as:
# y_{t} = \Mu F_t + \epsilon_t     (1)
#where y_{t} is (J+1)x1, \Mu is (J+1)xK, F_t is kx1

#Number of factors (K)
K = 2

#Group distribution of each factor (groups will be equally sized, with treated unit
#pertaining to first group)
group_distribution = list(
  "lambda1" = c(1,0),
  "lambda2" = c(0,1)
)


#Each factor is generated as a Gaussian AR(1)
#F_{kt} = \alpha + \rho F_{kt-1} + u_{t}  (2)

#Specify rho
rho = 0.5

#Intercept. Set it equal to mean*(1-rho) to define mean of process
alpha = 0*(1-rho)

#Specify variance of u_t. Set it to (1 - rho^2) will lead to var(\lambda^k_t) = 1
var_u = (1-rho^2)

# Specify variance of transitory shocks in Factor model equation 
var_epsilon = 1


#Number of post-treatment periods
T1 = 1

#Post-treatment effects (may change over time)
post_effect = 0


results = c()

# own to understand
J = 4
T_spec = 1
replication = 5


#Loop in T0 and J
for(J in c(4,10,50,100))
  for(T_spec in c(1,2))
    for(replication in 1:nreps){
      
      if(T_spec == 1)
        T0 = 2*J else T0 = J+5
        
        print(paste("J is ", J,". T0 is ",T0, ". T1 is ", T1,". Replication is ", replication, ".", sep=""))
        
        Factors = sapply(1:K, function(k){simulate_ar1(rho = rho, var_shock = var_u, T0 = T0+T1, intercept = alpha)}) 
        
        Mu = matrix(0, nrow = J+1, ncol = K)
        for(k in 1:K)
        {
          fac_distr = group_distribution[[k]]
          for(pos in 1:length(fac_distr))
            if(pos==1)
              Mu[1:(1 + J%/%length(fac_distr)),k] = fac_distr[pos] else 
                if(pos<length(fac_distr))
                  Mu[(2+(pos-1)*J%/%length(fac_distr)):(1 +pos*J%/%length(fac_distr)),k] = fac_distr[pos] else
                    Mu[(2+(pos-1)*J%/%length(fac_distr)):nrow(Mu),k] = fac_distr[pos]
                  
        }
        transitory_shocks = matrix(rnorm((J+1)*(T0+T1), sd = sqrt(var_epsilon)), nrow = (T0+T1), ncol = J+1)
        
        y = Factors%*%t(Mu) + transitory_shocks
        
        #Adding effect
        y[(T0+1):(T0+T1),1] = post_effect + y[(T0+1):(T0+T1),1] 
        
        y_before= y[1:T0,]
        y_after = y[(T0+1):(T0+T1),]
        
        
        #####1. Original SC
        synth_1 = synth_control_est(y_before, y_after)
        
        eff_1 = mean(synth_1$effects)
        
        weights_1 = apply(Mu, 2, function(x){
          sum(synth_1$w[as.logical(x[-1])])
        })
        names(weights_1) = paste("sc_sum_w_",names(group_distribution),sep="")
        
        ####2. OLS Unconstrained
        ols_unconstrained = lm(y_before[,1]~-1+y_before[,-1])$coefficients
        eff_unconstrained = mean(y_after%*%c(1,-ols_unconstrained))
        
        weights_unconstrained = apply(Mu, 2, function(x){
          sum(ols_unconstrained[as.logical(x[-1])])
        })
        names(weights_unconstrained) = paste("ols_unconstrained_sum_w_",names(group_distribution),sep="")
        
        
        ###3. OLS Adding up
        y_add_up = y_before[,1] - y_before[,2]
        X_add_up = y_before[,-c(1,2)] - y_before[,2]
        ols_add_up = lm(y_add_up~-1+X_add_up)$coefficients
        ols_add_up = c(1-sum(ols_add_up),ols_add_up)
        
        eff_add_up = mean(y_after%*%c(1,-ols_add_up))
        
        
        weights_add_up = apply(Mu, 2, function(x){
          sum(ols_add_up[as.logical(x[-1])])
        })
        names(weights_add_up) = paste("ols_add_up_sum_w_",names(group_distribution),sep="")
        
        
        
        results = rbind(results,
                        cbind("T0"=T0,"J" = J, "T1" = T1, "replication" = replication, "sc_eff" = eff_1, t(weights_1),
                              "ols_unconstrained_eff"=eff_unconstrained, t(weights_unconstrained),
                              "ols_add_up_eff"=eff_add_up, t(weights_add_up)))
        
        write.csv(results, "results.csv")
    }




