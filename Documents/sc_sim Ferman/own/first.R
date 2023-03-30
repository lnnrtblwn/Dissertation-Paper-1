#Loads required package, installing if needed
for(pack in c("nnls","quadprog"))
  if(!require(pack, character.only = T))
  {
    install.packages(pack)
    require(pack, character.only = T)
  }


#Function to simulate AR(1) for T0 periods. If is.na(start) == T, then initial value
#is drawn from stationary distribution. We discard the initial draw after that.
simulate_ar1 <- function(rho, var_shock, T0, intercept = 0, start = NA)
{
  y = ifelse(is.na(start), rnorm(1,mean=intercept/(1-rho),sd = sqrt(var_shock/(1-rho^2))), start) # generates RV with mean 0 and sd 1
  
  for(t in 1:T0)
    y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
  
  return(y[-1])
}

#Synthetic control estimator
synth_control_est <-function(y_before, y_after)
{
  y= y_before[,1]
  X =y_before[,-1]
  
  Dmat = t(X)%*%X
  dvec = t(X)%*%y
  
  Amat = t(rbind(rep(1,ncol(X)),diag(ncol(X))))
  bvec = c(1, rep(0,ncol(X)))
  
  synth_model = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  
  w = synth_model$solution
  effects = y_after%*%c(1,-w)
  
  return(list("w" = w, "effects" = effects, "value_min" = (2*synth_model$value +y%*%y) ))
  
}