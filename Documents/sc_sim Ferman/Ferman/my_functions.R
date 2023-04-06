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
  y = ifelse(is.na(start), rnorm(1,mean=intercept/(1-rho),sd = sqrt(var_shock/(1-rho^2))), start)
  
  for(t in 1:T0)
    y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
  
  return(y[-1])
}

#Synthetic control estimator
synth_control_est <-function(y_before, y_after)
{
  y= y_before[,1]
  X =y_before[,-1]
  
  # D is the square matrix in the middle of quadratic form multiplied by 2
  Dmat = t(X)%*%X
  
  # d is the column vector specifying the linear part in the objective function
  dvec = t(X)%*%y
  
  # first col specifies sum, remaining cols specify individual coefficients
  Amat = t(rbind(rep(1,ncol(X)),diag(ncol(X)),-1*diag(ncol(X))))
  
  # b is the column vector of decision variables. 
  # First element specifies coefficient sum
  # Second element specifies minimum value the coefficients can take on. Here 0
  # Third element specifies maximum value the coefficients can take on. Here 1
  bvec = c(1, rep(0, ncol(X)), rep(-1,ncol(X)))
  
  synth_model = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  
  w = synth_model$solution
  round(w,2)
  sum(round(w,2))
  
  # Difference between observation (y_after) and counterfactual (synthetic control)
  effects = y_after %*% c(1,-w)
  
  # w = SC-weights
  # effects = estimated treatment effect
  # value_min = some measure of precision. synth_model$value = scalar-solution of minimization, y = pre-treatment-values of treated unit 
  # don't understand this statistic entirely. 
  return(list("w" = w, "effects" = effects, "value_min" = (2*synth_model$value +y%*%y) ))
  
}
