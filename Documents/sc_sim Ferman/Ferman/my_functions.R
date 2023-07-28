for(pack in c("nnls","quadprog", "glmnet"))
  if(!require(pack, character.only = T))
  {
    install.packages(pack)
    require(pack, character.only = T)
  }
rm(pack)

#Function to simulate AR(1) for T0 periods. If is.na(start) == T, then initial value
#is drawn from stationary distribution. We discard the initial draw after that.

simulate_ar1 = function(rho, var_shock, T0, intercept = 0, start = NA) {
  y = ifelse(is.na(start), rnorm(1,mean=intercept/(1-rho),sd = sqrt(var_shock/(1-rho^2))), start)
  
  for(t in 1:T0)
    y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
  
  return(y[-1])
}

#Synthetic control estimator
synth_control_est <-function(y_before, y_after) {
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

regularized_ols1 = function(x, y, l1, l2) {
  
  Y = y %>%
    as.matrix()
  
  I = diag(ncol(x))
  U = matrix(data = 1,
             nrow = ncol(x),
             ncol = 1)
  
  X = x %>%
    as.matrix()
  
  beta = solve(t(X) %*% X + l1 * I + l2 * (U %*% t(U))) %*% ((t(X) %*% Y) + l2 * U)
  # solve(t(X) %*% X) %*% (t(X) %*% Y)
  y_hat = X %*% beta
  
  synth_out = list()
  
  synth_out$y_hat = y_hat
  synth_out$beta = beta
  
  precision = c(sqrt(mean((Y - y_hat) ^ 2)),
                mean(abs(Y - y_hat)))
  names(precision) = c("RMSPE", "MAPE")
  
  synth_out$precision = precision
  
  return(synth_out)
}
regularized_ols2 = function(x, y, l1, l2) {

  Y = y %>%
    as.matrix()

  I = diag(ncol(x) + 1)
  U = matrix(data = 1,
             nrow = ncol(x) + 1,
             ncol = 1)

  X = x %>%
    mutate(x0 = 1) %>%
    dplyr::select(x0, everything()) %>%
    as.matrix()

  I[1, 1] = 0
  U[1, 1] = 0

  beta = solve(t(X) %*% X + l1 * I + l2 * (U %*% t(U))) %*% ((t(X) %*% Y) + l2 * U)

  return(beta)
}

REGOLS_LASSO <- function(X, y, lambda1, lambda2, max_iter = 1000) {
  
  p <- ncol(X)
  beta_init <- rep(0, p)
  
  result <- optim(beta_init,
                  fn = loss, 
                  gr = gradient,
                  X = X,
                  y = y,
                  lambda1 = lambda1,
                  lambda2 = lambda2,
                  method = "L-BFGS-B",
                  control = list(maxit = max_iter))  # Maximum number of iterations
  
  return(result$par)  # Return the optimized coefficients
}
loss <- function(beta, X, y, lambda1, lambda2) {
  
  n = nrow(X)
  lambda3 = mean(c(lambda1 ,lambda2))
  
  mse <- (1/n) * sum((y - X %*% beta)^2) 
  
  # Lasso
  l1 <- exp(lambda1) * sum(abs(beta[-1]))
  
  # Ridge
  #l2 <- exp(lambda2) * sum(beta[-1]^2) 
  
  # Own
  l2 <- exp(lambda2) * (1 - sum(beta[-1]))^2 
  
  return(mse + l1 + l2)
}
gradient <- function(beta, X, y, lambda1, lambda2) {
  
  n = nrow(X)
  lambda3 = mean(c(lambda1 ,lambda2))
  
  mse <- (-2/n) * t(X) %*% (y - X %*% beta) 
  
  # Lasso
  l1 <- c(0, exp(lambda1) * sign(beta[-1]))
  
  # Ridge
  #l2 <- c(0, 2 * exp(lambda2) * beta[-1])
  
  # Own
  l2 <- c(0, -2 * exp(lambda2) * (1 - sum(beta[-1])))
  
  return(mse + l1 + l2)
}

simulation_factor = function(J, simu_type = 'Factor'){
  
  #### ab hier ----
  
  if (simu_type == 'Factor'){
  
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
  
  } else {
    est_coefs = VAR_est(J=J, p = p)
    stat_test = Stat_test(est_coefs)
    y = tail(VAR_simu(est_coefs),(T0+T1))
    
    while (stat_test > 0.99) {
      est_coefs = VAR_est(J=J, p = p)
      stat_test = Stat_test(est_coefs)
      y = tail(VAR_simu(est_coefs),(T0+T1))
    }
  }
  
  results = list()
  
  # Define Series specific intercepts for Factor Simu
  if (simu_type == 'Factor'){
  # curve(dnorm(x, 0, 10), from=-4, to=4)
  my_means = c(rnorm(1, mean = treat_inter, sd = 1), rnorm(J, mean = 0, sd = 1))
  
  
  
  # interval check: 1 if treated series within donor series, 0 if not.
  # results[["bound_check"]] = ifelse(my_means[1] > min(my_means[c(FALSE, Mu[-1,1] == 1)]) & 
  #                                my_means[1] < max(my_means[c(FALSE, Mu[-1,1] == 1)]),
  #                              1, 0)
  
  # Das macht an dieser Stelle doch gar keinen Sinn? Die Means werden doch erst unten addiert. Nach unten kopiert
  # results[["bound_check"]] = as.numeric(rank(colMeans(y))[1])
  
  
  for (i in 1:J) {
    y[,i] = y[,i] + my_means[i]
  }
  
  results[["bound_check"]] = as.numeric(rank(colMeans(y))[1])
  } else {
    # if (colMeans(y)[1] > max(colMeans(y)[-1])) {results[["bound_check"]] = 1
    # } else if (colMeans(y)[1] < min(colMeans(y)[-1])) {results[["bound_check"]] = -1
    # } else {results[["bound_check"]] = 0}
    results[["bound_check"]] = as.numeric(rank(colMeans(y))[1])
  }
  
  # Adding effect
  y[(T0+1):(T0+T1),1] = post_effect + y[(T0+1):(T0+T1),1] 
  
  y_pre = y[1:T0,1]
  y_post = y[(T0+1):(T0+T1),1]
  
  x_pre = y[1:T0, -1]
  x_post = y[(T0+1):(T0+T1), -1]
  
  # matplot(ts(y),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "Overview: Simulated Data",
  #         xlab = "Time",
  #         ylab = "Value")

  #rm(list = setdiff(ls(), c("y", "y_pre", "y_post", "x_pre", "x_post", "Mu", "T0", "T1")))
  
  # ESTIMATION
  
  # SC
  Dmat = t(x_pre) %*% x_pre
  dvec = t(x_pre) %*% y_pre
  Amat = t(rbind(rep(1, ncol(x_pre)), diag(ncol(x_pre)), -1*diag(ncol(x_pre))))
  bvec = c(1, rep(0, ncol(x_pre)), rep(-1,ncol(x_pre)))
  
  synth_model = tryCatch({
    quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)}, 
    error = function(e) {
    # Return a default value to indicate failure
    return(NA)
      })
  
  if (any(is.na(synth_model))) {
    y_sc_pre = rep(NA, T0)
    y_sc_post = rep(NA, T1)
  } else {
    w_sc = synth_model$solution
    y_sc_pre = x_pre %*% w_sc
    y_sc_post = x_post %*% w_sc}
    
  y_treat_sc = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_sc$y_hat = c(y_sc_pre, y_sc_post)
  
  # matplot(ts(y_treat_sc),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "SC Path",
  #         xlab = "Time",
  #         ylab = "Value")

  results_SC = c()
  
  results_SC["PRE_SC_RMSPE"] = sqrt(mean((y_pre - y_sc_pre)^2)) 
  results_SC["PRE_SC_BIAS"] = mean(y_sc_pre - y_pre)
  results_SC["PRE_SC_VAR"] = mean((y_sc_pre - mean(y_sc_pre))^2)
  
  results_SC["POST_SC_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_sc_post)^2))
  results_SC["POST_SC_BIAS"] = mean(y_sc_post - (y_post-post_effect))
  results_SC["POST_SC_VAR"] = mean((y_sc_post - mean(y_sc_post))^2)
  
  results[["SC"]] = results_SC
  
  
  # results[2] = sqrt(mean((y_pre - y_sc_pre)^2))
  # results[3] = sqrt(mean(((y_post-post_effect) - y_sc_post)^2))
  
  # OLS
  
  if(ncol(x_pre) >= nrow(x_pre)){
    w_ols = rep(NA, ncol(x_pre)+1)
  } else {
    w_ols = summary(lm(y_pre ~ x_pre))$coefficients[,1]
  }
  
  y_ols_pre = cbind(rep(1, (T0)), x_pre) %*% w_ols
  y_ols_post = cbind(rep(1, (T1)), x_post) %*% w_ols

  y_treat_ols = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_ols$y_hat = c(y_ols_pre, y_ols_post)

  # matplot(ts(y_treat_ols),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "OLS Path",
  #         xlab = "Time",
  #         ylab = "Value")

  results_OLS = c()
  
  results_OLS["PRE_OLS_RMSPE"] = sqrt(mean((y_pre - y_ols_pre)^2)) 
  results_OLS["PRE_OLS_BIAS"] = mean(y_ols_pre - y_pre)
  results_OLS["PRE_OLS_VAR"] = mean((y_ols_pre - mean(y_ols_pre))^2)
  
  results_OLS["POST_OLS_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_ols_post)^2))
  results_OLS["POST_OLS_BIAS"] = mean(y_ols_post - (y_post-post_effect))
  results_OLS["POST_OLS_VAR"] = mean((y_ols_post - mean(y_ols_post))^2)
  
  results[["OLS"]] = results_OLS
  
  # results[4] = sqrt(mean((y_pre - y_ols_pre)^2))
  # results[5] = sqrt(mean(((y_post-post_effect) - y_ols_post)^2))
  
  # Regularized OLS. Simple Cross-Validation with 2-folds
  
  # x, y as data frame, de-meaned s.t. constant is independent of regularization
  # l1, l2 as numerics
  # intercept as boolean
  
  # random hyperparameter grid search
  
  param_grid = expand.grid(
    l1 = 5^seq(1, 5, length.out = 50), 
    l2 = 10^seq(1, 7, length.out = 50)) %>% 
    sample_n(400)
  
  param_grid$coeff_sum = NA
  param_grid$RMSPE = NA
  param_grid$RMSFE = NA
  
  for (grid in 1:nrow(param_grid)) {
    
    current_params = param_grid[grid, ]
    
    # Prepare demeaned series
    my = mean(y_pre[1:(length(y_pre)*(CV_share))])
    y0 = y_pre[1:(length(y_pre)*(CV_share))] - my
    
    mx0 = colMeans(x_pre[1:(nrow(x_pre)*(CV_share)),])
    X0 = x_pre[1:(nrow(x_pre)*(CV_share)),] 
    for (i in seq_along(mx0)) {
      X0[, i] = X0[, i] - mx0[i]
    }
    
    mx1 = colMeans(x_pre[((length(y_pre) * (CV_share)) + 1):length(y_pre),])
    X1 = x_pre[((length(y_pre) * (CV_share)) + 1):length(y_pre),]
    for (i in seq_along(mx1)) {
      X1[, i] = X1[, i] - mx1[i]
    }
      
    model = regularized_ols1(
      x = X0 %>% 
        as.data.frame(),
      
      y = y0 %>% 
        as.data.frame(),
      
      l1 = current_params[["l1"]],
      l2 = current_params[["l2"]])
      
    # train performance
    
    y_regols_pre_CV1 = my + X0 %*% model$beta
    y_pre_CV1 = y_pre[1:(length(y_pre)*(CV_share))]
    param_grid$RMSPE[grid] = sqrt(mean((y_pre_CV1 - y_regols_pre_CV1)^2)) 
    
    # test performance
    
    y_regols_pre_CV2 = my + X1 %*% model$beta
    y_pre_CV2 = y_pre[((length(y_pre) * (CV_share)) + 1):length(y_pre)]
    param_grid$RMSFE[grid] = sqrt(mean((y_pre_CV2 - y_regols_pre_CV2)^2)) 

  }
  
  # here comes the second step
  
  param_grid[which.min(param_grid$RMSFE),]
  
  param_grid_2nd = param_grid[which.min(param_grid$RMSFE),] %>% 
    bind_rows(expand.grid(
      l1 = c(param_grid[which.min(param_grid$RMSFE),1],
             param_grid[which.min(param_grid$RMSFE),1] + 2^seq(1, 5, length.out = 5),
             param_grid[which.min(param_grid$RMSFE),1] - 2^seq(1, 5, length.out = 5)),
      l2 = c(param_grid[which.min(param_grid$RMSFE),2],
             param_grid[which.min(param_grid$RMSFE),2] + 5^seq(1, 10, length.out = 5),
             param_grid[which.min(param_grid$RMSFE),2] - 5^seq(1, 10, length.out = 5)))) %>% 
    filter(l1 > 0,
           l2 > 0)
  
  for (grid in 1:nrow(param_grid_2nd)) {
    
    current_params = param_grid_2nd[grid, ]
    
    # Prepare demeaned series
    my = mean(y_pre[1:(length(y_pre)*(CV_share))])
    y0 = y_pre[1:(length(y_pre)*(CV_share))] - my
    
    mx0 = colMeans(x_pre[1:(nrow(x_pre)*(CV_share)),])
    X0 = x_pre[1:(nrow(x_pre)*(CV_share)),] 
    for (i in seq_along(mx0)) {
      X0[, i] = X0[, i] - mx0[i]
    }
    
    mx1 = colMeans(x_pre[((length(y_pre) * (CV_share)) + 1):length(y_pre),])
    X1 = x_pre[((length(y_pre) * (CV_share)) + 1):length(y_pre),]
    for (i in seq_along(mx1)) {
      X1[, i] = X1[, i] - mx1[i]
    }
    
    model = regularized_ols1(
      x = X0 %>% 
        as.data.frame(),
      
      y = y0 %>% 
        as.data.frame(),
      
      l1 = current_params[["l1"]],
      l2 = current_params[["l2"]])
    
    # train performance
    
    y_regols_pre_CV1 = my + X0 %*% model$beta
    y_pre_CV1 = y_pre[1:(length(y_pre)*(CV_share))]
    param_grid_2nd$RMSPE[grid] = sqrt(mean((y_pre_CV1 - y_regols_pre_CV1)^2)) 
    
    # test performance
    
    y_regols_pre_CV2 = my + X1 %*% model$beta
    y_pre_CV2 = y_pre[((length(y_pre) * (CV_share)) + 1):length(y_pre)]
    param_grid_2nd$RMSFE[grid] = sqrt(mean((y_pre_CV2 - y_regols_pre_CV2)^2)) 
  }
  
  best_params = param_grid_2nd[which.min(param_grid_2nd$RMSFE),] 
  
  # now extract cv-parameter combination
  
  # Prepare demeaned series
  my = mean(y_pre)
  y1 = y_pre - my
  
  mx1 = colMeans(x_pre)
  X1 = x_pre
  for (i in seq_along(mx1)) {
    X1[, i] = X1[, i] - mx1[i]
  }
  
  mx1_post = colMeans(x_post)
  X1_post = x_post
  for (i in seq_along(mx1_post)) {
    X1_post[, i] = X1_post[, i] - mx1_post[i]
  }
  
  
  w_regols = regularized_ols1(
    x = X1 %>% 
      as.data.frame(),
    
    y = y1 %>% 
      as.data.frame(),
    
    l1 = best_params[["l1"]],
    l2 = best_params[["l2"]])$beta
  
  y_treat_regols = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  
  y_regols_pre = my + X1 %*% w_regols
  y_regols_post = my + X1_post %*% w_regols
  
  y_treat_regols$y_hat = c(y_regols_pre, 
                           y_regols_post)
    
  matplot(ts(y_treat_regols),
          type = "l",
          lty = 1,
          lwd = 2,
          main = "Regularized OLS Path (Jörg)",
          xlab = "Time",
          ylab = "Value")

  results_REGOLS1 = c()
  
  results_REGOLS1["PRE_REGOLS_RMSPE"] = sqrt(mean((y_pre - y_regols_pre)^2)) 
  results_REGOLS1["PRE_REGOLS_BIAS"] = mean(y_regols_pre- y_pre)
  results_REGOLS1["PRE_REGOLS_VAR"] = mean((y_regols_pre - mean(y_regols_pre))^2)
  
  results_REGOLS1["POST_REGOLS_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_regols_post)^2))
  results_REGOLS1["POST_REGOLS_BIAS"] = mean(y_regols_post - (y_post-post_effect))
  results_REGOLS1["POST_REGOLS_VAR"] = mean((y_regols_post - mean(y_regols_post))^2)
  
  results[["REGOLS1"]] = results_REGOLS1
  
  # REGOLS2
  
  # Regularized OLS. Simple Cross-Validation with 2-folds
  
  # x, y as data frame
  # l1, l2 as numerics
  # intercept as boolean
  
  # random hyperparameter grid search
  
  param_grid = expand.grid(
    l1 = 5^seq(1, 5, length.out = 50),
    l2 = 10^seq(1, 7, length.out = 50)) %>%
    as.data.frame() %>%
    sample_n(400) %>%
    bind_rows(c(l1 = 1,
                l2 = 1))
  
  #param_grid$coeff_sum = NA
  param_grid$RMSPE = NA
  #param_grid$MAPE = NA
  param_grid$RMSFE = NA
  #param_grid$MAFE = NA
  
  for (grid in 1:nrow(param_grid)) {
    
    current_params = param_grid[grid, ]
    
    model = tryCatch({
      regularized_ols2(
        x = x_pre[1:(nrow(x_pre)*(CV_share)),] %>% 
          as.data.frame(),
        y = y_pre[1:(length(y_pre)*(CV_share))] %>% 
          as.data.frame(),
        l1 = current_params[["l1"]],
        l2 = current_params[["l2"]])},
      
      error = function(e) {
        # Error handling
        # print("An error occurred:")
        # print(e$message)
        # Return value to indicate failure
        return(list(NA))
      })

      # train performance
      
      y_train = y_pre[1:(length(y_pre) * (CV_share))] %>%
        as.matrix()
      x_train = x_pre[1:(nrow(x_pre) * (CV_share)), ] %>%
        as.matrix()
      
      y_regols_pre_CV1 = as.matrix(cbind(1, x_train)) %*% model
      
      param_grid$RMSPE[grid] = sqrt(mean((y_train - y_regols_pre_CV1)^2))
      
      # test performance
      
      y_test = y_pre[(1+length(y_pre) * (CV_share)):length(y_pre)] %>%
        as.matrix()
      x_test = x_pre[(1+nrow(x_pre) * (CV_share)):nrow(x_pre), ] %>%
        as.matrix()
      
      y_regols_pre_CV2 = as.matrix(cbind(1, x_test)) %*% model
      
      param_grid$RMSFE[grid] = sqrt(mean((y_test - y_regols_pre_CV2) ^ 2))

  }
  
  # here comes the second step
  
  param_grid[which.min(param_grid$RMSFE),]
  
  param_grid_2nd = param_grid[which.min(param_grid$RMSFE),] %>%
    bind_rows(expand.grid(
      l1 = c(param_grid[which.min(param_grid$RMSFE),1],
             param_grid[which.min(param_grid$RMSFE),1] + 2^seq(1, 5, length.out = 5),
             param_grid[which.min(param_grid$RMSFE),1] - 2^seq(1, 5, length.out = 5)),
      l2 = c(param_grid[which.min(param_grid$RMSFE),2],
             param_grid[which.min(param_grid$RMSFE),2] + 5^seq(1, 10, length.out = 5),
             param_grid[which.min(param_grid$RMSFE),2] - 5^seq(1, 10, length.out = 5)))) %>%
    filter(l1 > 0,
           l2 > 0)
  
  for (grid in 1:nrow(param_grid_2nd)) {
    
    current_params = param_grid_2nd[grid, ]
    
    model = tryCatch({
      regularized_ols2(
        x = x_pre[1:(nrow(x_pre)*(CV_share)),] %>% 
          scale(center = FALSE, scale = FALSE) %>% 
          as.data.frame(),
        y = y_pre[1:(length(y_pre)*(CV_share))] %>% 
          scale(center = FALSE, scale = FALSE) %>% 
          as.data.frame(),
        l1 = current_params[["l1"]],
        l2 = current_params[["l2"]])},
      
      error = function(e) {
        return(list(NA))
      })
    
    # train performance
    
    y_train = y_pre[1:(length(y_pre) * (CV_share))] %>%
      as.matrix()
    x_train = x_pre[1:(nrow(x_pre) * (CV_share)), ] %>%
      as.matrix()
    
    y_regols_pre_CV1 = as.matrix(cbind(1, x_train)) %*% model
    
    param_grid_2nd$RMSPE[grid] = sqrt(mean((y_train - y_regols_pre_CV1)^2))
    
    # test performance
    
    y_test = y_pre[(1+length(y_pre) * (CV_share)):length(y_pre)] %>%
      as.matrix()
    x_test = x_pre[(1+nrow(x_pre) * (CV_share)):nrow(x_pre), ] %>%
      as.matrix()
    
    y_regols_pre_CV2 = as.matrix(cbind(1, x_test)) %*% model
    
    param_grid_2nd$RMSFE[grid] = sqrt(mean((y_test - y_regols_pre_CV2) ^ 2))
  }
  
  best_params = param_grid_2nd[which.min(param_grid_2nd$RMSFE),] 
  best_params_REGOLS = best_params
  
  # now extract cv-parameter combination
  
  w_regols = regularized_ols2(
    x = x_pre %>%
      as.data.frame(),
    y = y_pre %>%
      as.data.frame(),
    l1 = best_params[["l1"]],
    l2 = best_params[["l2"]])
  
  y_regols_pre = as.matrix(cbind(1, x_pre)) %*% w_regols
  y_regols_post = as.matrix(cbind(1, x_post)) %*% w_regols 

  y_treat_regols = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_regols$y_hat = c(y_regols_pre, y_regols_post)
  
  matplot(ts(y_treat_regols),
          type = "l",
          lty = 1,
          lwd = 2,
          main = "Regularized OLS Path (My)",
          xlab = "Time",
          ylab = "Value")
  
  results_REGOLS2 = c()
  
  results_REGOLS2["PRE_REGOLS_RMSPE"] = sqrt(mean((y_pre - y_regols_pre)^2)) 
  results_REGOLS2["PRE_REGOLS_BIAS"] = mean(y_regols_pre- y_pre)
  results_REGOLS2["PRE_REGOLS_VAR"] = mean((y_regols_pre - mean(y_regols_pre))^2)
  
  results_REGOLS2["POST_REGOLS_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_regols_post)^2))
  results_REGOLS2["POST_REGOLS_BIAS"] = mean(y_regols_post - (y_post-post_effect))
  results_REGOLS2["POST_REGOLS_VAR"] = mean((y_regols_post - mean(y_regols_post))^2)
  
  results[["REGOLS2"]] = results_REGOLS2
  
  # REGOLS_LASSO
  # 
  # param_grid = expand.grid(
  #   lambda1 = seq(-20, 10, length.out = 200), 
  #   lambda2 = seq(-20, 10, length.out = 200)) %>% 
  #   sample_n(500) %>% 
  #   mutate(RMSPE = NA,
  #          RMSFE = NA)
  # 
  # for (grid in 1:nrow(param_grid)) {
  #   
  #   current_params = param_grid[grid, ]
  #   
  #   model = REGOLS_LASSO(
  #     X = cbind(1, x_pre[1:(nrow(x_pre)*(CV_share)),]),
  #     y = y_pre[1:(length(y_pre)*(CV_share))],
  #     
  #     lambda1 = current_params[["lambda1"]],
  #     lambda2 = current_params[["lambda2"]])
  #   
  #   y_hat_train = cbind(1, x_pre[1:(nrow(x_pre)*(CV_share)),]) %*% model  
  #   y_hat_test = cbind(1, x_pre[((nrow(x_pre) * (CV_share)) + 1):nrow(x_pre), ]) %*% model
  #   
  #   # train performance
  #   
  #   param_grid[grid,3] = sqrt(mean((y_pre[1:(length(y_pre)*(CV_share))] - y_hat_train)^2))
  #   
  #   # test performance
  #   
  #   param_grid[grid,4] = sqrt(mean((y_pre[((length(y_pre) * (CV_share)) + 1):length(y_pre)] - y_hat_test)^2))
  #   
  # }
  # 
  # param_grid$sum = param_grid$lambda1 + param_grid$lambda2
  # 
  # param_grid_2nd = param_grid[which.min(param_grid$RMSFE),] %>% 
  #   bind_rows(expand.grid(
  #     lambda1 = c(param_grid[which.min(param_grid$RMSFE),1],
  #                 param_grid[which.min(param_grid$RMSFE),1] + seq(1, 1.5, length.out = 5),
  #                 param_grid[which.min(param_grid$RMSFE),1] - seq(1, 1.5, length.out = 5)),
  #     lambda2 = c(param_grid[which.min(param_grid$RMSFE),2],
  #                 param_grid[which.min(param_grid$RMSFE),2] + seq(1, 1.5, length.out = 5),
  #                 param_grid[which.min(param_grid$RMSFE),2] - seq(1, 1.5, length.out = 5)))) %>% 
  #   slice(-1)
  # 
  # 
  # for (grid in 1:nrow(param_grid_2nd)) {
  #   
  #   current_params = param_grid_2nd[grid, ]
  #   
  #   model = REGOLS_LASSO(
  #     X = cbind(1, x_pre[1:(nrow(x_pre)*(CV_share)),]),
  #     y = y_pre[1:(length(y_pre)*(CV_share))],
  #     
  #     lambda1 = current_params[["lambda1"]],
  #     lambda2 = current_params[["lambda2"]])
  #   
  #   y_hat_train = cbind(1, x_pre[1:(nrow(x_pre)*(CV_share)),]) %*% model  
  #   y_hat_test = cbind(1, x_pre[((nrow(x_pre) * (CV_share)) + 1):nrow(x_pre), ]) %*% model
  #   
  #   # train performance
  #   
  #   param_grid_2nd[grid,3] = sqrt(mean((y_pre[1:(length(y_pre)*(CV_share))] - y_hat_train)^2))
  #   
  #   # test performance
  #   
  #   param_grid_2nd[grid,4] = sqrt(mean((y_pre[((length(y_pre) * (CV_share)) + 1):length(y_pre)] - y_hat_test)^2))
  #   
  # }
  # 
  # best_params = param_grid_2nd[which.min(param_grid_2nd$RMSFE),] 
  # 
  # w_regols = REGOLS_LASSO(
  #   cbind(1, x_pre),
  #   y_pre,
  #   lambda1 = best_params[["lambda1"]],
  #   lambda2 = best_params[["lambda2"]])
  # 
  # y_regols_pre = as.matrix(cbind(1, x_pre)) %*% w_regols
  # y_regols_post = as.matrix(cbind(1, x_post)) %*% w_regols
  # 
  # y_treat_regols = as.data.frame(c(y_pre, y_post)) %>%
  #   rename(y = c(1))
  # y_treat_regols$y_hat = c(y_regols_pre, y_regols_post)
  # 
  # # matplot(ts(y_treat_regols),
  # #         type = "l",
  # #         lty = 1,
  # #         lwd = 2,
  # #         main = "Regularized OLS Path (LASSO)",
  # #         xlab = "Time",
  # #         ylab = "Value")
  # 
  # results_REGOLS_LASSO = c()
  # 
  # results_REGOLS_LASSO["PRE_NET_RMSPE"] = sqrt(mean((y_pre - y_regols_pre)^2)) 
  # results_REGOLS_LASSO["PRE_NET_BIAS"] = mean(y_regols_pre - y_pre)
  # results_REGOLS_LASSO["PRE_NET_VAR"] = mean((y_regols_pre - mean(y_regols_pre))^2)
  # 
  # results_REGOLS_LASSO["POST_NET_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_regols_post)^2))
  # results_REGOLS_LASSO["POST_NET_BIAS"] = mean(y_regols_post - (y_post-post_effect))
  # results_REGOLS_LASSO["POST_NET_VAR"] = mean((y_regols_post - mean(y_regols_post))^2)
  # 
  # results[["REGOLS_LASSO"]] = results_REGOLS_LASSO
  
  
  # GLMNET estimation
  
  my_alpha = seq(0,1, by = 0.1)
  cv_fit = data.frame(matrix(NA, nrow = length(my_alpha), ncol = 3)) %>%
    rename(
      lambda = c(1),
      alpha = c(2),
      CVM = c(3))
  i = 1
  for (a in my_alpha) {
    
    cvfit = cv.glmnet(x_pre, y_pre, alpha = a, type.measure = "mse")
    cv_fit$lambda[i] = cvfit$lambda.min
    cv_fit$alpha[i] = a
    cv_fit$CVM[i] = min(cvfit$cvm)
    i = i+1
  }
  
  best_params = cv_fit[cv_fit$CVM == min(cv_fit$CVM),]
  cvfit = cv.glmnet(x_pre, y_pre, alpha = best_params$alpha, type.measure = "mse")
  
  y_net_pre = predict(cvfit, newx = x_pre, s = best_params$lambda)
  y_net_post = predict(cvfit, newx = x_post, s = best_params$lambda)
  
  y_treat_net = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_net$y_hat = c(y_net_pre, y_net_post)
# 
#   matplot(ts(y_treat_net),
#           type = "l",
#           lty = 1,
#           lwd = 2,
#           main = "Elastic Net Path",
#           xlab = "Time",
#           ylab = "Value")

  results_NET = c()
  
  results_NET["PRE_NET_RMSPE"] = sqrt(mean((y_pre - y_net_pre)^2)) 
  results_NET["PRE_NET_BIAS"] = mean(y_net_pre - y_pre)
  results_NET["PRE_NET_VAR"] = mean((y_net_pre - mean(y_net_pre))^2)
  
  results_NET["POST_NET_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_net_post)^2))
  results_NET["POST_NET_BIAS"] = mean(y_net_post - (y_post-post_effect))
  results_NET["POST_NET_VAR"] = mean((y_net_post - mean(y_net_post))^2)
  
  results[["NET"]] = results_NET
  
  # results[8] = sqrt(mean((y_pre - y_net_pre)^2))
  # results[9] = sqrt(mean(((y_post-post_effect) - y_net_post)^2))
  
  # FACTOR
  
  V    = cov(x_pre)
  eig  = eigen(V)
  w    = eig$vectors[,1:K]
  fhat = x_pre %*% w
  
  V_post    = cov(x_post)
  eig_post  = eigen(V_post)
  w_post    = eig_post$vectors[,1:K]
  fhat_post = x_post %*% w
  
  model = lm(y_pre ~ fhat)
  yest = model$fitted.values
  
  y_factor_pre = yest
  # Verify both approaches yield the same: plot(as.matrix(cbind(1, x_pre)) %*% w_factor, yest)
 
  w_factor = as.numeric(lm(yest ~ x_pre)$coefficients)
    
  # if (ncol(x_pre) < nrow(x_pre)){} else {}
  # y_factor_post = as.matrix(cbind(1, x_post)) %*% w_factor
  
  y_factor_post = as.matrix(cbind(1, fhat_post)) %*% model$coefficients
      
  y_treat_factor = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_factor$y_hat = c(y_factor_pre, y_factor_post)
  
  # matplot(ts(y_treat_factor),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "Factor Path",
  #         xlab = "Time",
  #         ylab = "Value")

  results_FACTOR = c()
  
  results_FACTOR["PRE_FACTOR_RMSPE"] = sqrt(mean((y_pre - y_factor_pre)^2)) 
  results_FACTOR["PRE_FACTOR_BIAS"] = mean(y_factor_pre - y_pre)
  results_FACTOR["PRE_FACTOR_VAR"] = mean((y_factor_pre - mean(y_factor_pre))^2)
  
  results_FACTOR["POST_FACTOR_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_factor_post)^2))
  results_FACTOR["POST_FACTOR_BIAS"] = mean(y_factor_post - (y_post-post_effect))
  results_FACTOR["POST_FACTOR_VAR"] = mean((y_factor_post - mean(y_factor_post))^2)
  
  results[["FACTOR"]] = results_FACTOR

  
  # results[10] = sqrt(mean((y_pre - y_factor_pre)^2))
  # results[11] = sqrt(mean(((y_post-post_effect) - y_factor_post)^2))
  
  #### bis hier ----
  
  # UNIDYN1: Jörg Fall korrigiert

  w0 = matrix(1 / J, J, 1)
  # w0 = matrix(c(rep(1/((J-1)/2), (J-1)/2), rep(0, 1 + J/2)), 
  #             J, 1)
  S0 = 100000000000000000000
  S1 = 90000000000000000000
  p_uni = 3
  
  y0 = y_pre[(p_uni + 1):T0]
  
  iter = 0
  while ((S0 - S1) > 0.00001) {
    # Mean of donor series 
    z = x_pre %*% w0
    
    # Lagged to allow dynamic modelling
    lagz = z[(p_uni + 1):T0]
    
    # Matrix of lagged z's. First col (t), second col (t-1) and so on. 
    for (i in (1:p_uni)) {
      lagz = cbind(lagz, z[(p_uni + 1 - i):(T0 - i)])
    }
    colnames(lagz) = paste0("lag", 0:p_uni)
    
    # Regression of y_0 on z-matrix 
    outreg = lm(y0 ~ lagz)
    alphas = outreg$coefficients[2:(p_uni + 2)]
    
    # Coefficient of lag0 multiplied with current donor values
    x1 = alphas[1] * x_pre[(p_uni + 1):T0, ]
    
    for (i in (1:(p_uni-1))) {
      x1 = x1 + alphas[i+1] * x_pre[(p_uni+1-i):(T0-i),]
    }
    
    # for (i in (1:p_uni - 1)) {
    #   x1 = x1 + alphas[i + 1] * x_pre[(p_uni + 1 - i):(T0 - i), ]
    # }

    outreg = lm(y0 ~ x1)
    
    # De-Meaning
    y0m = y0 - mean(y0)
    x1m = sweep(x1, 2, colMeans(x1), `-`)
    
    # Regularization w/o CV
    lam1 = best_params_REGOLS[["l1"]]
    lam2 = best_params_REGOLS[["l2"]]
    
    w0 = regularized_ols2(as.data.frame(x1), 
                     as.data.frame(y0), 
                     lam1, lam2)[2:(J+1)]
    
    y_hat = x1 %*% w0
    
    A = t(x1m) %*% x1m + lam1 * diag(J) + lam2 * matrix(1, J, J)   #  hier hatte ich nicht mit den mittelwertbereinigten Größen gearbeitet.
    w0 = solve(A) %*% (t(x1m) %*% y0m + lam2 * matrix(1, J, 1))
    xb = x1m %*% w0
    y0hat = mean(y0) + xb
    # 
    # plot(y_hat, y0hat)
    # y_hat == y0hat
    
    S0 = S1
    S1 = sd(y0 - y_hat)
    iter = iter + 1
    if (iter == 20) {
      S0 = S1
    }
  }
  
  plot(y_hat, y0hat)
  
  y_unidyn_pre = y_hat + mean(y_pre)
  
  # building x1_post
  
  x1_post = alphas[1] * x_post

  for (i in (1:(p_uni-1))) {
    x1_post = x1_post + alphas[i + 1] * rbind(x_pre[(T0 + 1 - i):T0, ],
                                              x_post[1:(T1 - i), ])
  }
  
  # for (i in (1:p_uni - 1)) {
  # 
  #   if (i == 0) {
  #     x1_post_lag = x_post
  #   } else {
  #     x1_post_lag = rbind(x_pre[(T0 +1 - i):T0, ],
  #                              x_post[1:(T1 - i), ])
  #   }
  # 
  #   x1_post = x1_post + alphas[i + 1] * x1_post_lag
  # }
  

  y_unidyn_post = x1_post %*% w0 + mean(y_pre)
  
  # De-Meaning
  # x1m_post = sweep(x1_post, 2, colMeans(x1_post), `-`)
  # 
  # y_unidyn_post = mean(y0) + x1m_post %*% w0
  
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

  # kurz ggplot um besser exportieren zu können

  df_gg = y_treat_unidyn %>%
    gather() %>%
    mutate(id = c(1:(T0+T1), 1:(T0+T1)))

  plot = ggplot(df_gg) +
    aes(x = id, y = value, colour = key) +
    geom_line(linewidth = 1.0) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    labs(title = "Dynamic (Univar) Path",
         subtitle = paste0("rho_u = ", round(rho_u,4), ", Donors = ", J,", rho_factor = ", rho))+
    theme(
      plot.title = element_text(size = 15L,
                                hjust = 0.5),
      plot.subtitle = element_text(size = 13L,
                                   hjust = 0.5))


  results[["Plots"]] = plot

  #ts.plot(x_pre)
  #ts.plot(x_post)

  results_UNIDYN1 = c()

  results_UNIDYN1["PRE_UNIDYN_RMSPE"] = sqrt(mean((tail(y_pre,T0-p_uni) - y_unidyn_pre)^2))
  results_UNIDYN1["PRE_UNIDYN_BIAS"] = mean(y_unidyn_pre - tail(y_pre,T0-p_uni))
  results_UNIDYN1["PRE_UNIDYN_VAR"] = mean((y_unidyn_pre - mean(y_unidyn_pre))^2)

  results_UNIDYN1["POST_UNIDYN_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_unidyn_post)^2))
  results_UNIDYN1["POST_UNIDYN_BIAS"] = mean(y_unidyn_post - (y_post-post_effect))
  results_UNIDYN1["POST_UNIDYN_VAR"] = mean((y_unidyn_post - mean(y_unidyn_post))^2)

  #sqrt(mean(((y_post-post_effect) - y_unidyn_post)^2))

  results[["UNIDYN1"]] = results_UNIDYN1
  
  # UNIDYN2: originär vor Korrektur
  
  w0 = matrix(1 / J, J, 1)
  # w0 = matrix(c(rep(1/((J-1)/2), (J-1)/2), rep(0, 1 + J/2)), 
  #             J, 1)
  S0 = 100000000000000000000
  S1 = 90000000000000000000
  p_uni = 3
  
  y0 = y_pre[(p_uni + 1):T0]
  
  iter = 0
  while ((S0 - S1) > 0.00001) {
    # Mean of donor series 
    z = x_pre %*% w0
    
    # Lagged to allow dynamic modelling
    lagz = z[(p_uni + 1):T0]
    
    # Matrix of lagged z's. First col (t), second col (t-1) and so on. 
    for (i in (1:p_uni)) {
      lagz = cbind(lagz, z[(p_uni + 1 - i):(T0 - i)])
    }
    colnames(lagz) = paste0("lag", 0:p_uni)
    
    # Regression of y_0 on z-matrix 
    outreg = lm(y0 ~ lagz)
    alphas = outreg$coefficients[2:(p_uni + 2)]
    
    # Coefficient of lag0 multiplied with current donor values
    x1 = alphas[1] * x_pre[(p_uni + 1):T0, ]
    
    for (i in (1:(p_uni))) {
      x1 = x1 + alphas[i+1] * x_pre[(p_uni+1-i):(T0-i),]
    }
    
    # for (i in (1:p_uni - 1)) {
    #   x1 = x1 + alphas[i + 1] * x_pre[(p_uni + 1 - i):(T0 - i), ]
    # }
    
    outreg = lm(y0 ~ x1)
    
    # Regularization w/o CV
    lam1 = best_params_REGOLS[["l1"]]
    lam2 = best_params_REGOLS[["l2"]]
    
    # De-Meaning
    y0m = y0 - mean(y0)
    x1m = sweep(x1, 2, colMeans(x1), `-`)

    A = t(x1m) %*% x1m + lam1 * diag(J) + lam2 * matrix(1, J, J)   #  hier hatte ich nicht mit den mittelwertbereinigten Größen gearbeitet.
    w0 = solve(A) %*% (t(x1m) %*% y0m + lam2 * matrix(1, J, 1))
    xb = x1m %*% w0
    y_hat = mean(y0) + xb
    
    # w0 = regularized_ols2(as.data.frame(x1), 
    #                       as.data.frame(y0), 
    #                       lam1, lam2)
    # 
    # y_hat = as.matrix(cbind(1, x1)) %*% w0
    
    S0 = S1
    S1 = sd(y0 - y_hat)
    iter = iter + 1
    if (iter == 20) {
      S0 = S1
    }
  }
  
  y_unidyn_pre = y_hat
  
  # building x1_post
  
  x1_post = alphas[1] * x_post
  
  for (i in (1:(p_uni))) {
    x1_post = x1_post + alphas[i + 1] * rbind(x_pre[(T0 + 1 - i):T0, ],
                                              x_post[1:(T1 - i), ])
  }

  # y_unidyn_post = as.matrix(cbind(1, x1_post)) %*% w0 
  
  
  # De-Meaning
  x1m_post = sweep(x1_post, 2, colMeans(x1_post), `-`)

  y_unidyn_post = mean(y0) + x1m_post %*% w0
  
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
  
  results_UNIDYN2 = c()
  
  results_UNIDYN2["PRE_UNIDYN_RMSPE"] = sqrt(mean((tail(y_pre,T0-p_uni) - y_unidyn_pre)^2))
  results_UNIDYN2["PRE_UNIDYN_BIAS"] = mean(y_unidyn_pre - tail(y_pre,T0-p_uni))
  results_UNIDYN2["PRE_UNIDYN_VAR"] = mean((y_unidyn_pre - mean(y_unidyn_pre))^2)
  
  results_UNIDYN2["POST_UNIDYN_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_unidyn_post)^2))
  results_UNIDYN2["POST_UNIDYN_BIAS"] = mean(y_unidyn_post - (y_post-post_effect))
  results_UNIDYN2["POST_UNIDYN_VAR"] = mean((y_unidyn_post - mean(y_unidyn_post))^2)
  
  #sqrt(mean(((y_post-post_effect) - y_unidyn_post)^2))
  
  results[["UNIDYN2"]] = results_UNIDYN2
  
  # UNIDYN3
  p_uni = 3
  
  lam1 = best_params_REGOLS[["l1"]]
  lam2 = best_params_REGOLS[["l2"]]
  
  T = T0 + T1
  
  w0 = matrix(1 / J, J, 1)
  y0 = y_pre[(p_uni + 1):T0]
  
  z = matrix(NA, T0 - p_uni, J)
  zhat = matrix(NA, T - p_uni, J)
  
  for (k in (1:J)) {
    lagx = y[(p_uni + 1):T, k + 1]
    for (i in (1:(p_uni - 1))) {
      lagx = cbind(lagx, y[(p_uni + 1 - i):(T - i), k + 1])
    }
    outz = lm(y0 ~ lagx[1:(T0 - p_uni), ])
    z[, k] = outz$fitted.values
    ahat = as.matrix(outz$coefficients)
    zhat[, k] = ahat[1] + lagx %*% ahat[2:(p_uni + 1), 1]
  }
  
  w0 = regularized_ols2(as.data.frame(z), 
                        as.data.frame(y0), 
                        lam1, lam2)[2:(J+1)]
  
  y_hat = z %*% w0
  
  y0m = y0 - mean(y0)
  zhat = sweep(zhat, 2, colMeans(zhat), `-`)
  z = sweep(z, 2, colMeans(z), `-`)

  mat1 = matrix(1, J, J)
  A = t(z) %*% z + lam1 * diag(J) + lam2 * mat1
  w0 = solve(A) %*% (t(z) %*% y0m + lam2 * matrix(1, J, 1))
  xb = z %*% w0
  y_unidyn2_pre = mean(y0) + xb
  # 
  # #  computing ex-post estimates of Y0
  yhat = mean(y0) + zhat %*% w0
  y_unidyn2_post = yhat[(T0 - p_uni + 1):(T - p_uni)]

  y_treat_unidyn2 = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))

  y_treat_unidyn2$y_hat = c(rep(NA, p_uni), y_unidyn2_pre, y_unidyn2_post)
  
  matplot(ts(y_treat_unidyn2),
          type = "l",
          lty = 1,
          lwd = 2,
          main = paste0("UNIDYN3. \n","rho_u = ", round(rho_u,4)),
          xlab = "Time",
          ylab = "Value")
  
  results_UNIDYN3 = c()
  
  results_UNIDYN3["PRE_UNIDYN_RMSPE"] = sqrt(mean((tail(y_pre,T0-p_uni) - y_unidyn2_pre)^2))
  results_UNIDYN3["PRE_UNIDYN_BIAS"] = mean(y_unidyn2_pre - tail(y_pre,T0-p_uni))
  results_UNIDYN3["PRE_UNIDYN_VAR"] = mean((y_unidyn2_pre - mean(y_unidyn2_pre))^2)
  
  results_UNIDYN3["POST_UNIDYN_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_unidyn2_post)^2))
  results_UNIDYN3["POST_UNIDYN_BIAS"] = mean(y_unidyn2_post - (y_post-post_effect))
  results_UNIDYN3["POST_UNIDYN_VAR"] = mean((y_unidyn2_post - mean(y_unidyn2_post))^2)
  
  results[["UNIDYN3"]] = results_UNIDYN3
  
  return(results)
}

simulation_VAR <- function(J) {
  
  T0 = obs/2
  T1 = T0
  mu = c(rep(1, times = J))
  A = matrix(runif(J * J), nrow = J, ncol = J)
  A = t(A) %*% A
  diag(A) = diag(A)*J^(0.6)*my_var
  
  # Randomly shrink covariances with donor unit to prevent perfect fit
  my_sample = sample(c(2:nrow(A)), round(J/2,0))
  scaling = 0
  
  for (i in my_sample) {
    A[i, 1] = A[i, 1]*scaling
    A[1, i] = A[1, i]*scaling
  }
  
  # Generating, de-mean data and splitting it
  df = MASS::mvrnorm(n = obs, mu = mu, Sigma = A, tol = T) 
  
  df_pre = df[1:(obs/2),]
  df_post = df[((obs/2)+1):obs,]
  x_pre = df_pre[,-1]
  x_post = df_post[,-1]
  y_pre = df_pre[,1]
  y_post = df_post[,1] + c
  
  df = cbind(c(y_pre, y_post), rbind(x_pre, x_post)) 
  
  matplot(ts(df),
          type = "l",
          lty = 1,
          lwd = 2,
          main = "Overview: Simulated Data",
          xlab = "Time",
          ylab = "Value")
  
  # OLS estimation
  # summary(lm(y_pre ~ x_pre))
  w_ols = summary(lm(y_pre ~ x_pre))$coefficients[,1]
  y_ols_pre = cbind(rep(1, (obs/2)), x_pre) %*% w_ols
  y_ols_post = cbind(rep(1, (obs/2)), x_post) %*% w_ols

  results = c()
  results[1] = sqrt(mean((y_pre - y_ols_pre)^2))
  results[2] = sqrt(mean(((y_post-c) - y_ols_post)^2))
  
  # y_treat_OLS = as.data.frame(c(y_pre, y_post)) %>%
  #   rename(y = c(1))
  # y_treat_OLS$y_hat = c(y_ols_pre, y_ols_post)
  # 
  # matplot(ts(y_treat_OLS),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "OLS Path",
  #         xlab = "Time",
  #         ylab = "Value")
  
  # SC estimation
  Dmat = t(x_pre) %*% x_pre
  dvec = t(x_pre) %*% y_pre
  Amat = t(rbind(rep(1, ncol(x_pre)), diag(ncol(x_pre)), -1*diag(ncol(x_pre))))
  bvec = c(1, rep(0, ncol(x_pre)), rep(-1,ncol(x_pre)))
  synth_model = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  w_sc = synth_model$solution
  y_sc_pre = x_pre %*% w_sc
  y_sc_post = x_post %*% w_sc
  
  results[3] = sqrt(mean((y_pre - y_sc_pre)^2))
  results[4] = sqrt(mean(((y_post-c) - y_sc_post)^2))
  
  # y_treat_SC = as.data.frame(c(y_pre, y_post)) %>%
  #   rename(y = c(1))
  # y_treat_SC$y_hat = c(y_sc_pre, y_sc_post)
  # 
  # matplot(ts(y_treat_SC),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "SC Path",
  #         xlab = "Time",
  #         ylab = "Value")
  
  # Regularized OLS. Simple Cross-Validation with 2-folds
  
  # x, y as data frame
  # l1, l2 as numerics
  # intercept as boolean
  
  # hyperparameter grid search
  
  param_grid <- expand.grid(
    l1 = seq(0,max((obs)/4, J*4), by = 5), # unsure about specific numbers. I see that more reg. is needed, the larger T1 and J
    l2 = seq(0,max((obs)/4, J*4), by = 5), # commented out for sake of simplicity
    intercept = c(TRUE, FALSE))
  
  param_grid$coeff_sum = NA
  param_grid$RMSPE = NA
  #param_grid$MAPE = NA
  param_grid$RMSFE = NA
  #param_grid$MAFE = NA
  
  for (grid in 1:nrow(param_grid)) {
    
    current_params = param_grid[grid, ]
    
    model = regularized_ols(
      x = x_pre[1:(nrow(x_pre)*(CV_share)),] %>% 
        scale(center = FALSE, scale = FALSE) %>% 
        as.data.frame(),
      y = y_pre[1:(length(y_pre)*(CV_share))] %>% 
        scale(center = FALSE, scale = FALSE) %>% 
        as.data.frame(),
      l1 = current_params[["l1"]],
      l2 = current_params[["l2"]], # same as above. If commented in, have to compute n^2 instead of n loops
      #l2 = current_params[["l1"]],
      intercept = current_params[["intercept"]])
    
    param_grid$coeff_sum[grid] = sum(model$beta[rownames(model$beta) != "x0"])
    
    # test performance
    
    param_grid$RMSPE[grid] = model$precision["RMSPE"]
    
    # train performance
    
    y_test = y_pre[((length(y_pre)*(CV_share)) + 1):length(y_pre)] %>% 
      as.matrix()
    x_test = x_pre[((nrow(x_pre)*(CV_share))+1):nrow(x_pre),] %>%  
      scale(center = FALSE, scale = FALSE) %>% 
      as.matrix()
    
    if ("x0" %in% rownames(model$beta)) {
      y_hat = as.matrix(cbind(1, x_test)) %*% model$beta
    } else {
      y_hat = x_test %*% model$beta
    }
    
    param_grid$RMSFE[grid] = sqrt(mean((y_test - y_hat)^2))
    # param_grid$MAFE[i] =  mean(abs(y_test - y_hat))
    
    # svMisc::progress(grid, nrow(param_grid))
    
  } 
  
  best_params = param_grid[which.min(param_grid$RMSFE),]
  
  # now extract cv-parameter combination
  
  w_regols = regularized_ols(
    x = x_pre %>% 
      scale(center = FALSE, scale = FALSE) %>% 
      as.data.frame(),
    y = y_pre %>% 
      scale(center = FALSE, scale = FALSE) %>% 
      as.data.frame(),
    l1 = best_params[["l1"]],
    # l2 = current_params[["l2"]] # same as above. If commented in, have to compute n^2 instead of n loops
    l2 = best_params[["l1"]],
    intercept = best_params[["intercept"]])$beta
  
  if ("x0" %in% rownames(w_regsynth)) {
    y_regols_pre = as.matrix(cbind(rep(1, (T0)), x_pre)) %*% w_regols
    y_regols_post = as.matrix(cbind(rep(1, (T0)), x_post)) %*% w_regols 
  } else {
    y_hat = x_test %*% model$beta
    y_regols_pre = as.matrix(x_pre) %*% w_regols
    y_regols_post = as.matrix(x_post) %*% w_regols
  }
  
  y_treat_regols = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_regols$y_hat = c(y_regols_pre, y_regols_post)
  
  # matplot(ts(y_treat_regsynth),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "Regularized OLS Path",
  #         xlab = "Time",
  #         ylab = "Value")
  
  results[5] = sqrt(mean((y_pre - y_regsynth_pre)^2))
  results[6] = sqrt(mean(((y_post-c) - y_regsynth_post)^2))
  
  return(results)
}


