for(pack in c("nnls","quadprog", "glmnet"))
  if(!require(pack, character.only = T))
  {
    install.packages(pack)
    require(pack, character.only = T)
  }
rm(pack)

# Function to simulate AR(1)
simulate_ar1 = function(rho, var_shock, T0, intercept = 0, start = NA) {
  y = ifelse(is.na(start), rnorm(1,mean=intercept/(1-rho),sd = sqrt(var_shock/(1-rho^2))), start)
  
  for(t in 1:T0)
    y = c(y,intercept + rho*y[t]+rnorm(1,mean=0,sd = sqrt(var_shock)))
  
  return(y[-1])
}

# SC estimator
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

# REGOLS estimator
regularized_ols = function(x, y, l1, l2) {

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

regularized_ols_multi = function(x, y, l1, l2){
  
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
  
  # Intercept (not regularized)
  I[1, 1] = 0
  U[1, 1] = 0
  
  # Only shrink contemporaneous values to 1
  U[c((J+2):length(U))] = 0
  
  beta = solve(t(X) %*% X + l1 * I + l2 * (U %*% t(U))) %*% ((t(X) %*% Y) + l2 * U)
  
  return(beta)
}

# Functions for numerical optimization
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

# all estimators
simulation_factor = function(J, simu_type = 'Factor'){
  
  # ab hier ----
  
  results = list()
  
  if (simu_type == 'Factor'){
    
    # print TS-structure in results-frame
    if (exists("rho_u_left")){
      rho_u = runif(1, rho_u_left, rho_u_right)
    }
    
    results[["rho_factor"]] = rho_factor
    results[["rho_error"]] = rho_u
    
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
  
  Factors = sapply(1:K, function(k){simulate_ar1(rho = rho_factor, var_shock = var_factor, T0 = T0+T1, intercept = alpha)}) 
  
  # AR(1) shocks
  transitory_shocks = sapply(1:dim(Mu)[1], function(k){simulate_ar1(rho = rho_u, var_shock = var_u, T0 = T0+T1, intercept = 0)}) 
  
  y = Factors %*% t(Mu) + transitory_shocks
  y = y + (1:nrow(y))^c
  
  } else {
    rho_u = rho_factor = NA
    
    est_coefs = VAR_est(J=J, J_max=J_max, p = p)
    stat_test = Stat_test(est_coefs)
    
    while (stat_test > 0.99) {
      est_coefs = VAR_est(J=J, p = p)
      stat_test = Stat_test(est_coefs)
     
      
    }
     y = tail(VAR_simu(est_coefs),(T0+T1))
      
      # Subset the dataframe using the randomly selected columns
      y <- y[, 1:(J+1)]
  }
  
  
  
  if (simu_type == 'Factor'){
  my_means = c(rnorm(1, mean = treat_inter, sd = 1), rnorm(J, mean = 0, sd = 1))

  
  for (i in 1:J) {
    y[,i] = y[,i] + my_means[i]
  }
  
  results[["bound_check"]] = as.numeric(rank(colMeans(y))[1])
  } else {
    results[["bound_check"]] = as.numeric(rank(colMeans(y))[1])
  }
  
  # Adding a treatment effect
  y[(T0+1):(T0+T1),1] = post_effect + y[(T0+1):(T0+T1),1] 
  
  # getting rid of year-label in VAR-data
  if (simu_type == "VAR"){
    
    y_matrix = as.numeric(y[,1])
    
    for (i in 1:J) {
      y_matrix = cbind(y_matrix, as.numeric(y[,i+1]))
    }
    
    y = y_matrix
    rm(y_matrix)
    
  }
  
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
  
  y_sc_forecast = y_treat_sc %>% 
    slice(-c(1:T0)) %>% 
    mutate(y = y - post_effect)
  
  # ggplot(y_sc_forecast) +
  #   aes(x = y, y = y_hat) +
  #   geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  #   geom_smooth(span = 0.75, method = "lm") +
  #   theme_minimal()
  
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
  
  #results_SC["MZ_REG"] = list(summary(lm(y ~ y_hat, data = y_sc_forecast))$coefficients[,1])
  results_SC["MZ_REG"] = list(summary(lm(as.vector(y_post-post_effect) ~ as.vector(y_sc_post)))$coefficients[,1])
  
  
  
  results[["SC"]] = results_SC
  
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
  
  y_ols_forecast = y_treat_ols %>% 
    slice(-c(1:T0)) %>% 
    mutate(y = y - post_effect)
  
  # ggplot(y_ols_forecast) +
  #   aes(x = y_hat, y = y) +
  #   geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  #   geom_smooth(span = 0.75, method = "lm") +
  #   theme_minimal()

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
  
  #results_OLS["MZ_REG"] = list(summary(lm(y ~ y_hat, data = y_ols_forecast))$coefficients[,1])
  results_OLS["MZ_REG"] = list(summary(lm(as.vector(y_post-post_effect) ~ as.vector(y_ols_post)))$coefficients[,1])
  
  results[["OLS"]] = results_OLS
  
  # REGOLS
  
  # 2-fold CV with random hyperparameter grid search
  # x, y as data frame
  # l1, l2 as numerics
  
  param_grid = expand.grid(
    l1 = 5^seq(1, 5, length.out = 50),
    l2 = 10^seq(1, 7, length.out = 50)) %>%
    as.data.frame() %>%
    sample_n(400) %>%
    bind_rows(c(l1 = 1,
                l2 = 1))
  
  param_grid$RMSPE = NA
  param_grid$RMSFE = NA
  
  for (grid in 1:nrow(param_grid)) {
    
    current_params = param_grid[grid, ]
    
    model = tryCatch({
      regularized_ols(
        x = x_pre[1:(nrow(x_pre)*(CV_share)),] %>% 
          as.data.frame(),
        y = y_pre[1:(length(y_pre)*(CV_share))] %>% 
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
      
      param_grid$RMSPE[grid] = sqrt(mean((y_train - y_regols_pre_CV1)^2))
      
      # test performance
      
      y_test = y_pre[(1+length(y_pre) * (CV_share)):length(y_pre)] %>%
        as.matrix()
      x_test = x_pre[(1+nrow(x_pre) * (CV_share)):nrow(x_pre), ] %>%
        as.matrix()
      
      y_regols_pre_CV2 = as.matrix(cbind(1, x_test)) %*% model
      
      param_grid$RMSFE[grid] = sqrt(mean((y_test - y_regols_pre_CV2) ^ 2))

  }
  
  # Second step
  
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
      regularized_ols(
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
  
  best_params_REGOLS = param_grid_2nd[which.min(param_grid_2nd$RMSFE),] 
   
  # Extract best CV-Parameter Combination
  
  w_regols = regularized_ols(
    x = x_pre %>%
      as.data.frame(),
    y = y_pre %>%
      as.data.frame(),
    l1 = best_params_REGOLS[["l1"]],
    l2 = best_params_REGOLS[["l2"]])
  
  y_regols_pre = as.matrix(cbind(1, x_pre)) %*% w_regols
  y_regols_post = as.matrix(cbind(1, x_post)) %*% w_regols 

  y_treat_regols = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_regols$y_hat = c(y_regols_pre, y_regols_post)
  
  y_regols_forecast = y_treat_regols %>% 
    slice(-c(1:T0)) %>% 
    mutate(y = y - post_effect)
  
  # ggplot(y_regols_forecast) +
  #   aes(x = y_hat, y = y) +
  #   geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  #   geom_smooth(span = 0.75, method = "lm") +
  #   theme_minimal()
  
  # matplot(ts(y_treat_regols),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "Regularized OLS Path",
  #         xlab = "Time",
  #         ylab = "Value")
  
  # Export to ggplot
  
  df_gg = y_treat_regols %>%
    gather() %>%
    mutate(id = c(1:(T0+T1), 1:(T0+T1)))
  
  plot = ggplot(df_gg) +
    aes(x = id, y = value, colour = key) +
    geom_line(linewidth = 1.0) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    labs(title = "REGOLS Path",
         subtitle = paste0("rho_u = ", round(rho_u,4), ", Donors = ", J,", rho_factor = ", rho_factor))+
    theme(
      plot.title = element_text(size = 15L,
                                hjust = 0.5),
      plot.subtitle = element_text(size = 13L,
                                   hjust = 0.5))
  
  
  results[["Plots_REGOLS"]] = plot
  
  results_REGOLS = c()
  
  results_REGOLS["PRE_REGOLS_RMSPE"] = sqrt(mean((y_pre - y_regols_pre)^2)) 
  results_REGOLS["PRE_REGOLS_BIAS"] = mean(y_regols_pre- y_pre)
  results_REGOLS["PRE_REGOLS_VAR"] = mean((y_regols_pre - mean(y_regols_pre))^2)
  
  results_REGOLS["POST_REGOLS_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_regols_post)^2))
  results_REGOLS["POST_REGOLS_BIAS"] = mean(y_regols_post - (y_post-post_effect))
  results_REGOLS["POST_REGOLS_VAR"] = mean((y_regols_post - mean(y_regols_post))^2)
  
  results_REGOLS["l1"] = c(best_params_REGOLS[1])
  results_REGOLS["l2"] = c(best_params_REGOLS[2])
  
  #results_REGOLS["MZ_REG"] = list(summary(lm(y ~ y_hat, data = y_regols_forecast))$coefficients[,1])
  results_REGOLS["MZ_REG"] = list(summary(lm(as.vector(y_post-post_effect) ~ as.vector(y_regols_post)))$coefficients[,1])
  
  
  results[["REGOLS"]] = results_REGOLS
  
  # GLMNET estimation
  
  my_alpha = seq(0,1, by = 0.1)
  cv_fit = data.frame(matrix(NA, nrow = length(my_alpha), ncol = 3)) %>%
    rename(
      lambda = c(1),
      alpha = c(2),
      CVM = c(3))
  i = 1
  for (a in my_alpha) {
    
    cvfit = cv.glmnet(x_pre, y_pre, alpha = a, type.measure = "mse", nfolds = 3)
    cv_fit$lambda[i] = cvfit$lambda.min
    cv_fit$alpha[i] = a
    cv_fit$CVM[i] = min(cvfit$cvm)
    i = i+1
  }
  
  best_params = cv_fit[cv_fit$CVM == min(cv_fit$CVM),]
  cvfit = cv.glmnet(x_pre, y_pre, alpha = best_params$alpha, type.measure = "mse", nfolds = 3)
  
  y_net_pre = predict(cvfit, newx = x_pre, s = best_params$lambda)
  y_net_post = predict(cvfit, newx = x_post, s = best_params$lambda)
  
  y_treat_net = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_net$y_hat = c(y_net_pre, y_net_post)
  
  y_net_forecast = y_treat_net %>% 
    slice(-c(1:T0)) %>% 
    mutate(y = y - post_effect)
  
  # ggplot(y_net_forecast) +
  #   aes(x = y_hat, y = y) +
  #   geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  #   geom_smooth(span = 0.75, method = "lm") +
  #   theme_minimal()
  

  # matplot(ts(y_treat_net),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "Elastic Net Path",
  #         xlab = "Time",
  #         ylab = "Value")

  results_NET = c()
  
  results_NET["PRE_NET_RMSPE"] = sqrt(mean((y_pre - y_net_pre)^2)) 
  results_NET["PRE_NET_BIAS"] = mean(y_net_pre - y_pre)
  results_NET["PRE_NET_VAR"] = mean((y_net_pre - mean(y_net_pre))^2)
  
  results_NET["POST_NET_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_net_post)^2))
  results_NET["POST_NET_BIAS"] = mean(y_net_post - (y_post-post_effect))
  results_NET["POST_NET_VAR"] = mean((y_net_post - mean(y_net_post))^2)
  
  #results_NET["MZ_REG"] = list(summary(lm(y ~ y_hat, data = y_net_forecast))$coefficients[,1])
  results_NET["MZ_REG"] = list(summary(lm(as.vector(y_post-post_effect) ~ as.vector(y_net_post)))$coefficients[,1])
  results[["NET"]] = results_NET
  
  # bis hier ----
  
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
  
  y_factor_forecast = y_treat_factor %>% 
    # slice(-c(1:T0)) %>% 
    # mutate(y = y - post_effect) 
    slice(c(1:T0))  
  
  
  ggplot(y_factor_forecast) +
    aes(x = y_hat, y = y) +
    geom_point(shape = "circle", size = 1.5, colour = "#112446") +
    geom_smooth(span = 0.75, method = "lm") +
    theme_minimal()
  
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
  
  #results_FACTOR["MZ_REG"] = list(summary(lm(y ~ y_hat, data = y_factor_forecast))$coefficients[,1])
  results_FACTOR["MZ_REG"] = list(summary(lm(as.vector(y_post-post_effect) ~ as.vector(y_factor_post)))$coefficients[,1])
  list(summary(lm(as.vector(y_post - post_effect) ~ as.vector(y_factor_post)))$coefficients[, 1])
  list(summary(lm(as.vector(y_pre) ~ as.vector(y_factor_pre)))$coefficients[, 1])
  
  results[["FACTOR"]] = results_FACTOR
  
  # if (dynamic == "yes"){
  
  # UNIDYN

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
  
  # 2-fold CV with random hyperparameter grid search
  # x, y as data frame
  # l1, l2 as numerics
  
  param_grid = expand.grid(
    l1 = 5^seq(1, 5, length.out = 50),
    l2 = 10^seq(1, 7, length.out = 50)) %>%
    as.data.frame() %>%
    sample_n(400) %>%
    bind_rows(c(l1 = 1,
                l2 = 1))
  
  param_grid$RMSPE = NA
  param_grid$RMSFE = NA
  
  for (grid in 1:nrow(param_grid)) {
    
    current_params = param_grid[grid, ]
    
    model = tryCatch({
      regularized_ols(
        x = z[1:(nrow(z)*(CV_share)),] %>% 
          as.data.frame(),
        y = y0[1:(length(y0)*(CV_share))] %>% 
          as.data.frame(),
        l1 = current_params[["l1"]],
        l2 = current_params[["l2"]])},
      
      error = function(e) {
        return(list(NA))
      })
    
    # train performance
    
    y_train = y0[1:(length(y0) * (CV_share))] %>%
      as.matrix()
    x_train = z[1:(nrow(z) * (CV_share)), ] %>%
      as.matrix()
    
    y_regols_pre_CV1 = as.matrix(cbind(1, x_train)) %*% model
    
    param_grid$RMSPE[grid] = sqrt(mean((y_train - y_regols_pre_CV1)^2))
    
    # test performance
    
    y_test = y0[(1+length(y0) * (CV_share)):length(y0)] %>%
      as.matrix()
    x_test = z[(1+nrow(z) * (CV_share)):nrow(z), ] %>%
      as.matrix()
    
    y_regols_pre_CV2 = as.matrix(cbind(1, x_test)) %*% model
    
    param_grid$RMSFE[grid] = sqrt(mean((y_test - y_regols_pre_CV2) ^ 2))
    
  }
  
  # Second step
  
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
      regularized_ols(
        x = z[1:(nrow(z)*(CV_share)),] %>% 
          as.data.frame(),
        y = y0[1:(length(y0)*(CV_share))] %>% 
          as.data.frame(),
        l1 = current_params[["l1"]],
        l2 = current_params[["l2"]])},
      
      error = function(e) {
        return(list(NA))
      })
    
    # train performance
    
    y_train = y0[1:(length(y0) * (CV_share))] %>%
      as.matrix()
    x_train = z[1:(nrow(z) * (CV_share)), ] %>%
      as.matrix()
    
    y_regols_pre_CV1 = as.matrix(cbind(1, x_train)) %*% model
    
    param_grid_2nd$RMSPE[grid] = sqrt(mean((y_train - y_regols_pre_CV1)^2))
    
    # test performance
    
    y_test = y0[(1+length(y0) * (CV_share)):length(y0)] %>%
      as.matrix()
    x_test = z[(1+nrow(z) * (CV_share)):nrow(z), ] %>%
      as.matrix()
    
    y_regols_pre_CV2 = as.matrix(cbind(1, x_test)) %*% model
    
    param_grid_2nd$RMSFE[grid] = sqrt(mean((y_test - y_regols_pre_CV2) ^ 2))
  }
  
  
  best_params = param_grid_2nd[which.min(param_grid_2nd$RMSFE),] 
  
  # Check to see if CV pays off. It does
  # lam1 = best_params_REGOLS[["l1"]]
  # lam2 = best_params_REGOLS[["l2"]]
  # lam1 = lam2 = 0
  
  lam1 = best_params[["l1"]]
  lam2 = best_params[["l2"]]
  
  
  w0 = regularized_ols(as.data.frame(z),
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
  
  # computing ex-post estimates of Y0
  yhat = mean(y0) + zhat %*% w0
  y_unidyn2_post = yhat[(T0 - p_uni + 1):(T - p_uni)]

  y_treat_unidyn2 = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))

  y_treat_unidyn2$y_hat = c(rep(NA, p_uni), y_unidyn2_pre, y_unidyn2_post)
  
  y_unidyn_forecast = y_treat_unidyn2 %>% 
    slice(-c(1:T0)) %>% 
    mutate(y = y - post_effect)
  
  # ggplot(y_unidyn_forecast) +
  #   aes(x = y_hat, y = y) +
  #   geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  #   geom_smooth(span = 0.75, method = "lm") +
  #   theme_minimal()
  
  # matplot(ts(y_treat_unidyn2),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = paste0("UNIDYN. \n","rho_u = ", round(rho_u,4)),
  #         xlab = "Time",
  #         ylab = "Value")
  
  
  # Export to ggplot
  
  df_gg = y_treat_unidyn2 %>%
    gather() %>%
    mutate(id = c(1:(T0+T1), 1:(T0+T1)))
  
  plot = ggplot(df_gg) +
    aes(x = id, y = value, colour = key) +
    geom_line(linewidth = 1.0) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    labs(title = "UNIDYN2 Path",
         subtitle = paste0("rho_u = ", round(rho_u,4), ", Donors = ", J,", rho_factor = ", rho_factor))+
    theme(
      plot.title = element_text(size = 15L,
                                hjust = 0.5),
      plot.subtitle = element_text(size = 13L,
                                   hjust = 0.5))
  
  results[["Plots_UNIDYN"]] = plot
  
  results_UNIDYN = c()
  
  results_UNIDYN["PRE_UNIDYN_RMSPE"] = sqrt(mean((tail(y_pre,T0-p_uni) - y_unidyn2_pre)^2))
  results_UNIDYN["PRE_UNIDYN_BIAS"] = mean(y_unidyn2_pre - tail(y_pre,T0-p_uni))
  results_UNIDYN["PRE_UNIDYN_VAR"] = mean((y_unidyn2_pre - mean(y_unidyn2_pre))^2)
  
  results_UNIDYN["POST_UNIDYN_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_unidyn2_post)^2))
  results_UNIDYN["POST_UNIDYN_BIAS"] = mean(y_unidyn2_post - (y_post-post_effect))
  results_UNIDYN["POST_UNIDYN_VAR"] = mean((y_unidyn2_post - mean(y_unidyn2_post))^2)
  
  #results_UNIDYN["MZ_REG"] = list(summary(lm(y ~ y_hat, data = y_unidyn_forecast))$coefficients[,1])
  results_UNIDYN["MZ_REG"] = list(summary(lm(as.vector(y_post-post_effect) ~ as.vector(y_unidyn2_post)))$coefficients[,1])
  
  results[["UNIDYN"]] = results_UNIDYN
  
  # MULTIDYN1
  
  lagx = x_pre[(p_multi + 1):T0,]
  
  for (i in (1:p_multi)) {
    lagx = cbind(lagx, x_pre[(p_multi + 1 - i):(T0 - i),])
  }
  
  lagy = y_pre[(p_multi):(T0-1)]
  
  for (i in (1:(p_multi-1))) {
    lagy = cbind(lagy, y_pre[(p_multi - i):(T0 - 1- i)])
  }
  
  xfull = cbind(lagx, lagy)
  yfull = y_pre[(p_multi+1):(T0)]
  
  colnames(xfull) = c(paste0(paste0("lagx", 1:ncol(x_pre)), "_", rep(0:(p_multi), each= ncol(x_pre))),
                      paste0("lagy_", 1:(p_multi)))
  
  # Regularized regression
  
  param_grid = expand.grid(
    l1 = 5^seq(1, 5, length.out = 50),
    l2 = 10^seq(1, 7, length.out = 50)) %>%
    as.data.frame() %>%
    sample_n(400) %>%
    bind_rows(c(l1 = 1,
                l2 = 1))
  
  param_grid$RMSPE = NA
  param_grid$RMSFE = NA
  
  for (grid in 1:nrow(param_grid)) {
    
    current_params = param_grid[grid, ]
    
    model = tryCatch({
      regularized_ols(
        x = xfull[1:(nrow(xfull)*(CV_share)),] %>% 
          as.data.frame(),
        y = yfull[1:(length(yfull)*(CV_share))] %>% 
          as.data.frame(),
        l1 = current_params[["l1"]],
        l2 = current_params[["l2"]])},
      
      error = function(e) {
        return(list(NA))
      })
    
    # train performance
    
    y_train = yfull[1:(round(length(yfull) * (CV_share))+.1)] %>%
      as.matrix()
    x_train = xfull[1:(round(nrow(xfull) * (CV_share))+.1), ] %>%
      as.matrix()
    
    y_multidyn_pre_CV1 = as.matrix(cbind(1, x_train)) %*% model
    
    param_grid$RMSPE[grid] = sqrt(mean((y_train - y_multidyn_pre_CV1)^2))
    
    # test performance
    
    y_test = yfull[(1 + (round((length(yfull) * (CV_share))+.1))):length(yfull)] %>%
      as.matrix()
    x_test = xfull[(1 + (round((nrow(xfull) * (CV_share))+.1))):nrow(xfull), ] %>%
      as.matrix()
    
    y_multidyn_pre_CV2 = as.matrix(cbind(1, x_test)) %*% model
    
    param_grid$RMSFE[grid] = sqrt(mean((y_test - y_multidyn_pre_CV2)^2))
    
  }
  
  # Second step
  
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
      regularized_ols(
        x = xfull[1:(nrow(xfull)*(CV_share)),] %>% 
          as.data.frame(),
        y = yfull[1:(length(yfull)*(CV_share))] %>% 
          as.data.frame(),
        l1 = current_params[["l1"]],
        l2 = current_params[["l2"]])},
      
      error = function(e) {
        return(list(NA))
      })
    
    # train performance
    
    y_train = yfull[1:(round(length(yfull) * (CV_share))+.1)] %>%
      as.matrix()
    x_train = xfull[1:(round(nrow(xfull) * (CV_share))+.1), ] %>%
      as.matrix()
    
    y_multidyn_pre_CV1 = as.matrix(cbind(1, x_train)) %*% model
    
    param_grid_2nd$RMSPE[grid] = sqrt(mean((y_train - y_multidyn_pre_CV1)^2))
    
    # test performance
    
    y_test = yfull[(1 + (round((length(yfull) * (CV_share))+.1))):length(yfull)] %>%
      as.matrix()
    x_test = xfull[(1 + (round((nrow(xfull) * (CV_share))+.1))):nrow(xfull), ] %>%
      as.matrix()
    
    y_multidyn_pre_CV2 = as.matrix(cbind(1, x_test)) %*% model
    
    param_grid_2nd$RMSFE[grid] = sqrt(mean((y_test - y_multidyn_pre_CV2)^2))
  }
  
  best_params = param_grid_2nd[which.min(param_grid_2nd$RMSFE),] 

  
  w_multidyn1 = regularized_ols(
    x = xfull[1:(nrow(xfull) * (CV_share)), ] %>%
      as.data.frame(),
    y = yfull[1:(length(yfull) * (CV_share))] %>%
      as.data.frame(),
    l1 = best_params[["l1"]],
    l2 = best_params[["l2"]])
  
  y_multidyn1_pre = as.matrix(cbind(1, xfull)) %*% w_multidyn1
  
  # bulding x_full_post
  
  x_prepost = rbind(x_pre[((T0+1)-p_multi):T0,],
                    x_post[1:T1,])
  
  lagx_post = x_prepost[(p_multi + 1):nrow(x_prepost),]
  
  for (i in (1:p_multi)) {
    lagx_post = cbind(lagx_post, x_prepost[(p_multi + 1 - i):(nrow(x_prepost) - i),])
  }
  
  y_prepost = c(y_pre[((T0+1)-p_multi):T0],
                rep(NA, T1))
  
  lagy_post = y_prepost[p_multi:(length(y_prepost)-1)]
  
  for (i in (1:(p_multi-1))) {
    lagy_post = cbind(lagy_post, y_prepost[(p_multi - i):(length(y_prepost) - 1- i)])
  }
  
  xfull_post = cbind(lagx_post, lagy_post)
  
  y_multidyn1_post = rep(NA, T1)
  
  for (i in 1:(T1-1)) {
    y_multidyn1_post[i] = as.matrix(cbind(1, xfull_post))[i, ] %*% w_multidyn1
    
    # updating y in xfull_post
    
    xfull_post[i + 1, ncol(lagx_post)+1] = y_multidyn1_post[i]
    
    if (i + 2 <= T1 & p_multi >= 2){
      xfull_post[i + 2, ncol(lagx_post)+2] = y_multidyn1_post[i]
    } 
    if (i + 3 <= T1 & p_multi >= 3){
      xfull_post[i + 3, ncol(lagx_post)+3] = y_multidyn1_post[i]
    }  
    if (i + 4 <= T1 & p_multi >= 4){
      xfull_post[i + 4, ncol(lagx_post)+4] = y_multidyn1_post[i]
    }  
  }
  
  # last period
  y_multidyn1_post[T1] = as.matrix(cbind(1, xfull_post))[T1, ] %*% w_multidyn1
  
  
  y_treat_multidyn1 = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_multidyn1$y_hat = c(rep(NA, p_multi),
                             y_multidyn1_pre, 
                             y_multidyn1_post)
  
  y_multidyn1_forecast = y_treat_multidyn1 %>% 
    slice(-c(1:T0)) %>% 
    mutate(y = y - post_effect)
  
  # ggplot(y_multidyn1_forecast) +
  #   aes(x = y_hat, y = y) +
  #   geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  #   geom_smooth(span = 0.75, method = "lm") +
  #   theme_minimal()
  
  # matplot(ts(y_treat_multidyn1),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "MULTIDYN1",
  #         xlab = "Time",
  #         ylab = "Value")
  
  # Export to ggplot
  
  df_gg = y_treat_multidyn1 %>%
    gather() %>%
    mutate(id = c(1:(T0+T1), 1:(T0+T1)))
  
  plot = ggplot(df_gg) +
    aes(x = id, y = value, colour = key) +
    geom_line(linewidth = 1.0) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    labs(title = "MULTIDYN1 Path",
         subtitle = paste0("rho_u = ", round(rho_u,4), ", Donors = ", J,", rho_factor = ", rho_factor))+
    theme(
      plot.title = element_text(size = 15L,
                                hjust = 0.5),
      plot.subtitle = element_text(size = 13L,
                                   hjust = 0.5))
  
  results[["Plots_MULTIDYN1"]] = plot
  
  results_MULTIDYN1 = c()
  
  results_MULTIDYN1["PRE_MULTIDYN1_RMSPE"] = sqrt(mean((y_pre[(p_multi+1):T0] -  y_multidyn1_pre)^2)) 
  results_MULTIDYN1["PRE_MULTIDYN1_BIAS"] = mean(y_multidyn1_pre- y_pre[(p_multi+1):T0])
  results_MULTIDYN1["PRE_MULTIDYN1_VAR"] = mean((y_multidyn1_pre - mean(y_multidyn1_pre))^2)
  
  results_MULTIDYN1["POST_MULTIDYN1_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_multidyn1_post)^2))
  results_MULTIDYN1["POST_MULTIDYN1_BIAS"] = mean(y_multidyn1_post - (y_post-post_effect))
  results_MULTIDYN1["POST_MULTIDYN1_VAR"] = mean((y_multidyn1_post - mean(y_multidyn1_post))^2)
  
  #results_MULTIDYN1["MZ_REG"] = list(summary(lm(y ~ y_hat, data = y_multidyn1_forecast))$coefficients[,1])
  results_MULTIDYN1["MZ_REG"] = list(summary(lm(as.vector(y_post-post_effect) ~ as.vector(y_multidyn1_post)))$coefficients[,1])
  results[["MULTIDYN1"]] = results_MULTIDYN1
  
  # MULTIDYN2
  
  lagx = x_pre[(p_multi + 1):T0,]
  
  for (i in (1:p_multi)) {
    lagx = cbind(lagx, x_pre[(p_multi + 1 - i):(T0 - i),])
  }
  
  lagy = y_pre[(p_multi):(T0-1)]
  
  for (i in (1:(p_multi-1))) {
    lagy = cbind(lagy, y_pre[(p_multi - i):(T0 - 1- i)])
  }
  
  xfull = cbind(lagx, lagy)
  yfull = y_pre[(p_multi+1):(T0)]
  
  colnames(xfull) = c(paste0(paste0("lagx", 1:ncol(x_pre)), "_", rep(0:(p_multi), each= ncol(x_pre))),
                      paste0("lagy_", 1:(p_multi)))
  
  # Regularized regression: only contemporaneous values to 1
  
  param_grid = expand.grid(
    l1 = 5^seq(1, 5, length.out = 50),
    l2 = 10^seq(1, 7, length.out = 50)) %>%
    as.data.frame() %>%
    sample_n(400) %>%
    bind_rows(c(l1 = 1,
                l2 = 1))
  
  param_grid$RMSPE = NA
  param_grid$RMSFE = NA
  
  for (grid in 1:nrow(param_grid)) {
    
    current_params = param_grid[grid, ]
    
    model = tryCatch({
      regularized_ols_multi(
        x = xfull[1:(nrow(xfull)*(CV_share)),] %>% 
          as.data.frame(),
        y = yfull[1:(length(yfull)*(CV_share))] %>% 
          as.data.frame(),
        l1 = current_params[["l1"]],
        l2 = current_params[["l2"]])},
      
      error = function(e) {
        return(list(NA))
      })
    
    # train performance
    
    y_train = yfull[1:(round(length(yfull) * (CV_share))+.1)] %>%
      as.matrix()
    x_train = xfull[1:(round(nrow(xfull) * (CV_share))+.1), ] %>%
      as.matrix()
    
    y_multidyn_pre_CV1 = as.matrix(cbind(1, x_train)) %*% model
    
    param_grid$RMSPE[grid] = sqrt(mean((y_train - y_multidyn_pre_CV1)^2))
    
    # test performance
    
    y_test = yfull[(1 + (round((length(yfull) * (CV_share))+.1))):length(yfull)] %>%
      as.matrix()
    x_test = xfull[(1 + (round((nrow(xfull) * (CV_share))+.1))):nrow(xfull), ] %>%
      as.matrix()
    
    y_multidyn_pre_CV2 = as.matrix(cbind(1, x_test)) %*% model
    
    param_grid$RMSFE[grid] = sqrt(mean((y_test - y_multidyn_pre_CV2)^2))
    
  }
  
  # Second step
  
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
      regularized_ols_multi(
        x = xfull[1:(nrow(xfull)*(CV_share)),] %>% 
          as.data.frame(),
        y = yfull[1:(length(yfull)*(CV_share))] %>% 
          as.data.frame(),
        l1 = current_params[["l1"]],
        l2 = current_params[["l2"]])},
      
      error = function(e) {
        return(list(NA))
      })
    
    # train performance
    
    y_train = yfull[1:(round(length(yfull) * (CV_share))+.1)] %>%
      as.matrix()
    x_train = xfull[1:(round(nrow(xfull) * (CV_share))+.1), ] %>%
      as.matrix()
    
    y_multidyn_pre_CV1 = as.matrix(cbind(1, x_train)) %*% model
    
    param_grid_2nd$RMSPE[grid] = sqrt(mean((y_train - y_multidyn_pre_CV1)^2))
    
    # test performance
    
    y_test = yfull[(1 + (round((length(yfull) * (CV_share))+.1))):length(yfull)] %>%
      as.matrix()
    x_test = xfull[(1 + (round((nrow(xfull) * (CV_share))+.1))):nrow(xfull), ] %>%
      as.matrix()
    
    y_multidyn_pre_CV2 = as.matrix(cbind(1, x_test)) %*% model
    
    param_grid_2nd$RMSFE[grid] = sqrt(mean((y_test - y_multidyn_pre_CV2)^2))
  }
  
  best_params = param_grid_2nd[which.min(param_grid_2nd$RMSFE),] 
  
  
  w_multidyn2 = regularized_ols_multi(
    x = xfull[1:(nrow(xfull) * (CV_share)), ] %>%
      as.data.frame(),
    y = yfull[1:(length(yfull) * (CV_share))] %>%
      as.data.frame(),
    l1 = best_params[["l1"]],
    l2 = best_params[["l2"]])
  
  y_multidyn2_pre = as.matrix(cbind(1, xfull)) %*% w_multidyn2
  
  # bulding x_full_post
  
  x_prepost = rbind(x_pre[((T0+1)-p_multi):T0,],
                    x_post[1:T1,])
  
  lagx_post = x_prepost[(p_multi + 1):nrow(x_prepost),]
  
  for (i in (1:p_multi)) {
    lagx_post = cbind(lagx_post, x_prepost[(p_multi + 1 - i):(nrow(x_prepost) - i),])
  }
  
  y_prepost = c(y_pre[((T0+1)-p_multi):T0],
                rep(NA, T1))
  
  lagy_post = y_prepost[p_multi:(length(y_prepost)-1)]
  
  for (i in (1:(p_multi-1))) {
    lagy_post = cbind(lagy_post, y_prepost[(p_multi - i):(length(y_prepost) - 1- i)])
  }
  
  xfull_post = cbind(lagx_post, lagy_post)
  
  y_multidyn2_post = rep(NA, T1)
  
  for (i in 1:(T1-1)) {
    y_multidyn2_post[i] = as.matrix(cbind(1, xfull_post))[i, ] %*% w_multidyn2
    
    # updating y in xfull_post
    
    xfull_post[i + 1, ncol(lagx_post)+1] = y_multidyn2_post[i]
    
    if (i + 2 <= T1 & p_multi >= 2){
      xfull_post[i + 2, ncol(lagx_post)+2] = y_multidyn2_post[i]
    } 
    if (i + 3 <= T1 & p_multi >= 3){
      xfull_post[i + 3, ncol(lagx_post)+3] = y_multidyn2_post[i]
    }  
    if (i + 4 <= T1 & p_multi >= 4){
      xfull_post[i + 4, ncol(lagx_post)+4] = y_multidyn2_post[i]
    }  
  }
  
  # last period
  y_multidyn2_post[T1] = as.matrix(cbind(1, xfull_post))[T1, ] %*% w_multidyn2
  
  
  y_treat_multidyn2 = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_multidyn2$y_hat = c(rep(NA, p_multi),
                              y_multidyn2_pre, 
                              y_multidyn2_post)
  
  y_multidyn2_forecast = y_treat_multidyn2 %>% 
    slice(-c(1:T0)) %>% 
    mutate(y = y - post_effect)
  
  # ggplot(y_multidyn2_forecast) +
  #   aes(x = y_hat, y = y) +
  #   geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  #   geom_smooth(span = 0.75, method = "lm") +
  #   theme_minimal()
  
  
  # matplot(ts(y_treat_multidyn2),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "MULTIDYN2",
  #         xlab = "Time",
  #         ylab = "Value")
  
  # Export to ggplot
  
  df_gg = y_treat_multidyn2 %>%
    gather() %>%
    mutate(id = c(1:(T0+T1), 1:(T0+T1)))
  
  plot = ggplot(df_gg) +
    aes(x = id, y = value, colour = key) +
    geom_line(linewidth = 1.0) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    labs(title = "MULTIDYN2 Path",
         subtitle = paste0("rho_u = ", round(rho_u,4), ", Donors = ", J,", rho_factor = ", rho_factor))+
    theme(
      plot.title = element_text(size = 15L,
                                hjust = 0.5),
      plot.subtitle = element_text(size = 13L,
                                   hjust = 0.5))
  
  results[["Plots_MULTIDYN2"]] = plot
  
  results_MULTIDYN2 = c()
  
  results_MULTIDYN2["PRE_MULTIDYN2_RMSPE"] = sqrt(mean((y_pre[(p_multi+1):T0] -  y_multidyn2_pre)^2)) 
  results_MULTIDYN2["PRE_MULTIDYN2_BIAS"] = mean(y_multidyn2_pre- y_pre[(p_multi+1):T0])
  results_MULTIDYN2["PRE_MULTIDYN2_VAR"] = mean((y_multidyn2_pre - mean(y_multidyn2_pre))^2)
  
  results_MULTIDYN2["POST_MULTIDYN2_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_multidyn2_post)^2))
  results_MULTIDYN2["POST_MULTIDYN2_BIAS"] = mean(y_multidyn2_post - (y_post-post_effect))
  results_MULTIDYN2["POST_MULTIDYN2_VAR"] = mean((y_multidyn2_post - mean(y_multidyn2_post))^2)
  
  #results_MULTIDYN2["MZ_REG"] = list(summary(lm(y ~ y_hat, data = y_multidyn2_forecast))$coefficients[,1])
  results_MULTIDYN2["MZ_REG"] = list(summary(lm(as.vector(y_post-post_effect) ~ as.vector(y_multidyn2_post)))$coefficients[,1])
  
  results[["MULTIDYN2"]] = results_MULTIDYN2
  
  # MULTIDYN3
  
  lagx = x_pre[(p_multi + 1):T0,]
  
  for (i in (1:p_multi)) {
    lagx = cbind(lagx, x_pre[(p_multi + 1 - i):(T0 - i),])
  }
  
  lagy = y_pre[(p_multi):(T0-1)]
  
  for (i in (1:(p_multi-1))) {
    lagy = cbind(lagy, y_pre[(p_multi - i):(T0 - 1- i)])
  }
  
  xfull = cbind(lagx, lagy)
  yfull = y_pre[(p_multi+1):(T0)]
  
  colnames(xfull) = c(paste0(paste0("lagx", 1:ncol(x_pre)), "_", rep(0:(p_multi), each= ncol(x_pre))),
                      paste0("lagy_", 1:(p_multi)))
  
  # Elastic net regression
  
  my_alpha = seq(0,1, by = 0.1)
  cv_fit = data.frame(matrix(NA, nrow = length(my_alpha), ncol = 3)) %>%
    rename(
      lambda = c(1),
      alpha = c(2),
      CVM = c(3))
  i = 1
  for (a in my_alpha) {
    
    cvfit = cv.glmnet(xfull, yfull, alpha = a, type.measure = "mse", nfolds = 3)
    cv_fit$lambda[i] = cvfit$lambda.min
    cv_fit$alpha[i] = a
    cv_fit$CVM[i] = min(cvfit$cvm)
    i = i+1
  }
  
  best_params = cv_fit[cv_fit$CVM == min(cv_fit$CVM),]
  cvfit = cv.glmnet(xfull, yfull, alpha = best_params$alpha, type.measure = "mse", nfolds = 3)

  y_multidyn3_pre = predict(cvfit, newx = xfull, s = best_params$lambda)
  
  # bulding x_full_post
  
  x_prepost = rbind(x_pre[((T0+1)-p_multi):T0,],
                    x_post[1:T1,])
  
  lagx_post = x_prepost[(p_multi + 1):nrow(x_prepost),]
  
  for (i in (1:p_multi)) {
    lagx_post = cbind(lagx_post, x_prepost[(p_multi + 1 - i):(nrow(x_prepost) - i),])
  }
  
  y_prepost = c(y_pre[((T0+1)-p_multi):T0],
                rep(NA, T1))
  
  lagy_post = y_prepost[p_multi:(length(y_prepost)-1)]
  
  for (i in (1:(p_multi-1))) {
    lagy_post = cbind(lagy_post, y_prepost[(p_multi - i):(length(y_prepost) - 1- i)])
  }
  
  xfull_post = cbind(lagx_post, lagy_post)
  
  y_multidyn3_post = rep(NA, T1)
  
  for (i in 1:(T1-1)) {
    y_multidyn3_post[i] = predict(cvfit, newx = xfull_post[i,], s = best_params$lambda)
    
    # updating y in xfull_post
    
    xfull_post[i + 1, ncol(lagx_post)+1] = y_multidyn3_post[i]
    
    if (i + 2 <= T1 & p_multi >= 2){
      xfull_post[i + 2, ncol(lagx_post)+2] = y_multidyn3_post[i]
    } 
    if (i + 3 <= T1 & p_multi >= 3){
      xfull_post[i + 3, ncol(lagx_post)+3] = y_multidyn3_post[i]
    }  
    if (i + 4 <= T1 & p_multi >= 4){
      xfull_post[i + 4, ncol(lagx_post)+4] = y_multidyn3_post[i]
    }  
  }
  
  # last period
  y_multidyn3_post[T1] = predict(cvfit, newx = xfull_post[T1,], s = best_params$lambda)
  
  y_treat_multidyn3 = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_multidyn3$y_hat = c(rep(NA, p_multi),
                              y_multidyn3_pre, 
                              y_multidyn3_post)
  
  y_multidyn3_forecast = y_treat_multidyn3 %>% 
    slice(-c(1:T0)) %>% 
    mutate(y = y - post_effect)
  
  # ggplot(y_multidyn3_forecast) +
  #   aes(x = y_hat, y = y) +
  #   geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  #   geom_smooth(span = 0.75, method = "lm") +
  #   theme_minimal()
  
  # matplot(ts(y_treat_multidyn3),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "MULTIDYN3",
  #         xlab = "Time",
  #         ylab = "Value")
  
  # Export to ggplot
  
  df_gg = y_treat_multidyn3 %>%
    gather() %>%
    mutate(id = c(1:(T0+T1), 1:(T0+T1)))
  
  plot = ggplot(df_gg) +
    aes(x = id, y = value, colour = key) +
    geom_line(linewidth = 1.0) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    labs(title = "MULTIDYN3 Path",
         subtitle = paste0("rho_u = ", round(rho_u,4), ", Donors = ", J,", rho_factor = ", rho_factor))+
    theme(
      plot.title = element_text(size = 15L,
                                hjust = 0.5),
      plot.subtitle = element_text(size = 13L,
                                   hjust = 0.5))
  
  results[["Plots_MULTIDYN3"]] = plot
  
  results_MULTIDYN3 = c()
  
  results_MULTIDYN3["PRE_MULTIDYN3_RMSPE"] = sqrt(mean((y_pre[(p_multi+1):T0] -  y_multidyn3_pre)^2)) 
  results_MULTIDYN3["PRE_MULTIDYN3_BIAS"] = mean(y_multidyn3_pre- y_pre[(p_multi+1):T0])
  results_MULTIDYN3["PRE_MULTIDYN3_VAR"] = mean((y_multidyn3_pre - mean(y_multidyn3_pre))^2)
  
  results_MULTIDYN3["POST_MULTIDYN3_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_multidyn3_post)^2))
  results_MULTIDYN3["POST_MULTIDYN3_BIAS"] = mean(y_multidyn3_post - (y_post-post_effect))
  results_MULTIDYN3["POST_MULTIDYN3_VAR"] = mean((y_multidyn3_post - mean(y_multidyn3_post))^2)
  
  #results_MULTIDYN3["MZ_REG"] = list(summary(lm(y ~ y_hat, data = y_multidyn3_forecast))$coefficients[,1])
  results_MULTIDYN3["MZ_REG"] = list(summary(lm(as.vector(y_post-post_effect) ~ as.vector(y_multidyn3_post)))$coefficients[,1])
  results[["MULTIDYN3"]] = results_MULTIDYN3
  
  # VAR
  
  lagx = x_pre[(p_multi):(T0-1),]
  
  for (i in (2:p_multi)) {
    lagx = cbind(lagx, x_pre[(p_multi + 1 - i):(T0 - i),])
  }
  
  lagy = y_pre[(p_multi):(T0-1)]
  
  for (i in (1:(p_multi-1))) {
    lagy = cbind(lagy, y_pre[(p_multi - i):(T0 - 1- i)])
  }
  
  xfull = cbind(1, lagy, lagx)
  yfull = y_pre[(p_multi+1):(T0)]
  
  colnames(xfull) = c(paste0("c"),
                      paste0("lagy_", 1:(p_multi)),
                      paste0(paste0("lagx", 1:ncol(x_pre)), "_", rep(1:(p_multi), each= ncol(x_pre))))
  
  if (ncol(xfull) > nrow(xfull)) {
    y_VAR_pre = rep(NA, T0-p_multi)
    y_VAR_post = rep(NA, T1)
  } else {
    
    w_VAR = solve(t(xfull) %*% xfull) %*% (t(xfull) %*% yfull)
    y_VAR_pre = xfull %*% w_VAR
    
    # bulding x_full_post
    
    x_prepost = rbind(x_pre[((T0+1)-p_multi):T0,],
                      x_post[1:T1,])
    
    lagx_post = x_prepost[(p_multi):(nrow(x_prepost)-1),]
    
    for (i in (2:p_multi)) {
      lagx_post = cbind(lagx_post, x_prepost[(p_multi + 1 - i):(nrow(x_prepost) - i),])
    }
    
    y_prepost = c(y_pre[((T0+1)-p_multi):T0],
                  rep(NA, T1))
    
    lagy_post = y_prepost[p_multi:(length(y_prepost)-1)]
    
    for (i in (1:(p_multi-1))) {
      lagy_post = cbind(lagy_post, y_prepost[(p_multi - i):(length(y_prepost) - 1- i)])
    }
    
    xfull_post = cbind(lagy_post, lagx_post)
    
    y_VAR_post = rep(NA, T1)
    
    for (i in 1:(T1-1)) {
      y_VAR_post[i] = as.matrix(cbind(1, xfull_post))[i, ] %*% w_VAR
      
      # updating y in xfull_post
      
      xfull_post[i + 1, 1] = y_VAR_post[i]
      
      if (i + 2 <= T1 & p_multi >= 2){
        xfull_post[i + 2, 2] = y_VAR_post[i]
      } 
      if (i + 3 <= T1 & p_multi >= 3){
        xfull_post[i + 3, 3] = y_VAR_post[i]
      }  
      if (i + 4 <= T1 & p_multi >= 4){
        xfull_post[i + 4, 3] = y_VAR_post[i]
      }  
    }
    
    # last period
    y_VAR_post[T1] = as.matrix(cbind(1, xfull_post))[T1, ] %*% w_VAR
    
  }
  

  
  y_treat_VAR = as.data.frame(c(y_pre, y_post)) %>%
    rename(y = c(1))
  y_treat_VAR$y_hat = c(rep(NA, p_multi),
                        y_VAR_pre, 
                        y_VAR_post)
  

  y_var_forecast = y_treat_VAR %>% 
    # slice(-c(1:T0)) %>%
    # mutate(y = y - post_effect)
    slice(c(1:T0))
  
  ggplot(y_var_forecast) +
    aes(x = y_hat, y = y) +
    geom_point(shape = "circle", size = 1.5, colour = "#112446") +
    geom_smooth(span = 0.75, method = "lm") +
    theme_minimal()
  
  # matplot(ts(y_treat_VAR),
  #         type = "l",
  #         lty = 1,
  #         lwd = 2,
  #         main = "VAR",
  #         xlab = "Time",
  #         ylab = "Value")
  
  # Export to ggplot
  
  df_gg = y_treat_VAR %>%
    gather() %>%
    mutate(id = c(1:(T0+T1), 1:(T0+T1)))
  
  plot = ggplot(df_gg) +
    aes(x = id, y = value, colour = key) +
    geom_line(linewidth = 1.0) +
    scale_color_hue(direction = 1) +
    theme_minimal() +
    labs(title = "VAR Path",
         subtitle = paste0("rho_u = ", round(rho_u,4), ", Donors = ", J,", rho_factor = ", rho_factor))+
    theme(
      plot.title = element_text(size = 15L,
                                hjust = 0.5),
      plot.subtitle = element_text(size = 13L,
                                   hjust = 0.5))
  
  results[["Plots_VAR"]] = plot
  
  results_VAR = c()
  
  results_VAR["PRE_VAR_RMSPE"] = sqrt(mean((y_pre[(p_multi+1):T0] -  y_VAR_pre)^2)) 
  results_VAR["PRE_VAR_BIAS"] = mean(y_VAR_pre- y_pre[(p_multi+1):T0])
  results_VAR["PRE_VAR_VAR"] = mean((y_VAR_pre - mean(y_VAR_pre))^2)
  
  results_VAR["POST_VAR_RMSFE"] = sqrt(mean(((y_post-post_effect) - y_VAR_post)^2))
  results_VAR["POST_VAR_BIAS"] = mean(y_VAR_post - (y_post-post_effect))
  results_VAR["POST_VAR_VAR"] = mean((y_VAR_post - mean(y_VAR_post))^2)
  
  #results_VAR["MZ_REG"] = list(summary(lm(y ~ y_hat, data = y_var_forecast))$coefficients[,1])
  results_VAR["MZ_REG_post"] = list(summary(lm(as.vector(y_post-post_effect) ~ as.vector(y_VAR_post)))$coefficients[,1])
  results_VAR["MZ_REG_pre"] = list(summary(lm(as.vector(y_pre[(p_multi+1):T0]) ~ as.vector(y_VAR_pre)))$coefficients[,1])
  
  results[["VAR"]] = results_VAR
  #}
  
  return(results)
}
