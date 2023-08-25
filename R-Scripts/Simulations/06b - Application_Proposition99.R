library(Synth)
library(tidysynth)
library(tidyverse)

rm(list = ls())

if (Sys.info()[6] == "jctoe"){
  setwd("C:/Promotion/SC_Paper/") 
} else {
  setwd("~/Diss/Topics/Synthetic Control/") 
}

source("R-Scripts/Simulations/Functions/my_functions.R")

# 02 PROPOSITION 99 ----

smoking <- read.csv("R-Scripts/Simulations/Application Data/smoking_data.csv") 

## 02.1 Levels ----

# SC (original)

smoking_out = smoking %>%
  synthetic_control(outcome = cigsale, # outcome
                    unit = state, # unit index in the panel data
                    time = year, # time index in the panel data
                    i_unit = "California", # unit where the intervention occurred
                    i_time = 1988, # time period when the intervention occurred
                    generate_placebos = F) %>%
  generate_predictor(time_window = 1980:1988,
                     ln_income = mean(lnincome, na.rm = T),
                     ret_price = mean(retprice, na.rm = T),
                     youth = mean(age15to24, na.rm = T)) %>%
  generate_predictor(time_window = 1984:1988,
                     beer_sales = mean(beer, na.rm = T)) %>%
  generate_predictor(time_window = 1975,
                     cigsale_1975 = cigsale) %>%
  generate_predictor(time_window = 1980,
                     cigsale_1980 = cigsale) %>%
  generate_predictor(time_window = 1988,
                     cigsale_1988 = cigsale) %>%
  generate_weights(optimization_window = 1970:1988, # time to use in the optimization task
                   margin_ipop = .02,sigf_ipop = 7,bound_ipop = 6) %>% 
  generate_control()

weights = smoking_out$.unit_weights[[1]]

df_result = smoking_out %>% 
  grab_synthetic_control() %>% 
  rename(Time = time_unit,
         SC_original = synth_y)

# SC (no covariates)

# Too many donors for 19 pre-treatment observations

# # Way 1: only selects Donors with positive weight
# non_zero_weight = weights %>% 
#   filter(weight > 0.01) %>% 
#   select(unit) %>% 
#   pull()

# Way 2: only selects Donors with large correlation
x = smoking %>% 
  filter(year <= 1988) %>% 
  select(state, cigsale) %>% 
  group_by(state) %>% 
  mutate(row_id = row_number()) %>% 
  ungroup() %>% 
  spread(key = state, value = cigsale) %>% 
  select(-c(row_id))

top_cor = cor(x) %>% 
  as.data.frame() %>% 
  select(California) %>% 
  arrange(desc(California)) %>% 
  slice(c(2:16))

large_corr = rownames(top_cor)
rm(x)

x_pre = smoking_out$.original_data[[2]] %>% 
  select(year, state, cigsale) %>% 
  filter(state %in% large_corr) %>% 
  arrange(year, state) %>% 
  filter(year <= 1988) %>% 
  spread(year, cigsale) %>% 
  select(-state) %>% 
  t() %>% 
  as.matrix()

y_pre = smoking_out$.original_data[[1]] %>% 
  filter(year <= 1988) %>% 
  select(cigsale) %>% 
  as.matrix()

Dmat = t(x_pre) %*% x_pre
dvec = t(x_pre) %*% y_pre
Amat = t(rbind(rep(1, ncol(x_pre)), diag(ncol(x_pre)), -1*diag(ncol(x_pre))))
bvec = c(1, rep(0, ncol(x_pre)), rep(-1,ncol(x_pre)))
sc2 = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_sc2 = sc2$solution %>% 
  as.matrix()

x_full = smoking_out$.original_data[[2]] %>% 
  select(year, state, cigsale) %>% 
  filter(state %in% large_corr) %>% 
  arrange(year, state) %>% 
  spread(year, cigsale) %>% 
  select(-state) %>% 
  t() %>% 
  as.matrix()

y_pred = x_full %*% w_sc2
df_result$SC_no_covariates = as.numeric(y_pred)

# REGOLS

x_pre = smoking_out$.original_data[[2]] %>% 
  select(year, state, cigsale) %>% 
  arrange(year, state) %>% 
  filter(year <= 1988) %>% 
  spread(year, cigsale) %>% 
  select(-state) %>% 
  t() %>% 
  as.matrix()

x_full = smoking_out$.original_data[[2]] %>% 
  select(year, state, cigsale) %>% 
  arrange(year, state) %>% 
  spread(year, cigsale) %>% 
  select(-state) %>% 
  t() %>% 
  as.matrix()

param_grid = expand.grid(
  l1 = 5^seq(1, 5, length.out = 50),
  l2 = 10^seq(1, 7, length.out = 50)) %>%
  as.data.frame() %>%
  sample_n(400) %>%
  bind_rows(c(l1 = 1,
              l2 = 1))

param_grid$RMSPE = NA
param_grid$RMSFE = NA
CV_share = .5

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

y_pred = cbind(1, x_full) %*% w_regols
df_result$REGSC = as.numeric(y_pred)

# NET

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

y_pred = predict(cvfit, newx = x_full, s = best_params$lambda)
df_result$NET = as.numeric(y_pred)

# MULTIDYN

p_multi = 1
T0 = 19
T1 = 12

x_post = x_full[c(20:31),] 
# dim(x_full) = 31,38
# x_post: 1989-2000

lagx = x_pre[(p_multi + 1):T0,]

for (i in (1:p_multi)) {
  lagx = cbind(lagx, x_pre[(p_multi + 1 - i):(T0 - i),])
}

lagy = y_pre[(p_multi):(T0-1)]

if (p_multi > 1){
  for (i in (1:(p_multi-1))) {
    lagy = cbind(lagy, y_pre[(p_multi - i):(T0 - 1- i)])
  }
}

xfull_d = cbind(lagx, lagy)
yfull_d = y_pre[(p_multi+1):(T0)]

colnames(xfull_d) = c(paste0(paste0("lagx", 1:ncol(x_pre)), "_", rep(0:(p_multi), each= ncol(x_pre))),
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
  
  cvfit = cv.glmnet(xfull_d, yfull_d, alpha = a, type.measure = "mse", nfolds = 3)
  cv_fit$lambda[i] = cvfit$lambda.min
  cv_fit$alpha[i] = a
  cv_fit$CVM[i] = min(cvfit$cvm)
  i = i+1
}

best_params = cv_fit[cv_fit$CVM == min(cv_fit$CVM),]
cvfit = cv.glmnet(xfull_d, yfull_d, alpha = best_params$alpha, type.measure = "mse", nfolds = 3)

y_multidyn3_pre = predict(cvfit, newx = xfull_d, s = best_params$lambda)

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

if (p_multi > 1){
  for (i in (1:(p_multi-1))) {
    lagy_post = cbind(lagy_post, y_prepost[(p_multi - i):(length(y_prepost) - 1- i)])
  }
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

df_result$MULTIDYN = c(rep(NA, p_multi),
                       y_multidyn3_pre,
                       y_multidyn3_post)

df_result = df_result %>% 
  select(Time, real_y, everything()) %>% 
  rename(`SC (no covariates)` = SC_no_covariates,
         `SC (original)` = SC_original,
         `Original Series` = real_y)

df_result_long = df_result %>%
  gather(type, value, `Original Series`:MULTIDYN) 

df_result_long = df_result_long %>% 
  mutate(type = case_when(
    type == "Original Series" ~ 1,
    type == "SC (original)" ~ 2,
    type == "SC (no covariates)" ~ 3,
    type == "NET" ~ 4,
    type == "REGSC" ~ 5,
    type == "MULTIDYN" ~ 6)) %>% 
  mutate(type = dplyr::recode_factor(type,
                                     `1` = "Original Series",
                                     `2` = "SC (original)",
                                     `3` = "SC (no covariates)",
                                     `4` = "NET",
                                     `5` = "REGSC",
                                     `6` = "MULTIDYN"))

ggplot(df_result_long) +
  aes(x = Time, y = value, colour = type) +
  geom_line() +
  scale_color_hue(direction = 1) +
  theme_minimal()

rm(list = setdiff(ls(), c("df_result_long", "df_result", "smoking")))
source("R-Scripts/Simulations/Functions/my_functions.R")

## 02.2 Differences ----

d_smoking = smoking %>% 
  group_by(state) %>% 
  summarise(year = year ,
            g_cigsale = cigsale - lag(cigsale))

smoking = smoking %>% 
  left_join(d_smoking, by = c("year" = "year",
                             "state" = "state"), na_matches = "never") %>% 
  select(-cigsale) %>% 
  rename(cigsale = g_cigsale) %>% 
  select(state, year, cigsale, everything()) %>% 
  filter(year >= 1971)
rm(d_smoking)

# SC (original)

smoking_out = smoking %>%
  synthetic_control(outcome = cigsale, # outcome
                    unit = state, # unit index in the panel data
                    time = year, # time index in the panel data
                    i_unit = "California", # unit where the intervention occurred
                    i_time = 1988, # time period when the intervention occurred
                    generate_placebos = F) %>%
  generate_predictor(time_window = 1980:1988,
                     ln_income = mean(lnincome, na.rm = T),
                     ret_price = mean(retprice, na.rm = T),
                     youth = mean(age15to24, na.rm = T)) %>%
  generate_predictor(time_window = 1984:1988,
                     beer_sales = mean(beer, na.rm = T)) %>%
  generate_predictor(time_window = 1975,
                     cigsale_1975 = cigsale) %>%
  generate_predictor(time_window = 1980,
                     cigsale_1980 = cigsale) %>%
  generate_predictor(time_window = 1988,
                     cigsale_1988 = cigsale) %>%
  generate_weights(optimization_window = 1971:1988, # time to use in the optimization task
                   margin_ipop = .02,sigf_ipop = 7,bound_ipop = 6) %>% 
  generate_control()

weights = smoking_out$.unit_weights[[1]]

df_result_d = smoking_out %>% 
  grab_synthetic_control() %>% 
  rename(Time = time_unit,
         SC_original = synth_y)

# SC (no covariates)

# Too many donors for 19 pre-treatment observations

# # Way 1: only selects Donors with positive weight
# non_zero_weight = weights %>% 
#   filter(weight > 0.01) %>% 
#   select(unit) %>% 
#   pull()

# Way 2: only selects Donors with large correlation
x = smoking %>% 
  filter(year <= 1988) %>% 
  select(state, cigsale) %>% 
  group_by(state) %>% 
  mutate(row_id = row_number()) %>% 
  ungroup() %>% 
  spread(key = state, value = cigsale) %>% 
  select(-c(row_id))

top_cor = cor(x) %>% 
  as.data.frame() %>% 
  select(California) %>% 
  arrange(desc(California)) %>% 
  slice(c(2:16))

large_corr = rownames(top_cor)
rm(x)

x_pre = smoking_out$.original_data[[2]] %>% 
  select(year, state, cigsale) %>% 
  filter(state %in% large_corr) %>% 
  arrange(year, state) %>% 
  filter(year <= 1988) %>% 
  spread(year, cigsale) %>% 
  select(-state) %>% 
  t() %>% 
  as.matrix()

y_pre = smoking_out$.original_data[[1]] %>% 
  filter(year <= 1988) %>% 
  select(cigsale) %>% 
  as.matrix()

Dmat = t(x_pre) %*% x_pre
dvec = t(x_pre) %*% y_pre
Amat = t(rbind(rep(1, ncol(x_pre)), diag(ncol(x_pre)), -1*diag(ncol(x_pre))))
bvec = c(1, rep(0, ncol(x_pre)), rep(-1,ncol(x_pre)))
sc2 = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_sc2 = sc2$solution %>% 
  as.matrix()

x_full = smoking_out$.original_data[[2]] %>% 
  select(year, state, cigsale) %>% 
  filter(state %in% large_corr) %>% 
  arrange(year, state) %>% 
  spread(year, cigsale) %>% 
  select(-state) %>% 
  t() %>% 
  as.matrix()

y_pred = x_full %*% w_sc2
df_result_d$SC_no_covariates = as.numeric(y_pred)

# REGOLS

x_pre = smoking_out$.original_data[[2]] %>% 
  select(year, state, cigsale) %>% 
  arrange(year, state) %>% 
  filter(year <= 1988) %>% 
  spread(year, cigsale) %>% 
  select(-state) %>% 
  t() %>% 
  as.matrix()

x_full = smoking_out$.original_data[[2]] %>% 
  select(year, state, cigsale) %>% 
  arrange(year, state) %>% 
  spread(year, cigsale) %>% 
  select(-state) %>% 
  t() %>% 
  as.matrix()

param_grid = expand.grid(
  l1 = 5^seq(1, 5, length.out = 50),
  l2 = 10^seq(1, 7, length.out = 50)) %>%
  as.data.frame() %>%
  sample_n(400) %>%
  bind_rows(c(l1 = 1,
              l2 = 1))

param_grid$RMSPE = NA
param_grid$RMSFE = NA
CV_share = .5

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

y_pred = cbind(1, x_full) %*% w_regols
df_result_d$REGSC = as.numeric(y_pred)

# NET

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

y_pred = predict(cvfit, newx = x_full, s = best_params$lambda)
df_result_d$NET = as.numeric(y_pred)

# MULTIDYN

p_multi = 1
T0 = 18
T1 = 12

x_post = x_full[c(19:30),]
# dim(x_full) = 30,38
# x_post: 1989-2000

lagx = x_pre[(p_multi + 1):T0,]

for (i in (1:p_multi)) {
  lagx = cbind(lagx, x_pre[(p_multi + 1 - i):(T0 - i),])
}

lagy = y_pre[(p_multi):(T0-1)]

if (p_multi > 1){
  for (i in (1:(p_multi-1))) {
    lagy = cbind(lagy, y_pre[(p_multi - i):(T0 - 1- i)])
  }
}

xfull_d = cbind(lagx, lagy)
yfull_d = y_pre[(p_multi+1):(T0)]

colnames(xfull_d) = c(paste0(paste0("lagx", 1:ncol(x_pre)), "_", rep(0:(p_multi), each= ncol(x_pre))),
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
  
  cvfit = cv.glmnet(xfull_d, yfull_d, alpha = a, type.measure = "mse", nfolds = 3)
  cv_fit$lambda[i] = cvfit$lambda.min
  cv_fit$alpha[i] = a
  cv_fit$CVM[i] = min(cvfit$cvm)
  i = i+1
}

best_params = cv_fit[cv_fit$CVM == min(cv_fit$CVM),]
cvfit = cv.glmnet(xfull_d, yfull_d, alpha = best_params$alpha, type.measure = "mse", nfolds = 3)

y_multidyn3_pre = predict(cvfit, newx = xfull_d, s = best_params$lambda)

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

if (p_multi > 1){
  for (i in (1:(p_multi-1))) {
    lagy_post = cbind(lagy_post, y_prepost[(p_multi - i):(length(y_prepost) - 1- i)])
  }
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

df_result_d$MULTIDYN = c(rep(NA, p_multi),
                       y_multidyn3_pre,
                       y_multidyn3_post)

df_result_d = df_result_d %>% 
  select(Time, real_y, everything()) %>% 
  rename(`SC (no covariates)` = SC_no_covariates,
         `SC (original)` = SC_original,
         `Original Series` = real_y)

df_result_long_d = df_result_d %>%
  gather(type, value, `Original Series`:MULTIDYN) 

df_result_long_d = df_result_long_d %>% 
  mutate(type = case_when(
    type == "Original Series" ~ 1,
    type == "SC (original)" ~ 2,
    type == "SC (no covariates)" ~ 3,
    type == "NET" ~ 4,
    type == "REGSC" ~ 5,
    type == "MULTIDYN" ~ 6)) %>% 
  mutate(type = dplyr::recode_factor(type,
                                     `1` = "Original Series",
                                     `2` = "SC (original)",
                                     `3` = "SC (no covariates)",
                                     `4` = "NET",
                                     `5` = "REGSC",
                                     `6` = "MULTIDYN"))

ggplot(df_result_long_d) +
  aes(x = Time, y = value, colour = type) +
  geom_line() +
  scale_color_hue(direction = 1) +
  theme_minimal()

rm(list = setdiff(ls(), c("df_result_long_d", "df_result_d", "df_result", "df_result_long")))

## 02.3 Visualization ----

color_codes <- viridis::viridis_pal(option = "mako")(10)
custom_colors <- c(color_codes[2], color_codes[3],
                   color_codes[5], color_codes[6], 
                   color_codes[8], color_codes[9])

# Levels

p1 = ggplot(df_result_long) +
  aes(x = Time, y = value, colour = type, linetype = type) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c("dotted", rep("solid", times = 5)))+
  geom_vline(xintercept = c(1988), linetype= "dashed", color = "red", linewidth = .8)+
  labs(x = "Year", y = "Cigarette Sales (in packs) per Capita") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L)) + 
  guides(color = guide_legend(title = "Method", override.aes = list(linetype = c("dotted", rep("solid", times = 5)))),
         linetype  = "none")

png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F04.png', width = 7, height = 3.5, units = 'in', res = 1000)
p1
dev.off()

# Differences

p2 = ggplot(df_result_long_d) +
  aes(x = Time, y = value, colour = type, linetype = type) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c("dotted", rep("solid", times = 5)))+
  geom_vline(xintercept = c(1988), linetype= "dashed", color = "red", linewidth = .8)+
  labs(x = "Year", y = "First Diff. of Cigarette Sales (in packs) per Capita") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L)) + 
  guides(color = guide_legend(title = "Method", override.aes = list(linetype = c("dotted", rep("solid", times = 5)))),
         linetype  = "none")

png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F13.png', width = 7, height = 3.5, units = 'in', res = 1000)
p2
dev.off()

# Joint plot

my_coef = 8

df_joint = df_result_long %>% 
  left_join(df_result_long_d %>% 
              rename(value_d = value),
            by = c("Time", "type")) %>% 
  mutate(value_d_scale = value_d * my_coef) %>% 
  rename(year = Time)


p3 = ggplot(df_joint) +
  geom_line(aes(x = year, y = value, colour = type, linetype = type), linewidth = 1) + 
  geom_line(aes(x = year, y = value_d_scale, colour = type, linetype = type), linewidth = .75) +
  labs(x = "Year", y = "Cigarette Sales (in pakcs) per Capita") +
  scale_y_continuous(
    # limits = c(-4,12),
    # breaks = seq(-4, 12, by = 4),  
    # labels = seq(-4, 12, by = 4),
    sec.axis = sec_axis(~ . / my_coef, name = "First Diff. of Cigarette Sales (in packs) per Capita")) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c("dotted", rep("solid", times = 5)))+
  geom_vline(xintercept = c(1988), linetype= "dashed", color = "red", linewidth = .8)+
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L),
    legend.position="bottom") + 
  guides(color = guide_legend(title = "Method", override.aes = list(linetype = c("dotted", rep("solid", times = 5)))),
         linetype  = "none")

png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F14.png', width = 7, height = 5, units = 'in', res = 1000)
p3
dev.off()





