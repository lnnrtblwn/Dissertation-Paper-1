library(Synth)
library(tidysynth)
library(tidyverse)
library(glmnet)

rm(list = ls())

if (Sys.info()[6] == "jctoe"){
  setwd("C:/Promotion/SC_Paper/") 
} else {
  setwd("~/Diss/Topics/Synthetic Control/") 
}

source("R-Scripts/Simulations/Functions/my_functions.R")


# 03 GERMAN UNIFICATION ----

unification = foreign::read.dta("R-Scripts/Simulations/Application Data/repgermany.dta")

## 03.1 Levels ----

# SC (original)

dataprep_out <-
  dataprep(
    foo = unification,
    predictors    = c("gdp","trade","infrate"),
    dependent     = "gdp",
    unit.variable = "index",
    time.variable = "year",
    special.predictors = list(
      list("industry" ,1981:1990, c("mean")),
      list("schooling",c(1980,1985), c("mean")),
      list("invest80" ,1980, c("mean"))
    ),
    treatment.identifier = 7,
    controls.identifier = unique(unification$index)[-7],
    time.predictors.prior = 1981:1990,
    time.optimize.ssr = 1960:1989,
    unit.names.variable = "country",
    time.plot = 1960:2003)

synth_out <- synth(dataprep_out)

sc1 = round(synth_out$solution.w,5) %>% 
  as.data.frame() %>% 
  {. ->> sc1} %>% 
  mutate(index = as.numeric(rownames(sc1)))

x_full = unification[,c(1,3,4)] %>% 
  filter(!index %in% c(7)) %>% 
  arrange(year, index) %>% 
  spread(year, gdp) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

sc1 = sc1 %>% 
  select(-index) %>% 
  as.matrix()

y_pred = x_full %*% sc1

df_results = unification[,c(1,3,4)] %>% 
  filter(index == 7)

df_results$SC_original= as.numeric(y_pred)

# SC (no covariates)

x_pre = unification[,c(1,3,4)] %>% 
  filter(!index %in% c(7)) %>% 
  arrange(year, index) %>%
  filter(year <= 1989) %>% 
  spread(year, gdp) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

x_post = unification[,c(1,3,4)] %>% 
  filter(!index %in% c(7)) %>% 
  arrange(year, index) %>%
  filter(year > 1989) %>% 
  spread(year, gdp) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

y_pre = df_results %>% 
  filter(year <= 1989) %>% 
  select(gdp) %>% 
  as.matrix()

y_post = df_results %>% 
  filter(year > 1989) %>% 
  select(gdp) %>% 
  as.matrix()

Dmat = t(x_pre) %*% x_pre
dvec = t(x_pre) %*% y_pre
Amat = t(rbind(rep(1, ncol(x_pre)), diag(ncol(x_pre)), -1*diag(ncol(x_pre))))
bvec = c(1, rep(0, ncol(x_pre)), rep(-1,ncol(x_pre)))
sc2 = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_sc2 = sc2$solution %>% 
  as.matrix()

y_pred = x_full %*% w_sc2
df_results$SC_no_covariates = as.numeric(y_pred)

# REGSC

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

y_pred = cbind(1,x_full) %*% w_regols
df_results$REGSC = as.numeric(y_pred)

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

y_net_pre = predict(cvfit, newx = x_pre, s = best_params$lambda)
y_net_post = predict(cvfit, newx = x_post, s = best_params$lambda)

y_pred = predict(cvfit, newx = x_full, s = best_params$lambda)
df_results$NET = as.numeric(y_pred)

# MULTIDYN

p_multi = 1
T0 = 30
T1 = 14

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

df_results$MULTIDYN = c(rep(NA, p_multi),
                  y_multidyn3_pre,
                  y_multidyn3_post)

df_result_long = df_results %>%
  rename(`SC (original)` = SC_original,
         `SC (no covariates)` = SC_no_covariates,
         `Original Series` = gdp) %>% 
  select(-index) %>% 
  select(year, `Original Series`, everything()) %>% 
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
  aes(x = year, y = value, colour = type) +
  geom_line() +
  scale_color_hue(direction = 1) +
  theme_minimal()

rm(list = setdiff(ls(), c("df_result_long", "df_results", "unification")))
source("R-Scripts/Simulations/Functions/my_functions.R")

## 03.2 Differences ----

d_unification = unification %>% 
  group_by(country) %>% 
  summarise(year = year ,
            g_gdp = gdp - lag(gdp))

unification = unification %>% 
  left_join(d_unification, by = c("year" = "year",
                              "country" = "country"), na_matches = "never") %>% 
  select(-gdp) %>% 
  rename(gdp = g_gdp) %>% 
  select(country, year, gdp, everything()) %>% 
  filter(year >= 1961)
rm(d_unification)


# SC (original)
dataprep_out <-
  dataprep(
    foo = unification,
    predictors    = c("gdp","trade","infrate"),
    dependent     = "gdp",
    unit.variable = "index",
    time.variable = "year",
    special.predictors = list(
      list("industry" ,1981:1990, c("mean")),
      list("schooling",c(1980,1985), c("mean")),
      list("invest80" ,1980, c("mean"))
    ),
    treatment.identifier = 7,
    controls.identifier = unique(unification$index)[-7],
    time.predictors.prior = 1981:1990,
    time.optimize.ssr = 1961:1989,
    unit.names.variable = "country",
    time.plot = 1961:2003)

synth_out <- synth(dataprep_out)

sc1 = round(synth_out$solution.w,5) %>% 
  as.data.frame() %>% 
  {. ->> sc1} %>% 
  mutate(index = as.numeric(rownames(sc1)))

x_full = unification %>% 
  select(c(index, year, gdp)) %>% 
  filter(!index %in% c(7)) %>% 
  arrange(year, index) %>% 
  spread(year, gdp) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

sc1 = sc1 %>% 
  select(-index) %>% 
  as.matrix()

y_pred = x_full %*% sc1

df_results_d = unification %>% 
  select(c(index, year, gdp)) %>% 
  filter(index == 7)

df_results_d$SC_original= as.numeric(y_pred)

# SC (no covariates)

x_pre = unification %>% 
  select(c(index, year, gdp)) %>% 
  filter(!index %in% c(7)) %>% 
  arrange(year, index) %>%
  filter(year <= 1989) %>% 
  spread(year, gdp) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

x_post = unification %>% 
  select(c(index, year, gdp)) %>% 
  filter(!index %in% c(7)) %>% 
  arrange(year, index) %>%
  filter(year > 1989) %>% 
  spread(year, gdp) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

y_pre = df_results_d %>% 
  filter(year <= 1989) %>% 
  select(gdp) %>% 
  as.matrix()

y_post = df_results_d %>% 
  filter(year > 1989) %>% 
  select(gdp) %>% 
  as.matrix()

Dmat = t(x_pre) %*% x_pre
dvec = t(x_pre) %*% y_pre
Amat = t(rbind(rep(1, ncol(x_pre)), diag(ncol(x_pre)), -1*diag(ncol(x_pre))))
bvec = c(1, rep(0, ncol(x_pre)), rep(-1,ncol(x_pre)))
sc2 = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_sc2 = sc2$solution %>% 
  as.matrix()

y_pred = x_full %*% w_sc2
df_results_d$SC_no_covariates = as.numeric(y_pred)

# REGSC

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

y_pred = cbind(1,x_full) %*% w_regols
df_results_d$REGSC = as.numeric(y_pred)

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

y_net_pre = predict(cvfit, newx = x_pre, s = best_params$lambda)
y_net_post = predict(cvfit, newx = x_post, s = best_params$lambda)

y_pred = predict(cvfit, newx = x_full, s = best_params$lambda)
df_results_d$NET = as.numeric(y_pred)

# MULTIDYN

p_multi = 1
T0 = 29
T1 = 14

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

df_results_d$MULTIDYN = c(rep(NA, p_multi),
                        y_multidyn3_pre,
                        y_multidyn3_post)

df_result_long_d = df_results_d %>%
  rename(`SC (original)` = SC_original,
         `SC (no covariates)` = SC_no_covariates,
         `Original Series` = gdp) %>% 
  select(-index) %>% 
  select(year, `Original Series`, everything()) %>% 
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
  aes(x = year, y = value, colour = type) +
  geom_line() +
  scale_color_hue(direction = 1) +
  theme_minimal()

rm(list = setdiff(ls(), c("df_results", "df_results_d", "df_result_long", "df_result_long_d")))

## 03.3 Visualizations ----


color_codes <- viridis::viridis_pal(option = "mako")(10)
custom_colors <- c(color_codes[2], color_codes[3],
                   color_codes[5], color_codes[6], 
                   color_codes[8], color_codes[9])

# Levels

p1 = ggplot(df_result_long) +
  aes(x = year, y = value, colour = type, linetype = type) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c("dotted", rep("solid", times = 5)))+
  geom_vline(xintercept = c(1990), linetype= "dashed", color = "red", linewidth = .8)+
  labs(x = "Year", y = "GDP per Capita") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L)) + 
  guides(color = guide_legend(title = "Method", override.aes = list(linetype = c("dotted", rep("solid", times = 5)))),
         linetype  = "none")

# png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F05.png', width = 7, height = 3.5, units = 'in', res = 1000)
# p1
# dev.off()

# Differences

p2 = ggplot(df_result_long_d) +
  aes(x = year, y = value, colour = type, linetype = type) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c("dotted", rep("solid", times = 5)))+
  geom_vline(xintercept = c(1990), linetype= "dashed", color = "red", linewidth = .8)+
  labs(x = "Year", y = "First Differences of GDP per Capita") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L)) + 
  guides(color = guide_legend(title = "Method", override.aes = list(linetype = c("dotted", rep("solid", times = 5)))),
         linetype  = "none")


# png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F15.png', width = 7, height = 3.5, units = 'in', res = 1000)
# p2
# dev.off()

# Joint plot

my_coef = 8

df_joint = df_result_long %>% 
  left_join(df_result_long_d %>% 
              rename(value_d = value),
            by = c("year", "type")) %>% 
  mutate(value_d_scale = value_d * my_coef)


p3 = ggplot(df_joint) +
  geom_line(aes(x = year, y = value, colour = type, linetype = type), linewidth = 1) + 
  geom_line(aes(x = year, y = value_d_scale, colour = type, linetype = type), linewidth = .75) +
  labs(x = "Year", y = "GDP per Capita") +
  scale_y_continuous(
    # limits = c(-4,12),
    # breaks = seq(-4, 12, by = 4),  
    # labels = seq(-4, 12, by = 4),
    sec.axis = sec_axis(~ . / my_coef, name = "First Difference of GDP per Capita")) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c("dotted", rep("solid", times = 5)))+
  geom_vline(xintercept = c(1990), linetype= "dashed", color = "red", linewidth = .8)+
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L),
    legend.position="bottom") + 
  guides(color = guide_legend(title = "Method", override.aes = list(linetype = c("dotted", rep("solid", times = 5)))),
         linetype  = "none")

# png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F16.png', width = 7, height = 5, units = 'in', res = 1000)
# p3
# dev.off()






