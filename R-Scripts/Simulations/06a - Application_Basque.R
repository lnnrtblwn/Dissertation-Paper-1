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

# 01 BASQUE TERROR ----

basque <- read.csv("~/Diss/Topics/Synthetic Control/R-Scripts/Simulations/Application Data/basque_data.csv") %>% 
  select(-`Unnamed..0`)

## 01.1 Levels ----

# SC (original)

dataprep.out <- dataprep(
  foo = basque,
  predictors = c(
    "school.illit",
    "school.prim",
    "school.med",
    "school.high",
    "school.post.high",
    "invest"),
  predictors.op = "mean",
  time.predictors.prior = 1964:1969,
  special.predictors = list(
    list("gdpcap", 1960:1969 , "mean"),
    list("sec.agriculture", seq(1961, 1969, 2), "mean"),
    list("sec.energy", seq(1961, 1969, 2), "mean"),
    list("sec.industry", seq(1961, 1969, 2), "mean"),
    list("sec.construction", seq(1961, 1969, 2), "mean"),
    list("sec.services.venta", seq(1961, 1969, 2), "mean"),
    list("sec.services.nonventa", seq(1961, 1969, 2), "mean"),
    list("popdens", 1969, "mean")),
  dependent = "gdpcap",
  unit.variable = "regionno",
  unit.names.variable = "regionname",
  time.variable = "year",
  treatment.identifier = 17,
  controls.identifier = c(2:16, 18),
  time.optimize.ssr = 1960:1969,
  time.plot = 1955:1997)

synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS")

sc1 = round(synth.out$solution.w,5) %>% 
  as.data.frame() %>% 
  {. ->> sc1} %>% 
  mutate(regionno = as.numeric(rownames(sc1)))

x_full = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  spread(year, gdpcap) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

sc1 = sc1 %>% 
  select(-regionno) %>% 
  as.matrix()

y_pred = x_full %*% sc1

df_result = basque[,c(1,3,4)] %>% 
  filter(regionno == 17)

df_result$SC_original= as.numeric(y_pred)

# SC (constraint)

x_pre = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  filter(year <= 1970) %>% 
  spread(year, gdpcap) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

x_post = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  filter(year > 1970) %>% 
  spread(year, gdpcap) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

y_pre = df_result %>% 
  filter(year <= 1970) %>% 
  select(gdpcap) %>% 
  as.matrix()

y_post = df_result %>% 
  filter(year > 1970) %>% 
  select(gdpcap) %>% 
  as.matrix()

Dmat = t(x_pre) %*% x_pre
dvec = t(x_pre) %*% y_pre
Amat = t(rbind(rep(1, ncol(x_pre)), diag(ncol(x_pre)), -1*diag(ncol(x_pre))))
bvec = c(1, rep(0, ncol(x_pre)), rep(-1,ncol(x_pre)))
sc2 = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_sc2 = sc2$solution %>% 
  as.matrix()

y_pred = x_full %*% w_sc2
df_result$SC_no_covariates = as.numeric(y_pred)

# REGOLS

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
T0 = 16
T1 = 27

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

xfull_d= cbind(lagx, lagy)
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

df_result_long = df_result %>%
  rename(`SC (original)` = SC_original,
         `SC (no covariates)` = SC_no_covariates,
         `Original Series` = gdpcap) %>% 
  select(-regionno) %>% 
  select(year, `Original Series`, everything()) %>% 
  gather(type, value, `Original Series`:MULTIDYN) 

color_codes <- viridis::viridis_pal(option = "mako")(10)
custom_colors <- c(color_codes[2], color_codes[3],
                   color_codes[5], color_codes[6], 
                   color_codes[8], color_codes[9])

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

rm(list = setdiff(ls(), c("df_result_long", "df_result", "basque")))
source("R-Scripts/Simulations/Functions/my_functions.R")

## 01.2 Differences ----

d_basque = basque %>% 
  group_by(regionno) %>% 
  summarise(year = year ,
            g_gdp = gdpcap - lag(gdpcap))

basque = basque %>% 
  left_join(d_basque, by = c("year" = "year",
                             "regionno" = "regionno"), na_matches = "never") %>% 
  select(-gdpcap) %>% 
  rename(gdpcap = g_gdp) %>% 
  select(regionno, regionname, year, gdpcap, everything()) %>% 
  filter(year >= 1956)
rm(d_basque)

# SC (original)

dataprep.out <- dataprep(
  foo = basque,
  predictors = c(
    "school.illit",
    "school.prim",
    "school.med",
    "school.high",
    "school.post.high",
    "invest"),
  predictors.op = "mean",
  time.predictors.prior = 1964:1969,
  special.predictors = list(
    list("gdpcap", 1960:1969 , "mean"),
    list("sec.agriculture", seq(1961, 1969, 2), "mean"),
    list("sec.energy", seq(1961, 1969, 2), "mean"),
    list("sec.industry", seq(1961, 1969, 2), "mean"),
    list("sec.construction", seq(1961, 1969, 2), "mean"),
    list("sec.services.venta", seq(1961, 1969, 2), "mean"),
    list("sec.services.nonventa", seq(1961, 1969, 2), "mean"),
    list("popdens", 1969, "mean")),
  dependent = "gdpcap",
  unit.variable = "regionno",
  unit.names.variable = "regionname",
  time.variable = "year",
  treatment.identifier = 17,
  controls.identifier = c(2:16, 18),
  time.optimize.ssr = 1960:1969,
  time.plot = 1956:1997)

synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS")
sc1 = round(synth.out$solution.w,5) %>% 
  as.data.frame() %>% 
  {. ->> sc1} %>% 
  mutate(regionno = as.numeric(rownames(sc1)))

x_full = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  spread(year, gdpcap) %>%
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

sc1 = sc1 %>% 
  select(-regionno) %>% 
  as.matrix()

y_pred = x_full %*% sc1

df_result_d = basque[,c(1,3,4)] %>% 
  filter(regionno == 17)
df_result_d$SC_original= as.numeric(y_pred)

# SC (no covariates)

x_pre = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  filter(year <= 1970) %>% 
  spread(year, gdpcap) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

x_post = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  filter(year > 1970) %>% 
  spread(year, gdpcap) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

y_pre = df_result_d %>% 
  filter(year <= 1970) %>% 
  select(gdpcap) %>% 
  as.matrix()

check = cor(y_pre, x_pre) %>% 
  t() %>% 
  as.data.frame()

# remove 5 donors with smallest correlation to make computation work

x_pre_SC = x_pre[,-c(5,12,11)]

y_post = df_result_d %>% 
  filter(year > 1970) %>% 
  select(gdpcap) %>% 
  as.matrix()


Dmat = t(x_pre_SC) %*% x_pre_SC
dvec = t(x_pre_SC) %*% y_pre
Amat = t(rbind(rep(1, ncol(x_pre_SC)), diag(ncol(x_pre_SC)), -1*diag(ncol(x_pre_SC))))
bvec = c(1, rep(0, ncol(x_pre_SC)), rep(-1,ncol(x_pre_SC)))
sc2 = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_sc2 = sc2$solution %>% 
  as.matrix()

y_pred = x_full[,-c(5,12,11)] %*% w_sc2
df_result_d$SC_no_covariates = as.numeric(y_pred)
rm(x_pre_SC)

# REGOLS

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

w_regols = regularized_ols(
  x = x_pre %>%
    as.data.frame(),
  y = y_pre %>%
    as.data.frame(),
  l1 = best_params_REGOLS[["l1"]],
  l2 = best_params_REGOLS[["l2"]])

y_pred = cbind(1,x_full) %*% w_regols
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
T0 = 15
T1 = 27

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

# visualization

df_result_long_d = df_result_d %>%
  rename(`SC (original)` = SC_original,
         `SC (no covariates)` = SC_no_covariates,
         `Original Series` = gdpcap) %>% 
  select(-regionno) %>% 
  select(year, `Original Series`, everything()) %>% 
  gather(type, value, `Original Series`:MULTIDYN) 

color_codes <- viridis::viridis_pal(option = "mako")(10)
custom_colors <- c(color_codes[2], color_codes[3],
                   color_codes[5], color_codes[6], 
                   color_codes[8], color_codes[9])

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

rm(list = setdiff(ls(), c("df_result", "df_result_d", "df_result_long", "df_result_long_d")))

## 01.3 Visualizations ----


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
  geom_vline(xintercept = c(1970), linetype= "dashed", color = "red", linewidth = .8)+
  labs(x = "Year", y = "GDP per Capita") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L)) + 
  guides(color = guide_legend(title = "Method", override.aes = list(linetype = c("dotted", rep("solid", times = 5)))),
         linetype  = "none")

# png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F03.png', width = 7, height = 3.5, units = 'in', res = 1000)
# p1
# dev.off()

# Differences

p2 = ggplot(df_result_long_d) +
  aes(x = year, y = value, colour = type, linetype = type) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c("dotted", rep("solid", times = 5)))+
  geom_vline(xintercept = c(1970), linetype= "dashed", color = "red", linewidth = .8)+
  labs(x = "Year", y = "First Difference of GDP per Capita") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L)) + 
  guides(color = guide_legend(title = "Method", override.aes = list(linetype = c("dotted", rep("solid", times = 5)))),
         linetype  = "none")

# png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F10.png', width = 7, height = 3.5, units = 'in', res = 1000)
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
    limits = c(-4,12),
    breaks = seq(-4, 12, by = 4),  
    labels = seq(-4, 12, by = 4),
    sec.axis = sec_axis(~ . / my_coef, name = "First Difference of GDP per Capita")) +
  scale_color_manual(values = custom_colors) +
  scale_linetype_manual(values = c("dotted", rep("solid", times = 5)))+
  geom_vline(xintercept = c(1970), linetype= "dashed", color = "red", linewidth = .8)+
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L),
    legend.position="bottom") + 
  guides(color = guide_legend(title = "Method", override.aes = list(linetype = c("dotted", rep("solid", times = 5)))),
         linetype  = "none")

# png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F12.png', width = 7, height = 5, units = 'in', res = 1000)
# p3
# dev.off()





