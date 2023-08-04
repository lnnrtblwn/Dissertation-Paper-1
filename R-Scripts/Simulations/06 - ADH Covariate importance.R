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

# wie gut performt SC?

test = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  spread(year, gdpcap) %>% 
  t() %>% 
  as.matrix()

test = test[-1,]
sc1 = sc1 %>% 
  select(-regionno) %>% 
  as.matrix()

y_pred = test %*% sc1
df_basque = basque[,c(1,3,4)] %>% 
  filter(regionno == 17)
df_basque$SC_original= as.numeric(y_pred)

# now only constraint regression.

x_pre = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  filter(year <= 1970) %>% 
  spread(year, gdpcap) %>% 
  t() %>% 
  as.matrix()

x_post = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  filter(year > 1970) %>% 
  spread(year, gdpcap) %>% 
  t() %>% 
  as.matrix()

x_pre = x_pre[-1,]
x_post = x_post[-1,]

y_pre = df_basque %>% 
  filter(year <= 1970) %>% 
  select(gdpcap) %>% 
  as.matrix()

y_post = df_basque %>% 
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

y_pred = test %*% w_sc2
df_basque$SC_no_covariates = as.numeric(y_pred)

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

# Extract best CV-Parameter Combination

w_regols = regularized_ols(
  x = x_pre %>%
    as.data.frame(),
  y = y_pre %>%
    as.data.frame(),
  l1 = best_params_REGOLS[["l1"]],
  l2 = best_params_REGOLS[["l2"]])

y_pred = cbind(1,test) %*% w_regols
df_basque$REGOLS = as.numeric(y_pred)

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

y_pred = predict(cvfit, newx = test, s = best_params$lambda)
df_basque$NET = as.numeric(y_pred)

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

df_basque$MULTIDYN = c(rep(NA, p_multi),
                       y_multidyn3_pre,
                       y_multidyn3_post)

# visualization

df_basque_long = df_basque %>%
  select(-regionno) %>% 
  gather(type, value, gdpcap:MULTIDYN) 


ggplot(df_basque_long) +
  aes(x = year, y = value, colour = type) +
  geom_line(size = 0.8) +
  geom_vline(xintercept=c(1970), linetype="dotted", size = 1)+
  scale_color_hue(direction = 1) +
  labs(x = "Year", y = "GDP per Capita", colour = "Method") +
  theme_minimal() + 
  theme(legend.position="bottom")

sqrt(mean((df_basque$gdpcap[1:16] - df_basque$prediction_sc[1:16])^2))
sqrt(mean((df_basque$gdpcap[1:16] - df_basque$prediction_sc_constr[1:16])^2))
sqrt(mean((df_basque$gdpcap[1:16] - df_basque$prediction_REGOLS[1:16])^2))
sqrt(mean((df_basque$gdpcap[1:16] - df_basque$prediction_NET[1:16])^2))
sqrt(mean((df_basque$gdpcap[2:16] - df_basque$prediction_MULTIDYN[2:16])^2))

test = df_basque %>% 
  mutate(diff_sc = gdpcap - prediction_sc,
         diff_sc_constr = gdpcap - prediction_sc_constr) %>% 
  filter(year <= 1970)

sqrt(mean((test$diff_sc)^2))
sqrt(mean((test$diff_sc_constr)^2))

# 02 PROPOSITION 99 ----

data("smoking")


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

df = smoking_out %>% grab_synthetic_control()

# weights = as.matrix(weights$weight)
# 
# data = smoking_out$.original_data[[2]] %>% 
#   select(c(1:3)) %>% 
#   spread(year, cigsale) %>% 
#   select(-c(1)) %>% 
#   t() %>% 
#   as.matrix()
# 
# dim(data)
# 
# data %*% weights

# restricted calculation

non_zero_weight = weights %>% 
  filter(weight > 0.01) %>% 
  select(unit) %>% 
  pull()

x_pre = smoking_out$.original_data[[2]] %>% 
  select(year, state, cigsale) %>% 
  filter(state %in% non_zero_weight) %>% 
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

x = smoking_out$.original_data[[2]] %>% 
  select(year, state, cigsale) %>% 
  filter(state %in% non_zero_weight) %>% 
  arrange(year, state) %>% 
  spread(year, cigsale) %>% 
  select(-state) %>% 
  t() %>% 
  as.matrix()

y_pred = x %*% w_sc2
df$synth_constrained_y = as.numeric(y_pred)

# visualization

df_long = df %>%
  gather(type, value, real_y:synth_constrained_y) 

ggplot(df_long) +
  aes(x = time_unit, y = value, colour = type) +
  geom_line(size = 0.8) +
  geom_vline(xintercept=c(1988), linetype="dotted")+
  scale_color_hue(direction = 1) +
  theme_minimal()

sqrt(mean((df$real_y[1:18] - df$synth_y[1:18])^2))
sqrt(mean((df$real_y[1:18] - df$synth_constrained_y[1:18])^2))
