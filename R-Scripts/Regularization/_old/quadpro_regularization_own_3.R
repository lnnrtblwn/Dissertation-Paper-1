
# PACKAGE INSTALLATION ####

suppressPackageStartupMessages({
  list.of.packages = c("tidyverse", 
                       "caret")
  new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages, require, character.only = T,)
})

# DATA SIMULATION ####

gen_data = function(obs, k){
  n = k
  mu = c(rep(1, times = n))
  A = matrix(runif(n * n), nrow = n, ncol = n)
  A = t(A) %*% A
  diag(A) = my_var
  
  df = MASS::mvrnorm(n = obs,
                     mu = mu,
                     Sigma = A,
                     tol = T) %>% 
    as.data.frame()
  
  colnames(df) = c("y", paste0("x", 1:(k-1)))
  
  return(df)
}

obs = 3000
k = 5
my_var = 1.5*k/2

df = gen_data(obs, k)

x = df[,2:k] %>% 
  as.matrix %>% 
  scale(center=T,scale=T)

y = df[,1] %>% 
  as.matrix %>% 
  scale(center=T,scale=T)

# REGULARIZATION INSIDE CARET #####

getModelInfo(model = "glmnet", regex = FALSE)[[1]]

modelInfo = list(label = "Linear Regression",
                  library = NULL,
                  type = "Regression",
                  parameters = data.frame(parameter = "intercept",
                                          class = "logical",
                                          label = "intercept"),
                  grid = function (x, y, len = NULL, search = "grid") 
                    data.frame(intercept = TRUE),
                  loop = NULL,
                  fit = function (x, y, wts, param, lev, last, classProbs, ...) 
                  {
                    dat = if(is.data.frame(x)) x else as.data.frame(x)
                    dat$y <- y
                    lm(y ~ ., data = dat)$coefficients
                  },

                  predict = function (modelFit, newdata, submodels = NULL) 
                  {
                    if (!is.data.frame(newdata)) 
                      newdata <- as.data.frame(newdata, stringsAsFactors = TRUE)
                    predict(modelFit, newdata)
                  },
                  prob = NULL,
                  levels = NULL,
                  sort = function(x) x)

train(y ~ ., 
              data = df, 
              method = modelInfo,
              trControl = trainControl(method = "repeatedcv",
                                                repeats = 5))

test <- train(y ~ ., 
              data = df, 
              method = "lm",
              trControl = trainControl(method = "repeatedcv",
                                       repeats = 5))

