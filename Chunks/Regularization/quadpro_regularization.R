library(quadprog)
library(MASS)
library(tidyverse)


obs = 1000
my_var = 1.5
k = 10

n = k
mu = c(rep(1, times = n))
A = matrix(runif(n * n), nrow = n, ncol = n)
A = t(A) %*% A
diag(A) = n / my_var

# Generating, de-mean data and splitting it
df = MASS::mvrnorm(n = obs,
                   mu = mu,
                   Sigma = A,
                   tol = T) %>% 
  as.data.frame()
colnames(df) = c("y", paste0("x", 1:(k-1)))

summary(lm(y ~ ., data = df))

ggplot(df) +
  aes(x = y, y = x1) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  geom_smooth(span = 0.75, method = "lm") +
  theme_minimal()

# regularization.

df_m = as.matrix(df)
y = df_m[,1]
X = df_m[,-1]

Dmat = t(X) %*% X
dvec = t(X) %*% y
Amat = t(rbind(rep(1, ncol(X)), diag(ncol(X)), -1*diag(ncol(X))))
bvec = c(1, rep(0, ncol(X)), rep(-1,ncol(X)))
synth_model = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_sc = synth_model$solution
y_sc_pre = X_pre %*% w_sc

