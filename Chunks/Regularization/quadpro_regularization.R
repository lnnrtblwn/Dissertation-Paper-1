library(quadprog)
library(MASS)
library(tidyverse)


obs = 20
my_var = 1.5
k = 25

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

df = gen_data(obs, k)

summary(lm(y ~ ., data = df))

ggplot(df) +
  aes(x = x1, y = y) +
  geom_point(shape = "circle", size = 1.5, colour = "#112446") +
  geom_smooth(span = 0.75, method = "lm") +
  theme_minimal()

# regularization.

lasso.qp <- function(X,y,t) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  D <- t(X) %*% X
  P <- 2*D
  d <- 2 * t(y) %*% X
  A <- as.matrix(expand.grid(rep(list(c(-1,1)),p)))
  b0 <- rep(-t,2^p)
  beta.lasso <- solve.QP(Dmat=P,dvec=d,Amat=t(A),bvec=b0)
  return(beta.lasso$solution)
}

# Load in the data
X = as.matrix(df[,-1]) %>% 
  scale(center=T,scale=F)
y = as.matrix(df[,1]) %>% 
  scale(center=T,scale=F)

lasso.qp(X = X,
         y = y, 
         t = 100) %>% round(2)

glmnet(X, y)

####### GLM packages

data(QuickStartExample)
n = 100 # all data

x <- QuickStartExample$x[1:n,]
y <- QuickStartExample$y[1:n,]

summary(lm(y ~ x))
fit <- glmnet(x, y)
plot(fit,label = TRUE)
# number of nonzero coefficients (Df)
# the percent (of null) deviance explained (%dev)
# the value of λ
print(fit)

coef(fit, s = 1)

predict(fit, newx = x, s = c(0.1, 0.05))

cvfit <- cv.glmnet(x, y)
plot(cvfit)

