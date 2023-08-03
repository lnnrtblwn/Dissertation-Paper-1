library(glmnet)
library(tidyverse)
library(quadprog)

data("Boston", package = "MASS")

# Load in the data
X <- as.matrix(Boston[,c('rm','dis')]) %>% scale(center=T,scale=F)
y <- as.matrix(Boston[,'medv']) %>% scale(center=T,scale=F)

# Find the OLS estimates
ls.beta <- lm(y~-1+X)
t.max <- ls.beta  %>% coef %>% abs %>% sum

# Define the function
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

# Run the LASSO for 50% shrinkage
lasso.qp(X=X,y=y,t=0.5*t.max) %>% round(2)
