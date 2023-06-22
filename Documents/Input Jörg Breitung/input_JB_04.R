# Unrestricted ----

# Theoretical

A = matrix(c(1, 0.1, 0.4, 0.1, 1, 0.5, 0.4, 0.5, 1), nrow = 3)

A = matrix(c(1, 0.999, 0.999, 0.999, 1, 0.999, 0.999, 0.999, 1), nrow = 3)

w1 = solve(A[2:3, 2:3]) %*% A[2:3, 1]

1 - w1[1] * 1 - w1[2] * 1

# Empirical 

library(MASS)

y = mvrnorm(n = 10000, mu = c(1, 1, 1), Sigma = A)

w2 = lm(y[,1] ~ y[,-1])$coefficients
fit = lm(y[,1] ~ y[,-1])$fitted.values

round(var(c(y[,1] - fit)), 4)

# Restricted ----

library(quadprog)

Dmat = t(y[,-1]) %*% y[,-1]
dvec = t(y[,-1]) %*% y[,1]
Amat = t(rbind(rep(1,ncol(y[,-1])),diag(ncol(y[,-1])),-1*diag(ncol(y[,-1]))))
bvec = c(1, rep(0, ncol(y[,-1])), rep(-1,ncol(y[,-1])))

synth_model = solve.QP(Dmat, dvec, Amat, bvec, meq = 1)

w3 = synth_model$solution

fit = t(t(matrix(w3)) %*% t(y[,-1]))
var(c(y[,1] - fit))

# Manually computing the Variances:

# Var(Y_0 - aY_1 - bY_2) = Var(Y_0) + a^2 * Var(Y_1) + b^2 * Var(Y_2) - 2 * a * Cov(Y_0, Y_1) - 2 * b * Cov(Y_0, Y_2) + 2 * a * b * Cov(Y_1, Y_2)
# 0.2 vs. -0.1333
# 0.8 vs. 0.4667

1 + 0.2^2 + 0.8^2 - 
  2* 0.2* 0.1 -
  2* 0.8* 0.4 +
  2* 0.2* 0.8* 0.5

1 + 0.1333^2 + 0.4667^2 +
  2* 0.1333* 0.1 -
  2* 0.4667* 0.4 -
  2* 0.1333* 0.4667* 0.5


