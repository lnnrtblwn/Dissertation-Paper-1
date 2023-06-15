# Theoretical

A = matrix(c(1, 0.1, 0.4, 0.1, 1, 0.5, 0.4, 0.5, 1), nrow = 3)

w1 = solve(A[2:3, 2:3]) %*% A[2:3, 1]

1 - w1[1] * 1 - w1[2] * 1

# Empirical

library(MASS)

y = mvrnorm(n = 100000, mu = c(1, 1, 1), Sigma = A)

w2 = lm(y[,1] ~ y[,-1])$coefficients
