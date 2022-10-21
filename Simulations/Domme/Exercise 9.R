# Stochastic Models and Processes

# Exercise 9

simSizePower = function(m, n, mux = 1, muy, sigma2x = 1, sigma2y = 1,
                        nreps = 10000, alpha = 0.05) {
  x = replicate(nreps, rnorm(n = m, mean = mux, sd = sqrt(sigma2x)))
  y = replicate(nreps, rnorm(n = n, mean = muy, sd = sqrt(sigma2y)))
  
  xbar = colMeans(x)
  s2x = apply(X = x, MAR = 2, FUN = var)
  ybar = colMeans(y)
  s2y = apply(X = y, MAR = 2, FUN = var)
  Te = sqrt(m * n * (m + n - 2) / (m + n)) * (xbar - ybar) / 
    sqrt((m - 1) * s2x + (n - 1) * s2y)
  mean(abs(Te) > qt(p = 1 - alpha / 2, df = m + n - 2))
}

simSizePower(m = 10, n = 10, muy = 1)
simSizePower(m = 20, n = 20, muy = 1)
simSizePower(m = 20, n = 30, muy = 1)

simSizePower(m = 10, n = 10, muy = 2)
simSizePower(m = 20, n = 20, muy = 2)
simSizePower(m = 20, n = 30, muy = 2)

simSizePower(m = 10, n = 10, muy = 3)
simSizePower(m = 20, n = 20, muy = 3)
simSizePower(m = 20, n = 30, muy = 3)