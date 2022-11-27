# Manual regression analysis

library(tidyverse)

height <- rnorm(200,175,sd = 15)
age = rnorm(200, 25, sd = 5)

weight <- height - 80 + 1.02 * (height)^0.01 * + 0.8*age + rnorm(200,0,sd = 15) * rnorm(200,1.1,sd = 0.2) 

df <- data.frame(height, age, weight)

ggplot(df, aes(height, weight)) +
  geom_point() +
  geom_smooth(method = "lm")

reg = lm(weight ~ height + age, data=df)

# Manually
X <- as.matrix(cbind(1,df$height, df$age))
y <- as.matrix(df$weight)

#  beta = ((X'X)^(-1))X'y
beta = solve(t(X) %*% X) %*% t(X) %*% y

# residuals
res = as.matrix(y - beta[1] - beta[2] * X[, 2] - beta[3] * X[, 3])

# define the number of observations (n) and the number of parameters (k)
n = nrow(df)
k = ncol(X)

# calculate the Variance-Covariance matrix (VCV)
#  VCV = (1/(n-k))res'res(X'X)^(-1)
VCV = 1 / (n - k) * as.numeric(t(res) %*% res) * solve(t(X) %*% X)

# calculate standard errors (se) of coefficients
se = sqrt(diag(VCV))

# calculate the p-values
p_value = rbind(2 * pt(abs(beta[1] / se[1]), df = n - k,
                       lower.tail = FALSE),
                2 * pt(abs(beta[2] / se[2]), df = n - k,
                       lower.tail = FALSE),
                2 * pt(abs(beta[3] / se[3]), df = n - k,
                       lower.tail = FALSE))

# combine all necessary information
output = as.data.frame(cbind(c("(Intercept)", "height", "age"),
                             beta, se, p_value))
names(output) <- c("Coefficients:", "Estimate",
                   "Std. Error", "Pr(>|t|)")

output
summary(reg)

# Nice, now the same for VAR-data