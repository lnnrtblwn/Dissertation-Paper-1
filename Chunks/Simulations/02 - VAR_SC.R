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

set.seed(123)

# Generate sample
t = 200 
k = 3 
p = 1 

# Generate coefficient matrices. Adjust these values in simulation
A = matrix(c(-.3, .6, -.4,
              .5, .3, .4,
              .2, .3, .1), k) # Coefficient matrix of lag 1

# Generate series
series = matrix(0, k, t + 2*p) # Raw series with zeros
for (i in (p + 1):(t + 2 * p)) {
  series[, i] <- A %*% series[, i - 1] + rnorm(k, 0, .5) # Generate series with e ~ N(0,0.5)
}

series <- ts(t(series[, -(1:p)])) # Convert to time series format
names <- c("V1", "V2", "V3") # Rename variables

plot.ts(series) # Plot the series



# Estimation. Similar to OLS 
# y = constant + current donors values (d) + all lagged values (a).

d = as.data.frame(series[,-1]) %>% # take 1:t or 2:(t+1)? System is singular for 1:t
  slice(2:(t+1)) %>%
  as.matrix()

a = dplyr::lag(as.data.frame(series)) %>% 
  slice(2:(t + 1)) %>%
  as.matrix()

X = as.matrix(cbind(1, d, a))
y = as.matrix(series[c(2:201),1])

df = data.frame(y, X[,-1])

var_sc = lm(y ~ ., data = df)
summary(var_sc)

df$predicted = var_sc$fitted.values

df_plot = df %>% 
  mutate(time = 1:n()) %>% 
  select(time, y, predicted) %>% 
  gather(type, value, y:predicted)

ggplot(df_plot) +
  aes(x = time, y = value, colour = type) +
  geom_line(size = 1) +
  scale_color_hue(direction = 1) +
  theme_minimal()

  




