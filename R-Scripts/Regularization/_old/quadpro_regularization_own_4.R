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
  as.matrix

y = df[,1] %>% 
  as.matrix 

I = diag(k-1)
U = matrix(data = 1, nrow = k-1, ncol = 1)

x = scale(x)
y = scale(y)

l1 = 1
l2 = 1

beta = solve(t(x) %*% x + l1 * I + l2 * (U %*% t(U))) %*% ((t(x) %*% y) + l2 * U)
sum(solve(t(x) %*% x + l1 * I + l2 * (U %*% t(U))) %*% ((t(x) %*% y) + l2 * U))

library(glmnet)


# Ridge stimmt überein
solve(t(x) %*% x + l1 * I) %*% ((t(x) %*% y))
glmnet(x, y, alpha = 0, lambda = l1, standardize = F)$beta






