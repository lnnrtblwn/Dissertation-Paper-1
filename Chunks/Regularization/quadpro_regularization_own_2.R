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

# Applying Jörgs solution

obs = 3000
k = 5
my_var = 1.5*k/2

df = gen_data(obs, k)

X = df[,2:k] %>% 
  as.matrix %>% 
  scale(center=T,scale=T)

y = df[,1] %>% 
  as.matrix %>% 
  scale(center=T,scale=T)

# My approach

I = diag(k-1)
U = matrix(data = 1, nrow = k-1, ncol = 1)

l1 = 1
l2 = 1

solve(t(X) %*% X + l1 * I + l2 * (U %*% t(U))) %*% ((t(X) %*% y) + l2 * U)
round(sum(solve(t(X) %*% X + l1 * I + l2 * (U %*% t(U))) %*% ((t(X) %*% y) + l2 * U)), 3)

cv = 


