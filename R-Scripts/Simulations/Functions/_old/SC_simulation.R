# Generate SC Series
SC_simulation = function(observations, n_variables) {
  
  params = tibble(beta = rep(NA, n_variables),
                  value = rep(NA, n_variables))
  
  df = tibble(y1 = rep(NA, observations))
  
  for (i in 1:n_variables) {
    df = df %>%
      bind_cols(rnorm(observations, mean = 1, sd = 0.5))
    
    params$beta[i] = paste("beta", i + 1, sep  = "")
    params$value[i] = runif(1)
    
    rm(i)
  }
  
  params = params %>% 
    mutate(value_sum = value / sum(value))
  
  df = df %>%
    rename_at(vars(1:ncol(.)), ~ paste("y", 1:length(.), sep = "")) %>% 
    mutate(y1 = rowSums(as.matrix(df[,-1]) %*% params$value_sum) + rnorm(observations, 0, 0.5))
  
  print(params)
  return(df)
}