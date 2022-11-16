library(tidyverse)

N = 250
delta = 0.8
variance = 1
treat = 0

# simulate random walk with drift plus treatment after half of time period
# weiterüberlegen und hier eine Funktion bauen mit der MCMC gemacht werden kann
RW1 <- function(N, delta, variance, treat) {
  z = cumsum(rnorm(n=N/2, mean=0, sd=sqrt(variance)))
  t = 1:(N/2)
  x1 = t*delta + z
  
  z = cumsum(rnorm(n=N/2, mean=0, sd=sqrt(variance)))
  t = ((N/2)+1) :N
  x2 = t*delta + z + treat
  x = c(x1, x2)
  return(x)}

RW2 <- function(N, delta, variance) {
  z = cumsum(rnorm(n=N, mean=0, sd=sqrt(variance)))
  t = 1:N
  x = t*delta + z
  return(x)}

df_treat = RW1(250,0.5,1, 25) %>% 
  as.data.frame() %>% 
  mutate(time = 1:250) %>% 
  rename(treated = c(1)) %>% 
  select(time, treated)

ggplot(df_treat) +
  aes(x = time, y = treated) +
  geom_line(size = 1, colour = "#112446") +
  geom_vline(xintercept = 125, linetype="dotted", colour = "red", size = 1) + 
  theme_minimal()

# including covariates

df = df_treat %>% 
  bind_cols(replicate(5, RW2(250, 0.5, 1))) %>% 
  rename_with(.fn = ~paste0("x", 1:5), 
              .cols = c(3:7)) %>% 
  select(time, everything()) %>% 
  gather(type, value, treated:x5)

ggplot(df) +
  aes(x = time, y = value, colour = type, linetype = type) +
  geom_line(size = 1) +
  scale_color_hue(direction = 1) +
  scale_linetype_manual(values=c("solid", "dotted","dotted","dotted","dotted","dotted")) +
  theme_minimal()


