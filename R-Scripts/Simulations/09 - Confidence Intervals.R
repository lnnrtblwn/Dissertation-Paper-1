library(Synth)
library(tidysynth)
library(tidyverse)

rm(list = ls())

if (Sys.info()[6] == "jctoe"){
  setwd("C:/Promotion/SC_Paper/") 
} else {
  setwd("~/Diss/Topics/Synthetic Control/") 
}

source("R-Scripts/Simulations/Functions/my_functions.R")

basque <- read.csv("~/Diss/Topics/Synthetic Control/R-Scripts/Simulations/Application Data/basque_data.csv") %>% 
  select(-`Unnamed..0`)

df_result = basque[,c(1,3,4)] %>% 
  filter(regionno == 17)

# SC (constraint)

x_pre = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  filter(year <= 1970) %>% 
  spread(year, gdpcap) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

x_post = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  filter(year > 1970) %>% 
  spread(year, gdpcap) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

y_pre = df_result %>% 
  filter(year <= 1970) %>% 
  select(gdpcap) %>% 
  as.matrix()

y_post = df_result %>% 
  filter(year > 1970) %>% 
  select(gdpcap) %>% 
  as.matrix()

Dmat = t(x_pre) %*% x_pre
dvec = t(x_pre) %*% y_pre
Amat = t(rbind(rep(1, ncol(x_pre)), diag(ncol(x_pre)), -1*diag(ncol(x_pre))))
bvec = c(1, rep(0, ncol(x_pre)), rep(-1,ncol(x_pre)))
sc2 = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_sc2 = sc2$solution %>% 
  as.matrix()

x_full = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  spread(year, gdpcap) %>% 
  select(-c(1)) %>% 
  t() %>% 
  as.matrix()

y_pred = x_full %*% w_sc2

df_result$SC_no_covariates = as.numeric(y_pred)

# compute bootstrap CI

data = cbind(y_pre, x_pre)
custom_model = function(df){
  
  x = df[,-1]
  y = df[,1]
  
  Dmat = t(x) %*% x
  dvec = t(x) %*% y
  Amat = t(rbind(rep(1, ncol(x)), diag(ncol(x)), -1*diag(ncol(x))))
  bvec = c(1, rep(0, ncol(x)), rep(-1,ncol(x)))
  sc2 = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
    
  w_sc2 = sc2$solution %>% 
      as.matrix()
    
  y_hat <- x %*% w_sc2
  residuals <- y - y_hat
  residual_se <- sqrt(sum(residuals^2) / (length(y) - ncol(x)))
  
  return(residual_se)

}

n_reps = 1500
result = matrix(NA, nrow = n_reps, ncol = 1)

block_size <- 2
num_ids <- nrow(data) - block_size + 1
ids <- c()
blocks <- c()

for (block in 1:(num_ids - block_size + 1)) {
  ids <- c(ids, seq(block, block + block_size - 1))
  blocks <- c(blocks, rep(block, block_size))
}
block <- data.frame(ids = ids, blocks = blocks)

for (i in 1:n_reps) {
  boot_block = sample(x = c(1:n_distinct(blocks)), size = 30, replace = T) 
  boot_sample = c()
  
  for (j in 1:30) {
    block_ids = block %>%
      filter(blocks == boot_block[j]) %>%
      pull(ids)
    
    boot_sample <- c(boot_sample, block_ids)
  }
  
  boot_result = tryCatch({
    custom_model(data[boot_sample, ])
  },
  
  error = function(e) {
    return(list(NA))
  })
  boot_result

  result[i,1] = as.numeric(boot_result)
}

hist(result[, 1],  
     main = "Histogram of Bootstrapped SE",
     xlab = "SE",
     ylab = "Frequency",
     col = "lightblue",
     border = "black")
sum(is.na(result))
sd = mean(result, na.rm = T)

t = qt(0.975, df = 16)

df_result = df_result %>% 
  mutate(CI_helper = c(rep(0, times = nrow(x_pre)), 1:nrow(x_post))) %>% 
  mutate(CI_lo = ifelse(CI_helper == 0, SC_no_covariates - t * sd, SC_no_covariates - t * sd * CI_helper^(1/2)),
         CI_hi = ifelse(CI_helper == 0, SC_no_covariates + t * sd, SC_no_covariates + t * sd * CI_helper^(1/2))) %>% 
  select(-CI_helper)

# visualize

df_result_long = df_result %>%
  rename(`SC (no covariates)` = SC_no_covariates,
         `Original Series` = gdpcap) %>% 
  select(-regionno) %>% 
  select(year, `Original Series`, everything()) %>% 
  gather(type, value, `Original Series`:CI_hi) 

color_codes <- viridis::viridis_pal(option = "mako")(10)
custom_colors <- c(color_codes[2], color_codes[3],
                   color_codes[5], color_codes[6], 
                   color_codes[8], color_codes[9])

df_result_long = df_result_long %>% 
  mutate(type = case_when(
    type == "Original Series" ~ 1,
    type == "SC (no covariates)" ~ 2,
    type == "CI_lo" ~ 3,
    type == "CI_hi" ~ 4)) %>% 
  mutate(type = dplyr::recode_factor(type,
                                     `1` = "Original Series",
                                     `2` = "SC (no covariates)",
                                     `3` = "lower 95%-CI",
                                     `4` = "upper 95%-CI"))

ggplot(df_result_long) +
  aes(x = year, y = value, colour = type, linetype = type) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(custom_colors[1:2], rep(custom_colors[3], times = 2))) +
  scale_linetype_manual(values = c("solid", "dotted", rep("dashed", times = 2))) +
  geom_vline(xintercept = c(1970), linetype= "dashed", color = "red", linewidth = .8)+
  labs(x = "Year", y = "GDP per Capita") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L)) + 
  guides(color = guide_legend(title = "Method", override.aes = list(linetype = c("solid", "dotted", rep("dashed", times = 2)))),
         linetype  = "none")



