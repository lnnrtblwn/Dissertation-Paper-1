library(tidyverse)

load("~/Diss/Topics/Synthetic Control/Chunks/Simulations/MC_Outputs/MC_meta_k20.Rda")
# load("~/Diss/Topics/Synthetic Control/Chunks/Simulations/MC_Outputs/MC_meta_k05.Rda")
# load("~/Diss/Topics/Synthetic Control/Chunks/Simulations/MC_Outputs/MC_meta_k50.Rda")
#load("~/Diss/Topics/Synthetic Control/Chunks/Simulations/MC_Outputs/MC_meta_k50_Vsh20.Rda")

# Analyzing the rsults
df_meta = data.frame(matrix(NA, nrow = 0, ncol = 9)) %>% 
  rename(VAR_Share = c(1),
         RMSPE_SC = c(2),
         RMSPE_OLS = c(3),
         RMSPE_VAR = c(4),
         RMSPE_VAR_SC = c(5),
         RMSFE_SC = c(6),
         RMSFE_OLS = c(7),
         RMSFE_VAR = c(8),
         RMSFE_VAR_SC = c(9))

share_help = seq(0,1, by = 1/20)
iter =  50

for (i in 1:length(share_help)) {
  name = paste("data")
  assign(
    name,
    data.frame(matrix(NA, nrow = iter, ncol = 9)) %>%
      rename(VAR_Share = c(1),
             RMSPE_SC = c(2),
             RMSPE_OLS = c(3),
             RMSPE_VAR = c(4),
             RMSPE_VAR_SC = c(5),
             RMSFE_SC = c(6),
             RMSFE_OLS = c(7),
             RMSFE_VAR = c(8),
             RMSFE_VAR_SC = c(9)))
  
  data$VAR_Share = MCMC_meta[[i]]$VAR_Share
  
  data$RMSPE_SC = MCMC_meta[[i]]$RMSPE_SC
  data$RMSPE_OLS = MCMC_meta[[i]]$RMSPE_OLS
  data$RMSPE_VAR = MCMC_meta[[i]]$RMSPE_VAR
  data$RMSPE_VAR_SC = MCMC_meta[[i]]$RMSPE_VAR_SC
  
  data$RMSFE_SC = MCMC_meta[[i]]$RMSFE_SC
  data$RMSFE_OLS = MCMC_meta[[i]]$RMSFE_OLS
  data$RMSFE_VAR = MCMC_meta[[i]]$RMSFE_VAR
  data$RMSFE_VAR_SC = MCMC_meta[[i]]$RMSFE_VAR_SC
  
  df_meta = df_meta %>%
    bind_rows(data)
  rm(data)
}

# summary(df_meta %>% 
#           filter(VAR_Share == 0))
# 
# df_meta = df_meta %>%
#   filter(VAR_Share == 0) %>% 
#   gather(type, value, RMSPE_SC:RMSFE_VAR_SC, -VAR_Share) %>%
#   mutate(VAR_Share = as.factor(round(VAR_Share,2)))
# 
# p1 = df_meta %>%
#   filter(type %in% c("RMSPE_SC", "RMSPE_OLS", "RMSPE_VAR", "RMSPE_VAR_SC")) %>%
#   ggplot() +
#   aes(x = value, fill = type) +
#   geom_density(adjust = 1L, alpha = 0.5) +
#   scale_fill_viridis_d(option = "viridis", direction = 1) +
#   labs(title = "RMSPE")+
#   theme_minimal()
# 
# p2 = df_meta %>%
#   filter(type %in% c("RMSFE_SC", "RMSFE_OLS", "RMSFE_VAR", "RMSFE_VAR_SC")) %>%
#   ggplot() +
#   aes(x = value, fill = type) +
#   geom_density(adjust = 1L, alpha = 0.5) +
#   scale_fill_viridis_d(option = "viridis", direction = 1) +
#   labs(title = "RMSFE")+
#   theme_minimal()
# 
#  ggpubr::ggarrange(p1, p2, ncol=2)
 
df_meta = df_meta %>%
  gather(type, value, RMSFE_SC:RMSFE_VAR_SC, -VAR_Share) %>%
  mutate(VAR_Share = as.factor(round(VAR_Share,2)))

ggplot(df_meta) +
  aes(x = VAR_Share, y = value, fill = type, colour = type) +
  geom_boxplot() +
  scale_fill_viridis_d(option = "viridis", direction = 1) +
  scale_color_viridis_d(option = "viridis", direction = 1) +
  theme_minimal() +
  facet_wrap(vars(type))

df_meta = df_meta %>%
  group_by(VAR_Share) %>%
  mutate(RMSFE_Min = ifelse(value == min(value, na.rm = T), 1, 0)) %>%
  ungroup()


table(df_meta$type ,df_meta$RMSFE_Min)

test = df_meta %>% 
  filter(RMSFE_Min == 1) %>% 
  select(VAR_Share, type, value) %>% 
  arrange(VAR_Share)
