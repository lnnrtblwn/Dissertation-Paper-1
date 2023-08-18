library(tidyverse)
library(gridExtra)
library(xtable)
library(readxl)
library(magrittr)

rm(list = ls())
options(scipen=999)

if (Sys.info()[6] == "jctoe"){
  setwd("C:/Promotion/SC_Paper/Documents/sc_sim Ferman/Ferman") 
} else {
  setwd("C:/Users/lbolwin/Documents/Diss/Topics/Synthetic Control/R-Scripts/Simulations/Results/Factor/20230817") 
}

# PREPARATION ----

files = list.files()
files = files[str_detect(files, "xlsx")]

df = read_excel(files[1]) %>% 
  mutate(Specification = str_replace(str_sub(files[1],-11,-6), "^_+", "")) %>% 
  select(Specification, everything())

for (i in 2:9) {
  df %<>% 
    bind_rows(read_excel(files[i]) %>% 
                mutate(Specification = str_replace(str_sub(files[i],-11,-6), "^_+", "")) %>% 
                select(Specification, everything()))
}

df = df %>% 
  mutate(MZ_SC = ifelse(SC_MZ_REG_pval > 0.05, 1, 0),
         MZ_OLS = ifelse(OLS_MZ_REG_pval > 0.05, 1, 0),
         MZ_REGOLS = ifelse(REGOLS_MZ_REG_pval > 0.05, 1, 0),
         MZ_NET = ifelse(NET_MZ_REG_pval > 0.05, 1, 0),
         MZ_FACTOR = ifelse(FACTOR_MZ_REG_pval > 0.05, 1, 0))

results_mean = df %>% 
  group_by(Donors, Specification) %>% 
  summarise_at(.vars = dplyr::vars(PRE_SC_RMSPE:MZ_FACTOR),
               .funs = list(mean = ~ mean(., na.rm = TRUE))) %>% 
  ungroup() %>% 
  mutate(Specification = case_when(
    Specification == "20_30" ~ 1,
    Specification == "20_20" ~ 2,
    Specification == "20_10" ~ 3,
    Specification == "50_30" ~ 4,
    Specification == "50_20" ~ 5,
    Specification == "50_10" ~ 6,
    Specification == "100_30" ~ 7,
    Specification == "100_20" ~ 8,
    Specification == "100_10" ~ 9)) %>% 
  mutate(Specification = as.factor(Specification)) %>% 
  mutate(Specification = recode_factor(Specification,
                              `1` = "20_30",
                              `2` = "20_20",
                              `3` = "20_10",
                              `4` = "50_30",
                              `5` = "50_20",
                              `6` = "50_10",
                              `7` = "100_30",
                              `8` = "100_20",
                              `9` = "100_10")) %>% 
  arrange(Specification)

names(results_mean) = str_remove(names(results_mean), "_mean$")

t_FACTOR = results_mean %>% 
  select(Specification, Donors, POST_FACTOR_RMSFE, POST_FACTOR_BIAS, MZ_FACTOR, POST_FACTOR_VAR) %>% 
  gather(Type, FACTOR, POST_FACTOR_RMSFE:POST_FACTOR_VAR, -c(Specification, Donors)) %>% 
  arrange(Specification)

t_SC = results_mean %>% 
  select(Specification, Donors, POST_SC_RMSFE, POST_SC_BIAS, MZ_SC, POST_SC_VAR) %>% 
  gather(Type, SC, POST_SC_RMSFE:POST_SC_VAR, -c(Specification, Donors)) %>% 
  arrange(Specification)

t_REGOLS = results_mean %>% 
  select(Specification, Donors, POST_REGOLS_RMSFE, POST_REGOLS_BIAS, MZ_REGOLS, POST_REGOLS_VAR) %>% 
  gather(Type, REGOLS, POST_REGOLS_RMSFE:POST_REGOLS_VAR, -c(Specification, Donors)) %>% 
  arrange(Specification)

t_NET = results_mean %>% 
  select(Specification, Donors, POST_NET_RMSFE, POST_NET_BIAS, MZ_NET, POST_NET_VAR) %>% 
  gather(Type, NET, POST_NET_RMSFE:POST_NET_VAR, -c(Specification, Donors)) %>% 
  arrange(Specification)

t_OLS = results_mean %>% 
  select(Specification, Donors, POST_OLS_RMSFE, POST_OLS_BIAS, MZ_OLS, POST_OLS_VAR) %>% 
  gather(Type, OLS, POST_OLS_RMSFE:POST_OLS_VAR, -c(Specification, Donors)) %>% 
  arrange(Specification) 

# LATEX-Export ----

my_bracket = function(x, i, j, bracket) {
  if (bracket == "[") {
    x[i, j] = sprintf("[%s]", x[i, j])
  } else if (bracket == "(") {
    x[i, j] = sprintf("(%s)", x[i, j])
  } else if (bracket == "{") {
    x[i, j] = sprintf("{%s}", x[i, j])
  }
  return(x)
}


# T1 ----
t1 = t_FACTOR %>% 
  filter(Donors == 5) %>% 
  {. ->> t1} %>% 
  mutate(T0 = unlist(str_split(t1$Specification, "_"))[seq(1,72, by = 2)],
         T1 = unlist(str_split(t1$Specification, "_"))[seq(2,73, by = 2)]) %>% 
  select(T0, T1, FACTOR) %>% 
  bind_cols(t_SC %>% 
              filter(Donors == 5) %>% 
              select(SC),
            t_REGOLS %>% 
              filter(Donors == 5) %>% 
              select(REGOLS),
            t_NET %>% 
              filter(Donors == 5) %>% 
              select(NET),
            t_OLS %>% 
              filter(Donors == 5) %>% 
              select(OLS)) %>% 
  mutate_at(vars(FACTOR:OLS), funs(format(round(., 4), nsmall = 4))) %>% 
  mutate_at(vars(FACTOR:OLS), funs(as.character(.)))
  
# Brackets for Var, parenthesis for Bias
for (i in seq(4,36, by = 4)) {
  for (j in 3:7) {
    t1 = my_bracket(t1, i, j, "{")
  }
}

for (i in seq(3,35, by = 4)) {
  for (j in 3:7) {
    t1 = my_bracket(t1, i, j, "[")
  }
}

for (i in seq(2,34, by = 4)) {
  for (j in 3:7) {
    t1 = my_bracket(t1, i, j, "(")
  }
}

# Remove unneeded T0/T1-entries
t1$T0 = NA
t1$T0[1] = "20"
t1$T0[13] = "50"
t1$T0[25] = "100"

t1$T1 = NA
t1$T1[1] = "30"
t1$T1[5] = "20"
t1$T1[9] = "10"
t1$T1[13] = "30"
t1$T1[17] = "20"
t1$T1[21] = "10"
t1$T1[25] = "30"
t1$T1[29] = "20"
t1$T1[33] = "10"

# T2 ----
t2 = t_FACTOR %>% 
  filter(Donors == 10) %>% 
  {. ->> t2} %>% 
  mutate(T0 = unlist(str_split(t2$Specification, "_"))[seq(1,72, by = 2)],
         T1 = unlist(str_split(t2$Specification, "_"))[seq(2,73, by = 2)]) %>% 
  select(T0, T1, FACTOR) %>% 
  bind_cols(t_SC %>% 
              filter(Donors == 10) %>% 
              select(SC),
            t_REGOLS %>% 
              filter(Donors == 10) %>% 
              select(REGOLS),
            t_NET %>% 
              filter(Donors == 10) %>% 
              select(NET),
            t_OLS %>% 
              filter(Donors == 10) %>% 
              select(OLS)) %>% 
  mutate_at(vars(FACTOR:OLS), funs(format(round(., 4), nsmall = 4))) %>% 
  mutate_at(vars(FACTOR:OLS), funs(as.character(.)))

# Brackets for Var, parenthesis for Bias
for (i in seq(4,36, by = 4)) {
  for (j in 3:7) {
    t2 = my_bracket(t2, i, j, "{")
  }
}

for (i in seq(3,35, by = 4)) {
  for (j in 3:7) {
    t2 = my_bracket(t2, i, j, "[")
  }
}

for (i in seq(2,34, by = 4)) {
  for (j in 3:7) {
    t2 = my_bracket(t2, i, j, "(")
  }
}

# Remove unneeded T0/T1-entries
t2$T0 = NA
t2$T0[1] = "20"
t2$T0[13] = "50"
t2$T0[25] = "100"

t2$T1 = NA
t2$T1[1] = "30"
t2$T1[5] = "20"
t2$T1[9] = "10"
t2$T1[13] = "30"
t2$T1[17] = "20"
t2$T1[21] = "10"
t2$T1[25] = "30"
t2$T1[29] = "20"
t2$T1[33] = "10"

# T3 ----
t3 = t_FACTOR %>% 
  filter(Donors == 15) %>% 
  {. ->> t3} %>% 
  mutate(T0 = unlist(str_split(t3$Specification, "_"))[seq(1,72, by = 2)],
         T1 = unlist(str_split(t3$Specification, "_"))[seq(2,73, by = 2)]) %>% 
  select(T0, T1, FACTOR) %>% 
  bind_cols(t_SC %>% 
              filter(Donors == 15) %>% 
              select(SC),
            t_REGOLS %>% 
              filter(Donors == 15) %>% 
              select(REGOLS),
            t_NET %>% 
              filter(Donors == 15) %>% 
              select(NET),
            t_OLS %>% 
              filter(Donors == 15) %>% 
              select(OLS)) %>% 
  mutate_at(vars(FACTOR:OLS), funs(format(round(., 4), nsmall = 4))) %>% 
  mutate_at(vars(FACTOR:OLS), funs(as.character(.)))

# Brackets for Var, parenthesis for Bias
for (i in seq(4,36, by = 4)) {
  for (j in 3:7) {
    t3 = my_bracket(t3, i, j, "{")
  }
}

for (i in seq(3,35, by = 4)) {
  for (j in 3:7) {
    t3 = my_bracket(t3, i, j, "[")
  }
}

for (i in seq(2,34, by = 4)) {
  for (j in 3:7) {
    t3 = my_bracket(t3, i, j, "(")
  }
}

# Remove unneeded T0/T1-entries
t3$T0 = NA
t3$T0[1] = "20"
t3$T0[13] = "50"
t3$T0[25] = "100"

t3$T1 = NA
t3$T1[1] = "30"
t3$T1[5] = "20"
t3$T1[9] = "10"
t3$T1[13] = "30"
t3$T1[17] = "20"
t3$T1[21] = "10"
t3$T1[25] = "30"
t3$T1[29] = "20"
t3$T1[33] = "10"

# T4 ----
t4 = t_FACTOR %>% 
  filter(Donors == 20) %>% 
  {. ->> t4} %>% 
  mutate(T0 = unlist(str_split(t4$Specification, "_"))[seq(1,72, by = 2)],
         T1 = unlist(str_split(t4$Specification, "_"))[seq(2,73, by = 2)]) %>% 
  select(T0, T1, FACTOR) %>% 
  bind_cols(t_SC %>% 
              filter(Donors == 20) %>% 
              select(SC),
            t_REGOLS %>% 
              filter(Donors == 20) %>% 
              select(REGOLS),
            t_NET %>% 
              filter(Donors == 20) %>% 
              select(NET),
            t_OLS %>% 
              filter(Donors == 20) %>% 
              select(OLS)) %>% 
  mutate_at(vars(FACTOR:OLS), funs(format(round(., 4), nsmall = 4))) %>% 
  mutate_at(vars(FACTOR:OLS), funs(as.character(.)))

# Brackets for Var, parenthesis for Bias
for (i in seq(4,36, by = 4)) {
  for (j in 3:7) {
    t4 = my_bracket(t4, i, j, "{")
  }
}

for (i in seq(3,35, by = 4)) {
  for (j in 3:7) {
    t4 = my_bracket(t4, i, j, "[")
  }
}

for (i in seq(2,34, by = 4)) {
  for (j in 3:7) {
    t4 = my_bracket(t4, i, j, "(")
  }
}

# Remove unneeded T0/T1-entries
t4$T0 = NA
t4$T0[1] = "20"
t4$T0[13] = "50"
t4$T0[25] = "100"

t4$T1 = NA
t4$T1[1] = "30"
t4$T1[5] = "20"
t4$T1[9] = "10"
t4$T1[13] = "30"
t4$T1[17] = "20"
t4$T1[21] = "10"
t4$T1[25] = "30"
t4$T1[29] = "20"
t4$T1[33] = "10"

# T5 ----
t5 = t_FACTOR %>% 
  filter(Donors == 25) %>% 
  {. ->> t5} %>% 
  mutate(T0 = unlist(str_split(t5$Specification, "_"))[seq(1,72, by = 2)],
         T1 = unlist(str_split(t5$Specification, "_"))[seq(2,73, by = 2)]) %>% 
  select(T0, T1, FACTOR) %>% 
  bind_cols(t_SC %>% 
              filter(Donors == 25) %>% 
              select(SC),
            t_REGOLS %>% 
              filter(Donors == 25) %>% 
              select(REGOLS),
            t_NET %>% 
              filter(Donors == 25) %>% 
              select(NET),
            t_OLS %>% 
              filter(Donors == 25) %>% 
              select(OLS)) %>% 
  mutate_at(vars(FACTOR:OLS), funs(format(round(., 4), nsmall = 4))) %>% 
  mutate_at(vars(FACTOR:OLS), funs(as.character(.)))

# Brackets for Var, parenthesis for Bias
for (i in seq(4,36, by = 4)) {
  for (j in 3:7) {
    t5 = my_bracket(t5, i, j, "{")
  }
}

for (i in seq(3,35, by = 4)) {
  for (j in 3:7) {
    t5 = my_bracket(t5, i, j, "[")
  }
}

for (i in seq(2,34, by = 4)) {
  for (j in 3:7) {
    t5 = my_bracket(t5, i, j, "(")
  }
}

# Remove unneeded T0/T1-entries
t5$T0 = NA
t5$T0[1] = "20"
t5$T0[13] = "50"
t5$T0[25] = "100"

t5$T1 = NA
t5$T1[1] = "30"
t5$T1[5] = "20"
t5$T1[9] = "10"
t5$T1[13] = "30"
t5$T1[17] = "20"
t5$T1[21] = "10"
t5$T1[25] = "30"
t5$T1[29] = "20"
t5$T1[33] = "10"

# T6 ----
t6 = t_FACTOR %>% 
  filter(Donors == 30) %>% 
  {. ->> t6} %>% 
  mutate(T0 = unlist(str_split(t6$Specification, "_"))[seq(1,72, by = 2)],
         T1 = unlist(str_split(t6$Specification, "_"))[seq(2,73, by = 2)]) %>% 
  select(T0, T1, FACTOR) %>% 
  bind_cols(t_SC %>% 
              filter(Donors == 30) %>% 
              select(SC),
            t_REGOLS %>% 
              filter(Donors == 30) %>% 
              select(REGOLS),
            t_NET %>% 
              filter(Donors == 30) %>% 
              select(NET),
            t_OLS %>% 
              filter(Donors == 30) %>% 
              select(OLS)) %>% 
  mutate_at(vars(FACTOR:OLS), funs(format(round(., 4), nsmall = 4))) %>% 
  mutate_at(vars(FACTOR:OLS), funs(as.character(.)))

# Brackets for Var, parenthesis for Bias
for (i in seq(4,36, by = 4)) {
  for (j in 3:7) {
    t6 = my_bracket(t6, i, j, "{")
  }
}

for (i in seq(3,35, by = 4)) {
  for (j in 3:7) {
    t6 = my_bracket(t6, i, j, "[")
  }
}

for (i in seq(2,34, by = 4)) {
  for (j in 3:7) {
    t6 = my_bracket(t6, i, j, "(")
  }
}

# Remove unneeded T0/T1-entries
t6$T0 = NA
t6$T0[1] = "20"
t6$T0[13] = "50"
t6$T0[25] = "100"

t6$T1 = NA
t6$T1[1] = "30"
t6$T1[5] = "20"
t6$T1[9] = "10"
t6$T1[13] = "30"
t6$T1[17] = "20"
t6$T1[21] = "10"
t6$T1[25] = "30"
t6$T1[29] = "20"
t6$T1[33] = "10"

# SAVE ----

t1 = xtable(t1, include.rownames = FALSE)
t1 = print.xtable(t1, include.rownames = FALSE)
writeLines(t1, "Latex_Export/t1.tex")

t2 = xtable(t2, include.rownames = FALSE)
t2 = print.xtable(t2, include.rownames = FALSE)
writeLines(t2, "Latex_Export/t2.tex")

t3 = xtable(t3, include.rownames = FALSE)
t3 = print.xtable(t3, include.rownames = FALSE)
writeLines(t3, "Latex_Export/t3.tex")

t4 = xtable(t4, include.rownames = FALSE)
t4 = print.xtable(t4, include.rownames = FALSE)
writeLines(t4, "Latex_Export/t4.tex")

t5 = xtable(t5, include.rownames = FALSE)
t5 = print.xtable(t5, include.rownames = FALSE)
writeLines(t5, "Latex_Export/t5.tex")

t6 = xtable(t6, include.rownames = FALSE)
t6 = print.xtable(t6, include.rownames = FALSE)
writeLines(t6, "Latex_Export/t6.tex")

# FIGURES ----

# color_codes <- viridis::viridis_pal(option = "mako")(10)
# custom_colors <- c(color_codes[2], color_codes[3],
#                    color_codes[5], color_codes[6], 
#                    color_codes[8])

## 01 Bias of SC ----

df_no = read_excel("Bias/Factor_results_50_20_null.xlsx") %>% 
  filter(Donors == 5) %>% 
  select(POST_SC_BIAS, POST_REGOLS_BIAS) %>% 
  mutate(Intercept = "Treatment = Donor = N(0,1)")

df_neg = read_excel("Bias/Factor_results_50_20_neg.xlsx") %>% 
  select(POST_SC_BIAS, POST_REGOLS_BIAS) %>% 
  mutate(Intercept = "Treatment = N(1,1), Donor = N(0,1)")

df_pos = read_excel("Bias/Factor_results_50_20_pos.xlsx") %>% 
  select(POST_SC_BIAS, POST_REGOLS_BIAS) %>% 
  mutate(Intercept = "Treatment = N(-1,1), Donor = N(0,1)")

df = df_no %>% 
  bind_rows(df_neg,
            df_pos) %>% 
  mutate(Intercept = factor(Intercept))


p1 = ggplot(df) +
  aes(x = POST_SC_BIAS,
      fill = Intercept,
      colour = Intercept) +
  geom_density(adjust = 1L, alpha = 0.3) +
  scale_fill_hue(direction = 1) +
  #scale_x_continuous(breaks = seq(-3,3,1), limits = c(-3,3))+
  scale_color_hue(direction = 1) +
  theme_minimal()


p2 = ggplot(df) +
  aes(x = POST_REGOLS_BIAS,
    fill = Intercept,
    colour = Intercept) +
  geom_density(adjust = 1L, alpha = 0.3) +
  scale_fill_hue(direction = 1) +
  #scale_x_continuous(breaks = seq(-3,3,1), limits = c(-3,3))+
  scale_color_hue(direction = 1) +
  theme_minimal()

library(ggpubr)
ggarrange(p1, p2, nrow = 1, ncol = 2,
          common.legend = TRUE, legend = "bottom")



## 02 Dot-Plots ----

df_RMSFE = df %>% 
  mutate(Specification = str_extract(Specification, "^[^_]+")) %>% 
  group_by(Donors, Specification) %>% 
  summarise_at(.vars = dplyr::vars(PRE_SC_RMSPE:FACTOR_MZ_REG_pval),
               .funs = mean) %>% 
  ungroup() %>% 
  select(c("Donors", "Specification", "POST_SC_RMSFE", "POST_OLS_RMSFE", "POST_REGOLS_RMSFE", "POST_NET_RMSFE", "POST_FACTOR_RMSFE")) %>% 
  gather(type, value, POST_SC_RMSFE:POST_FACTOR_RMSFE, -c(Donors, Specification)) %>% 
  mutate(type = case_when(
    type == "POST_FACTOR_RMSFE" ~ 1,
    type == "POST_SC_RMSFE" ~ 2,
    type == "POST_REGOLS_RMSFE" ~ 3,
    type == "POST_NET_RMSFE" ~ 4,
    type == "POST_OLS_RMSFE" ~ 5)) %>% 
  mutate(type = dplyr::recode_factor(type,
                                     `1` = "FACTOR",
                                     `2` = "SC",
                                     `3` = "REGSC",
                                     `4` = "NET",
                                     `5` = "OLS"))

custom_colors = c("#DEF5E5", "#54C9AD", "#348AA6", "#40498E", "black")

p1 = ggplot(df_RMSFE %>%
         filter(Specification == "20"))+
  aes(x = Donors, y = value) + 
  geom_point(aes(fill = type), 
             colour = "black",
             pch = 21, size = 2.5, alpha = .95, stroke = .8) +
  #scale_fill_manual(values = custom_colors) +
  scale_fill_viridis_d(option = "mako")+
  labs(x = "Number of Donors", y = "RMSE-Average", fill = "Method")+
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L))

p2 = ggplot(df_RMSFE %>%
              filter(Specification == "50"))+
  aes(x = Donors, y = value) + 
  geom_point(aes(fill = type), 
             colour = "black",
             pch = 21, size = 2.5, alpha = .95, stroke = .8) +
  #scale_fill_manual(values = custom_colors) +
  scale_fill_viridis_d(option = "mako")+
  labs(x = "Number of Donors", y = "RMSE-Average", fill = "Method") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L))

p3 = ggplot(df_RMSFE %>%
              filter(Specification == "100"))+
  aes(x = Donors, y = value) + 
  geom_point(aes(fill = type), 
             colour = "black",
             pch = 21, size = 2.5, alpha = .95, stroke = .8) +
  #scale_fill_manual(values = custom_colors) +
  scale_fill_viridis_d(option = "mako")+
  labs(x = "Number of Donors", y = "RMSE-Average", fill = "Method") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L))

png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F06.png', width = 7, height = 3.5, units = 'in', res = 1000)
p1
dev.off()

png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F07.png', width = 7, height = 3.5, units = 'in', res = 1000)
p2
dev.off()

png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F08.png', width = 7, height = 3.5, units = 'in', res = 1000)
p3
dev.off()


## 03 Densities ----

df_density = df %>% 
  select(POST_REGOLS_BIAS, REGOLS_MZ_REG_pval)

p4 = ggplot(df_density) +
  aes(x = POST_REGOLS_BIAS, y = REGOLS_MZ_REG_pval) +
  geom_point(size = .5, colour = "black", shape=21) +
  labs(x = "Bias", y = "MZ p-Values") +
  theme_minimal()+
  theme(
    axis.title.y = element_text(size = 12L),
    axis.title.x = element_text(size = 12L),
    plot.caption = element_text(size = 10L))

png('~/Diss/Topics/Synthetic Control/Latex/Paper/images/F09.png', width = 7, height = 3.5, units = 'in', res = 1000)
p4
dev.off()

