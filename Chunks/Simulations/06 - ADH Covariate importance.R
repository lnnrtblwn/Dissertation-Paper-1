library(Synth)
library(tidyverse)

data("basque")

basque[85:89, 1:4]

dataprep.out <- dataprep(
  foo=basque,
  predictors=c("school.illit","school.prim","school.med","school.high", "school.post.high", "invest"),
  predictors.op="mean",
  time.predictors.prior=1964:1969,
  special.predictors=list(
    list("gdpcap", 1960:1969 , "mean"),
    list("sec.agriculture", seq(1961, 1969, 2), "mean"),
    list("sec.energy", seq(1961, 1969, 2), "mean"),
    list("sec.industry", seq(1961, 1969, 2), "mean"),
    list("sec.construction", seq(1961, 1969, 2), "mean"),
    list("sec.services.venta", seq(1961, 1969, 2), "mean"),
    list("sec.services.nonventa", seq(1961, 1969, 2), "mean"),
    list("popdens", 1969, "mean")),
  dependent="gdpcap",
  unit.variable="regionno",
  unit.names.variable="regionname",
  time.variable="year",
  treatment.identifier=17,
  controls.identifier=c(2:16,18),
  time.optimize.ssr=1960:1969,
  time.plot=1955:1997)

synth.out = synth(data.prep.obj = dataprep.out, method = "BFGS")
sc1 = round(synth.out$solution.w,5) %>% 
  as.data.frame() %>% 
  {. ->> sc1} %>% 
  mutate(regionno = as.numeric(rownames(sc1)))

# wie gut performt SC?

test = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  spread(year, gdpcap) %>% 
  t() %>% 
  as.matrix()

test = test[-1,]
sc1 = sc1 %>% 
  select(-regionno) %>% 
  as.matrix()

y_pred = test %*% sc1
df_basque = basque[,c(1,3,4)] %>% 
  filter(regionno == 17)
df_basque$prediction_sc = as.numeric(y_pred)

# now only constraint regression.

x_pre = basque[,c(1,3,4)] %>% 
  filter(!regionno %in% c(1,17)) %>% 
  arrange(year, regionno) %>% 
  filter(year <= 1970) %>% 
  spread(year, gdpcap) %>% 
  t() %>% 
  as.matrix()

x_pre = x_pre[-1,]
y_pre = df_basque %>% 
  filter(year <= 1970) %>% 
  select(gdpcap) %>% 
  as.matrix()


Dmat = t(x_pre) %*% x_pre
dvec = t(x_pre) %*% y_pre
Amat = t(rbind(rep(1, ncol(x_pre)), diag(ncol(x_pre)), -1*diag(ncol(x_pre))))
bvec = c(1, rep(0, ncol(x_pre)), rep(-1,ncol(x_pre)))
sc2 = quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
w_sc2 = sc2$solution %>% 
  as.matrix()

y_pred = test %*% w_sc2
df_basque$prediction_sc_constr = as.numeric(y_pred)

# visualization

df_basque_long = df_basque %>%
  select(-regionno) %>% 
  gather(type, value, gdpcap:prediction_sc_constr) 


ggplot(df_basque_long) +
  aes(x = year, y = value, colour = type) +
  geom_line(size = 0.8) +
  geom_vline(xintercept=c(1970), linetype="dotted")+
  scale_color_hue(direction = 1) +
  theme_minimal()

sqrt(mean((df_basque$gdpcap[1:16] - df_basque$prediction_sc[1:16])^2))
sqrt(mean((df_basque$gdpcap[1:16] - df_basque$prediction_sc_constr[1:16])^2))

test = df_basque %>% 
  mutate(diff_sc = gdpcap - prediction_sc,
         diff_sc_constr = gdpcap - prediction_sc_constr) %>% 
  filter(year <= 1970)

sqrt(mean((test$diff_sc)^2))
sqrt(mean((test$diff_sc_constr)^2))

