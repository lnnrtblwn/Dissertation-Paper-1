


model <- gen_dfm(x = as.ts(df_spread_resid), p = 4, n = 3,
                 iterations = 5000, burnin = 1000)

object <- add_priors(model,
                     lambda = list(v_i = .01),
                     a = list(v_i = .01),
                     sigma_v = list(shape = 5, rate = 4),
                     sigma_u = list(shape = 5, rate = 4))

object = draw_posterior(object)

plot(object)




install.packages("bayesdfa")

library(bayesdfa)

set.seed(1)

sim_dat <- sim_dfa(
  num_trends = 2, 
  num_years = 20,
  num_ts = 15
)

matplot(t(sim_dat$y_sim), type = "l")

matplot(t(data.frame(diff(as.matrix(sim_dat$y_sim)))), type = "l")

as.data.frame(sim_dat$y_sim) %>%
  mutate_each(funs(. - lag(.))) %>%
  na.omit() 


x <- sim_dfa(num_trends = 2)
names(x)
matplot(t(x$y_sim), type = "l")
matplot(t(x$x), type = "l")

set.seed(42)
x <- sim_dfa(extreme_value = -4, extreme_loc = 10)
matplot(t(x$x), type = "l")
abline(v = 10)
matplot(t(x$pred), type = "l")
abline(v = 10)

set.seed(42)
x <- sim_dfa()
matplot(t(x$x), type = "l")
abline(v = 10)
matplot(t(x$pred), type = "l")
abline(v = 10)