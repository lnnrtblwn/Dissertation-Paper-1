# Install devtools if noy already installed
install.packages("devtools", repos='http://cran.us.r-project.org')

# Install augsynth from github
devtools::install_github("ebenmichael/augsynth")

## augsynth: The Augmented Synthetic Control Method
# https://github.com/ebenmichael/augsynth/blob/master/vignettes/singlesynth-vignette.md

library(magrittr)
library(dplyr)
library(augsynth)
data(kansas)

kansas %>% 
  select(year, qtr, year_qtr, state, treated, gdp, lngdpcapita) %>% 
  filter(state == "Kansas" & year_qtr >= 2012 & year_qtr < 2013) 

syn <- augsynth(lngdpcapita ~ treated, fips, year_qtr, kansas,
                progfunc = "None", scm = T)

summary(syn)

summary(syn, stat_func = function(x) -sum(x))

summary(syn, stat_func = function(x) abs(sum(x)))

plot(syn)
plot(syn, inf_type = "jackknife+")

# Augmenting synth with an outcome model

asyn <- augsynth(lngdpcapita ~ treated, fips, year_qtr, kansas,
                 progfunc = "Ridge", scm = T)

plot(asyn, cv = T)

summary(asyn)
plot(asyn)

# Adding covariates

covsyn <- augsynth(lngdpcapita ~ treated | lngdpcapita + log(revstatecapita) +
                     log(revlocalcapita) + log(avgwklywagecapita) +
                     estabscapita + emplvlcapita,
                   fips, year_qtr, kansas,
                   progfunc = "ridge", scm = T)

summary(covsyn)

plot(covsyn)

# Covariates plus residual adjustment

covsyn_resid <- augsynth(lngdpcapita ~ treated | lngdpcapita + log(revstatecapita) +
                           log(revlocalcapita) + log(avgwklywagecapita) +
                           estabscapita + emplvlcapita,
                         fips, year_qtr, kansas,
                         progfunc = "ridge", scm = T, lambda = asyn$lambda,
                         residualize = T)
summary(covsyn_resid)

plot(covsyn_resid)

# Augment synth with different outcome models

desyn <- augsynth(lngdpcapita ~ treated,
                  fips, year_qtr, kansas,
                  progfunc = "none", scm = T,
                  fixedeff = T)

summary(desyn)

plot(desyn)

# We can incorporate other outcome models by changing the progfunc. Several outcome models are available, including, 
# fitting the factor model directly with gsynth, general elastic net regression, bayesian structural time series estimation 
# with CausalImpact, and matrix completion with MCPanel. For each outcome model you can supply an optional set of parameters, 
# see documentation for details.
