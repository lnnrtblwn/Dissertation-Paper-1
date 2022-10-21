photon.data = read.csv("C:/Users/lbolwin/OneDrive - IW/Desktop/Diss/Kurse/Stochastic Processes/PS/R/photonen.txt", sep="")
  
photon = photon.data[, 1] / 60 # in minutes 

N_t = length(photon)

plot(c(0, photon), 0:N_t, type = "s", xlab = "time", ylab = "number of photons",
     main = "Poisson process of photons")
# type = "s" stands for step function (not visible here due to large data set)

# How to estimate the intensity parameter:
# The waiting times are i.i.d. Exp(lambda)-distributed.
# So, we have a sample of size N_t of an Exp(lambda)-distribution.
# The Maximum-Likelihood estimator for this is 1 divided by the arithmetic mean.
# The sum of the waiting times is the last arrival time photon[N_t] and the
# sample size is N_t, so the Maximum-Likelihood estimator for lambda is:
lambda = N_t / photon[N_t]
lambda


# N_t is Poi(lambda * t)-distributed, the expected value of N_t is lambda * t.
# Add the expected value of N_t to the plot.
lines(x = c(0, photon[N_t]), y = c(0, N_t), lty = 2)
abline(coef = c(0, lambda), lty = 2, xlim = c(0, photon[N_t]))


# Take a look at a small part of the process, e.g. until the fifth arrival.
arrivals = 20
plot(c(0, photon[1:arrivals]), 0:arrivals, xlab = "time", 
     ylab = "number of photons", 
     main = paste("Poisson process of photons up to arrival", arrivals))
plot(c(0, photon[1:arrivals]), 0:arrivals, xlab = "time", 
     ylab = "number of photons", type = "s",
     main = paste("Poisson process of photons up to arrival", arrivals))
# Here you can see the step function .
# (Actually, the vertical lines do not belong to the function.)


# Write own function to plot a step function:
my.steps = function(x, y, xlab, ylab, main) {
  plot(x = x, y = y, xlab = xlab, ylab = ylab, type = "n", main = main)
  for (i in 1:(length(x) - 1)) {
    lines(x = c(x[i], x[i + 1]), y = rep(i - 1, 2), type = "l", lty = 1)
    lines(x = rep(x[i + 1], 2), y = c(i - 1, i), type = "l", lty = 2)
  }
}

my.steps(c(0, photon[1:arrivals]), 0:arrivals, xlab = "time", 
         ylab = "number of photons", 
         main = paste("Poisson process of photons up to arrival", arrivals))
