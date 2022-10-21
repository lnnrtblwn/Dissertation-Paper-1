# Problem 23

dax = read.csv("C:/Users/lbolwin/OneDrive - IW/Desktop/Diss/Kurse/Stochastic Processes/PS/R/DAX.csv")
dax$Date = as.Date(dax$Date, "%Y-%m-%d")

n = dim(dax)[1]
price = dax$Close[n:1]

# a)

changes.abs = price[2:n] - price[1:(n - 1)] # absolute daily changes
changes.abs = diff(price) # shorter

changes.rel = price[2:n] / price[1:(n - 1)] # relative daily changes


plot(dax$Date[n:2], changes.abs, type = "l", xlab = "date", ylab = "absolute daily change")

plot(dax$Date[n:2], changes.rel, type = "l", xlab = "date", ylab = "relative daily change")

# Both types of changes seem independent but not stationary (volatility cluster)
# You can also check the autocorrelation function of the series 
# (cf. a lecture covering time series analysis)
acf(changes.abs)
acf(changes.rel)

# b)

mean.abs = rep(NA, n - 1)
var.abs = rep(NA, n - 1)

mean.rel = rep(NA, n - 1)
var.rel = rep(NA, n - 1)

for (i in 30:(n - 1)) {
  mean.abs[i] = mean(changes.abs[(i - 29):i])
  var.abs[i] = var(changes.abs[(i - 29):i])
  
  mean.rel[i] = mean(changes.rel[(i - 29):i])
  var.rel[i] = var(changes.rel[(i - 29):i])
}

plot(dax$Date[n:2], mean.abs, type = "l", xlab = "date", ylab = "mean")
plot(dax$Date[n:2], var.abs, type = "l", xlab = "date", ylab = "variance")

plot(dax$Date[n:2], mean.rel, type = "l", xlab = "date", ylab = "mean")
plot(dax$Date[n:2], var.rel, type = "l", xlab = "date", ylab = "variance")

setwd("C:/Users/lbolwin/OneDrive - IW/Desktop/Diss/Kurse/Stochastic Processes/Lecture/R")

pdf("absoluteDailyChanges.pdf", width = 10)
par(mfrow = c(3, 1), las = 1, mar = c(5.1, 5.1, 4.1, 2.1))
plot(dax$Date[n:2], changes.abs, type = "l", xlab = "date", ylab = "")
mtext("absolute daily change", side = 2, line = 4, las = 0, cex = 0.65)
plot(dax$Date[n:2], mean.abs, type = "l", xlab = "date", ylab = "")
mtext("mean", side = 2, line = 4, las = 0, cex = 0.65)
plot(dax$Date[n:2], var.abs, type = "l", xlab = "date", ylab = "")
mtext("variance", side = 2, line = 4, las = 0, cex = 0.65)
dev.off()

pdf("relativeDailyChanges.pdf", width = 10)
par(mfrow = c(3, 1), las = 1, mar = c(5.1, 5.1, 4.1, 2.1))
plot(dax$Date[n:2], changes.rel, type = "l", xlab = "date", ylab = "")
mtext("relative daily change", side = 2, line = 4, las = 0, cex = 0.65)
plot(dax$Date[n:2], mean.rel, type = "l", xlab = "date", ylab = "")
mtext("mean", side = 2, line = 4, las = 0, cex = 0.65)
plot(dax$Date[n:2], var.rel, type = "l", xlab = "date", ylab = "")
mtext("variance", side = 2, line = 4, las = 0, cex = 0.65)
dev.off()
