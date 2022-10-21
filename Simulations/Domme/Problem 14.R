# Problem 14

qwilcox(p = 0.05, m = 10, n = 10) + (10 * (10 +1)) / 2
pwilcox(80:85 - 55, m = 10, n = 10)
# 82 is the critical value
# (see section 'Note' on the help page of wilcox.test)

createPlot = function(dist, reps = 1000) {
  if (dist == "norm") obs = replicate(reps, rnorm(n = 20))
  if (dist == "t2") obs = replicate(reps, rt(n = 20, df = 2))
  if (dist == "t3") obs = replicate(reps, rt(n = 20, df = 3))
  if (dist == "exp") obs = replicate(reps, rexp(n = 20))
  
  xbar = colMeans(obs[1:10, ])
  s2x = apply(X = obs[1:10, ], MAR = 2, FUN = var)
  ybar = colMeans(obs[11:20, ])
  s2y = apply(X = obs[11:20, ], MAR = 2, FUN = var)
  
  rejectionsW = numeric(21)
  rejectionsWapprox = numeric(21)
  rejectionsT = numeric(21)
  for (i in 1:21) {
    W_N = apply(X = obs, MAR = 2, function(z) sum(rank(z)[1:10]))
    rejectionsW[i] = mean(W_N <= 82)
    Wapprox = (W_N - 105) / sqrt(175)
    rejectionsWapprox[i] = mean(Wapprox <= qnorm(0.05))
    Te = sqrt(90) * (xbar - ybar) / sqrt(9 * s2x + 9 * s2y)
    rejectionsT[i] = mean(Te <= qt(0.05, df = 18))
    obs[11:20, ] = obs[11:20, ] + 0.1
    ybar = ybar + 0.1
  }
  plot(seq(0, 2, by = 0.1), rejectionsW, type = "l", ylim = c(0, 1), xlim = c(0, 2),
       xlab = expression(delta), ylab = "rejection frequency")
  points(seq(0, 2, by = 0.1), rejectionsWapprox, type = "l", lty = 2, col = "red")
  points(seq(0, 2, by = 0.1), rejectionsT, type = "l", lty = 2, col = "blue")
  legend("bottomright", lty = c(1, 2, 2), col = c("black", "red", "blue"),
         legend = c("Wilcoxon test", "approx. Wilcoxon test", "t test"))
  return(rbind(rejectionsW, rejectionsWapprox, rejectionsT))
}


createPlot("norm")
createPlot("t2")
createPlot("t3")
createPlot("exp")

