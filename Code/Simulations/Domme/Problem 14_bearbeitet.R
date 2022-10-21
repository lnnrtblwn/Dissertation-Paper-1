# Problem 14

x = 14
y = 14

E_W = (x*(x+y+1))/2
V_W = (x*y*(x+y+1))/12

t_1 = (x*y)/(x+y)/ (1/(x+y-2))
cv_wx = (qwilcox(p = 0.05, m = x, n = y) + (x * (x +1)) / 2)-1

qwilcox(p = 0.05, m = 10, n = 10) + (10 * (10 +1)) / 2
pwilcox(80:85 - 55, m = 10, n = 10)
# 82 is the critical value
# (see section 'Note' on the help page of wilcox.test)

my.createPlot = function(dist, reps = 100) {
  if (dist == "norm") obs = replicate(reps, rnorm(n = x+y))
  if (dist == "t2") obs = replicate(reps, rt(n = x+y, df = 2))
  if (dist == "t3") obs = replicate(reps, rt(n = x+y, df = 3))
  if (dist == "exp") obs = replicate(reps, rexp(n = x+y))
  
  xbar = colMeans(obs[1:(x+y)/2, ])
  s2x = apply(X = obs[1:(x+y)/2, ], MAR = 2, FUN = var)
  ybar = colMeans(obs[(1+(x+y)/2):(x+y), ])
  s2y = apply(X = obs[(1+(x+y)/2):(x+y), ], MAR = 2, FUN = var)
  
  rejectionsW = numeric((x+y)+1)
  rejectionsWapprox = numeric((x+y)+1)
  rejectionsT = numeric((x+y)+1)
  for (i in 1:(x+y)+1) {
    W_N = apply(X = obs, MAR = 2, function(z) sum(rank(z)[1:x]))
    rejectionsW[i] = mean(W_N <= cv_wx)
    Wapprox = (W_N - E_W) / sqrt(V_W)
    rejectionsWapprox[i] = mean(Wapprox <= qnorm(0.05))
    Te = sqrt(t_1) * (xbar - ybar) / sqrt((x-1) * s2x + (y-1) * s2y)
    rejectionsT[i] = mean(Te <= qt(0.05, df = 18))
    obs[(x+1):(2*x), ] = obs[(x+1):(2*x), ] + 0.1
    ybar = ybar + 0.1
  }
  plot(seq(0, 2, by = 1/x), rejectionsW, type = "l", ylim = c(0, 1), xlim = c(0, 2),
       xlab = expression(delta), ylab = "rejection frequency")
  points(seq(0, 2, by = 1/x), rejectionsWapprox, type = "l", lty = 2, col = "red")
  points(seq(0, 2, by = 1/x), rejectionsT, type = "l", lty = 2, col = "blue")
  legend("bottomright", lty = c(1, 2, 2), col = c("black", "red", "blue"),
         legend = c("Wilcoxon test", "approx. Wilcoxon test", "t test"))
  return(rbind(rejectionsW, rejectionsWapprox, rejectionsT))
}


my.createPlot("norm")
my.createPlot("t2")
my.createPlot("t3")
my.createPlot("exp")

