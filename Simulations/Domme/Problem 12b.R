# Problem 12 b)

my.wilcox.test = function(x, y) {
  W_N = sum(rank(c(x, y))[1:3])
  if (W_N <= 6 | W_N >= 21) return (1)
  return (0)
}

sim = function(reps = 10000, s2 = 1) {
  testresult = numeric(reps)
  for (i in 1:reps) {
    testresult[i] = my.wilcox.test(x = rnorm(3, sd = sqrt(s2)), y = rnorm(5))
  }
  mean(testresult)
}

# easier
sim2 = function(reps = 10000, s2 = 1) {
  x = replicate(n = reps, expr = rnorm(3, sd = sqrt(s2)))
  y = replicate(n = reps, expr = rnorm(5))
  
  W_N = apply(X = rbind(x, y), MAR = 2, FUN = function(z) sum(rank(z)[1:3]))
  testresult = W_N <= 6 | W_N >= 21
  mean(testresult)
}

sim2(s2 = 1)
sim2(s2 = 5)
sim2(s2 = 10)
sim2(s2 = 100)
sim2(s2 = 1000)
sim2(s2 = 10000)
sim2(s2 = 100000)
## Considerable size distortions