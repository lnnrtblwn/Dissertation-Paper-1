# a)  NEUES PROGRAMM zum ausprobieren:

lam1 = 0
lam2 = 0
T = T0 + T1

w0 = matrix(1 / J, J, 1)
y0 = y_pre[(p_uni + 1):T0]

z = matrix(0, T0 - p_uni, J)
zhat = matrix(0, T - p_uni, J)
for (k in (1:J)) {
  lagx = y[(p_uni + 1):T, k + 1]
  for (i in (1:(p_uni - 1))) {
    lagx = cbind(lagx, y[(p_uni + 1 - i):(T - i), k + 1])
  }
  outz = lm(y0 ~ lagx[1:(T0 - p_uni), ])
  z[, k] = outz$fitted.values
  ahat = as.matrix(outz$coefficients)
  zhat[, k] = ahat[1] + lagx %*% ahat[2:(p_uni + 1), 1]
}
y0m = demean(y0)
zhat = zhat - mean(y0)
z = demean(z)
mat1 = matrix(1, J, J)
A = t(z) %*% z + lam1 * diag(J) + lam2 * mat1
w0 = solve(A) %*% (t(z) %*% y0m + lam2 * matrix(1, J, 1))
xb = z %*% w0
y_unidyn_pre = mean(y0) + xb

#  computing ex-post estimates of Y0
yhat = mean(y0) + zhat %*% w0
y_unidyn_post = yhat[(T0 - p_uni + 1):(T - p_uni)]

# b)  Korrigiertes Programm für UNIDYN:
p_uni = 2
lam1 = 1
lam2 = 1

w0 = matrix(1 / J, J, 1)
S0 = 100000000000000000000
S1 = 90000000000000000000
y0 = y_pre[(p_uni + 1):T0]

iter = 0
while ((S0 - S1) > 0.00001) {
  z = x_pre %*% w0
  lagz = z[(p_uni + 1):T0]
  for (i in (1:p_uni)) {
    lagz = cbind(lagz, z[(p_uni + 1 - i):(T0 - i)])
  }
  outreg = lm(y0 ~ lagz)
  alphas = outreg$coefficients[2:(p_uni + 2)]
  
  x1 = alphas[1] * x_pre[(p_uni + 1):T0, ]
  for (i in (1:p_uni - 1)) {
    # hier war ein kleiner Fehler
    x1 = x1 + alphas[i + 1] * x_pre[(p_uni + 1 - i):(T0 - i), ]
  }
  outreg = lm(y0 ~ x1)
  y0m = demean(y0)
  x1m = demean(x1)
  A = t(x1m) %*% x1m + lam1 * diag(J) + lam2 * matrix(1, J, J)                             #  hier hatte ich nicht mit den mittelwertbereinigten Größen gearbeitet.
  w0 = solve(A) %*% (t(x1m) %*% y0m + lam2 * matrix(1, J, 1))
  xb = x1m %*% w0
  y0hat = mean(y0) + xb
  S0 = S1
  S1 = sd(y0 - y0hat)
  iter = iter + 1
  if (iter == 20) {
    S0 = S1
  }
}