

I = diag(ncol(x))
U = matrix(data = 1,  nrow = ncol(x) + 0, ncol = 1)

x = x_pre[1:(nrow(x_pre) * (CV_share)), ] %>%
  scale(center = FALSE, scale = FALSE) %>%
  as.data.frame()

y = y_pre[1:(length(y_pre) * (CV_share))] %>%
  scale(center = FALSE, scale = FALSE) %>%
  as.data.frame()

Y = y %>%
  as.matrix()

X = x %>%
  mutate(x0 = 1) %>%
  dplyr::select(x0, everything()) %>%
  as.matrix()

l1 = l2 = 1

solve(t(X) %*% X + l1 * I + l2 * (U %*% t(U))) %*% ((t(X) %*% Y) + l2 * U)
summary(lm(Y ~ -1+X))$coefficients

my = mean(Y)
y0 = Y - my
mx = colMeans(X)
X0 = sweep(X, 2, colMeans(X), `-`)[,-1]
#X0 = X[,-1]

beta = summary(lm(y0 ~ -1+X0))$coefficients[,1]

pred = my + X0 %*% beta

plot(pred, Y)

# now without demeaning and with constant

beta1 = summary(lm(Y ~ -1 + X))$coefficients[,1] 

pred1 = X %*% beta1

plot(pred, pred1)

Y = y_pre[1:(length(y_pre)*(CV_share))]
X = x_pre[1:(nrow(x_pre)*(CV_share)),]
