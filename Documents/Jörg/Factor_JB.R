

# Als (N+1)-Vektor w wird zunächst die Konstante und dann die N Gewichte für die N Donors ausgewiesen, 
# die man dann für das SC verwenden kann.

N=20
T0=20
r=2
b=matrix(c(1,0.5),2,1)

#  generate series: y=treated  x=donors
f=matrix(rnorm(T0*r),T0,r)
lam=matrix(rnorm(N*r),r,N)
y=f%*%b+rnorm(T)
x=f%*%lam+matrix(rnorm(N*T0),T0,N)

#  program for estimating the weights based on r factors and T0 prior observations
#  INPUT: T0 observations for Y (treatment) and X (donors)
#  OUTPUT:yest
V    = cov(x)
eig  = eigen(V)
w    = eig$vectors[,1:r]
fhat = x%*%w
yest = lm(y~fhat)$fitted.values
w    = lm(yest~x)$coefficients

