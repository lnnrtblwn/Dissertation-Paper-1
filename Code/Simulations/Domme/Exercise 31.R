## Exercise 31

## a

n <- 1000
increments <- rnorm(n,0,1/sqrt(n))
bmotion <- cumsum(increments)
#bbridge <- bmotion-(1:n)/n*bmotion[n]
plot(bmotion,type="l",xaxt="n")
#plot(bbridge,type="l",xaxt="n")
axis(1,c(1,n),c(0,1))

## b

m <- 1000
n <- 1000
counter <- rep(0,m)
for(i in 1:m)
{
increments <- rnorm(n,0,1/sqrt(n))
bmotion <- cumsum(increments)
bbridge <- bmotion-(1:n)/n*bmotion[n]
counter[i] <- max(abs(bbridge))
}
sort(counter)[0.95*m]

## Implemented distribution function in R:
x <- 1.358
distributionf <- .Call(stats:::C_pKS2, p = x, tol = 10^(-6))