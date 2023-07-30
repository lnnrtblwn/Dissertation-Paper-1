OLSdist = function(y,T0,lam1,lam2,flag) {
  
dimy=dim(y)
J=dimy[2]
T=dimy[1]
T1=T-T0
x=y[,2:J]
J=J-1
y0=y[,1]
y0m=demean(y0[1:T0])
xm=demean(x[1:T0,])
mx=x[1,]-xm[1,]
mat1=matrix(1,J,J)
vec1=matrix(1,J,1)
y_OLSdist_post=matrix(1,T1,1)

for (i in 1:T1 ) {
   dist=matrix(1,1,T0)
     x0=x[(T0+i),]
    if (flag==1) {
    for (j in 1:T0) {
       d=x0-x[j,]
       dist[j]=1/sqrt(sum(d^2))
     }
   }
   wmat=vec1%*%dist
   xdist=xm*t(wmat)
   ydist=y0m*t(dist)
   A = t(xdist)%*%xdist + lam1*diag(J) + lam2*mat1
   w0 = solve(A)%*%(t(xdist)%*%ydist +lam2*vec1)
   y_OLSdist_post[i] = mean(y0)+(x0-mx)%*%w0
}
return(y_OLSdist_post)
}