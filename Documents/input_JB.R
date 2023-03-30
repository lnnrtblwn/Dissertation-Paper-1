# Beispiel von JB

my = matrix(c(1,1,1), byrow = T, nrow= 3)
SIG = matrix(c(1, 0.1, 0.4, 0.1, 1, 0.5, 0.4, 0.5, 1), byrow = T, nrow = 3)

r = chol(SIG)

t(r) %*% r == SIG # only R'R = SIG, not RR'. Why?

# R
Q = solve(r)
q = matrix(Q[1,], nrow = 1)

t(q) %*% my # dimension error. Error in notation?
q %*% my
# t(q) %*% t(my)

w = c(-1*q[,2]/q[,1], -1*q[,3]/q[,1]) # wrong result

# t(R)
Q = solve(t(r))
q = matrix(Q[1,], nrow = 1)

t(q) %*% my # dimension error. Error in notation?
q %*% my
# t(q) %*% t(my)

w = c(-q[,2]/q[,1], -q[,3]/q[,2]) # wrong result
