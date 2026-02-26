library(mvtnorm)

# Simulation issues with a non-invertible matrix
mat1 <- diag(2)
mat2 <- diag(c(0,0.25))
matsim <- kronecker(matrix(1,2,2), mat1) + kronecker(diag(1,2,2), mat2)
# checks - not invertible (duh)
chol(matsim)

rmvnorm(1,sigma=matsim)

eigen(matsim)
svd(matsim)

rmvnorm(1,sigma=matsim, method="chol")

rmvnorm(1,sigma=matsim, method="eigen")

rmvnorm(1,sigma=matsim, method="svd")

#################################################################
# Simulation times with rmvnorm
# Speed - svd faster than eigen, similar to chol
t1<-Sys.time()
x<-rmvnorm(10000,sigma=matsim, method="eigen")
t2<-Sys.time()
cat("Computation time, eigen:",t2-t1,"\n")

t1<-Sys.time()
x<-rmvnorm(10000,sigma=matsim, method="chol")
t2<-Sys.time()
cat("Computation time, eigen:",t2-t1,"\n")

t1<-Sys.time()
x<-rmvnorm(10000,sigma=matsim, method="svd")
t2<-Sys.time()
cat("Computation time, svd:",t2-t1,"\n")

#################################################################
# For 200 subjects, 4 random effects (2 IIV+2 occasions)
t1<-Sys.time()
for(i in 1:100) x<-rmvnorm(200,sigma=matsim, method="eigen")
t2<-Sys.time()
cat("Computation time, eigen:",t2-t1,"\n")

t1<-Sys.time()
for(i in 1:100) x<-rmvnorm(200,sigma=matsim, method="svd")
t2<-Sys.time()
cat("Computation time, svd:",t2-t1,"\n")

#################################################################
# Simulation times - comparing rmvtnorm with svd and first computing svd
# 3 times faster when using decomposition before
t1<-Sys.time()
x<-NULL
for(i in 1:100) x<-rbind(x,rmvnorm(200,sigma=matsim, method="svd"))
t2<-Sys.time()
cat("Computation time, svd:",t2-t1,"\n")
eta.rmv<-x

t1<-Sys.time()
cmat<-svd(matsim)
rmat <- t(cmat$v %*% (t(cmat$u) * sqrt(pmax(cmat$d, 0))))
x<-NULL
for(i in 1:100) x<-rbind(x,matrix(rnorm(800),ncol=4)%*%rmat)
t2<-Sys.time()
cat("Computation time, svd:",t2-t1,"\n")
eta.rno<-x

# Checking matrices - both similar to matsim
cov(eta.rmv)
cov(eta.rno)
matsim

#################################################################
# Simulation times - comparing svd and chol
# Comparing with an invertible matrix

# using rmvnorm directly
# 40 times slower :-(
mat3<-diag(c(1,0.09,0.25))
mat3[2,3]<-mat3[3,2]<-0.4*sqrt(mat3[2,2]*mat3[3,3])

t1<-Sys.time()
for(i in 1:1000) x<-rmvnorm(3,sigma=mat3, method="svd")
t2<-Sys.time()
cat("Computation time, svd:",t2-t1,"\n")

t1<-Sys.time()
cmat3<-chol(mat3)
for(i in 1:1000) x<-matrix(rnorm(9),ncol=3)%*%cmat3
t2<-Sys.time()
cat("Computation time, inverted using chol:",t2-t1,"\n")

# By hand inverting with svd before - 2.5 times slower
# probably because svd slower than chol
cmat3<-svd(mat3)
rmat <- t(cmat3$v %*% (t(cmat3$u) * sqrt(pmax(cmat3$d, 0))))
t1<-Sys.time()
for(i in 1:1000) x<-matrix(rnorm(9),ncol=3)%*%rmat
t2<-Sys.time()
cat("Computation time, inverted using svd by hand:",t2-t1,"\n")

#################################################################
# Chol versus svd (svd slower by 2, and eigen by a factor 17)
t1<-Sys.time()
for(i in 1:10000) cmat3<-chol(mat3)
t2<-Sys.time()
cat("Computation time for chol:",t2-t1,"\n")

t1<-Sys.time()
for(i in 1:10000) cmat3<-svd(mat3)
t2<-Sys.time()
cat("Computation time for svd:",t2-t1,"\n")

t1<-Sys.time()
for(i in 1:10000) cmat3<-eigen(mat3)
t2<-Sys.time()
cat("Computation time for eigen:",t2-t1,"\n")

#################################################################
# Applied to saemix



#################################################################
