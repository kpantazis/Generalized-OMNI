#' This code reproduces Figure 3 from the paper "The Importance of Being Correlated:
#' Implications of Dependence in Joint Spectral Inference across Multiple Networks".
#' This code plots the estimated differences between SBMs for Classical OMNI matrix
#' For Weighted OMNI example (w = c(1,1,10)) simply uncomment all once.
#' Requires ase and make_omni functions
#' 
#' 
#' @references The Importance of Being Correlated:Implications of
#'  Dependence in Joint Spectral Inference across Multiple Networks
#'
#' @author Konstantinos Pantazis <kpantazi@umd.edu>
#' 
set.seed(1234)
## Number of vertices
n <- 300

## Number of graphs
m <- 3

## Block probability matrix
L <- matrix(c(0.7,0.2,0.2,0.7),nrow=2)

## ASE with embedding dimension d=2
d <- 2
X <- ase(L,d)$X

## Differences between (scaled) estimates
# diffs2 <- list()
diffs <- list()
for (i in 1:(choose(m,2))){
  # diffs2[[i]] <- matrix(0,n,d)
  diffs[[i]] <- matrix(0,n,d)
}

## nMC=100 replicates
for(l in 1:100){
  ## Sampling SBMs with K=2 communities, probability matrix L and probability vector (1/2,1/2)
  A<-list()
  g1<-sample_sbm(n,L,c(floor(n/2),ceiling(n/2)))
  A[[1]] <- g1[]
  g2<-sample_sbm(n,L,c(floor(n/2),ceiling(n/2)))
  A[[2]] <- g2[]
  g3<-sample_sbm(n,L,c(floor(n/2),ceiling(n/2)))
  A[[3]] <- g3[]

  ## Weight vector // The particular one constructs the Classical OMNI matrix
  # w2<-c(1,1,10)
  w <- c(1,1,1)
  
  
  ## Weight Omnibus Matrix
  # M2<-make_Omni(A,w2)
  M <- make_Omni(A,w)
  
  
  
  ## OMNI embedding
  # X2omni <- ase(M2,d)$X
  Xomni <- ase(M,d)$X
  
  
  ## Unobseved latent positions (not aligned)
  Z<-rbind(t(outer(X[1,],rep(1,floor(n/2)))),t(outer(X[2,],rep(1,ceiling(n/2)))))
  Z <- rbind(Z,Z,Z)
  
  
  ## Procrustes alignment
  # W2 <- procrustes(X2omni,Z)
  # X2omni<-X2omni%*%(W2$W2)
  W<-procrustes(Xomni,Z)
  Xomni<-Xomni%*%(W$W)

  ## Scaled differences of estimated latent positions of vertex 1
  # diffs2[[1]][l,]<-sqrt(n)*(X2omni[1,]-X2omni[301,])
  # diffs2[[2]][l,]<-sqrt(n)*(X2omni[1,]-X2omni[601,])
  # diffs2[[3]][l,]<-sqrt(n)*(X2omni[301,]-X2omni[601,])
  diffs[[1]][l,]<-sqrt(n)*(Xomni[1,]-Xomni[301,])
  diffs[[2]][l,]<-sqrt(n)*(Xomni[1,]-Xomni[601,])
  diffs[[3]][l,]<-sqrt(n)*(Xomni[301,]-Xomni[601,])

  #print(l)
}

#############################################################
#
#############################################################
#
#############################################################
## Computation of the covariance matrix
x1<-X[1,]
x2<-X[2,]
t11<-x1%*%t(x1)
t12<-x1%*%t(x2)
t21<-x2%*%t(x1)
t22<-x2%*%t(x2)
p11<-0.7
p12<-0.2
p22<-0.7
del<-0.5*t11+0.5*t22
sig<- 0.5*(p11-p11^2)*t11 +0.5*(p12-p12^2)*t22  
sig<-ginv(del)%*%sig%*%ginv(del)

## Coefficients of the covariance matrices among pairs of estimated latent positions
  #Classical OMNI example
  s12<- 1/2
  s13<-1/2
  s23<-1/2

## For general Weighted OMNI matrices call the function "coeff_weighted"
## from coeff_cov_Weighted_OMNI.R file

## Plots
par(mfrow=c(2,2),mai = c(0.3, 0.3, 0.2, 0.2))
plot(diffs[[1]],main = "12", col='blue',xlim=c(-3,3),ylim=c(-3,3),xlab='',ylab='')
# points(diffs2[[1]],col='red')
ellipse(c(0,0), s12*sig, draw=TRUE, newplot = FALSE,col='blue',cex=0.1)
legend("topleft", legend=c("classical omni"),
       col=c( "blue"), lty=c(1,1), cex=0.7,bty="n")


plot(diffs[[2]],main = "13",col='blue',xlim=c(-3,3),ylim=c(-3,3),xlab='',ylab='')
# points(diffs2[[2]],col='red')
ellipse(c(0,0), s13*sig, draw=TRUE, newplot = FALSE,col='blue',cex=0.1)
legend("topleft", legend=c("classical omni"),
       col=c( "blue"), lty=c(1,1), cex=0.7,bty="n")


plot(diffs[[3]],main = "23",col='blue',xlim=c(-3,3),ylim=c(-3,3),xlab='',ylab='')
# points(diffs2[[3]],col='red')
ellipse(c(0,0), s23*sig, draw=TRUE, newplot = FALSE,col='blue',cex=0.1)
legend("topleft", legend=c("classical omni"),
       col=c( "blue"), lty=c(1,1), cex=0.7,bty="n")

