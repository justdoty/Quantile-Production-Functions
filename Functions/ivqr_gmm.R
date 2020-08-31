#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#https://faculty.missouri.edu/~kaplandm/
source('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Functions/gmmq.R')
# source('gmmq.R')
IDENTITY.FLAG <- FALSE
ZZ.FLAG <- FALSE

ivqr.gmm <- function(tau, mX, mlX, mZ, fitphi, fitlagphi, h=0, max.time, upper, lower, 
                     structure=c('iid','ts'), LRV.kernel=c('QS','Bartlett','uniform'), 
                     Lambda=function(theta, mX, mlX, fitphi, fitlagphi, tau){fitphi-mX%*%theta[1:(ncol(mX))]-theta[length(theta)]*(fitlagphi-mlX%*%theta[1:(ncol(mX))])}, 
                     theta.init=0 ) {
    
    n <- dim(mX)[1]
    dZ <- ncol(mZ) #number of instruments
    dX <- ncol(mX) #number of regressors

    Itilde <- Itilde.KS17
    Itilde.deriv <- Itilde.deriv.KS17
    
    # Add bandwidth as second argument to G() and G'() functions.
    Gfn <- function(v,h){      
      Itilde.KS17(v/h)    
    }
    Gpfn <- function(v,h){      
      Itilde.deriv.KS17(v/h)    
    }
    
    if(length(theta.init)!=length(lower)){stop("lower bound doesn't have right dimension.")}
    if(length(theta.init)!=length(upper)){stop("upper bound doesn't have right dimension.")}
    
    
    LRV.hat <- function(theta, h) {LRV.est.fn(tau=tau, mX=mX, mlX=mlX, mZ=mZ, fitphi=fitphi, fitlagphi=fitlagphi, 
                             Lambda=Lambda, theta, Itilde=Itilde.KS17, h, 
                              structure=structure, LRV.kernel=LRV.kernel)  } 
    ivqr.obj <- function(theta, h) {
      if (IDENTITY.FLAG==TRUE){
        W.hat <- diag(dZ)  
        } else if (ZZ.FLAG==TRUE) {
          W.hat <- solve(t(mZ)%*%mZ/n)
        } else {  
          W.hat <- solve(LRV.hat(theta, h))  
        }
        L <- matrix(Lambda(theta, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau), ncol=1)
        gni <- mZ*repmat((Gfn(-L,h)-tau),1,dZ) 
        g.bar <- as.matrix(colMeans(mZ*repmat((Gfn(-L,h)-tau),1,dZ))) 
        obj.fn <- t(g.bar)%*% W.hat %*%g.bar
        return(obj.fn)
    }
    obj.fn <- function(theta) {
      return(n*ivqr.obj(theta,  h=h))
    }
 
     
     M.hat <- function(theta) {
       L <- matrix(Lambda(theta, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau), ncol=1)
       gni <- mZ*repmat((Gfn(-L,h)-tau),1,dZ) 
       g.bar <- as.matrix(colMeans(mZ*repmat((Gfn(-L,h)-tau),1,dZ))) 
       return(g.bar)        
     }
     
     
     ivqr.gmm1 <- GenSA(par=theta.init, fn=obj.fn, lower=lower, upper=upper, control=list(max.time=max.time))
     theta.hat <- ivqr.gmm1$par
     
     J.stat <- obj.fn(theta.hat)
 
     return(list(theta=theta.hat, obj.fn=obj.fn, h=h, M.hat=M.hat, J.stat=J.stat))
}