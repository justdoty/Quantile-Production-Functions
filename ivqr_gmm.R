setwd('/Users/justindoty/Documents/Research/Structural_Estimation/Production/Heterogeneity_in_Firms/R_Code')
#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#https://faculty.missouri.edu/~kaplandm/

IDENTITY.FLAG <- TRUE
ZZ.FLAG <- FALSE

source("gmmq.R")

ivqr.gmm <- function(tau, X, Z, h=0, dB, max.time, upper, lower, 
                     structure=c('iid','ts'), LRV.kernel=c('QS','Bartlett','uniform'), 
                     Lambda=function(x,b){x[,1]-x[,-1]%*%b }, Lambda.derivative=function(x,b){-x[,-1]}, 
                     b.init=0 ) {
    
    n <- dim(X)[1]
    dZ <- ncol(Z) #number of instruments
    dX <- ncol(X[,-1]) #number of regressors

    Itilde <- Itilde.KS17
    Itilde.deriv <- Itilde.deriv.KS17
    
    # Add bandwidth as second argument to G() and G'() functions.
    Gfn <- function(v,h){      
      Itilde.KS17(v/h)    
    }
    Gpfn <- function(v,h){      
      Itilde.deriv.KS17(v/h)    
    }
    
    if(length(b.init)!=length(lower)){stop("lower bound doesn't have right dimension.")}
    if(length(b.init)!=length(upper)){stop("upper bound doesn't have right dimension.")}
    
    
    LRV.hat <- LRV.est.fn(tau=tau, X=X, Z=Z, 
                          Lambda=Lambda, beta.hat=b.init, Itilde=Itilde.KS17, h=h, 
                          structure=structure, LRV.kernel=LRV.kernel) 
    if (IDENTITY.FLAG) {
      W.hat <- diag(dZ) 
    } else if (ZZ.FLAG) {
      W.hat <- solve(t(Z)%*%Z/n)
    } else {   
      W.hat <- solve(LRV.hat) 
    }
      ivqr.obj <- function(b, h, W.hat) {
          L <- matrix(Lambda(x=X, b=b), ncol=1)

          gni <- Z*repmat((Gfn(-L,h)-tau),1,dZ) 
          g.bar <- as.matrix(colMeans(Z*repmat((Gfn(-L,h)-tau),1,dZ))) 

          obj.fn <- t(g.bar)%*% W.hat %*%g.bar
          return(obj.fn)
      }

      obj.fn <- function(b) n*ivqr.obj(b=b,  h=h, W.hat=W.hat)
 
     
       M.hat <- function(b) {
       L <- matrix(Lambda(x=X, b=b), ncol=1)
       gni <- Z*repmat((Gfn(-L,h)-tau),1,dZ) 
       g.bar <- as.matrix(colMeans(Z*repmat((Gfn(-L,h)-tau),1,dZ))) 
       return(g.bar)        
       }
     
     
     ivqr.gmm1 <- GenSA(par=b.init, fn=obj.fn,
                        lower=lower, upper=upper,
                  control=list(max.time=max.time))
     b <- ivqr.gmm1$par
     
     J.stat <- obj.fn(b)
 
     return(list(b=b, obj.fn=obj.fn, h=h, W.hat=W.hat, M.hat=M.hat, J.stat=J.stat))
  }


#EOF