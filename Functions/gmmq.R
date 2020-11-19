# source('PFQR/FUN/gmmq_aux.R')
# source('PFQR/FUN/QLP.R')
gmmq <- function(mY, mX, mlX, mZ, fitphi, fitlagphi, h, tau, b.init, VERBOSE){
	n <- nrow(mY)
	obj.fn <- function(b, h){ 
	    L <- Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau)
	    return(colMeans((mZ*array(data=Itilde(-L/h)-tau,dim=dim(mZ)))))
  	}
    x <- seq(0,1,by=.001)
    y <- sapply(x, function(b) obj.fn(b,h=.1))
    save(x,y,file='PFQR/plot.Rdata')
  	jac.fn <- function(b, h) { 
      L <- Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau)
      Lp <- Lambda.derivative(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau)
      return((t(mZ)%*%(Lp*array(data=Itilde.deriv(-L/h),dim=dim(mZ))*(-1/h)))/n)
    }
    #Compute estimator: if h=0, find smallest h newtonsys can use; else use user's h (if big enough)
	 if (h==0) {
	   h.init <- h.cur <- 0.001
	   hfac <- 100
	 } else {
	   h.init <- h.cur <- h
	   hfac <- 10
	 }
	 h.lo <- h
   	 h.hi <- NA
  	 MAXITER <- 400;  H.MAX <- 1e10
     if (VERBOSE==TRUE) {cat(sprintf('n=%d, tau=%g, dB=%d, desired h=%g\n',n,tau,length(b.init),h))}
     last.error <- NULL
  	 while (h.cur<=H.MAX && h.cur>.Machine$double.eps*1e2 && (is.na(h.hi) || (h.hi/h.lo)>1.4)){
      print(h.cur)
  	 	if (VERBOSE==TRUE) {cat(sprintf("h.cur=%g: ",h.cur))}
  	 	J <- function(b) jac.fn(b, h.cur)
      print(b.init)
  	 	soln <- tryCatch(GenSA(par=b.init, fn=function(b) obj.fn(b, h.cur), lower=0, upper=1, control=list(max.time=5, maxit=MAXITER)), 
  	 		warning=function(w) NA, error=function(e) {list(NA,e)})
      print(soln)
  	 	exitOK <- !(is.na(soln$par) || soln$counts==MAXITER)
  	 	if (!exitOK) {
  	  		last.error <- soln[2]
  	 		if (VERBOSE) {
  	    		if (is.na(soln[1])) cat(sprintf("no solution, soln=NA\n")) else cat(sprintf("no solution, soln$niter=%d\n",soln$niter))
  	  		}
  	  		h.lo <- h.cur
  	  		if (is.na(h.hi)){ #keep going up till find h.cur that actually works
  	   	    	h.cur <- h.cur*hfac
  	  		} else { #try arithmetic [not geometric] mean of h.lo and h.hi
  	   			h.cur <- (h.lo+h.hi)/2 
  	  		}	
  		} else { #newtonsys found a solution
  	  		if (VERBOSE) cat(sprintf("Solution found!\n"))
  	  			h.hi <- h.cur; b.best <- soln$zero; h.lo <- h
  	  		if (h.cur==h){ #found solution w/ user-requested h
  	    		break
		  	  } else if (h.hi/h.init>2.5 || h==0) { 
		  	    h.cur <- (h.lo+h.hi)*(2/3)
		  	  } else h.cur <- h.init
  		} 
  	}
  	if (h.cur>H.MAX){
    warning("(above) error from pracma::newtonsys")
    stop(last.error)
  	}
  	return(list(b=b.best,h=h.hi))
}	