# source('PFQR/FUN/gmmq_aux.R')
# source('PFQR/FUN/QLP.R')
# source('PFQR/FUN/gmmq.R')
# source('PFQR/FUN/1step.R')
twostep <- function(mY, mX, mlX, mZ, fitphi, fitlagphi, W, h, tau, b.init, gbartrue, upper, lower, maxtime){
	n <- nrow(mY)
	ivqr.obj <- function(b, h, weight.mtx) {
    	L <- matrix(Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau), ncol=1)
        gni <- mZ*repmat((Gfn(-L, h)-tau),1,ncol(mZ)) 
        g.bar <- as.matrix(colMeans(mZ*repmat((Gfn(-L,h)-tau),1,ncol(mZ))))-gbartrue 
        obj.fn <- t(g.bar)%*%weight.mtx%*%g.bar
        return(obj.fn)
    }
    obj.fn <- function(b){
    	return(n*ivqr.obj(b=b,  h=h, weight.mtx=W))
    } 

    G.hat <- function(b) {
       L <- matrix(Lambda(b=b, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau), ncol=1)
       gni <- mZ*repmat((Gfn(-L, h)-tau),1, ncol(mZ))-gbartrue 
       g.bar <- as.matrix(colMeans(mZ*repmat((Gfn(-L, h)-tau),1, ncol(mZ)))) 
       return(g.bar)        
   }

   ivqr.gmm1 <- GenSA(par=b.init, fn=obj.fn, lower=lower, upper=upper, control=list(max.time=maxtime))
   b <- ivqr.gmm1$par
   J <- obj.fn(b=b)
   G <- G.hat(b=b)
   return(list(b=b, h=h, W=W, J=J, G=G))
}