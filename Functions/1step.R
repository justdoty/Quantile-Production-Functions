source('PFQR/FUN/gmmq_aux.R')
source('PFQR/FUN/QLP.R')
source('PFQR/FUN/gmmq.R')
onestep <- function(mY, mX, mlX, mZ, fitphi, fitlagphi, h, tau, b.init, gbartrue, VERBOSE){
	if (h==0){
		ivqr.est <- gmmq(mY=mY, mX=mX, mlX, mZ=mZ[1:length(b.init)], fitphi=fitphi, fitlagphi=fitlagphi, h=0, tau=tau, b.init=b.init, VERBOSE=VERBOSE)
	} else {
		ivqr.est <- gmmq(mY=mY, mX=mX, mlX, mZ=mZ[1:length(b.init)], fitphi=fitphi, fitlagphi=fitlagphi, h=h, tau=tau, b.init=b.init, VERBOSE=VERBOSE)
	}
  	b.gmmq <- ivqr.est$b
  	h.gmmq <- ivqr.est$h
  	n <- nrow(mY)
  	L <- matrix(Lambda(b=b.gmmq, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau), ncol=1)
  	Ld <- Lambda.derivative(b=b.gmmq, mY=mY, mX=mX, mlX=mlX, fitphi=fitphi, fitlagphi=fitlagphi, tau=tau)
  	gni <- sweep(mZ*repmat((Gfn(-L,h.gmmq)-tau),1, ncol(mZ)), MARGIN=2, gbartrue, `-`)
  	g.bar <- colMeans(gni)
  	h.use <- max(sort(abs(L))[floor(n^(4/5))], h.gmmq)
  	jac <- t(mZ)%*%(Ld*repmat(Gpfn(-L,h.use),1,length(b.init)))/(-h.use)
  	W <- solve(t(gni)%*%gni/n)
  	b_onestep <- b.gmmq-solve(t(jac)%*%W%*%jac)%*%t(jac)%*%W%*%g.bar 
  	return(list(b=b_onestep,h=h.gmmq, G=g.bar, W=W))
}