require(dplyr)
###################################################################################################
#This file contains auxillary files for data subsetting and other miscellaneous tasks
#####################################################################################################

#####################################################################################################
#This function lags the data and returns both the lagged variables as well as its contemporaneous values
#####################################################################################################
lagdata <- function(idvar, X){
  condata <- data.frame(idvar, X) %>% group_by(idvar) %>% slice(-1)
  lagdata <- data.frame(idvar, X) %>% group_by(idvar) %>% slice(-n())
  data <- cbind(data.frame(condata), data.frame(lagdata[,-1]))
  return(data)
}
########################################################################################################
#boot resampling on IDs: bootstrapping on individuals (prodest.R)
#######################################################################################################
block.boot.resample <- function(idvar, R, seed){
  set.seed(seed)
  unique.ids <- unique(idvar) # find the unique values of panels in order to reshape the data
  panel.time.indices <- apply(unique.ids, 1, function(x) {return(list(which(idvar == x)))}) # find the time indices for each panel
  seq.indices <- 1:length(unique.ids) # the panel.time.indices list is indexed with sequential numbers: we mimic it
  boot.panel.id <- replicate(R, sample(seq.indices, replace = TRUE)) # generate a matrix of new IDs - R times
  new.indices <- list() # generate the matrix of the new indices
  ind <- 1:length(unique.ids)
  for (r in 1:R){ # for each boot rep we generate a vector of indices with rownames equal to a new - and fake - ID
    new.indices[[r]] <- cbind(unlist(mapply(function(x,y) {
      names(panel.time.indices[[x]][[1]]) <- rep(y,length(panel.time.indices[[x]][[1]]))
      return(list(panel.time.indices[[x]][[1]]))
    }, boot.panel.id[,r], ind))) # return a fake ID (sequential number) as row name and the index referring to the true ID
  }
  return(new.indices)
}
###########################################################################################################
###########################################################################################################
###########################################################################################################
lprq <- function(Y, X, Z, h, m, k, tau){
  xx <- apply(X, 2, function(x) seq(min(x), max(x), length=m))
  bx <- xx[1,]
  for (i in 1:m){
    x1 <- do.call(cbind, lapply(1:k, function(k) sweep(X, MARGIN=2, xx[i,], `-`)^k))
    hx <- sapply(1:k, function(k) h^k)
    wx <- sweep(x2, MARGIN=2, hx, `/`)
    r <- rq(Y, cbind(Z, x2), weights=wx, tau=tau)
    bx[i] <- r$coef[2]
  }
  b <- mean(bx)
  return(b)
}