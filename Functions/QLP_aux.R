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
##################################################################################################
# Smoothing function Itilde() and its derviative Itilde'(); Vector input supported
##################################################################################################
Itilde.KS17 <- function(u) {  ifelse(u >= 1, 1, ifelse(u > -1, 1/2 + (105/64)*(u-(5/3)*u^3+(7/5)*u^5 -(3/7)*u^7), 0)) }
Itilde.deriv.KS17 <- function(u) { ifelse(u > -1 & u < 1, (105/64)*(1-5*u^2+7*u^4-3*u^6), 0) }
Gfn <- function(v,h){      
  Itilde.KS17(v/h)    
    }
Gpfn <- function(v,h) {
  Itilde.deriv.KS17(v/h)
}
###########################################################################################################
###########################################################################################################
###########################################################################################################