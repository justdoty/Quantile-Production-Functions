setwd('/Users/justindoty/Documents/Research/Structural_Estimation/Production/Heterogeneity_in_Firms/R_Code')
set.seed(123456)
#This code is a modified version of Smoothed GMM for Quantile Models, de Castro, Galvao,
#Kaplan, and Liu (2018). See David Kaplan's website for more details
#https://faculty.missouri.edu/~kaplandm/

#Some data preparation follows prodest.R (Gabrielle Rovigatti)
source("gmmq.R")
source("ivqr.bw.R")
source("LP.R")
#Required for 1st step estimation
require(quantreg)
#Required for QGMM estimation variations
require(GenSA)
require(pracma)
#Optional for parallel computing (optional)
require(snow)
require(dplyr)
#Required for IVQR.BW (for two-step estimation)
require(MASS)
#Required for various codes to check ACF and LP estimates and bootstrapping
require(prodest)
QLP <- function(tau, va, state, free, proxy, idvar, timevar, h=0, b.init=NULL, R=20){
  #Make all data arguments into matrices
  va <- as.matrix(va)
  state <- as.matrix(state)
  free <- as.matrix(free)
  proxy <- as.matrix(proxy)
  idvar <- as.matrix(idvar)
  timevar<- as.matrix(timevar)
  snum <- ncol(state)
  fnum <- ncol(free)

  polyframe <- data.frame(state, proxy) # vars to be used in polynomial approximation
  mod <- model.matrix( ~.^2-1, data = polyframe) # generate the polynomial elements - this drops NAs
  mod <- mod[match(rownames(polyframe),rownames(mod)),] # replace NAs if there was any
  regvars <- cbind(free, mod, state^2, proxy^2)

  #Generate state lags
  lagstate <- state
  for (lag in 1:ncol(state)){
    lagstate[,lag] <- lagPanel(idvar=idvar, timevar=timevar, value=state[,lag])
  }
  #Generate free lags
  lagfree <- free
  for (lag in 1:ncol(free)){
    lagfree[,lag] <- lagPanel(idvar=idvar, timevar=timevar, value=free[,lag])
  }
  #Generate the matrix of data
  data <- suppressWarnings(as.matrix(data.frame(state=state, lagState=lagstate, free=free, lagFree=lagfree, va=va, idvar=idvar, timevar=timevar, regvars=regvars)))

  noboot <- finalQLP(tau=tau, h=h, ind=TRUE, data=data, fnum=fnum, snum=snum, b.init=b.init, boot=FALSE, gbar=TRUE)
  beta <- noboot[[1]]
  truegbar <- noboot[[2]]
  boot.indices <- block.boot.resample(idvar, R)
  boot.betas <- matrix(NA, R, (fnum+snum))

  for (i in 1:R){
    print(i)
    boot.betas[i,] <- finalQLP(tau=tau, h=h, ind=boot.indices[[i]], data=data, fnum=fnum, snum=snum, b.init=b.init, boot=TRUE, gbar=truegbar)
  }
  #Calculate standard deviations
  boot.beta <- apply(boot.betas, 2, mean)
  boot.errors <- apply(boot.betas, 2, sd)
  boot.CI<- apply(boot.betas, 2, function(x) quantile(x, c(0.025, 0.975)))
  #This returns the non-bootstrapped estimates, the matrix of bootstrapped estimates
  #The bootstrapped standard errors and the bootstrapped CI
  return(list(beta, boot.betas, boot.beta, boot.errors, boot.CI))
}
  #Function to estimate and to bootstrap QLP
finalQLP <- function(tau, h, ind, data, fnum, snum, b.init, boot, gbar){
  if (sum(as.numeric(ind))==length(ind)){ #if the ind variable is not always TRUE
    newid <- data[ind, 'idvar', drop = FALSE]
  } else {
    newid <- as.matrix(as.numeric(rownames(ind)))
    ind <- as.matrix(ind)
  }
  #change the index according to bootstrapped indices
  data <- data[ind,] 
  first.stage <- rq(data[,'va', drop = FALSE] ~ data[, grepl('regvars', colnames(data)), drop = FALSE], na.action = na.exclude, tau=tau)
  free <- data[, grepl('free', colnames(data), ignore.case=FALSE), drop = FALSE]
  phi <- fitted(first.stage)
  beta.free <-  as.numeric(coef(first.stage)[2:(1+fnum)])
  #If not specified, starting points are the first stage+normal noise
  if (is.null(b.init)){
    b.init <- coef(first.stage)[(2 + fnum):(1 + fnum + snum)]+rnorm((snum), 0, 0.01)
  } 
  #Clean Phi from the effects of free variables
  phi <- phi - (free%*%beta.free)
  newtime <- data[,'timevar', drop = FALSE]
  rownames(phi) <- NULL
  rownames(newtime) <- NULL
  #Lag Fitted Values  
  lag.phi <- lagPanel(value=phi, idvar=newid, timevar=newtime)
  #Clean the output from the effect of free variables
  va <- data[,'va', drop = FALSE] - (free%*%beta.free)
  state <- data[, grepl('state', colnames(data)), drop = FALSE]
  lagState <- data[, grepl('lagState', colnames(data)), drop = FALSE]
  lagFree <- data[, grepl('lagFree', colnames(data)), drop = FALSE]
  tmp.data <- na.omit(data.frame(state, lagState, lagFree, phi, lag.phi, va))
  lagFree <- tmp.data[, grepl('lagFree', colnames(tmp.data)), drop = FALSE]
  lagState <- tmp.data[, grepl('lagState', colnames(tmp.data)), drop = FALSE]
  Z <- as.matrix(cbind(tmp.data$state, lagState, lagFree))
  W <- solve(tau*(1-tau)*crossprod(Z)/nrow(Z))

  try.state <- GenSA(par=b.init, fn=goQLP, mZ=Z, mW=W, mX=tmp.data$state, mlX=tmp.data$lagState, 
  vphi=tmp.data$phi, vlag.phi=tmp.data$lag.phi, vres=tmp.data$va, tau=tau, h=h, gbar=gbar, lower=0, upper=1, control=list(max.time=5))
  beta.state <- as.numeric(try.state$par)

  if (gbar==TRUE){
    gbartrue <- g.bar(tau=tau, h=h, vtheta=beta.state, mZ=Z, mW=W, mX=tmp.data$state, mlX=tmp.data$lagState, 
      vphi=tmp.data$phi, vlag.phi=tmp.data$lag.phi, vres=tmp.data$va)
    return(list(c(beta.state, beta.free), gbartrue))
  } else {
    return(c(beta.free, beta.state))
  }
}

#Smoothing Kernels
  Itilde <- Itilde.KS17
  Itilde.deriv <- Itilde.deriv.KS17
#Add bandwidth as second argument to G() and G'() functions.
  Gfn <- function(v,h){      
    Itilde.KS17(v/h)    
    }
  Gpfn <- function(v,h){      
    Itilde.deriv.KS17(v/h)    
    }

g.bar <- function(tau, h, vtheta, mZ, mW, mX, mlX, vphi, vlag.phi, vres){
  vtheta <- as.matrix(as.numeric(vtheta))
  Omega <- vphi-mX%*%vtheta
  Omega_lag <- vlag.phi-mlX%*%vtheta
  Omega_lag_pol <- cbind(1, Omega_lag, Omega_lag^2, Omega_lag^3)
  g_b <- fitted(rq(Omega~Omega_lag+Omega_lag^2+Omega_lag^3, tau=tau))
  Lambda <- vres-(mX%*%vtheta)-g_b
  g.bar <- as.matrix(colMeans(mZ*repmat((Gfn(-Lambda, h)-tau), 1, ncol(mZ))))
  return(g.bar)
}

#QLP Objective Function
goQLP <- function(tau, h, vtheta, mZ, mW, mX, mlX, vphi, vlag.phi, vres, gbar){
  if (gbar==TRUE){
    g.bar <- g.bar(tau=tau, h=h, vtheta, mZ=mZ, mW=mW, mX=mX, mlX=mlX, vphi=vphi, vlag.phi=vlag.phi, vres=vres)
    crit <- nrow(mZ)*t(g.bar)%*%mW%*%g.bar
    return(crit)
  } else {
    g.bar.boot <- g.bar(tau=tau, h=h, vtheta, mZ=mZ, mW=mW, mX=mX, mlX=mlX, vphi=vphi, vlag.phi=vlag.phi, vres=vres)-gbar
    return(nrow(mZ)*t(g.bar.boot)%*%mW%*%g.bar.boot)
  }
}   

chile_panel <- read.csv('chile_panel.csv')
chile <- na.omit(subset(chile_panel, ciiu_3d==381))

results <- QLP(tau=0.5, va=chile$lnva, state=chile$lnk, free=cbind(chile$lnw, chile$lnb), proxy=chile$proxy_e, idvar=chile$id, timevar=chile$year, h=1e-6, b.init=NULL, R=20)
# for (q in 1:length(tau)){
#   results <- QLP(tau=tau[q], va=chile$lnva, state=chile$lnk, free=cbind(chile$lnw, chile$lnb), proxy=chile$proxy_e, id=chile$id, time=chile$year, h=1e-6)
#   print(results)
# }
#Load Chilean dataset
# chile_panel <- read.csv('chile_panel.csv')
# #Choose which industry to select
# industries <- c(311, 321, 381)
# #Industries that look best: 311, 384 (capital), 381
# #Industries that dont work at all: 312
# #Vector of quantiles
# tau <- seq(0.1, 0.9, by=0.05)
# #Store results for coefficients and standard errors
# # LP_coefficients <- replicate(length(industries), list(array(0, dim=4)))
# # LP_se <- replicate(length(industries), list(array(0, dim=4)))
# LP_coefficients <- replicate(length(industries), list(array(0, dim=3)))
# LP_se <- replicate(length(industries), list(array(0, dim=3)))
# coefficients <- replicate(length(industries), list(array(0, dim=c(length(tau), 4))))
# se <- replicate(length(industries), list(array(0, dim=c(length(tau), 4))))
# #Estimate for different industries
# #TO DO: Include Returns to Scale and test for returns to scale
# for (ISIC in 1:length(industries)){
#   print(industries[ISIC])
#   chile <- subset(chile_panel, ciiu_3d==industries[ISIC])
#   # LP_Results <- LP(va=chile$lnva, state=chile$lnk, free=cbind(chile$lnw, chile$lnb), proxy=chile$proxy_e, id=chile$id, time=chile$year, b.init=NULL)
#   # LP_coefficients[[ISIC]] <- LP_Results[[1]]
#   # LP_se[[ISIC]] <- LP_Results[[2]]
#   # print(cbind(LP_Results[[1]], LP_Results[[2]]))
#   LP_Results <- prodestLP(Y=chile$lnva, sX=chile$lnk, fX=cbind(chile$lnw, chile$lnb), pX=chile$proxy_e, idvar=chile$id, timevar=chile$year, opt='solnp')@Estimates
#   LP_coefficients[[ISIC]] <- as.numeric(LP_Results$pars)[c(3,1,2)]
#   LP_se[[ISIC]] <- as.numeric(LP_Results$std.errors)[c(3,1,2)]
#   print(cbind(LP_coefficients[[ISIC]], LP_se[[ISIC]]))
#   for (q in 1:length(tau)){
#     print(tau[q])
#     results <- QLP(tau=tau[q], va=chile$lnva, state=chile$lnk, free=cbind(chile$lnw, chile$lnb), proxy=chile$proxy_e, id=chile$id, time=chile$year, h=1e-6, b.init=NULL)
#     coefficients[[ISIC]][q,] <- results[[1]]
#     se[[ISIC]][q,] <- results[[2]]
#     print(cbind(results[[1]], results[[2]]))
#   }
# }

# #Make an estimates table for Quantile GMM
# estimates <- data.frame(cbind(rep(tau, length(industries)), cbind(do.call(rbind, coefficients), do.call(rbind, se))[,c(rbind(c(1:4), 4+(1:4)))]))
# colnames(estimates) <- c('Tau','K',"se_K", 'Lw', "se_Lw", 'Lb', "se_Lb",'Rho', "se_Rho")
# #Make an estimates table for LP
# # estimates_LP <- data.frame(cbind(do.call(rbind, LP_coefficients), do.call(rbind, LP_se))[,c(rbind(c(1:4), 4+(1:4)))])
# # colnames(estimates_LP) <- c('K', 'se_K', 'Lw', 'se_Lw', 'Lb', 'se_Lb', "Rho", "se_Rho")

# estimates_LP <- data.frame(cbind(do.call(rbind, LP_coefficients), do.call(rbind, LP_se))[,c(rbind(c(1:3), 3+(1:3)))])
# colnames(estimates_LP) <- c('K', 'se_K', 'Lw', 'se_Lw', 'Lb', 'se_Lb')
# #Prepare estimates for table in paper/presentation
# require(xtable)
# tau_table <- c(0.25, 0.5, 0.75)

# #Table Labels
# ISIC_labels <- array(NA, length(tau_table)*length(industries)); ISIC_labels[seq(1, length(tau_table)*length(industries), by=length(tau_table))] <- industries
# ISIC_labels[is.na(ISIC_labels)] <- ""

# estimates_table <- subset(cbind(ISIC_labels, estimates[rep(tau, length(industries))%in%tau_table, ]), select=-c(Rho, se_Rho))
# colnames(estimates_table) <- c("Industry (ISIC code)", "$\\tau$", "Coef.", 's.e.', "Coef.",'s.e.', "Coef.", 's.e.')

# estimates_table <- xtable(estimates_table, digits=c(0,0,2,3,4,3,4,3,4))
# align(estimates_table) <- rep('c', 9)
# addtorow <- list()
# addtorow$pos <- list(-1)
# addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Skilled Labor} & \\multicolumn{2}{c}{Unskilled Labor} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}'
# print(estimates_table, hline.after=c(0,nrow(estimates_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x)

# ############################Coefficicent Plots######################################
# require(ggplot2)
# require(cowplot)
# K_plot <- list(); Lw_plot <- list(); Lb_plot <- list()
# for (p in 1:length(industries)){
#   #Capital
#   ky <- split(estimates$K, ceiling(seq_along(estimates$K)/length(tau)))[[p]]
#   ksd <- split(estimates$se_K, ceiling(seq_along(estimates$se_K)/length(tau)))[[p]]
#   K_data <- data.frame(x=tau, y=ky, z=estimates_LP$K[p], lower=(ky-qnorm(0.95)*ksd), upper=(ky+qnorm(0.95)*ksd), lower_LP=(estimates_LP$K[p]-qnorm(0.95)*estimates_LP$se_K[p]), upper_LP=(estimates_LP$K[p]+qnorm(0.95)*estimates_LP$se_K[p]))
#   K_plot[[p]] <- ggplot(K_data, aes(x=x, y=y)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) +  geom_hline(yintercept=K_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(K_data$lower_LP, K_data$upper_LP), linetype='dashed', color='red') + expand_limits(y=c(0,.8))
#   #Skilled Labor
#   lwy <- split(estimates$Lw, ceiling(seq_along(estimates$Lw)/length(tau)))[[p]]
#   lwsd <- split(estimates$se_Lw, ceiling(seq_along(estimates$se_Lw)/length(tau)))[[p]]
#   Lw_data <- data.frame(x=tau, y=lwy, z=estimates_LP$Lw[p], lower=(lwy-qnorm(0.95)*lwsd), upper=(lwy+qnorm(0.95)*lwsd), lower_LP=(estimates_LP$Lw[p]-qnorm(0.95)*estimates_LP$se_Lw[p]), upper_LP=(estimates_LP$Lw[p]+qnorm(0.95)*estimates_LP$se_Lw[p]))
#   Lw_plot[[p]] <- ggplot(Lw_data, aes(x=x, y=y)) + xlab(expression('percentile-'*tau)) + ylab("Skilled Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) +  geom_hline(yintercept=Lw_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(Lw_data$lower_LP, Lw_data$upper_LP), linetype='dashed', color='red') + expand_limits(y=c(0,.8))
#   #Unskilled Labor
#   lby <- split(estimates$Lb, ceiling(seq_along(estimates$Lb)/length(tau)))[[p]]
#   lbsd <- split(estimates$se_Lb, ceiling(seq_along(estimates$se_Lb)/length(tau)))[[p]]
#   Lb_data <- data.frame(x=tau, y=lby, z=estimates_LP$Lb[p], lower=(lby-qnorm(0.95)*lbsd), upper=(lby+qnorm(0.95)*lbsd), lower_LP=(estimates_LP$Lb[p]-qnorm(0.95)*estimates_LP$se_Lb[p]), upper_LP=(estimates_LP$Lb[p]+qnorm(0.95)*estimates_LP$se_Lb[p]))
#   Lb_plot[[p]] <- ggplot(Lb_data, aes(x=x, y=y)) + xlab(expression('percentile-'*tau)) + ylab("Unskilled Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) +  geom_hline(yintercept=Lb_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(Lb_data$lower_LP, Lb_data$upper_LP), linetype='dashed', color='red') + expand_limits(y=c(0,.8))
# }
# #Combine plots across industries
# #Industry ISIC Code Plot Labels
# Industry_1 <- ggdraw() + draw_label("Food Products (311)", fontface='plain')
# Industry_2 <- ggdraw() + draw_label("Textiles (321)", fontface='plain')
# Industry_3 <- ggdraw() + draw_label("Fabriacted Metals (381)", fontface='plain')

# #Coefficient Plot for 1st industry
# L_Plot_1 <- plot_grid(Lw_plot[[1]], Lb_plot[[1]], rel_heights=c(0.1, 1), ncol=2)
# K_Plot_1 <- plot_grid(NULL, K_plot[[1]], NULL, ncol=3, rel_widths=c(0.25,0.5,0.25))
# Coefficient_Plot_1 <- plot_grid(Industry_1, L_Plot_1, K_Plot_1, ncol=1, rel_heights=c(0.1, 1, 1))
# save_plot("Coefficient_Plot_1.png", Coefficient_Plot_1, base_height = 8, base_width = 9)

# #Coefficient Plot for 2nd industry
# L_Plot_2 <- plot_grid(Lw_plot[[2]], Lb_plot[[2]], rel_heights=c(0.1, 1), ncol=2)
# K_Plot_2 <- plot_grid(NULL, K_plot[[2]], NULL, ncol=3, rel_widths=c(0.25,0.5,0.25))
# Coefficient_Plot_2 <- plot_grid(Industry_2, L_Plot_2, K_Plot_2, ncol=1, rel_heights=c(0.1, 1, 1))
# save_plot("Coefficient_Plot_2.png", Coefficient_Plot_2, base_height = 8, base_width = 9)

# #Coefficient Plot for 3rd industry
# L_Plot_3 <- plot_grid(Lw_plot[[3]], Lb_plot[[3]], rel_heights=c(0.1, 1), ncol=2)
# K_Plot_3 <- plot_grid(NULL, K_plot[[3]], NULL, ncol=3, rel_widths=c(0.25,0.5,0.25))
# Coefficient_Plot_3 <- plot_grid(Industry_3, L_Plot_3, K_Plot_3, ncol=1, rel_heights=c(0.1, 1, 1))
# save_plot("Coefficient_Plot_3.png", Coefficient_Plot_3, base_height = 8, base_width = 9)




