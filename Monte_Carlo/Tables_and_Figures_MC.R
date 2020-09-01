setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo')
require(xtable)
require(ggplot2)
require(cowplot)
require(abind)

############################Load R Environments###############################
load("simulation_LP.Rdata")
# load("simulation_ACF.Rdata")
#####################Compute Estimates and Standard Errors##################
tau_table <- c(0.1, 0.25, 0.5, 0.75, 0.9)
dB <- 2
beta <- alpha

#Need to remove spurious minumum problem 0 capital coefficient and 1 for labor coefficient
# tol <- 0.05
# spurious_ACF <- list()
# for (i in 1:length(DGPs)){
# 	k0 <- resmat_ACF[,,i][,1][(abs(resmat_ACF[,,i][,1])>tol)|(abs(resmat_ACF[,,i][,2]-1)>tol)]
# 	l0 <- resmat_ACF[,,i][,2][(abs(resmat_ACF[,,i][,1])>tol)|(abs(resmat_ACF[,,i][,2]-1)>tol)]
# 	spurious_ACF[[i]] <- cbind(k0, l0)
# }
#Coefficients
# QACF_coef <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
QLP_coef <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
# ACF_coef <- replicate(length(DGPs), list(array(0, dim=2)))
LP_coef <- replicate(length(DGPs), list(array(0, dim=2)))
#Standard Deviations
# QACF_se <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
QLP_se <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
# ACF_se <- replicate(length(DGPs), list(array(0, dim=2)))
LP_se <- replicate(length(DGPs), list(array(0, dim=2)))
#Lower Confidence Bound
# QACF_Lower <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
QLP_Lower <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
# ACF_Lower <- replicate(length(DGPs), list(array(0, dim=2)))
LP_Lower <- replicate(length(DGPs), list(array(0, dim=2)))
#Upper Confidence Bound
# QACF_Upper <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
QLP_Upper <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
# ACF_Upper <- replicate(length(DGPs), list(array(0, dim=2)))
LP_Upper <- replicate(length(DGPs), list(array(0, dim=2)))

#Calculations
for (d in 1:length(DGPs)){
	#ACF and LP Coefficients
	# ACF_coef[[d]] <- round(apply(spurious_ACF[[d]], 2, mean), digits=3)
	LP_coef[[d]] <- round(apply(resmat_LP[,,d], 2, mean), digits=3)
	#ACF and LP Standard Deviations
	# ACF_se[[d]] <- round(apply(spurious_ACF[[d]], 2, sd), digits=4)
	LP_se[[d]] <- round(apply(resmat_LP[,,d], 2, sd), digits=4)
	#ACF and LP Lower Confidence Band (90%)
	# ACF_Lower[[d]] <- round(apply(spurious_ACF[[d]], 2, function(x) quantile(x, 0.05)), digits=3)
	LP_Lower[[d]] <- round(apply(resmat_LP[,,d], 2, function(x) quantile(x, 0.05)), digits=3)
	#ACF and LP Upper Confidence Band (90%)
	# ACF_Upper[[d]] <- round(apply(spurious_ACF[[d]], 2, function(x) quantile(x, 0.95)), digits=3)
	LP_Upper[[d]] <- round(apply(resmat_LP[,,d], 2, function(x) quantile(x, 0.95)), digits=3)
	for (q in 1:length(tau)){
	  #QACF and QLP Coefficients
	  # QACF_coef[[d]][q,] <- round(apply(resmat_ACFQ[,,,d][,,q], 2, mean), digits=3)
	  QLP_coef[[d]][q,] <- round(apply(resmat_LPQ[,,,d][,,q], 2, mean), digits=3)
	  #QACF and QLP Standard Deviations
	  # QACF_se[[d]][q,] <- round(apply(resmat_ACFQ[,,,d][,,q], 2, sd), digits=4)
	  QLP_se[[d]][q,] <- round(apply(resmat_LPQ[,,,d][,,q], 2, sd), digits=4)
	  #QACF and QLP Lower Confidence Band (90%)
	  # QACF_Lower[[d]][q,] <- round(apply(resmat_ACFQ[,,,d][,,q], 2, function(x) quantile(x, 0.05)), digits=3)
	  QLP_Lower[[d]][q,] <- round(apply(resmat_LPQ[,,,d][,,q], 2, function(x) quantile(x, 0.05)), digits=3)
	  #QACF and QLP Upper Confidence Band (90%)
	  # QACF_Upper[[d]][q,] <- round(apply(resmat_ACFQ[,,,d][,,q], 2, function(x) quantile(x, 0.95)), digits=3)
	  QLP_Upper[[d]][q,] <- round(apply(resmat_LPQ[,,,d][,,q], 2, function(x) quantile(x, 0.95)), digits=3)
	}
}
#####################Create a table for all estimates#######################
#For QACF and QLP
# QACF_estimates <- data.frame(cbind(rep(tau, length(DGPs)), cbind(do.call(rbind, QACF_coef), do.call(rbind, QACF_se))[,c(rbind(c(1:dB), dB+(1:dB)))]))
# colnames(QACF_estimates) <- c('Tau','K',"sd_K", 'L', "sd_L", 'Rho', "sd_Rho")
# QACF_Confidence <- data.frame(cbind(rep(tau, length(DGPs)), cbind(do.call(rbind, QACF_Lower), do.call(rbind, QACF_Upper))[,c(rbind(c(1:3), 3+(1:3)))]))
# colnames(QACF_Confidence) <- c('Tau', 'Lower_K', 'Upper_K', 'Lower_L', 'Upper_L', 'Lower_Rho', 'Upper_Rho')
QLP_estimates <- data.frame(cbind(rep(tau, length(DGPs)), cbind(do.call(rbind, QLP_coef), do.call(rbind, QLP_se))[,c(rbind(c(1:dB), dB+(1:dB)))]))
colnames(QLP_estimates) <- c('Tau','K',"sd_K", 'L', "sd_L")
QLP_Confidence <- data.frame(cbind(rep(tau, length(DGPs)), cbind(do.call(rbind, QLP_Lower), do.call(rbind, QLP_Upper))[,c(rbind(c(1:2), 2+(1:2)))]))
colnames(QLP_Confidence) <- c('Tau', 'Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
#For ACF and LP
# ACF_estimates <- data.frame(cbind(do.call(rbind, ACF_coef), do.call(rbind, ACF_se))[,c(rbind(c(1:2), 2+(1:2)))])
# colnames(ACF_estimates) <- c('K', 'sd_K', 'L', 'sd_L')
# ACF_Confidence <- data.frame(cbind(do.call(rbind, ACF_Lower), do.call(rbind, ACF_Upper))[,c(rbind(c(1:2), 2+(1:2)))])
# colnames(ACF_Confidence) <- c('Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
LP_estimates <- data.frame(cbind(do.call(rbind, LP_coef), do.call(rbind, LP_se))[,c(rbind(c(1:2), 2+(1:2)))])
colnames(LP_estimates) <- c('K', 'sd_K', 'L', 'sd_L')
LP_Confidence <- data.frame(cbind(do.call(rbind, LP_Lower), do.call(rbind, LP_Upper))[,c(rbind(c(1:2), 2+(1:2)))])
colnames(LP_Confidence) <- c('Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
#####################Compute Bias and MSE###################################
#Bias
#For QACF and QLP: Note that I exclude Rho as it is not of primary interest
# QACF_Bias <- round(beta-QACF_estimates[,c("K", "L")], digits=4)
QLP_Bias <- round(beta-QLP_estimates[,c("K", "L")], digits=4)
#For ACF and LP
# ACF_Bias <- round(beta-apply(ACF_estimates[, c("K", "L")], 2, function(x) rep(x, each=length(tau))), digits=4)
LP_Bias <- round(beta-apply(LP_estimates[, c("K", "L")], 2, function(x) rep(x, each=length(tau))), digits=4)
#MSE
#For QACF and QLP
# QACF_MSE <- round(QACF_estimates[, c("sd_K", "sd_L")]^2+QACF_Bias^2, digits=4)
QLP_MSE <- round(QLP_estimates[, c("sd_K", "sd_L")]^2+QLP_Bias^2, digits=4)
#For ACF and LP
# ACF_MSE <- round(apply(ACF_estimates[, c("sd_K", "sd_L")], 2, function(x) rep(x, each=length(tau)))^2+ACF_Bias^2, digits=4)
LP_MSE <- round(apply(LP_estimates[, c("sd_K", "sd_L")], 2, function(x) rep(x, each=length(tau)))^2+LP_Bias^2, digits=4)
#Create a table for Bias and MSE for QACF/ACF and QLP/LP
# ACF_Bias_MSE <- data.frame(cbind(tau, cbind(QACF_Bias, QACF_MSE)[,c(rbind(c(1:2), 2+(1:2)))], cbind(ACF_Bias, ACF_MSE)[, c(rbind(c(1:2), 2+(1:2)))]))
# colnames(ACF_Bias_MSE) <- c('Tau', 'QACF_Bias_K', 'QACF_MSE_K', 'QACF_Bias_L', 'QACF_MSE_L', 'ACF_Bias_K', 'ACF_MSE_K', 'ACF_Bias_L', 'ACF_MSE_L')
LP_Bias_MSE <- data.frame(cbind(tau, cbind(QLP_Bias, QLP_MSE)[,c(rbind(c(1:2), 2+(1:2)))], cbind(LP_Bias, LP_MSE)[, c(rbind(c(1:2), 2+(1:2)))]))
colnames(LP_Bias_MSE) <- c('Tau', 'QLP_Bias_K', 'QLP_MSE_K', 'QLP_Bias_L', 'QLP_MSE_L', 'LP_Bias_K', 'LP_MSE_K', 'LP_Bias_L', 'LP_MSE_L')
#####################Prepare Results for Xtable###############################
#################For Coefficients and Standard Errors########################
DGP_labels <- array(NA, length(tau_table)*length(DGPs)); DGP_labels[seq(1, length(tau_table)*length(DGPs), by=length(tau_table))] <- c(1:length(DGPs))
DGP_labels[is.na(DGP_labels)] <- ""
#For QACF and ACF
# QACF_estimates_table <- subset(cbind(DGP_labels, QACF_estimates[rep(tau, length(DGPs))%in%tau_table, ]), select=-c(Rho, sd_Rho))
# colnames(QACF_estimates_table) <- c("DGP", "$\\tau$","Coef.", 'Std. Dev.', "Coef.",'Std. Dev.')
# QACF_estimates_table <- xtable(QACF_estimates_table, digits=c(0,0,2,3,4,3,4))
# align(QACF_estimates_table) <- rep('c', 7)
# addtorow <- list()
# addtorow$pos <- list(-1)
# addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{$\\beta_{k}$}  & \\multicolumn{2}{c}{$\\beta_{l}$} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6}'
# print(QACF_estimates_table, hline.after=c(0,nrow(QACF_estimates_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/QACF_Sim_Estimates.txt")
#For QLP and LP
QLP_estimates_table <- subset(cbind(DGP_labels, QLP_estimates[rep(tau, length(DGPs))%in%tau_table, ]))
colnames(QLP_estimates_table) <- c("DGP", "$\\tau$","Coef.", 'Std. Dev.', "Coef.",'Std. Dev.')
QLP_estimates_table <- xtable(QLP_estimates_table, digits=c(0,0,2,3,4,3,4))
align(QLP_estimates_table) <- rep('c', 7)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{$\\beta_{k}$}  & \\multicolumn{2}{c}{$\\beta_{l}$} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6}'
print(QLP_estimates_table, hline.after=c(0,nrow(QLP_estimates_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/QLP_Sim_Estimates.txt")
#######################For Bias and MSE########################################
#For QACF and ACF
# ACF_Bias_MSE_table <- cbind(DGP_labels, ACF_Bias_MSE[rep(tau, length(DGPs))%in%tau_table, ])
# colnames(ACF_Bias_MSE_table) <- c("DGP", "$\\tau$","Bias", 'MSE', "Bias",'MSE',"Bias",'MSE', 'Bias', 'MSE')
# ACF_Bias_MSE_table <- xtable(ACF_Bias_MSE_table, digits=c(0,0,2,4,4,4,4,4,4,4,4))
# align(ACF_Bias_MSE_table) <- rep('c', 11)
# addtorow_ACF_Bias_MSE <- list()
# addtorow_ACF_Bias_MSE$pos <- list(-1)
# addtorow_ACF_Bias_MSE$command <- '\\hline\\hline & & \\multicolumn{4}{c}{QACF} & \\multicolumn{4}{c}{ACF} \\\\ \\cmidrule(lr){3-6} \\cmidrule(lr){7-10} \\\\& & \\multicolumn{2}{c}{$\\beta_{k}$}  & \\multicolumn{2}{c}{$\\beta_{l}$} &\\multicolumn{2}{c}{$\\beta_{k}$} & \\multicolumn{2}{c}{$\\beta_{l}$}  \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}'
# print(ACF_Bias_MSE_table, hline.after=c(0,nrow(ACF_Bias_MSE_table)), add.to.row=addtorow_ACF_Bias_MSE, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/QACF_Sim_Bias_MSE.txt")
#For QLP and LP
LP_Bias_MSE_table <- cbind(DGP_labels, LP_Bias_MSE[rep(tau, length(DGPs))%in%tau_table, ])
colnames(LP_Bias_MSE_table) <- c("DGP", "$\\tau$","Bias", 'MSE', "Bias",'MSE',"Bias",'MSE', 'Bias', 'MSE')
LP_Bias_MSE_table <- xtable(LP_Bias_MSE_table, digits=c(0,0,2,4,4,4,4,4,4,4,4))
align(LP_Bias_MSE_table) <- rep('c', 11)
addtorow_LP_Bias_MSE <- list()
addtorow_LP_Bias_MSE$pos <- list(-1)
addtorow_LP_Bias_MSE$command <- '\\hline\\hline & & \\multicolumn{4}{c}{QLP} & \\multicolumn{4}{c}{LP} \\\\ \\cmidrule(lr){3-6} \\cmidrule(lr){7-10} \\\\& & \\multicolumn{2}{c}{$\\beta_{k}$}  & \\multicolumn{2}{c}{$\\beta_{l}$} &\\multicolumn{2}{c}{$\\beta_{k}$} & \\multicolumn{2}{c}{$\\beta_{l}$}  \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}'
print(LP_Bias_MSE_table, hline.after=c(0,nrow(LP_Bias_MSE_table)), add.to.row=addtorow_LP_Bias_MSE, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/QLP_Sim_Bias_MSE.txt")
############################# Box Plots #############################################
tau_box <- c(0.1, 0.25, 0.5, 0.75, 0.9)
# ACF_K_Box_Plot <- list(); ACF_L_Box_Plot <- list()
LP_K_Box_Plot <- list(); LP_L_Box_Plot <- list()
for (p in 1:length(DGPs)){
	# #Capital QACF
	# ACF_K_Box_Data <- data.frame(x=factor(rep(tau_box, each=dim(resmat_ACFQ)[1])), y=c(resmat_ACFQ[,,,p][,,match(tau_box, tau)][,1,]))
	# ACF_K_Box_Plot[[p]] <- ggplot(ACF_K_Box_Data, aes(x=x, y=y)) +geom_boxplot() + xlab(expression('percentile-'*tau)) + ylab(expression(beta["k"]))
	# #Labor QACF
	# ACF_L_Box_Data <- data.frame(x=factor(rep(tau_box, each=dim(resmat_ACFQ)[1])), y=c(resmat_ACFQ[,,,p][,,match(tau_box, tau)][,2,]))
	# ACF_L_Box_Plot[[p]] <- ggplot(ACF_L_Box_Data, aes(x=x, y=y)) +geom_boxplot() + xlab(expression('percentile-'*tau)) + ylab(expression(beta["l"])) 
	# #Capital QLP
	LP_K_Box_Data <- data.frame(x=factor(rep(tau_box, each=dim(resmat_LPQ)[1])), y=c(resmat_LPQ[,,,p][,,match(tau_box, tau)][,1,]))
	LP_K_Box_Plot[[p]] <- ggplot(LP_K_Box_Data, aes(x=x, y=y)) +geom_boxplot() + xlab(expression('percentile-'*tau)) + ylab(expression(beta["k"]))
	#Labor QLP
	LP_L_Box_Data <- data.frame(x=factor(rep(tau_box, each=dim(resmat_LPQ)[1])), y=c(resmat_LPQ[,,,p][,,match(tau_box, tau)][,2,]))
	LP_L_Box_Plot[[p]] <- ggplot(LP_L_Box_Data, aes(x=x, y=y)) +geom_boxplot() + xlab(expression('percentile-'*tau)) + ylab(expression(beta["l"])) 
}
####################Combine plots to grid and save##############################
DGP1 <- ggdraw() + draw_label("DGP 1", fontface='plain')
DGP2 <- ggdraw() + draw_label("DGP 2", fontface='plain')
#For QACF
# QACF_Box_Plot_1 <- plot_grid(DGP1, plot_grid(ACF_K_Box_Plot[[1]], ACF_L_Box_Plot[[1]]), ncol=1, rel_heights=c(0.1, 1))
# QACF_Box_Plot_2 <- plot_grid(DGP2, plot_grid(ACF_K_Box_Plot[[2]], ACF_L_Box_Plot[[2]]), ncol=1, rel_heights=c(0.1, 1))
# QACF_Box_Plot <- plot_grid(QACF_Box_Plot_1, QACF_Box_Plot_2, nrow=2, align='v')
# save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/QACF_Box_Plot.png", QACF_Box_Plot, base_height = 8, base_width = 9)
# #For QLP
QLP_Box_Plot_1 <- plot_grid(DGP1, plot_grid(LP_K_Box_Plot[[1]], LP_L_Box_Plot[[1]]), ncol=1, rel_heights=c(0.1, 1))
QLP_Box_Plot_2 <- plot_grid(DGP2, plot_grid(LP_K_Box_Plot[[2]], LP_L_Box_Plot[[2]]), ncol=1, rel_heights=c(0.1, 1))
QLP_Box_Plot <- plot_grid(QLP_Box_Plot_1, QLP_Box_Plot_2, nrow=2, align='v')
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/QLP_Box_Plot.png", QLP_Box_Plot, base_height = 8, base_width = 9)
############################Coefficicent Plots######################################
# ACF_K_plot <- list(); ACF_L_plot <- list()
LP_K_plot <- list(); LP_L_plot <- list()
for (p in 1:length(DGPs)){
	# #ACF Capital
	# acfky <- split(QACF_estimates$K, ceiling(seq_along(QACF_estimates$K)/length(tau)))[[p]]
	# acfksd <- split(QACF_estimates$sd_K, ceiling(seq_along(QACF_estimates$sd_K)/length(tau)))[[p]]
	# acfklower <- split(QACF_Confidence$Lower_K, ceiling(seq_along(QACF_Confidence$Lower_K)/length(tau)))[[p]]
	# acfkupper <- split(QACF_Confidence$Upper_K, ceiling(seq_along(QACF_Confidence$Upper_K)/length(tau)))[[p]]
	# ACF_K_data <- data.frame(x=tau, y=acfky, z=ACF_estimates$K[p], lower=acfklower, upper=acfkupper, lower_ACF=ACF_Confidence$Lower_K[p], upper_ACF=ACF_Confidence$Upper_K[p])
	# ACF_K_plot[[p]] <- ggplot(ACF_K_data, aes(x=x, y=y)) + xlab(expression('percentile-'*tau)) + ylab(expression(beta["k"])) + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) +  geom_hline(yintercept=ACF_K_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(ACF_K_data$lower_ACF, ACF_K_data$upper_ACF), linetype='dashed', color='red')
	# #LP Capital
	lpky <- split(QLP_estimates$K, ceiling(seq_along(QLP_estimates$K)/length(tau)))[[p]]
	lpksd <- split(QLP_estimates$sd_K, ceiling(seq_along(QLP_estimates$sd_K)/length(tau)))[[p]]
	lpklower <- split(QLP_Confidence$Lower_K, ceiling(seq_along(QLP_Confidence$Lower_K)/length(tau)))[[p]]
	lpkupper <- split(QLP_Confidence$Upper_K, ceiling(seq_along(QLP_Confidence$Upper_K)/length(tau)))[[p]]
	LP_K_data <- data.frame(x=tau, y=lpky, z=LP_estimates$K[p], lower=lpklower, upper=lpkupper, lower_LP=LP_Confidence$Lower_K[p], upper_LP=LP_Confidence$Upper_K[p])
	LP_K_plot[[p]] <- ggplot(LP_K_data, aes(x=x, y=y)) + xlab(expression('percentile-'*tau)) + ylab(expression(beta["k"])) + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) +  geom_hline(yintercept=LP_K_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(LP_K_data$lower_LP, LP_K_data$upper_LP), linetype='dashed', color='red')
	#ACF Labor
	# acfly <- split(QACF_estimates$L, ceiling(seq_along(QACF_estimates$L)/length(tau)))[[p]]
	# acflsd <- split(QACF_estimates$sd_L, ceiling(seq_along(QACF_estimates$sd_L)/length(tau)))[[p]]
	# acfllower <- split(QACF_Confidence$Lower_L, ceiling(seq_along(QACF_Confidence$Lower_L)/length(tau)))[[p]]
	# acflupper <- split(QACF_Confidence$Upper_L, ceiling(seq_along(QACF_Confidence$Upper_L)/length(tau)))[[p]]
	# ACF_L_data <- data.frame(x=tau, y=acfly, z=ACF_estimates$L[p], lower=acfllower, upper=acflupper, lower_ACF=ACF_Confidence$Lower_L[p], upper_ACF=ACF_Confidence$Upper_L[p])
	# ACF_L_plot[[p]] <- ggplot(ACF_L_data, aes(x=x, y=y)) + xlab(expression('percentile-'*tau)) + ylab(expression(beta["l"])) + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) +  geom_hline(yintercept=ACF_L_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(ACF_L_data$lower_ACF, ACF_L_data$upper_ACF), linetype='dashed', color='red')
	# #LP Labor
	lply <- split(QLP_estimates$L, ceiling(seq_along(QLP_estimates$L)/length(tau)))[[p]]
	lplsd <- split(QLP_estimates$sd_L, ceiling(seq_along(QLP_estimates$sd_L)/length(tau)))[[p]]
	lpllower <- split(QLP_Confidence$Lower_L, ceiling(seq_along(QLP_Confidence$Lower_L)/length(tau)))[[p]]
	lplupper <- split(QLP_Confidence$Upper_L, ceiling(seq_along(QLP_Confidence$Upper_L)/length(tau)))[[p]]
	LP_L_data <- data.frame(x=tau, y=lply, z=LP_estimates$L[p], lower=lpllower, upper=lplupper, lower_LP=LP_Confidence$Lower_L[p], upper_LP=LP_Confidence$Upper_L[p])
	LP_L_plot[[p]] <- ggplot(LP_L_data, aes(x=x, y=y)) + xlab(expression('percentile-'*tau)) + ylab(expression(beta["l"])) + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) +  geom_hline(yintercept=LP_L_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(LP_L_data$lower_LP, LP_L_data$upper_LP), linetype='dashed', color='red')
}
#Combine plots to grid and save
#ACF
# ACF_Coefficient_Plot_1 <- plot_grid(DGP1, plot_grid(ACF_K_plot[[1]], ACF_L_plot[[1]]), ncol=1, rel_heights=c(0.1, 1))
# ACF_Coefficient_Plot_2 <- plot_grid(DGP2, plot_grid(ACF_K_plot[[2]], ACF_L_plot[[2]]), ncol=1, rel_heights=c(0.1, 1))
# ACF_Coefficient_Plot <- plot_grid(ACF_Coefficient_Plot_1, ACF_Coefficient_Plot_2, nrow=2, align='v')
# save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/ACF_Coefficient_Plot.png", ACF_Coefficient_Plot, base_height = 8, base_width = 9)
# #LP
LP_Coefficient_Plot_1 <- plot_grid(DGP1, plot_grid(LP_K_plot[[1]], LP_L_plot[[1]]), ncol=1, rel_heights=c(0.1, 1))
LP_Coefficient_Plot_2 <- plot_grid(DGP2, plot_grid(LP_K_plot[[2]], LP_L_plot[[2]]), ncol=1, rel_heights=c(0.1, 1))
LP_Coefficient_Plot <- plot_grid(LP_Coefficient_Plot_1, LP_Coefficient_Plot_2, nrow=2, align='v')
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/LP_Coefficient_Plot.png", LP_Coefficient_Plot, base_height = 8, base_width = 9)
######################Mean Square Error plots##########################################
# ACF_K_MSE_plot <- list(); ACF_L_MSE_plot <- list()
LP_K_MSE_plot <- list(); LP_L_MSE_plot <- list()
for (p in 1:length(DGPs)){
	#ACF Capital
	# qacfky <- split(ACF_Bias_MSE$QACF_MSE_K, ceiling(seq_along(ACF_Bias_MSE$QACF_MSE_K)/length(tau)))[[p]]
	# acfky <- split(ACF_Bias_MSE$ACF_MSE_K, ceiling(seq_along(ACF_Bias_MSE$ACF_MSE_K)/length(tau)))[[p]]
	# QACF_K_MSE_Data <- data.frame(x=tau, y=qacfky)
	# ACF_K_MSE_Data <- data.frame(x=tau, y=acfky)
	# ACF_K_MSE_plot[[p]] <- ggplot(QACF_K_MSE_Data, aes(x=x, y=y)) + geom_line() + geom_line(data=ACF_K_MSE_Data, color='red', linetype='dotdash') + xlab(expression('percentile-'*tau)) + ylab(expression('MSE'*(beta["k"])))  
	# #LP Capital
	qlpky <- split(LP_Bias_MSE$QLP_MSE_K, ceiling(seq_along(LP_Bias_MSE$QLP_MSE_K)/length(tau)))[[p]]
	lpky <- split(LP_Bias_MSE$LP_MSE_K, ceiling(seq_along(LP_Bias_MSE$LP_MSE_K)/length(tau)))[[p]]
	QLP_K_MSE_Data <- data.frame(x=tau, y=qlpky)
	LP_K_MSE_Data <- data.frame(x=tau, y=lpky)
	LP_K_MSE_plot[[p]] <- ggplot(QLP_K_MSE_Data, aes(x=x, y=y)) + geom_line() + geom_line(data=LP_K_MSE_Data, color='red', linetype='dotdash') + xlab(expression('percentile-'*tau)) + ylab(expression('MSE'*(beta["k"])))
	# #ACF Labor
	# qacfly <- split(ACF_Bias_MSE$QACF_MSE_L, ceiling(seq_along(ACF_Bias_MSE$QACF_MSE_L)/length(tau)))[[p]]
	# acfly <- split(ACF_Bias_MSE$ACF_MSE_L, ceiling(seq_along(ACF_Bias_MSE$ACF_MSE_L)/length(tau)))[[p]]
	# QACF_L_MSE_Data <- data.frame(x=tau, y=qacfly)
	# ACF_L_MSE_Data <- data.frame(x=tau, y=acfly)
	# ACF_L_MSE_plot[[p]] <- ggplot(QACF_L_MSE_Data, aes(x=x, y=y)) + geom_line() + geom_line(data=ACF_L_MSE_Data, color='red', linetype='dotdash') + xlab(expression('percentile-'*tau)) + ylab(expression('MSE'*(beta["l"])))
	# #LP Labor
	qlply <- split(LP_Bias_MSE$QLP_MSE_L, ceiling(seq_along(LP_Bias_MSE$QLP_MSE_L)/length(tau)))[[p]]
	lply <- split(LP_Bias_MSE$LP_MSE_L, ceiling(seq_along(LP_Bias_MSE$LP_MSE_L)/length(tau)))[[p]]
	QLP_L_MSE_Data <- data.frame(x=tau, y=qlply)
	LP_L_MSE_Data <- data.frame(x=tau, y=lply)
	LP_L_MSE_plot[[p]] <- ggplot(QLP_L_MSE_Data, aes(x=x, y=y)) + geom_line() + geom_line(data=LP_L_MSE_Data, color='red', linetype='dotdash') + xlab(expression('percentile-'*tau)) + ylab(expression('MSE'*(beta["l"])))
}
#Combine plots to grid and save
# ACF_MSE_Plot_1 <- plot_grid(DGP1, plot_grid(ACF_K_MSE_plot[[1]], ACF_L_MSE_plot[[1]]), ncol=1, rel_heights=c(0.1, 1))
# ACF_MSE_Plot_2 <- plot_grid(DGP2, plot_grid(ACF_K_MSE_plot[[2]], ACF_L_MSE_plot[[2]]), ncol=1, rel_heights=c(0.1, 1))
# ACF_MSE_Plot <- plot_grid(ACF_MSE_Plot_1, ACF_MSE_Plot_2, nrow=2, align='v')
# save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/ACF_MSE_Plot.png", ACF_MSE_Plot, base_height = 8, base_width = 9)
# #LP
LP_MSE_Plot_1 <- plot_grid(DGP1, plot_grid(LP_K_MSE_plot[[1]], LP_L_MSE_plot[[1]]), ncol=1, rel_heights=c(0.1, 1))
LP_MSE_Plot_2 <- plot_grid(DGP2, plot_grid(LP_K_MSE_plot[[2]], LP_L_MSE_plot[[2]]), ncol=1, rel_heights=c(0.1, 1))
LP_MSE_Plot <- plot_grid(LP_MSE_Plot_1, LP_MSE_Plot_2, nrow=2, align='v')
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/LP_MSE_Plot.png", LP_MSE_Plot, base_height = 8, base_width = 9)










