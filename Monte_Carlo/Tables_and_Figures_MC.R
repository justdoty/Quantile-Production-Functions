setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo')
require(xtable)
require(ggplot2)
require(cowplot)
require(abind)
############################ For LP Simulations###################################################
load("simulation_LP.Rdata")
#####################Compute Estimates and Standard Errors##################
#Store QLP Estimates and Standard Deviations
QLP_coef <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
QLP_se <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
#Store QR-LP Estimates and Standard Deviations
QRLP_coef <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
QRLP_se <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
#Calculations
for (d in 1:length(DGPs)){
	for (q in 1:length(tau)){
	  QLP_coef[[d]][q,] <- round(apply(resmat_LPQ[,,,d][,,q], 2, mean), digits=3)
	  QLP_se[[d]][q,] <- round(apply(resmat_LPQ[,,,d][,,q], 2, sd), digits=4)
	  QRLP_coef[[d]][q,] <- round(apply(resmat_QR[,,,d][,,q], 2, mean), digits=3)
	  QRLP_se[[d]][q,] <- round(apply(resmat_QR[,,,d][,,q], 2, sd), digits=4)
	}
}
#####################Create a table for all estimates#######################
QLP_estimates <- data.frame(cbind(rep(tau, length(DGPs)), cbind(do.call(rbind, QLP_coef), do.call(rbind, QLP_se))[,c(rbind(c(1:dB), dB+(1:dB)))]))
colnames(QLP_estimates) <- c('Tau','K',"sd_K", 'L', "sd_L")
QRLPhat <- data.frame(cbind(rep(tau, length(DGPs)), cbind(do.call(rbind, QRLP_coef), do.call(rbind, QRLP_se))[,c(rbind(c(1:dB), dB+(1:dB)))]))
colnames(QRLPhat) <- c('Tau','QRLP_K',"QRLPsd_K", 'QRLP_L', "QRLPsd_L")
#####################Compute Bias and MSE###################################
#Bias
QLP_Bias <- round(beta-QLP_estimates[,c("K", "L")], digits=4)
QR_LP_Bias <- round(beta-QRLPhat[,c("QRLP_K", "QRLP_L")], digits=4)
#MSE
QLP_MSE <- round(QLP_estimates[, c("sd_K", "sd_L")]^2+QLP_Bias^2, digits=4)
QR_LP_MSE <- round(QRLPhat[, c("QRLPsd_K", "QRLPsd_L")]^2+QR_Bias^2, digits=4)
#Combined
QLP_Bias_MSE <- data.frame(cbind(tau, cbind(QLP_Bias, QLP_MSE)[,c(rbind(c(1:dB), dB+(1:dB)))], cbind(QR_LP_Bias, QR_LP_MSE)[, c(rbind(c(1:dB), dB+(1:dB)))]))
colnames(QLP_Bias_MSE) <- c('Tau', 'QBias_K', 'QMSE_K', 'QBias_L', 'QMSE_L', 'QRBias_K', 'QRMSE_K', 'QRBias_L', 'QRMSE_L')
############################ For ACF Simulations###################################################
load("simulation_ACF.Rdata")
#####################Compute Estimates and Standard Errors##################
#Store QACF Estimates and Standard Deviations
QACF_coef <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
QACF_se <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
#Store QR-ACF Estimates and Standard Deviations
QRACF_coef <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
QRACF_se <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
#Calculations
for (d in 1:length(DGPs)){
	for (q in 1:length(tau)){
	  QACF_coef[[d]][q,] <- round(apply(resmat_ACFQ[,,,d][,,q], 2, mean), digits=3)
	  QACF_se[[d]][q,] <- round(apply(resmat_ACFQ[,,,d][,,q], 2, sd), digits=4)
	  QRACF_coef[[d]][q,] <- round(apply(resmat_QR[,,,d][,,q], 2, mean), digits=3)
	  QRACF_se[[d]][q,] <- round(apply(resmat_QR[,,,d][,,q], 2, sd), digits=4)
	}
}
#####################Create a table for all estimates#######################
QACF_estimates <- data.frame(cbind(rep(tau, length(DGPs)), cbind(do.call(rbind, QACF_coef), do.call(rbind, QACF_se))[,c(rbind(c(1:dB), dB+(1:dB)))]))
colnames(QACF_estimates) <- c('Tau','K',"sd_K", 'L', "sd_L")
QRACFhat <- data.frame(cbind(rep(tau, length(DGPs)), cbind(do.call(rbind, QRACF_coef), do.call(rbind, QRACF_se))[,c(rbind(c(1:dB), dB+(1:dB)))]))
colnames(QRACFhat) <- c('Tau','QRACF_K',"QRACFsd_K", 'QRACF_L', "QRACFsd_L")
#####################Compute Bias and MSE###################################
#Bias
QACF_Bias <- round(beta-QACF_estimates[,c("K", "L")], digits=4)
QR_ACF_Bias <- round(beta-QRACFhat[,c("QRACF_K", "QRACF_L")], digits=4)
#MSE
QACF_MSE <- round(QACF_estimates[, c("sd_K", "sd_L")]^2+QACF_Bias^2, digits=4)
QR_ACF_MSE <- round(QRACFhat[, c("QRACFsd_K", "QRACFsd_L")]^2+QR_Bias^2, digits=4)
#Combined
QACF_Bias_MSE <- data.frame(cbind(tau, cbind(QACF_Bias, QACF_MSE)[,c(rbind(c(1:dB), dB+(1:dB)))], cbind(QR_ACF_Bias, QR_ACF_MSE)[, c(rbind(c(1:dB), dB+(1:dB)))]))
colnames(QACF_Bias_MSE) <- c('Tau', 'QBias_K', 'QMSE_K', 'QBias_L', 'QMSE_L', 'QRBias_K', 'QRMSE_K', 'QRBias_L', 'QRMSE_L')
##############################################################################
#####################Prepare Results for Xtable###############################
##############################################################################
LP_DGP_labels <- array(NA, length(tau)*length(DGPs)); LP_DGP_labels[seq(1, length(tau)*length(DGPs), by=length(tau))] <- paste("LP", c(1:length(DGPs)))
ACF_DGP_labels <- array(NA, length(tau)*length(DGPs)); ACF_DGP_labels[seq(1, length(tau)*length(DGPs), by=length(tau))] <- paste("ACF", c(1:length(DGPs)))
DGP_labels <- c(LP_DGP_labels, ACF_DGP_labels)
DGP_labels[is.na(DGP_labels)] <- ""
#######################For Bias and MSE########################################
Bias_MSE <- rbind(QLP_Bias_MSE, QACF_Bias_MSE)
Bias_MSE_table <- cbind(DGP_labels, Bias_MSE)
colnames(Bias_MSE_table) <- c("DGP", "$\\tau$","Bias", 'MSE', "Bias",'MSE', "Bias", 'MSE', "Bias",'MSE')
Bias_MSE_table <- xtable(Bias_MSE_table, digits=c(0,0,2,4,4,4,4,4,4,4,4), type="latex")
align(Bias_MSE_table) <- rep('c', 11)
addtorow_Bias_MSE <- list()
addtorow_Bias_MSE$pos <- list(-1)
addtorow_Bias_MSE$command <- '\\hline\\hline & & & \\multicolumn{2}{c}{DS} & & & \\multicolumn{2}{c}{QR} \\\\ \\cmidrule(lr){3-6} \\cmidrule(lr){7-10} \\\\ & & \\multicolumn{2}{c}{Capital} & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Capital} & \\multicolumn{2}{c}{Labor} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}'
#For copy and pasting directly into latex file
print(Bias_MSE_table, hline.after=c(0, length(tau)*length(DGPs), nrow(Bias_MSE_table)), add.to.row=addtorow_Bias_MSE, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#For saving to file, more convenient, awful predefined formatting
print(Bias_MSE_table, hline.after=c(0, length(tau)*length(DGPs), nrow(Bias_MSE_table)), add.to.row=addtorow_Bias_MSE, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/Bias_MSE.tex")




