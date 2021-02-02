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
#Store QLP Estimates and Standard Deviations
QLP_coef <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
QLP_se <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
#Store LP Estimates and Standard Deviations
LP_coef <- replicate(length(DGPs), list(array(0, dim=2)))
LP_se <- replicate(length(DGPs), list(array(0, dim=2)))
#Store Qdif Estimates and Standard Deviations
Qdif_coef <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
Qdif_se <- replicate(length(DGPs), list(array(0, dim=c(length(tau), dB))))
#Calculations
for (d in 1:length(DGPs)){
	LP_coef[[d]] <- round(apply(resmat_LP[,,d], 2, mean), digits=3)
	LP_se[[d]] <- round(apply(resmat_LP[,,d], 2, sd), digits=4)
	for (q in 1:length(tau)){
	  QLP_coef[[d]][q,] <- round(apply(resmat_LPQ[,,,d][,,q], 2, mean), digits=3)
	  QLP_se[[d]][q,] <- round(apply(resmat_LPQ[,,,d][,,q], 2, sd), digits=4)
	  Qdif_coef[[d]][q,] <- round(apply(resmat_Qdif[,,,d][,,q], 2, mean), digits=3)
	  Qdif_se[[d]][q,] <- round(apply(resmat_Qdif[,,,d][,,q], 2, sd), digits=4)
	}
}
#####################Create a table for all estimates#######################
QLP_estimates <- data.frame(cbind(rep(tau, length(DGPs)), cbind(do.call(rbind, QLP_coef), do.call(rbind, QLP_se))[,c(rbind(c(1:dB), dB+(1:dB)))]))
colnames(QLP_estimates) <- c('Tau','K',"sd_K", 'L', "sd_L")
LP_estimates <- data.frame(cbind(do.call(rbind, LP_coef), do.call(rbind, LP_se))[,c(rbind(c(1:dB), dB+(1:dB)))])
colnames(LP_estimates) <- c('K', 'sd_K', 'L', 'sd_L')
Qdifhat <- data.frame(cbind(rep(tau, length(DGPs)), cbind(do.call(rbind, Qdif_coef), do.call(rbind, Qdif_se))[,c(rbind(c(1:dB), dB+(1:dB)))]))
colnames(Qdifhat) <- c('Tau','QdifK',"Qdifsd_K", 'QdifL', "Qdifsd_L")
#####################Compute Bias and MSE###################################
#Bias
QLP_Bias <- round(beta-QLP_estimates[,c("K", "L")], digits=4)
LP_Bias <- round(beta-apply(LP_estimates[, c("K", "L")], 2, function(x) rep(x, each=length(tau))), digits=4)
QLP_MSE <- round(QLP_estimates[, c("sd_K", "sd_L")]^2+QLP_Bias^2, digits=4)
LP_MSE <- round(apply(LP_estimates[, c("sd_K", "sd_L")], 2, function(x) rep(x, each=length(tau)))^2+LP_Bias^2, digits=4)
Bias_MSE <- data.frame(cbind(tau, cbind(QLP_Bias, QLP_MSE)[,c(rbind(c(1:dB), dB+(1:dB)))], cbind(LP_Bias, LP_MSE)[, c(rbind(c(1:dB), dB+(1:dB)))]))
colnames(Bias_MSE) <- c('Tau', 'QLP_Bias_K', 'QLP_MSE_K', 'QLP_Bias_L', 'QLP_MSE_L', 'LP_Bias_K', 'LP_MSE_K', 'LP_Bias_L', 'LP_MSE_L')
#####################Prepare Results for Xtable###############################
DGP_labels <- array(NA, length(tau_table)*length(DGPs)); DGP_labels[seq(1, length(tau_table)*length(DGPs), by=length(tau_table))] <- c(1:length(DGPs))
DGP_labels[is.na(DGP_labels)] <- ""
#######################For Bias and MSE########################################
Bias_MSE_table <- cbind(DGP_labels, Bias_MSE[rep(tau, length(DGPs))%in%tau_table, ])
colnames(Bias_MSE_table) <- c("DGP", "$\\tau$","Bias", 'MSE', "Bias",'MSE',"Bias",'MSE', 'Bias', 'MSE')
Bias_MSE_table <- xtable(Bias_MSE_table, digits=c(0,0,2,4,4,4,4,4,4,4,4), type="latex", caption="Bias and MSE")
align(Bias_MSE_table) <- rep('c', 11)
addtorow_Bias_MSE <- list()
addtorow_Bias_MSE$pos <- list(-1)
addtorow_Bias_MSE$command <- '\\hline\\hline & & \\multicolumn{4}{c}{QLP} & \\multicolumn{4}{c}{LP} \\\\ \\cmidrule(lr){3-6} \\cmidrule(lr){7-10} \\\\& & \\multicolumn{2}{c}{$\\beta_{k}$}  & \\multicolumn{2}{c}{$\\beta_{l}$} &\\multicolumn{2}{c}{$\\beta_{k}$} & \\multicolumn{2}{c}{$\\beta_{l}$}  \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}'
print(Bias_MSE_table, hline.after=c(0,nrow(Bias_MSE_table)), add.to.row=addtorow_Bias_MSE, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, caption.placement="top", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/Bias_MSE.tex")
############################# Plots #############################################
tau_box <- c(0.1, 0.25, 0.5, 0.75, 0.9)
KboxPlot <- list(); LboxPlot <- list()
KMSEplot <- list(); LMSEplot <- list()
QdifKplot <- list(); QdifLplot <- list()
DGP1 <- ggdraw() + draw_label("DGP 1", fontface='plain')
DGP2 <- ggdraw() + draw_label("DGP 2", fontface='plain')
for (p in 1:length(DGPs)){
	#MSE Plots
	MSEdata <- data.frame(tau=tau, QK=split(Bias_MSE$QLP_MSE_K, ceiling(seq_along(Bias_MSE$QLP_MSE_K)/length(tau)))[[p]], QL=split(Bias_MSE$QLP_MSE_L, ceiling(seq_along(Bias_MSE$QLP_MSE_L)/length(tau)))[[p]],
		K=split(Bias_MSE$LP_MSE_K, ceiling(seq_along(Bias_MSE$LP_MSE_K)/length(tau)))[[p]], L=split(Bias_MSE$LP_MSE_L, ceiling(seq_along(Bias_MSE$LP_MSE_L)/length(tau)))[[p]])
	KMSEplot[[p]] <- ggplot(MSEdata, aes(x=tau)) + geom_line(aes(y=QK)) + geom_line(aes(y=K), color='red', linetype='dotdash') + xlab(expression('percentile-'*tau)) + ylab(expression('MSE'*(beta["k"])))
	LMSEplot[[p]] <- ggplot(MSEdata, aes(x=tau)) + geom_line(aes(y=QL)) + geom_line(aes(y=L), color='red', linetype='dotdash') + xlab(expression('percentile-'*tau)) + ylab(expression('MSE'*(beta["l"])))
	#Qdif Plots
	Qdifdata <- data.frame(tau=tau, K=split(Qdifhat$QdifK, ceiling(seq_along(Qdifhat$QdifK)/length(tau)))[[p]], KLB=split(Qdifhat$QdifK-Qdifhat$Qdifsd_K*qnorm(1-.05/2), ceiling(seq_along(Qdifhat$QdifK)/length(tau)))[[p]],
		KUB=split(Qdifhat$QdifK+Qdifhat$Qdifsd_K*qnorm(1-.05/2), ceiling(seq_along(Qdifhat$QdifK)/length(tau)))[[p]], L=split(Qdifhat$QdifL, ceiling(seq_along(Qdifhat$QdifL)/length(tau)))[[p]],
		LLB=split(Qdifhat$QdifL-Qdifhat$Qdifsd_L*qnorm(1-.05/2), ceiling(seq_along(Qdifhat$QdifL)/length(tau)))[[p]], LUB=split(Qdifhat$QdifL+Qdifhat$Qdifsd_L*qnorm(1-.05/2), ceiling(seq_along(Qdifhat$QdifL)/length(tau)))[[p]])
	QdifKplot[[p]] <- ggplot(Qdifdata, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_point(aes(y=K)) + geom_errorbar(aes(ymin=KLB, ymax=KUB)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(Qdifdata$LLB), max(Qdifdata$KUB)))
	QdifLplot[[p]] <- ggplot(Qdifdata, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_point(aes(y=L)) + geom_errorbar(aes(ymin=LLB, ymax=LUB)) + geom_hline(yintercept=0, linetype='dashed', color='red')+ coord_cartesian(ylim =c(min(Qdifdata$LLB), max(Qdifdata$KUB)))
	#Box Plots
	BoxData <- data.frame(tau=factor(rep(tau_box, each=dim(resmat_LPQ)[1])), K=c(resmat_LPQ[,,,p][,,match(tau_box, tau)][,1,]), L=c(resmat_LPQ[,,,p][,,match(tau_box, tau)][,2,]))
	KboxPlot[[p]] <- ggplot(BoxData, aes(x=tau, y=K)) +geom_boxplot() + xlab(expression('percentile-'*tau)) + ylab(expression(beta["k"]))
	LboxPlot[[p]] <- ggplot(BoxData, aes(x=tau, y=L)) +geom_boxplot() + xlab(expression('percentile-'*tau)) + ylab(expression(beta["l"])) 
}
####################Combine plots to grid and save##############################
#MSE Plots
MSErow1 <- plot_grid(DGP1, plot_grid(KMSEplot[[1]], LMSEplot[[1]]), ncol=1, rel_heights=c(0.1, 1))
MSErow2 <- plot_grid(DGP2, plot_grid(KMSEplot[[2]], KMSEplot[[2]]), ncol=1, rel_heights=c(0.1, 1))
MSE_Plot <- plot_grid(MSErow1, MSErow2, nrow=2, align='v')
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/MSE_Plot.png", MSE_Plot, base_height = 8, base_width = 9)
#Qdif Plots
Qdifrow1 <- plot_grid(DGP1, plot_grid(QdifKplot[[1]], QdifLplot[[1]]), ncol=1, rel_heights=c(0.1, 1))
Qdifrow2 <- plot_grid(DGP2, plot_grid(QdifKplot[[2]], QdifLplot[[2]]), ncol=1, rel_heights=c(0.1, 1))
Qdif_Plot <- plot_grid(Qdifrow1, Qdifrow2, nrow=2, align='v')
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/Qdif_Plot.png", Qdif_Plot, base_height = 8, base_width = 9)
#Box Plots
BoxRow1 <- plot_grid(DGP1, plot_grid(KboxPlot[[1]], LboxPlot[[1]]), ncol=1, rel_heights=c(0.1, 1))
BoxRow2 <- plot_grid(DGP2, plot_grid(KboxPlot[[2]], LboxPlot[[2]]), ncol=1, rel_heights=c(0.1, 1))
Box_Plot <- plot_grid(BoxRow1, BoxRow2, nrow=2, align='v')
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Monte_Carlo/Box_Plot.png", Box_Plot, base_height = 8, base_width = 9)









