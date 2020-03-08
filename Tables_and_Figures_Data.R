setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code')
require(abind)
R <- 500
tau <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.75, 0.8, 0.9)
industries <- c(311, 381, 321, 331)
#Load the results from QLP, output is a list of estimates over quantiles for each industry
QLP_results <- list()
QLP_true <- list()
load("QLP_Estimation1.RData")
for (i in 1:length(tau)){
	load(sprintf("QLP_Estimation%s.RData", i))
	QLP_results[[i]] <- results
  QLP_true[[i]] <- true.beta

}
#Load the results from LP, output is a list of estimates over industries
LP_results <- list()
LP_true <- list()
for (i in 1:length(industries)){
	load(sprintf("LP_Estimation%s.RData", i))
	LP_results[[i]] <- results
  LP_true[[i]] <- true.beta.LP
}
#Load the results from OLS, output is a list of estimates over industries
OLS_results <- load('OLS_Estimates.Rdata')

#Prepare estimates for QLP: Outputs a list of quantile estimates over industry (columns)
#and estimates (rows)
#Bootstrapped QLP Coefficient Estimates
# QLP_Coef <- array(0, dim=c(3, length(industries), length(tau)))
QLP_Coef <- lapply(QLP_true, function(x) t(x))
QLP_SE <- array(0, dim=c(3, length(industries), length(tau)))
QLP_RtS <- lapply(QLP_Coef, rowSums)
QLP_RtS_SE <- array(0, dim=c(length(tau), length(industries)))
QLP_Lower <- array(0, dim=c(3, length(industries), length(tau)))
QLP_Upper <- array(0, dim=c(3, length(industries), length(tau)))
for (i in 1:length(tau)){
	for (j in 1:length(industries)){
		# QLP_Coef[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, mean)
		QLP_SE[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, sd)
    QLP_RtS_SE[i,j] <- sd(apply(QLP_results[[i]][,,j], 1, sum))
		QLP_Lower[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, function(x) quantile(x, 0.05))
		QLP_Upper[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, function(x) quantile(x, 0.95))
	}
}
#A little bit of reformating to obtain estimates over industries instead of quantiles
#Bootstrapped QLP Coefficient Estimates
# QLP_Coef <- aperm(QLP_Coef, c(3, 1, 2))
QLP_Coef <- array(as.numeric(unlist(QLP_Coef)), dim=c(4, 3, 10))
QLP_Coef <- aperm(QLP_Coef, dim=c(3, 1, 2))
QLP_SE <- aperm(QLP_SE, c(3, 1, 2))
QLP_RtS <- c(t(array(as.numeric(unlist(QLP_RtS)), dim=c(4, 10))))
QLP_RtS_SE <- c(QLP_RtS_SE)
QLP_Lower <- aperm(QLP_Lower, c(3, 1, 2))
QLP_Upper <- aperm(QLP_Upper, c(3, 1, 2))
#Prepare estimates for LP
#Bootstrapped LP Coefficient Estimates
# LP_Coef <- array(0, dim=c(3, length(industries)))
LP_Coef <- array(as.numeric(unlist(LP_true)), dim=c(3,4))
LP_SE <- array(0, dim=c(3, length(industries)))
LP_Lower <- array(0, dim=c(3, length(industries)))
LP_Upper <- array(0, dim=c(3, length(industries)))
for (i in 1:length(industries)){
	# LP_Coef[,i] <- colMeans(LP_results[[i]]) 
	LP_SE[,i] <- apply(LP_results[[i]], 2, sd) 
	LP_Lower[,i] <- apply(LP_results[[i]], 2, function(x) quantile(x, 0.05))
	LP_Upper[,i] <- apply(LP_results[[i]], 2, function(x) quantile(x, 0.95)) 
}
#Prepare estimates for OLS
OLS_Coef <- lapply(lm.soln, function(x) coefficients(x))
OLS_Coef <- lapply(OLS_Coef, function(x) as.numeric(x[2:4]))
OLS_CI <- lapply(lm.soln, function(x) confint(x, c('chile$lnk', 'chile$lnw', 'chile$lnb'), level=0.9))
OLS_Lower <- lapply(OLS_CI, function(x) as.numeric(x[,1]))
OLS_Upper <- lapply(OLS_CI, function(x) as.numeric(x[,2]))

#I prefer to list by industry
#For QLP
QLP_Coef <- lapply(seq(dim(QLP_Coef)[3]), function(x) QLP_Coef[ , , x])
QLP_SE <- lapply(seq(dim(QLP_SE)[3]), function(x) QLP_SE[ , , x])
QLP_Upper <- lapply(seq(dim(QLP_Upper)[3]), function(x) QLP_Upper[ , , x])
QLP_Lower <- lapply(seq(dim(QLP_Lower)[3]), function(x) QLP_Lower[ , , x])
#For LP
LP_Coef <- lapply(seq(dim(LP_Coef)[2]), function(x) LP_Coef[,x])
LP_SE <- lapply(seq(dim(LP_SE)[2]), function(x) LP_SE[,x])
LP_Upper <- lapply(seq(dim(LP_Upper)[2]), function(x) LP_Upper[,x])
LP_Lower <- lapply(seq(dim(LP_Lower)[2]), function(x) LP_Lower[,x])


#Make an estimates table for Quantile GMM
estimates <- data.frame(cbind(rep(tau, length(industries)), cbind(do.call(rbind, QLP_Coef), do.call(rbind, QLP_SE))[,c(rbind(c(1:3), 3+(1:3)))]))
estimates <- cbind(estimates, QLP_RtS, QLP_RtS_SE)
colnames(estimates) <- c('Tau','K',"se_K", 'Lw', "se_Lw", 'Lb', "se_Lb", 'RtS', 'RtS_SE')
#Make a Confidence Interval Table for Quantile GMM
QLP_CI <- data.frame(cbind(rep(tau, length(industries)), cbind(do.call(rbind, QLP_Lower), do.call(rbind, QLP_Upper))[,c(rbind(c(1:3), 3+(1:3)))]))
colnames(QLP_CI) <- c('Tau', 'Lower_K', 'Upper_K', 'Lower_Lw', 'Upper_Lw', 'Lower_Lb', 'Upper_Lb')
#Make an estimates table for LP
estimates_LP <- data.frame(cbind(do.call(rbind, LP_Coef), do.call(rbind, LP_SE))[,c(rbind(c(1:3), 3+(1:3)))])
colnames(estimates_LP) <- c('K', 'se_K', 'Lw', 'se_Lw', 'Lb', 'se_Lb')
#Make a Confidence Interval Table for Quantile GMM
LP_CI <- data.frame(cbind(do.call(rbind, LP_Lower), do.call(rbind, LP_Upper))[,c(rbind(c(1:3), 3+(1:3)))])
colnames(LP_CI) <- c('Lower_K', 'Upper_K', 'Lower_Lw', 'Upper_Lw', 'Lower_Lb', 'Upper_Lb')
#Make an estimates table for LP
estimates_OLS <- data.frame(do.call(rbind, OLS_Coef))
colnames(estimates_OLS) <- c('K', 'Lw', 'Lb')
#Make a Confidence Interval Table for Quantile GMM
OLS_CI <- data.frame(cbind(do.call(rbind, OLS_Lower), do.call(rbind, OLS_Upper))[,c(rbind(c(1:3), 3+(1:3)))])
colnames(OLS_CI) <- c('Lower_K', 'Upper_K', 'Lower_Lw', 'Upper_Lw', 'Lower_Lb', 'Upper_Lb')
#Prepare estimates for table in paper/presentation
require(xtable)
tau_table <- c(0.1, 0.25, 0.5, 0.75)

#Table Labels
ISIC_labels <- array(NA, length(tau_table)*length(industries)); ISIC_labels[seq(1, length(tau_table)*length(industries), by=length(tau_table))] <- industries
ISIC_labels[is.na(ISIC_labels)] <- ""

estimates_table <- cbind(ISIC_labels, estimates[rep(tau, length(industries))%in%tau_table, ])
colnames(estimates_table) <- c("Industry (ISIC code)", "$\\tau$", "Coef.", 's.e.', "Coef.",'s.e.', "Coef.", 's.e.', "Coef.", "s.e")

estimates_table <- xtable(estimates_table, digits=c(0,0,2,3,4,3,4,3,4,3,4))
align(estimates_table) <- rep('c', 11)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Skilled Labor} & \\multicolumn{2}{c}{Unskilled Labor} & \\multicolumn{2}{c}{Returns to Scale} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}'
print(estimates_table, hline.after=c(0,nrow(estimates_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x)

############################Coefficicent Plots######################################
require(ggplot2)
require(cowplot)
#The true values of coefficients from LP stata code
kcoef <- c(.16, .31, .125, .17)
kse <- c(0.0364, 0.0495, 0.077, 0.0421)
wcoef <- c(.259, .404, .478, .389)
wse <- c(0.0113, 0.0217, 0.0243, 0.0224)
bcoef <- c(.444, .455, .383, .409)
bse <- c(0.0121, 0.0212, 0.0276, 0.0262)
K_plot <- list(); Lw_plot <- list(); Lb_plot <- list()
for (p in 1:length(industries)){
  #Capital
  ky <- split(estimates$K, ceiling(seq_along(estimates$K)/length(tau)))[[p]]
  klow <- split(QLP_CI$Lower_K, ceiling(seq_along(QLP_CI$Lower_K)/length(tau)))[[p]]
  klow_LP <- OLS_CI$Lower_K[p]
  # klow_LP <- kcoef[p]-qnorm(0.95)*kse[p]
  kup <- split(QLP_CI$Upper_K, ceiling(seq_along(QLP_CI$Upper_K)/length(tau)))[[p]]
  kup_LP <- OLS_CI$Upper_K[p]
  # kup_LP <- kcoef[p]+qnorm(0.95)*kse[p]
  K_data <- data.frame(x=tau, y=ky, z=estimates_OLS$K[p], lower=klow, upper=kup, lower_LP=klow_LP, upper_LP=kup_LP)
  K_plot[[p]] <- ggplot(K_data, aes(x=x, y=y)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) +  geom_hline(yintercept=K_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(K_data$lower_LP, K_data$upper_LP), linetype='dashed', color='red')
  #Skilled Labor
  lwy <- split(estimates$Lw, ceiling(seq_along(estimates$Lw)/length(tau)))[[p]]
  lwlow <- split(QLP_CI$Lower_Lw, ceiling(seq_along(QLP_CI$Lower_Lw)/length(tau)))[[p]]
  lwlow_LP <- OLS_CI$Lower_Lw[p]
  # lwlow_LP <- wcoef[p]-qnorm(0.95)*wse[p]
  lwup <- split(QLP_CI$Upper_Lw, ceiling(seq_along(QLP_CI$Upper_Lw)/length(tau)))[[p]]
  lwup_LP <- OLS_CI$Upper_Lw[p]
  # lwup_LP <- wcoef[p]+qnorm(0.95)*wse[p]
  Lw_data <- data.frame(x=tau, y=lwy, z=estimates_OLS$Lw[p], lower=lwlow, upper=lwup, lower_LP=lwlow_LP, upper_LP=lwup_LP)
  Lw_plot[[p]] <- ggplot(Lw_data, aes(x=x, y=y)) + xlab(expression('percentile-'*tau)) + ylab("Skilled Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) +  geom_hline(yintercept=Lw_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(Lw_data$lower_LP, Lw_data$upper_LP), linetype='dashed', color='red')
  #Unskilled Labor
  lby <- split(estimates$Lb, ceiling(seq_along(estimates$Lb)/length(tau)))[[p]]
  lblow <- split(QLP_CI$Lower_Lb, ceiling(seq_along(QLP_CI$Lower_Lb)/length(tau)))[[p]]
  lblow_LP <- OLS_CI$Lower_Lb[p]
  # lblow_LP <- bcoef[p]-qnorm(0.95)*bse[p]
  lbup <- split(QLP_CI$Upper_Lb, ceiling(seq_along(QLP_CI$Upper_Lb)/length(tau)))[[p]]
  lbup_LP <- OLS_CI$Upper_Lb[p]
  # lbup_LP <- bcoef[p]+qnorm(0.95)*bse[p]
  Lb_data <- data.frame(x=tau, y=lby, z=estimates_OLS$Lb[p], lower=lblow, upper=lbup, lower_LP=lblow_LP, upper_LP=lbup_LP)
  Lb_plot[[p]] <- ggplot(Lb_data, aes(x=x, y=y)) + xlab(expression('percentile-'*tau)) + ylab("Unskilled Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) +  geom_hline(yintercept=Lb_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(Lb_data$lower_LP, Lb_data$upper_LP), linetype='dashed', color='red')
}
#Combine plots across industries
#Industry ISIC Code Plot Labels
Industry_1 <- ggdraw() + draw_label("Food Products (311)", fontface='plain')
Industry_2 <- ggdraw() + draw_label("Fabricated Metals (381)", fontface='plain')
Industry_3 <- ggdraw() + draw_label("Textiles (321)", fontface='plain')
Industry_4 <- ggdraw() + draw_label("Wood Products (331)", fontface='plain')

#Coefficient Plot for 1st industry
L_Plot_1 <- plot_grid(Lw_plot[[1]], Lb_plot[[1]], rel_heights=c(0.1, 1), ncol=2)
K_Plot_1 <- plot_grid(NULL, K_plot[[1]], NULL, ncol=3, rel_widths=c(0.25,0.5,0.25))
Coefficient_Plot_1 <- plot_grid(Industry_1, L_Plot_1, K_Plot_1, ncol=1, rel_heights=c(0.1, 1, 1))
save_plot("Coefficient_Plot_1.png", Coefficient_Plot_1, base_height = 8, base_width = 9)

#Coefficient Plot for 2nd industry
L_Plot_2 <- plot_grid(Lw_plot[[2]], Lb_plot[[2]], rel_heights=c(0.1, 1), ncol=2)
K_Plot_2 <- plot_grid(NULL, K_plot[[2]], NULL, ncol=3, rel_widths=c(0.25,0.5,0.25))
Coefficient_Plot_2 <- plot_grid(Industry_2, L_Plot_2, K_Plot_2, ncol=1, rel_heights=c(0.1, 1, 1))
save_plot("Coefficient_Plot_2.png", Coefficient_Plot_2, base_height = 8, base_width = 9)

#Coefficient Plot for 3rd industry
L_Plot_3 <- plot_grid(Lw_plot[[3]], Lb_plot[[3]], rel_heights=c(0.1, 1), ncol=2)
K_Plot_3 <- plot_grid(NULL, K_plot[[3]], NULL, ncol=3, rel_widths=c(0.25,0.5,0.25))
Coefficient_Plot_3 <- plot_grid(Industry_3, L_Plot_3, K_Plot_3, ncol=1, rel_heights=c(0.1, 1, 1))
save_plot("Coefficient_Plot_3.png", Coefficient_Plot_3, base_height = 8, base_width = 9)

#Coefficient Plot for 3rd industry
L_Plot_4 <- plot_grid(Lw_plot[[4]], Lb_plot[[4]], rel_heights=c(0.1, 1), ncol=2)
K_Plot_4 <- plot_grid(NULL, K_plot[[4]], NULL, ncol=3, rel_widths=c(0.25,0.5,0.25))
Coefficient_Plot_4 <- plot_grid(Industry_4, L_Plot_4, K_Plot_4, ncol=1, rel_heights=c(0.1, 1, 1))
save_plot("Coefficient_Plot_4.png", Coefficient_Plot_4, base_height = 8, base_width = 9)

















