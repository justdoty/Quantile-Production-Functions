library(stringr)
library(dplyr)
library(xtable)
#Load US dataset
USdata <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/USdata.csv') %>% 
  select(id, year, Y, K, L, M, naics3) %>% transmute(id=id, year=year, Y=log(Y), K=log(K), L=log(L), M=log(M), naics3=as.character(naics3), naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
#Industries as listed in Estimation_US.R file
NAICS <- c("31", "32", "33", "All")
H <- c(100, 10, 1, 0.01, 0.001)
############################################################################################################
#################################Load and prepare data frames for estimates#################################
############################################################################################################
#Vector of quantiles
tau <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.75, 0.8, 0.9)
#Number of parameters
dZ <- 2
R <- 500
alpha <- .1
require(abind)
#Load the results from QLP, output is a list of estimates over quantiles for each industry
QLP_results <- list()
QLP_true <- list()
for (i in 1:length(tau)){
	load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Environments/QLP_US_HTEST_Q%s.RData", i))
	QLP_results[[i]] <- results
  QLP_true[[i]] <- true.beta

}
#Prepare estimates for QLP: Outputs a list of quantile estimates over industry (columns)
#and estimates (rows)
#Bootstrapped QLP Coefficient Estimates
QLP_Boot <- array(0, dim=c(2, length(H), length(tau)))
#Matrix to store Bootstrap Bias estimates
QLP_Bias <- array(0, dim=c(2, length(H), length(tau)))
#Matrix to store logical whether bias correction is needed
QLP_Test <- array(0, dim=c(2, length(H), length(tau)))
#Estimates corrected for bias
QLP_BC <- array(0, dim=c(2, length(H), length(tau)))
#Bootstrapped Standard Errors
QLP_SE <- array(0, dim=c(2, length(H), length(tau)))
#Estimates of Returns to Scale
# QLP_RtS <- lapply(QLP_results, colSums)
QLP_RtS <- array(0, dim=c(length(tau), length(H)))
#Estimates of Standard Error Returns to Scale
QLP_RtS_SE <- array(0, dim=c(length(tau), length(H)))
#Lower and Upper bounds of CI
QLP_Lower <- array(0, dim=c(2, length(H), length(tau)))
QLP_Upper <- array(0, dim=c(2, length(H), length(tau)))
for (i in 1:length(tau)){
	for (j in 1:length(H)){
		QLP_Boot[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, mean)
    QLP_Bias[,,i][,j] <- QLP_Boot[,,i][,j]-QLP_true[[i]][,j]
		QLP_SE[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, sd)
    BC <- QLP_true[[i]][,j]-QLP_Bias[,,i][,j]
    QLP_Test[,,i][,j] <- QLP_Bias[,,i][,j]<0.25*QLP_SE[,,i][,j]
    for (z in 1:dZ){
        QLP_BC[,,i][z,j] <- if (QLP_Test[,,i][z,j]==TRUE) {QLP_true[[i]][z,j]} else {BC[z]} 
        #Bias centered student T confidence intervals
        # QLP_Lower[,,i][z,j] <- if (QLP_Test[,,i][z,j]==TRUE) {QLP_true[[i]][z,j]+qt(alpha/2, R)*QLP_SE[,,i][z,j]} else {QLP_true[[i]][z,j]+qt(alpha/2, R)*QLP_SE[,,i][z,j]}
        # QLP_Upper[,,i][z,j] <- if (QLP_Test[,,i][z,j]==TRUE) {QLP_true[[i]][z,j]+qt(1-alpha/2, R)*QLP_SE[,,i][z,j]} else {QLP_true[[i]][z,j]+qt(1-alpha/2, R)*QLP_SE[,,i][z,j]}
        #Bias Corrected Percentile Method
        # z0 <- qnorm(mean(QLP_results[[i]][,,j][,z]<QLP_true[[i]][z,j]))
        # QLP_Lower[,,i][z,j] <- if (QLP_Test[,,i][z,j]==TRUE) {quantile(QLP_results[[i]][,,j][,z], alpha/2)} else {quantile(QLP_results[[i]][,,j][,z], pnorm(2*z0+qnorm(alpha/2)))}
        # QLP_Upper[,,i][z,j] <- if (QLP_Test[,,i][z,j]==TRUE) {quantile(QLP_results[[i]][,,j][,z], 1-alpha/2)} else {quantile(QLP_results[[i]][,,j][,z], pnorm(2*z0-qnorm(alpha/2)))}
      }
    #Uncorrected Percentile Method
    QLP_Lower[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, function(x) quantile(x, alpha/2))
    QLP_Upper[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, function(x) quantile(x, 1-alpha/2))
    QLP_RtS[i,j] <- sum(QLP_Boot[,,i][,j])
    QLP_RtS_SE[i,j] <- sd(apply(QLP_results[[i]][,,j], 1, sum))
	}
}
#A little bit of reformating to obtain estimates over industries instead of quantiles
#Bootstrapped QLP Coefficient Estimates
QLP_Boot <- aperm(QLP_Boot, c(3, 1, 2))
QLP_Bias <- aperm(QLP_Bias, c(3, 1, 2))
QLP_BC <- aperm(QLP_BC, c(3, 1, 2))
QLP_Test <- aperm(QLP_Test, c(3, 1, 2))
# QLP_Coef <- aperm(array(as.numeric(unlist(QLP_Coef)), dim=c(dZ, length(NAICS), length(tau))), c(3, 1, 2))
QLP_SE <- aperm(QLP_SE, c(3, 1, 2))
QLP_RtS <- c(QLP_RtS)
QLP_RtS_SE <- c(QLP_RtS_SE)
QLP_Lower <- aperm(QLP_Lower, c(3, 1, 2))
QLP_Upper <- aperm(QLP_Upper, c(3, 1, 2))
#Listed by Industry
#For QLP
QLP_Boot <- lapply(seq(dim(QLP_Boot)[3]), function(x) QLP_Boot[ , , x])
QLP_BC <- lapply(seq(dim(QLP_BC)[3]), function(x) QLP_BC[ , , x])
QLP_SE <- lapply(seq(dim(QLP_SE)[3]), function(x) QLP_SE[ , , x])
QLP_Upper <- lapply(seq(dim(QLP_Upper)[3]), function(x) QLP_Upper[ , , x])
QLP_Lower <- lapply(seq(dim(QLP_Lower)[3]), function(x) QLP_Lower[ , , x])


#Make an estimates table for Quantile GMM
estimates <- data.frame(cbind(rep(tau, length(H)), cbind(do.call(rbind, QLP_Boot), do.call(rbind, QLP_SE))[,c(rbind(c(1:2), 2+(1:2)))]))
estimates <- cbind(estimates, QLP_RtS, QLP_RtS_SE)
colnames(estimates) <- c('Tau','K',"se_K", 'L', "se_L", 'RtS', 'RtS_SE')
#Make a Confidence Interval Table for Quantile GMM
QLP_CI <- data.frame(cbind(rep(tau, length(H)), cbind(do.call(rbind, QLP_Lower), do.call(rbind, QLP_Upper))[,c(rbind(c(1:2), 2+(1:2)))]))
colnames(QLP_CI) <- c('Tau', 'Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
#Prepare estimates for table in paper/presentation
tau_table <- c(0.1, 0.25, 0.5, 0.75)
#Table Labels
H_labels <- array(NA, length(tau_table)*length(H)); H_labels[seq(1, length(tau_table)*length(H), by=length(tau_table))] <- H
H_labels[is.na(H_labels)] <- ""

estimates_table <- cbind(H_labels, estimates[rep(tau, length(H))%in%tau_table, ])
colnames(estimates_table) <- c("Bandwidth", "$\\tau$", "Coef.", 's.e.', "Coef.", 's.e.', "Coef.", "s.e")

estimates_table <- xtable(estimates_table, digits=c(0,0,2,3,4,3,4,3,4), type="latex", caption="Coefficient Estimates and Standard Errors for US Manufacturing Firms")
align(estimates_table) <- rep('c', 9)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Returns to Scale} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}'
print(estimates_table, hline.after=c(0,nrow(estimates_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, caption.placement="top", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Estimates/US_HTEST_Estimates.tex")

############################Coefficicent Plots######################################
require(ggplot2)
require(cowplot)
require(reshape2)
#Bandwidth Plot Labels
QLP_K_plot <- list(); QLP_L_plot <- list()
for (p in 1:length(H)){
  #Capital Plot for QLP/LP######################################################
  #Capital Coefficient Values for QLP
  qlpk <- split(estimates$K, ceiling(seq_along(estimates$K)/length(tau)))[[p]]
  #Lower bound CI for QLP
  qlpklow <- split(QLP_CI$Lower_K, ceiling(seq_along(QLP_CI$Lower_K)/length(tau)))[[p]]
  #Upper bound CI for QLP
  qlpkup <- split(QLP_CI$Upper_K, ceiling(seq_along(QLP_CI$Upper_K)/length(tau)))[[p]]
  #Plotting Data
  QLPK_data <- data.frame(x=tau, y=qlpk, lower=qlpklow, upper=qlpkup)
  #Capital Plots
  QLP_K_plot[[p]] <- ggplot(QLPK_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70")+ geom_line(aes(y=y))
  #Labor Plot for QLP/LP#########################################################
  #Labor Coefficient Values for QLP
  qlpl <- split(estimates$L, ceiling(seq_along(estimates$L)/length(tau)))[[p]]
  #Lower bound CI for QLP
  qlpllow <- split(QLP_CI$Lower_L, ceiling(seq_along(QLP_CI$Lower_L)/length(tau)))[[p]]
  #Upper bound CI for QLP
  qlplup <- split(QLP_CI$Upper_L, ceiling(seq_along(QLP_CI$Upper_L)/length(tau)))[[p]]
  #Plotting Data
  QLPL_data <- data.frame(x=tau, y=qlpl, lower=qlpllow, upper=qlplup)
  #Labor Plots
  QLP_L_plot[[p]] <- ggplot(QLPL_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) 
###############################Combine Plots ##############################################
  ############################################################################################
  H_plots <- ggdraw() + draw_label(paste("H", H[p], sep=" "), fontface="plain", size=22) + theme(plot.title = element_text(hjust = 0.5))
  Coef_Plot <- plot_grid(QLP_K_plot[[p]], QLP_L_plot[[p]])
  Coef_Plot <- plot_grid(H_plots, Coef_Plot, ncol=1, align="h", rel_heights = c(0.3, 1, 1))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/Coef_Plot_H_", H[p], ".png", sep=""), Coef_Plot, base_height=8, base_width=7)
}






