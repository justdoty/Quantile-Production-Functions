library(stringr)
library(dplyr)
library(xtable)
#Load US dataset
USdata <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/USdata.csv') %>% 
  select(id, year, Y, K, L, M, naics3) %>% transmute(id=id, year=year, Y=log(Y), K=log(K), L=log(L), M=log(M), naics3=as.character(naics3), naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
#Industries as listed in Estimation_US.R file
NAICS <- c("31", "32", "33", "All")
########################################################################################################
##########################################Summary Statistics############################################
########################################################################################################
#Create table for all relevant summary statistics
sumNAICS <- group_by(USdata, naics2) %>% summarise_at(c("Y", "K", "L", "M"), list(Q1=~quantile(., 0.25), med=median, Q3=~quantile(.,0.75), mean=mean, sd=sd), na.rm=TRUE) 
sumALL <- cbind("All", summarise_at(USdata, c("Y", "K", "L", "M"), list(Q1=~quantile(., 0.25), med=median, Q3=~quantile(.,0.75), mean=mean, sd=sd), na.rm=TRUE))
colnames(sumALL)[1] <- "naics2"
sizeNAICS <- group_by(USdata, naics2) %>% summarise(Firms=length(unique(id)), Total=n())
sizeALL <- c("All", sum(sizeNAICS$Firms), sum(sizeNAICS$Total))
size <- rbind(sizeNAICS, sizeALL)
sumstat <- round(matrix(as.numeric(as.matrix(rbind(sumNAICS, sumALL))[,-1]), nrow=16, ncol=5), 2)
#Some pretty formatting
NAICS_labels <- array(NA, 4*length(NAICS)); NAICS_labels[seq(1, 4*length(NAICS), by=4)] <- paste(NAICS, paste("(N=", size$Total, ")", sep=""))
NAICS_labels[is.na(NAICS_labels)] <- ""
summary_table <- cbind(NAICS_labels, rep(c("Output", "Capital", "Labor", "Materials"), 4), sumstat)
colnames(summary_table) <- c("Industry (NAICS code)", " ", "1st Qu.", 'Median', "3rd Qu.", 'Mean', "sd")

summary_table <- xtable(summary_table, digits=c(2,2,0,4,4,2,2,2), type="latex", caption="Summary Statistics for the US Compustat Manufacturing Industries")
align(summary_table) <- rep('c', 8)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline '
print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, caption.placement="top", table.placement="H", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Estimates/US_Summary.tex")
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
	load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Environments/QLP_US_Q%s.RData", i))
	QLP_results[[i]] <- results
  QLP_true[[i]] <- true.beta

}
#Load the results from LP, output is a list of estimates over industries
LP_results <- list()
LP_true <- list()
for (i in 1:length(NAICS)){
	load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Environments/LP_US_NAICS_%s.RData", i))
	LP_results[[i]] <- results
  LP_true[[i]] <- true.beta.LP
}
#Load the results from OLS and QR
QR_OLS_results <- load('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Environments/US_QR_OLS.Rdata')
#Prepare estimates for QLP: Outputs a list of quantile estimates over industry (columns)
#and estimates (rows)
#Bootstrapped QLP Coefficient Estimates
QLP_Boot <- array(0, dim=c(2, length(NAICS), length(tau)))
#Matrix to store Bootstrap Bias estimates
QLP_Bias <- array(0, dim=c(2, length(NAICS), length(tau)))
#Matrix to store logical whether bias correction is needed
QLP_Test <- array(0, dim=c(2, length(NAICS), length(tau)))
#Estimates corrected for bias
QLP_BC <- array(0, dim=c(2, length(NAICS), length(tau)))
#Bootstrapped Standard Errors
QLP_SE <- array(0, dim=c(2, length(NAICS), length(tau)))
#Estimates of Returns to Scale
# QLP_RtS <- lapply(QLP_results, colSums)
QLP_RtS <- array(0, dim=c(length(tau), length(NAICS)))
#Estimates of Standard Error Returns to Scale
QLP_RtS_SE <- array(0, dim=c(length(tau), length(NAICS)))
#Lower and Upper bounds of CI
QLP_Lower <- array(0, dim=c(2, length(NAICS), length(tau)))
QLP_Upper <- array(0, dim=c(2, length(NAICS), length(tau)))
for (i in 1:length(tau)){
	for (j in 1:length(NAICS)){
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
#Prepare estimates for LP
#Bootstrapped LP Coefficient Estimates
LP_Coef <- array(0, dim=c(2, length(NAICS)))
# LP_Coef <- array(as.numeric(unlist(LP_true)), dim=c(dZ,length(NAICS)))
LP_SE <- array(0, dim=c(2, length(NAICS)))
LP_Lower <- array(0, dim=c(2, length(NAICS)))
LP_Upper <- array(0, dim=c(2, length(NAICS)))
for (i in 1:length(NAICS)){
	LP_Coef[,i] <- colMeans(LP_results[[i]]) 
	LP_SE[,i] <- apply(LP_results[[i]], 2, sd) 
	LP_Lower[,i] <- apply(LP_results[[i]], 2, function(x) quantile(x, 0.05))
	LP_Upper[,i] <- apply(LP_results[[i]], 2, function(x) quantile(x, 0.95)) 
}
#Listed by Industry
#For QLP
QLP_Boot <- lapply(seq(dim(QLP_Boot)[3]), function(x) QLP_Boot[ , , x])
QLP_BC <- lapply(seq(dim(QLP_BC)[3]), function(x) QLP_BC[ , , x])
QLP_SE <- lapply(seq(dim(QLP_SE)[3]), function(x) QLP_SE[ , , x])
QLP_Upper <- lapply(seq(dim(QLP_Upper)[3]), function(x) QLP_Upper[ , , x])
QLP_Lower <- lapply(seq(dim(QLP_Lower)[3]), function(x) QLP_Lower[ , , x])

#For LP
LP_Coef <- lapply(seq(dim(LP_Coef)[2]), function(x) LP_Coef[,x])
LP_SE <- lapply(seq(dim(LP_SE)[2]), function(x) LP_SE[,x])
LP_Upper <- lapply(seq(dim(LP_Upper)[2]), function(x) LP_Upper[,x])
LP_Lower <- lapply(seq(dim(LP_Lower)[2]), function(x) LP_Lower[,x])
#Prepare estimates for OLS and QR
OLS_Coef <- lm.coef
OLS_Lower <- lm.CI[,seq(1, ncol(lm.CI), by=2)]
OLS_Upper <- lm.CI[,seq(2, ncol(lm.CI), by=2)]
QR_Coef <- lapply(seq(dim(qr.coef)[3]), function(x) qr.coef[ , , x])
QR_Lower <- lapply(seq(dim(qr.CI)[3]), function(x) qr.CI[ , , x][,seq(1, ncol(qr.CI), by=2)])
QR_Upper <- lapply(seq(dim(qr.CI)[3]), function(x) qr.CI[ , , x][,seq(2, ncol(qr.CI), by=2)])

#Make an estimates table for Quantile GMM
estimates <- data.frame(cbind(rep(tau, length(NAICS)), cbind(do.call(rbind, QLP_Boot), do.call(rbind, QLP_SE))[,c(rbind(c(1:2), 2+(1:2)))]))
estimates <- cbind(estimates, QLP_RtS, QLP_RtS_SE)
colnames(estimates) <- c('Tau','K',"se_K", 'L', "se_L", 'RtS', 'RtS_SE')
#Make a Confidence Interval Table for Quantile GMM
QLP_CI <- data.frame(cbind(rep(tau, length(NAICS)), cbind(do.call(rbind, QLP_Lower), do.call(rbind, QLP_Upper))[,c(rbind(c(1:2), 2+(1:2)))]))
colnames(QLP_CI) <- c('Tau', 'Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
#Make an estimates table for LP
estimates_LP <- data.frame(cbind(do.call(rbind, LP_Coef), do.call(rbind, LP_SE))[,c(rbind(c(1:2), 2+(1:2)))])
colnames(estimates_LP) <- c('K', 'se_K', 'L', 'se_L')
#Make a Confidence Interval Table for LP
LP_CI <- data.frame(cbind(do.call(rbind, LP_Lower), do.call(rbind, LP_Upper))[,c(rbind(c(1:2), 2+(1:2)))])
colnames(LP_CI) <- c('Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
#Make an estimates table for OLS
estimates_OLS <- data.frame(OLS_Coef)
colnames(estimates_OLS) <- c('K', 'L')
#Make a Confidence Interval Table for OLS
OLS_CI <- data.frame(cbind(OLS_Lower, OLS_Upper)[,c(rbind(c(1:2), 2+(1:2)))])
colnames(OLS_CI) <- c('Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
#Make an estimates table for QR
estimates_QR <- data.frame(cbind(rep(tau, length(NAICS)), do.call(rbind, QR_Coef)))
colnames(estimates_QR) <- c("Tau", "K", "L")
#Make a Confidence Interval Table for Quantile Regression
QR_CI <- data.frame(cbind(rep(tau, length(NAICS)), cbind(do.call(rbind, QR_Lower), do.call(rbind, QR_Upper))[,c(rbind(c(1:2), 2+(1:2)))]))
colnames(QR_CI) <- c('Tau', 'Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
#Prepare estimates for table in paper/presentation
tau_table <- c(0.1, 0.25, 0.5, 0.75)

#Table Labels
NAICS_labels <- array(NA, length(tau_table)*length(NAICS)); NAICS_labels[seq(1, length(tau_table)*length(NAICS), by=length(tau_table))] <- NAICS
NAICS_labels[is.na(NAICS_labels)] <- ""

estimates_table <- cbind(NAICS_labels, estimates[rep(tau, length(NAICS))%in%tau_table, ])
colnames(estimates_table) <- c("Industry (NAICS code)", "$\\tau$", "Coef.", 's.e.', "Coef.", 's.e.', "Coef.", "s.e")

estimates_table <- xtable(estimates_table, digits=c(0,0,2,3,4,3,4,3,4), type="latex", caption="Coefficient Estimates and Standard Errors for US Manufacturing Firms")
align(estimates_table) <- rep('c', 9)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Returns to Scale} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}'
print(estimates_table, hline.after=c(0,nrow(estimates_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, caption.placement="top", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Estimates/US_Estimates.tex")

############################Coefficicent Plots######################################
require(ggplot2)
require(cowplot)
require(reshape2)
#Industry NAICS Code Plot Labels
QLP_K_plot <- list(); QLP_L_plot <- list()
QR_K_plot <- list(); QR_L_plot <- list()
for (p in 1:length(NAICS)){
  #Capital Plot for QLP/LP######################################################
  #Capital Coefficient Values for QLP
  qlpk <- split(estimates$K, ceiling(seq_along(estimates$K)/length(tau)))[[p]]
  #Lower bound CI for QLP
  qlpklow <- split(QLP_CI$Lower_K, ceiling(seq_along(QLP_CI$Lower_K)/length(tau)))[[p]]
  #Lower bound CI for LP
  lpklow <- LP_CI$Lower_K[p]
  #Upper bound CI for QLP
  qlpkup <- split(QLP_CI$Upper_K, ceiling(seq_along(QLP_CI$Upper_K)/length(tau)))[[p]]
  #Upper bound CI for LP
  lpkup <- LP_CI$Upper_K[p]
  #Plotting Data
  QLPK_data <- data.frame(x=tau, y=qlpk, z=estimates_LP$K[p], lower=qlpklow, upper=qlpkup, lower_LP=lpklow, upper_LP=lpkup)
  #Capital Plot for QR/OLS######################################################
  qrk <- split(estimates_QR$K, ceiling(seq_along(estimates_QR$K)/length(tau)))[[p]]
  #Lower bound CI for QR
  qrklow <- split(QR_CI$Lower_K, ceiling(seq_along(QR_CI$Lower_K)/length(tau)))[[p]]
  #Lower bound CI for OLS
  lmklow <- OLS_CI$Lower_K[p]
  #Upper bound CI for QR
  qrkup <- split(QR_CI$Upper_K, ceiling(seq_along(QR_CI$Upper_K)/length(tau)))[[p]]
  #Upper bound CI for OLS
  lmkup <- OLS_CI$Upper_K[p]
  #Plotting Data
  QRK_data <- data.frame(x=tau, y=qrk, z=estimates_OLS$K[p], lower=qrklow, upper=qrkup, lower_OLS=lmklow, upper_OLS=lmkup)
  #Capital Plots
  QLP_K_plot[[p]] <- ggplot(QLPK_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QLPK_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QLPK_data$lower_LP, QLPK_data$upper_LP), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(qlpklow, qrklow), max(qlpkup, qrkup)))
  QR_K_plot[[p]] <- ggplot(QRK_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QRK_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QRK_data$lower_OLS, QRK_data$upper_OLS), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(qlpklow, qrklow), max(qlpkup, qrkup)))
  #Labor Plot for QLP/LP#########################################################
  #Labor Coefficient Values for QLP
  qlpl <- split(estimates$L, ceiling(seq_along(estimates$L)/length(tau)))[[p]]
  #Lower bound CI for QLP
  qlpllow <- split(QLP_CI$Lower_L, ceiling(seq_along(QLP_CI$Lower_L)/length(tau)))[[p]]
  #Lower bound CI for LP
  lpllow <- LP_CI$Lower_L[p]
  #Upper bound CI for QLP
  qlplup <- split(QLP_CI$Upper_L, ceiling(seq_along(QLP_CI$Upper_L)/length(tau)))[[p]]
  #Upper bound CI for LP
  lplup <- LP_CI$Upper_L[p]
  #Plotting Data
  QLPL_data <- data.frame(x=tau, y=qlpl, z=estimates_LP$L[p], lower=qlpllow, upper=qlplup, lower_LP=lpllow, upper_LP=lplup)
  #Labor Plot for QR/OLS######################################################
  qrl <- split(estimates_QR$L, ceiling(seq_along(estimates_QR$L)/length(tau)))[[p]]
  #Lower bound CI for QR
  qrllow <- split(QR_CI$Lower_L, ceiling(seq_along(QR_CI$Lower_L)/length(tau)))[[p]]
  #Lower bound CI for OLS
  lmllow <- OLS_CI$Lower_L[p]
  #Upper bound CI for QR
  qrlup <- split(QR_CI$Upper_L, ceiling(seq_along(QR_CI$Upper_L)/length(tau)))[[p]]
  #Upper bound CI for OLS
  lmlup <- OLS_CI$Upper_L[p]
  #Plotting Data
  QRL_data <- data.frame(x=tau, y=qrl, z=estimates_OLS$L[p], lower=qrllow, upper=qrlup, lower_OLS=lmllow, upper_OLS=lmlup)
  #Labor Plots
  QLP_L_plot[[p]] <- ggplot(QLPL_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QLPL_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QLPL_data$lower_LP, QLPL_data$upper_LP), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(qlpllow, qrllow), max(qlplup, qrlup)))
  QR_L_plot[[p]] <- ggplot(QRL_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QRL_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QRL_data$lower_OLS, QRL_data$upper_OLS), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(qlpllow, qrllow), max(qlplup, qrlup)))
###############################Combine Plots ##############################################
  ############################################################################################
  NAICS_plots <- ggdraw() + draw_label(paste("NAICS", NAICS[p], sep=" "), fontface="plain", size=22) + theme(plot.title = element_text(hjust = 0.5))
  Lrow <- plot_grid(QLP_L_plot[[p]], QR_L_plot[[p]])
  Krow <- plot_grid(QLP_K_plot[[p]], QR_K_plot[[p]])
  Coef_Plot <- plot_grid(NAICS_plots, Lrow, Krow, ncol=1, align="h", rel_heights = c(0.3, 1, 1))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/Coef_Plot_NAICS_", NAICS[p], ".png", sep=""), Coef_Plot, base_height=8, base_width=7)
}
##################################Prepare Plots over Time##############################
###################Coefficients over Time #######################################
tau_t <- c(0.1, 0.3, 0.5, 0.7, 0.9)
T <- 5
dZ <- 2
time <- seq(min(USdata$year), max(USdata$year), by=T)
QLPT_Coef <- array(0, dim=c(length(tau_t), dZ, length(time)))
QLPT_True <- array(0, dim=c(length(tau_t), dZ, length(time)))
for (t in 1:length(time)){
  for (q in 1:length(tau_t)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Environments/QLPT_US_Q%s.RData", q))
    QLPT_Coef[,,t][q,] <- colMeans(results_T[,,t])
    QLPT_True[,,t][q,] <- true.beta_T[,t]
  }
}
pcolour <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
KT <- data.frame(cbind(time, t(QLPT_True[,1,][,])))
colnames(KT) <- c("Year", paste("Q", tau_t, sep=" "))
KT <- melt(KT, "Year")
KTplot <- ggplot(KT, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Capital") + scale_colour_manual(name=expression(tau), labels=tau_t, values = pcolour)
LT <- data.frame(cbind(time, t(QLPT_True[,2,][,])))
colnames(LT) <- c("Year", paste("Q", tau_t, sep=""))
LT <- melt(LT, "Year")
LTplot <- ggplot(LT, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Labor") + scale_colour_manual(name=expression(tau), labels=tau_t, values = pcolour)
Plot_Title <- ggdraw() + draw_label("Output Elasticities Over Time", fontface="plain", size=22) 
Time_Plot <- plot_grid(Plot_Title, plot_grid(LTplot, KTplot), ncol=1, rel_heights = c(0.3, 1))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/Time_Plot.png", Time_Plot, base_height=6, base_width=10)
###############TFP Over Time#######################################################
estimates_TFP <- estimates[tau%in%tau_t,]
All_NAICS_QLP <- data.frame(estimates_TFP[(nrow(estimates_TFP)-(length(tau_t)-1)):nrow(estimates_TFP), c(2,4)])
All_NAICS_LP <- estimates_LP[nrow(estimates_LP), c(1,3)]
LP_TFP <- exp(USdata$Y-cbind(USdata$K, USdata$L)%*%as.numeric(All_NAICS_LP))
QLP_TFP <- data.frame(cbind(USdata$id, USdata$year, apply(All_NAICS_QLP, 1, function(x) exp(USdata$Y-cbind(USdata$K, USdata$L)%*%as.numeric(x))), LP_TFP))
colnames(QLP_TFP) <- c("id", "Year", paste("Q", tau_t, sep=""), "LP")
TFP_Data <- group_by(QLP_TFP, Year) %>% summarise_at(c(paste("Q", tau_t, sep=""), "LP"), mean, na.rm=TRUE) %>% mutate_at(vars(-Year), function(x) x/x[1L]*100)

TFP_Plot_Title <- ggdraw() + draw_label("Productivity Over Time", fontface="plain", size=22)
TFP <- melt(TFP_Data, "Year")
TFP_Plot <- ggplot(TFP, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("") + ggtitle("Productivity Over Time") + theme(plot.title=element_text(size=22, face="plain"))+ scale_colour_manual(name="", labels=c(tau_t, "LP"), values=c(pcolour, "red"))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/TFP_Plot.png", TFP_Plot, base_height=8, base_width=10)








