library(stringr)
library(dplyr)
library(xtable)
library(reshape2)
#Load COL dataset
COLdata <- read.csv("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Colombia/COLdata.csv")
#Industries as listed in QLP_COL.R file
ISIC <- c("311", "322", "381", "All")
ISIC_des <- c("Food Products", "Apparel", "Fabricated Metal Products", "All Manufacturing")
########################################################################################################
##########################################Summary Statistics############################################
########################################################################################################
#Create table for all relevant summary statistics
sumISIC <- filter(COLdata, isic3==311|isic3==322|isic3==381) %>% group_by(isic3) %>% summarise_at(c("lny", "lnk", "lnl", "lnm"), list(Q1=~quantile(., 0.25), med=median, Q3=~quantile(.,0.75), mean=mean, sd=sd), na.rm=TRUE) 
sumALL <- cbind("All", summarise_at(COLdata, c("lny", "lnk", "lnl", "lnm"), list(Q1=~quantile(., 0.25), med=median, Q3=~quantile(.,0.75), mean=mean, sd=sd), na.rm=TRUE))
colnames(sumALL)[1] <- "isic3"
sizeISIC <- filter(COLdata, isic3==311|isic3==322|isic3==381) %>% group_by(isic3) %>% summarise(Firms=length(unique(id)), Total=n())
sizeALL <- c("All", length(unique(COLdata$id)), nrow(COLdata))
size <- rbind(sizeISIC, sizeALL)
sumstat <- round(matrix(as.numeric(as.matrix(rbind(sumISIC, sumALL))[,-1]), nrow=16, ncol=5), 2)
#Some pretty formatting
ISIC_labels <- array(NA, 4*length(ISIC)); ISIC_labels[seq(1, 4*length(ISIC), by=4)] <- paste(ISIC, paste("(Total=", size$Total, ")", sep=""))
ISIC_labels[is.na(ISIC_labels)] <- ""
summary_table <- cbind(ISIC_labels, rep(c("Output", "Capital", "Labor", "Materials"), 4), sumstat)
colnames(summary_table) <- c("Industry (ISIC code)", " ", "1st Qu.", 'Median', "3rd Qu.", 'Mean', "sd")

summary_table <- xtable(summary_table, digits=c(2,2,0,4,4,2,2,2), type="latex", caption="Summary Statistics (in logs) for Colombia Manufacturing Data")
align(summary_table) <- rep('c', 8)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline '
print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H", caption.placement="top", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Estimates/COL_Summary.tex")
############################################################################################################
#################################Load and prepare data frames for estimates#################################
############################################################################################################
#Vector of quantiles
tau <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.75, 0.8, 0.9)
#Number of parameters
dZ <- 2
require(abind)
#Load the results from QLP, output is a list of estimates over quantiles for each industry
QLP_results <- list()
QLP_true <- list()
for (i in 1:length(tau)){
  load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/QLP_COL_Q%s.RData", i))
  QLP_results[[i]] <- results
  QLP_true[[i]] <- true.beta

}
#Load the results from LP, output is a list of estimates over industries
LP_results <- list()
LP_true <- list()
for (i in 1:length(ISIC)){
  load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/LP_COL_ISIC_%s.RData", i))
  LP_results[[i]] <- results
  LP_true[[i]] <- true.beta.LP
}
#Load the results from OLS and QR
QR_OLS_results <- load('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/COL_QR_OLS.Rdata')
#Prepare estimates for QLP: Outputs a list of quantile estimates over industry (columns)
#and estimates (rows)
#True QLP Coefficient Estimates
QLP_Coef <- array(0, dim=c(2, length(ISIC), length(tau)))
#Bootstrapped QLP Coefficient Estimates
QLP_Boot <- array(0, dim=c(2, length(ISIC), length(tau)))
#Bootstrapped Standard Errors
QLP_SE <- array(0, dim=c(2, length(ISIC), length(tau)))
#Estimates of Returns to Scale
QLP_RtS <- array(0, dim=c(length(tau), length(ISIC)))
#Estimates of Standard Error Returns to Scale
QLP_RtS_SE <- array(0, dim=c(length(tau), length(ISIC)))
#Lower and Upper bounds of CI
QLP_Lower <- array(0, dim=c(2, length(ISIC), length(tau)))
QLP_Upper <- array(0, dim=c(2, length(ISIC), length(tau)))
for (i in 1:length(tau)){
  for (j in 1:length(ISIC)){
    QLP_Coef[,,i][,j] <- QLP_true[[i]][,j]
    QLP_Boot[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, mean)
    QLP_SE[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, sd)
    #Percentile Method
    QLP_Lower[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, function(x) quantile(x, alpha/2))
    QLP_Upper[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, function(x) quantile(x, 1-alpha/2))
    QLP_RtS[i,j] <- sum(QLP_Coef[,,i][,j])
    QLP_RtS_SE[i,j] <- sd(apply(QLP_results[[i]][,,j], 1, sum))
  }
}
#A little bit of reformating to obtain estimates over industries instead of quantiles
#Bootstrapped QLP Coefficient Estimates
QLP_Coef <- aperm(QLP_Coef, c(3, 1, 2))
QLP_Boot <- aperm(QLP_Boot, c(3, 1, 2))
QLP_SE <- aperm(QLP_SE, c(3, 1, 2))
QLP_RtS <- c(QLP_RtS)
QLP_RtS_SE <- c(QLP_RtS_SE)
QLP_Lower <- aperm(QLP_Lower, c(3, 1, 2))
QLP_Upper <- aperm(QLP_Upper, c(3, 1, 2))
#Prepare estimates for LP
#Bootstrapped LP Coefficient Estimates
LP_Coef <- array(0, dim=c(2, length(ISIC)))
LP_Boot <- array(0, dim=c(2, length(ISIC)))
LP_SE <- array(0, dim=c(2, length(ISIC)))
LP_Lower <- array(0, dim=c(2, length(ISIC)))
LP_Upper <- array(0, dim=c(2, length(ISIC)))
for (i in 1:length(ISIC)){
  LP_Coef[,i] <- LP_true[[i]]
  LP_Boot[,i] <- colMeans(LP_results[[i]]) 
  LP_SE[,i] <- apply(LP_results[[i]], 2, sd) 
  LP_Lower[,i] <- apply(LP_results[[i]], 2, function(x) quantile(x, 0.05))
  LP_Upper[,i] <- apply(LP_results[[i]], 2, function(x) quantile(x, 0.95)) 
}
#Listed by Industry
#For QLP
QLP_Coef <- lapply(seq(dim(QLP_Coef)[3]), function(x) QLP_Coef[ , , x])
QLP_Boot <- lapply(seq(dim(QLP_Boot)[3]), function(x) QLP_Boot[ , , x])
QLP_BC <- lapply(seq(dim(QLP_BC)[3]), function(x) QLP_BC[ , , x])
QLP_SE <- lapply(seq(dim(QLP_SE)[3]), function(x) QLP_SE[ , , x])
QLP_Upper <- lapply(seq(dim(QLP_Upper)[3]), function(x) QLP_Upper[ , , x])
QLP_Lower <- lapply(seq(dim(QLP_Lower)[3]), function(x) QLP_Lower[ , , x])

#For LP
LP_Coef <- lapply(seq(dim(LP_Coef)[2]), function(x) LP_Coef[,x])
LP_Boot <- lapply(seq(dim(LP_Boot)[2]), function(x) LP_Boot[,x])
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
estimates <- data.frame(cbind(rep(tau, length(ISIC)), cbind(do.call(rbind, QLP_Coef), do.call(rbind, QLP_SE))[,c(rbind(c(1:2), 2+(1:2)))]))
estimates <- cbind(estimates, QLP_RtS, QLP_RtS_SE)
colnames(estimates) <- c('Tau','K',"se_K", 'L', "se_L", 'RtS', 'RtS_SE')
#Make a Confidence Interval Table for Quantile GMM
QLP_CI <- data.frame(cbind(rep(tau, length(ISIC)), cbind(do.call(rbind, QLP_Lower), do.call(rbind, QLP_Upper))[,c(rbind(c(1:2), 2+(1:2)))]))
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
estimates_QR <- data.frame(cbind(rep(tau, length(ISIC)), do.call(rbind, QR_Coef)))
colnames(estimates_QR) <- c("Tau", "K", "L")
#Make a Confidence Interval Table for Quantile Regression
QR_CI <- data.frame(cbind(rep(tau, length(ISIC)), cbind(do.call(rbind, QR_Lower), do.call(rbind, QR_Upper))[,c(rbind(c(1:2), 2+(1:2)))]))
colnames(QR_CI) <- c('Tau', 'Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
#Prepare estimates for table in paper/presentation
tau_table <- c(0.1, 0.25, 0.5, 0.75)

#Table Labels
ISIC[length(ISIC)] <- "All"
ISIC_relabel <- ISIC %>% str_replace_all("[|]", ",") %>% str_replace_all("[\\.^]", "")
ISIC_labels <- array(NA, length(tau_table)*length(ISIC)); ISIC_labels[seq(1, length(tau_table)*length(ISIC), by=length(tau_table))] <- ISIC_relabel
ISIC_labels[is.na(ISIC_labels)] <- ""

estimates_table <- cbind(ISIC_labels, estimates[rep(tau, length(ISIC))%in%tau_table, ])
colnames(estimates_table) <- c("Industry (ISIC code)", "$\\tau$", "Coef.", 's.e.', "Coef.", 's.e.', "Coef.", "s.e")

estimates_table <- xtable(estimates_table, digits=c(0,0,2,3,4,3,4,3,4), type="latex", "Coefficient Estimates and Standard Errors for Colombia Manufacturing Firms")
align(estimates_table) <- rep('c', 9)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Returns to Scale} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}'
print(estimates_table, hline.after=c(0,nrow(estimates_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, caption.placement="top", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Estimates/COL_Estimates.tex")

############################Coefficicent Plots######################################
require(ggplot2)
require(cowplot)
#Industry ISIC Code Plot Labels
QLP_K_plot <- list(); QLP_L_plot <- list()
QR_K_plot <- list(); QR_L_plot <- list()
for (p in 1:length(ISIC)){
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
  QLP_K_plot[[p]] <- ggplot(QLPK_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QLPK_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QLPK_data$lower_LP, QLPK_data$upper_LP), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(qlpklow, lpklow, qrklow, lmklow), max(qlpkup, lpkup, qrkup, lmkup)))
  QR_K_plot[[p]] <- ggplot(QRK_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QRK_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QRK_data$lower_OLS, QRK_data$upper_OLS), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(qlpklow, lpklow, qrklow, lmklow), max(qlpkup, lpkup, qrkup, lmkup)))
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
  QLP_L_plot[[p]] <- ggplot(QLPL_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QLPL_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QLPL_data$lower_LP, QLPL_data$upper_LP), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(qlpllow, lpllow, qrllow, lmllow), max(qlplup, lplup, qrlup, lmlup)))
  QR_L_plot[[p]] <- ggplot(QRL_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QRL_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QRL_data$lower_OLS, QRL_data$upper_OLS), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(qlpllow, lpllow, qrllow, lmllow), max(qlplup, lplup, qrlup, lmlup)))
  ###############################Combine Plots ##############################################
  ############################################################################################
  ISIC_plots <- ggdraw() + draw_label(ISIC_des[p], fontface="plain", size=22) + theme(plot.title = element_text(hjust = 0.5))
  Lrow <- plot_grid(QLP_L_plot[[p]], QR_L_plot[[p]])
  Krow <- plot_grid(QLP_K_plot[[p]], QR_K_plot[[p]])
  Coef_Plot <- plot_grid(ISIC_plots, Lrow, Krow, ncol=1, align="h", rel_heights = c(0.3, 1, 1))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/Coef_Plot_ISIC_", ISIC_relabel[p], ".png", sep=""), Coef_Plot, base_height=8, base_width=7)
}
##################################Prepare Plots over Time##############################
###################Coefficients over Time #######################################
tau_t <- c(0.1, 0.3, 0.5, 0.7, 0.9)
T <- 2
dZ <- 2
time <- seq(1978, 1991, 2)
QLPT_Coef <- array(0, dim=c(length(tau_t), dZ, length(time)))
QLPT_True <- array(0, dim=c(length(tau_t), dZ, length(time)))
for (t in 1:length(time)){
  for (q in 1:length(tau_t)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/QLPT_COL_Q%s.RData", q))
    QLPT_Coef[,,t][q,] <- colMeans(results_T[,,t])
    QLPT_True[,,t][q,] <- true.beta_T[,t]
  }
}
pcolour <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
KT <- data.frame(cbind(time, t(QLPT_True[,1,][,])))
colnames(KT) <- c("Year", paste("Q", tau_t, sep=" "))
KT <- melt(KT, "Year")
KTplot <- ggplot(KT, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Capital") + scale_colour_manual(name=expression(tau), labels=tau_t, values = pcolour) + scale_x_continuous(breaks = time)
LT <- data.frame(cbind(time, t(QLPT_True[,2,][,])))
colnames(LT) <- c("Year", paste("Q", tau_t, sep=""))
LT <- melt(LT, "Year")
LTplot <- ggplot(LT, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Labor") + scale_colour_manual(name=expression(tau), labels=tau_t, values = pcolour) + scale_x_continuous(breaks = time)
Plot_Title <- ggdraw() + draw_label("Output Elasticities Over Time", fontface="plain", size=22) 
Time_Plot <- plot_grid(Plot_Title, plot_grid(LTplot, KTplot), ncol=1, rel_heights = c(0.3, 1))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/Time_Plot.png", Time_Plot, base_height=6, base_width=10)
###############TFP Over Time#######################################################
estimates_TFP <- estimates[tau%in%tau_t,]
All_ISIC_QLP <- data.frame(estimates_TFP[(nrow(estimates_TFP)-(length(tau_t)-1)):nrow(estimates_TFP), c(2,4)])
All_ISIC_LP <- estimates_LP[nrow(estimates_LP), c(1,3)]
LP_TFP <- exp(COLdata$lnva-cbind(COLdata$lnk, COLdata$lnl)%*%as.numeric(All_ISIC_LP))
QLP_TFP <- data.frame(cbind(COLdata$id, COLdata$year, apply(All_ISIC_QLP, 1, function(x) exp(COLdata$lnva-cbind(COLdata$lnk, COLdata$lnl)%*%as.numeric(x))), LP_TFP))
colnames(QLP_TFP) <- c("id", "Year", paste("Q", tau_t, sep=""), "LP")
TFP_Data <- group_by(QLP_TFP, Year) %>% summarise_at(c(paste("Q", tau_t, sep=""), "LP"), mean, na.rm=TRUE) %>% mutate_at(vars(-Year), function(x) x/x[1L]*100)
TFP_Data$Year <- seq(1978,1991,1)
TFP_Plot_Title <- ggdraw() + draw_label("Productivity Over Time", fontface="plain", size=22)
TFP <- melt(TFP_Data, "Year")
TFP_Plot <- ggplot(TFP, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("") + ggtitle("Productivity Over Time") + theme(plot.title=element_text(size=22, face="plain"))+ scale_colour_manual(name="", labels=c(tau_t, "LP"), values=c(pcolour, "red"))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/TFP_Plot.png", TFP_Plot, base_height=8, base_width=10)

















