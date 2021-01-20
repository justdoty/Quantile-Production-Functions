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
alpha <- .1
#Vector of quantiles of firm-size
tauvec <- seq(5, 95, length.out=19)/100
#Vector of quantiles of TFP
tfptau <- c(0.1, 0.25, 0.5, 0.9)
#Number of parameters
dZ <- 2
###############################################################################
#Store QLP Results
##############################################################################
#Store QLP Estimates and Standard Deviations
QLP_betahat <- array(0, c(length(tauvec), dZ, length(ISIC)))
QLP_betaSE <- array(0, c(length(tauvec), dZ, length(ISIC)))
#Store QR Estimates and Standard Deviations
QR_betahat <- array(0, c(length(tauvec), dZ, length(ISIC)))
QR_betaSE <- array(0, c(length(tauvec), dZ, length(ISIC)))
#Store QDIF Estimates and Standard Deviations
QDIF_hat <- array(0, c(length(tauvec), dZ, length(ISIC)))
QDIF_SE <- array(0, c(length(tauvec), dZ, length(ISIC)))
#Store QTFP Estimates and Standard Deviations
QTFP_hat <- array(0, c(length(tauvec), length(tfptau), length(ISIC)))
QTFP_SE <- array(0, c(length(tauvec), length(tfptau), length(ISIC)))
#Store QLP RTS Estimates and Standard Deviations
QLP_RTS <- array(0, c(length(tauvec), length(ISIC)))
QLP_RTS_SE <- array(0, c(length(tauvec), length(ISIC)))
#Store QLP Capital Intensity Estimates and Standard Deviations
QLP_IN <- array(0, c(length(tauvec), length(ISIC)))
QLP_IN_SE <- array(0, c(length(tauvec), length(ISIC)))
#############################################################################
#Store LP Results
#############################################################################
#Store LP Estimates and Standard Devitaions
LP_betahat <- array(0, c(length(ISIC), dZ))
LP_betaSE <- array(0, c(length(ISIC), dZ))
#Store LP QTFP Estimates and Standard Deviations
LP_QTFP_hat <- array(0, c(length(ISIC), length(tfptau)))
LP_QTFP_SE <- array(0, c(length(ISIC), length(tfptau)))
#Store LP RTS Standard Deviations
LP_RTS_SE <- array(0, c(length(ISIC), 1))
#Store LP Capital Intensity Standard Deviations
LP_IN_SE <- array(0, c(length(ISIC), 1))
##############################################################################
#Load LP and QLP Results
#############################################################################@
for (i in 1:length(ISIC)){
  load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/LP_COL_ISIC_%s.RData", i))
  #LP Estimates and Standard Deviations
  LP_betahat[i,] <- betahat
  LP_betaSE[i,] <- apply(betaboot, 2, sd)
  #LP QTFP Estimates and Standard Deviations
  LP_QTFP_hat[i,] <- QTFPhat
  LP_QTFP_SE[i,] <- apply(QTFPboot, 2, sd)
  #LP RTS Standard Deviations
  LP_RTS_SE[i,] <- sd(apply(betaboot, 1, sum))
  #LP Capital Intensity Standard Deviations
  LP_IN_SE[i,] <- sd(apply(betaboot, 1, function(x) x[1]/x[2]))
  for (j in 1:length(tauvec)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/QLP_Boot_COL_Q%s.RData", j))
    #QLP Estimates and Standard Deviations
    QLP_betahat[,,i][j,] <- betahat[,i]
    QLP_betaSE[,,i][j,] <- apply(betaboot[,,i], 2, sd)
    #QR Estimates and Standard Deviations
    QR_betahat[,,i][j,] <- qrhat[,i]
    QR_betaSE[,,i][j,] <- apply(qrboot[,,i], 2, sd)
    #QLP-QR Estimates and Standard Deviations
    QDIF_hat[,,i][j,] <- qdifhat[,i]
    QDIF_SE[,,i][j,] <- apply(qdifboot[,,i], 2, sd)
    #QLP QTFP Estimates and Standard Deviations
    QTFP_hat[,,i][j,] <- QTFPhat[,i]
    QTFP_SE[,,i][j,] <- apply(QTFPboot[,,i], 2, sd)
    #QLP RTS Estimates and Standard Deviations
    QLP_RTS[j,i] <- sum(betahat[,i])
    QLP_RTS_SE[j,i] <- sd(apply(betaboot[,,i], 1, sum))
    #QLP Capital Intensity Estimates and Standard Deviations
    QLP_IN[j,i] <- betahat[,i][1]/betahat[,i][2]
    QLP_IN_SE[j,i] <- sd(apply(betaboot[,,i], 1, function(x) x[1]/x[2]))
  }
}
#LP RTS Estimates
LP_RTS <- as.matrix(apply(LP_betahat, 1, sum))
#LP Capital Intensity Estimates
LP_IN <- as.matrix(apply(LP_betahat, 1, function(x) x[1]/x[2]))
#Make an estimates table for QLP Beta Estimates and Standard Deviations
QLP_betatable <- data.frame(cbind(rep(tauvec, length(ISIC)), cbind(do.call(rbind, lapply(seq(dim(QLP_betahat)[3]), function(x) QLP_betahat[ , , x])), do.call(rbind, lapply(seq(dim(QLP_betaSE)[3]), function(x) QLP_betaSE[ , , x])))[,c(rbind(c(1:dZ), dZ+(1:dZ)))]))
QLPestimates <- cbind(QLP_betatable, c(QLP_RTS), c(QLP_RTS_SE), c(QLP_IN), c(QLP_IN_SE))
colnames(QLPestimates) <- c('Tau','K',"se_K", 'L', "se_L", 'RTS', 'RTS_SE', 'In', 'In_SE')
#Make an estimates table for LP Beta Estimates and Standard Deviations
LP_betatable <- data.frame(cbind(LP_betahat, LP_betaSE)[,c(rbind(c(1:dZ), dZ+(1:dZ)))])
LPestimates <- cbind(LP_betatable, LP_RTS, LP_RTS_SE, LP_IN, LP_IN_SE)
colnames(LPestimates) <- c('K',"se_K", 'L', "se_L", 'RTS', 'RTS_SE', 'In', 'In_SE')
#Prepare estimates for table in paper/presentation
tau_table <- c(0.1, 0.25, 0.5, 0.9)

#Table Labels
ISIC_labels <- array(NA, length(tau_table)*length(ISIC)); ISIC_labels[seq(1, length(tau_table)*length(ISIC), by=length(tau_table))] <- ISIC
ISIC_labels[is.na(ISIC_labels)] <- ""

QLP_Table <- cbind(ISIC_labels, QLPestimates[rep(tauvec, length(ISIC))%in%tau_table, 1:(dZ*2+3)])
colnames(QLP_Table) <- c("Industry (ISIC code)", "$\\tau$", rep(c("Coef.", "s.e"), dZ+1))
QLP_Table_X <- xtable(QLP_Table, digits=c(0,0,2,rep(c(3,4), dZ+1)), type="latex", caption="Coefficient Estimates and Standard Errors for Colombian Manufacturing Firms")
align(QLP_Table_X) <- rep('c', 5+dZ*2)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Returns to Scale}\\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}'
print(QLP_Table, hline.after=c(0,nrow(QLP_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, caption.placement="top", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/COL/Estimates/COL_Beta_Estimates.tex")
##############################################################################################
############################Coefficicent Plots######################################
################################################################################################
require(ggplot2)
require(cowplot)
require(reshape2)
#Industry ISIC Code Plot Labels
QLP_Kplot <- list(); QLP_Lplot <- list(); QLP_RTSplot <- list()
QDIF_Kplot <- list(); QDIF_Lplot <- list(); QTFP_plot <- list()
pcolour <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
for (p in 1:length(ISIC)){
  #Plotting data for QLP
  QLPplotcoef <- apply(QLPestimates[c("K", "L", "RTS")], 2, function(x) split(x, ceiling(seq_along(x)/length(tauvec)))[[p]])
  QLPplotsd <- apply(QLPestimates[c("se_K", "se_L", "RTS_SE")], 2, function(x) split(x, ceiling(seq_along(x)/length(tauvec)))[[p]])
  QLPplotCI <- data.frame(LB=QLPplotcoef-QLPplotsd*qnorm(1-alpha/2), UB=QLPplotcoef+QLPplotsd*qnorm(1-alpha/2))
  QLPplotdat <- data.frame(tau=tauvec, QLPplotcoef, QLPplotsd, QLPplotCI)
  #Plotting data for LP
  LPplotcoef <- LPestimates[c("K", "L", "RTS")][p,]
  LPplotsd <- LPestimates[c("se_K", "se_L", "RTS_SE")][p,]
  LPplotCI <- data.frame(LB=LPplotcoef-LPplotsd*qnorm(1-alpha/2), UB=LPplotcoef+LPplotsd*qnorm(1-alpha/2))
  LPplotdat <- data.frame(LPplotcoef, LPplotsd, LPplotCI)
  #Plots
  QLP_Kplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=LB.K, ymax=UB.K), fill="grey70") + geom_line(aes(y=K)) + geom_hline(yintercept=LPplotdat$K, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.K, LPplotdat$UB.K), linetype='dashed', color='red')
  QLP_Lplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=LB.L, ymax=UB.L), fill="grey70") + geom_line(aes(y=L)) + geom_hline(yintercept=LPplotdat$L, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.L, LPplotdat$UB.L), linetype='dashed', color='red')
  QLP_RTSplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Returns to Scale") + geom_ribbon(aes(ymin=LB.RTS, ymax=UB.RTS), fill="grey70") + geom_line(aes(y=RTS)) + geom_hline(yintercept=LPplotdat$RTS, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.RTS, LPplotdat$UB.RTS), linetype='dashed', color='red')
  #QTFP Plots
  taufac <- as.factor(rep(tfptau, each=length(tauvec)))
  QLP_QTFPdat <- data.frame(tau=rep(tauvec, length(tfptau)), tfptau=taufac, qtfp=c(QTFP_hat[,,p]))
  QTFP_plot[[p]] <- ggplot(QLP_QTFPdat, aes(x=tau, y=qtfp, group=tfptau)) + xlab(expression('percentile-'*tau)) + ylab("")+ ggtitle(paste("ISIC", ISIC[p], sep=" "))+ geom_line(aes(colour=tfptau)) + scale_colour_manual(name=expression(tau), labels=tfptau, values = pcolour)
  #Plotting data for QDIF Plots
  QDIF_dat <- data.frame(tau=tauvec, coef=QDIF_hat[,,p], LB=QDIF_hat[,,p]-QDIF_SE[,,p]*qnorm(1-.05/2), UB=QDIF_hat[,,p]+QDIF_SE[,,p]*qnorm(1-.05/2))
  #QDIF Plots
  QDIF_Kplot[[p]] <- ggplot(QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_point(aes(y=coef.1)) + geom_errorbar(aes(ymin=LB.1, ymax=UB.1)) + geom_hline(yintercept=0, linetype='dashed', color='red')
  QDIF_Lplot[[p]] <- ggplot(QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_point(aes(y=coef.2)) + geom_errorbar(aes(ymin=LB.2, ymax=UB.2)) + geom_hline(yintercept=0, linetype='dashed', color='red')
  #Combine Plots #######################################################
  #QLP Coefficient Plots
  ISIC_plots <- ggdraw() + draw_label(paste("ISIC", ISIC[p], sep=" "), fontface="plain", size=22) + theme(plot.title = element_text(hjust = 0.5))
  coef_row1 <- plot_grid(QLP_Kplot[[p]], QLP_Lplot[[p]])
  coef_row2 <- plot_grid(QDIF_Kplot[[p]], QDIF_Lplot[[p]])
  Coef_Plot <- plot_grid(ISIC_plots, coef_row1, coef_row2, ncol=1, align="h", rel_heights = c(0.3, 1, 1))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/Coef_Plot_ISIC_", ISIC[p], ".png", sep=""), Coef_Plot, base_height=8, base_width=7)
}
#Combine the QTFP plots over industries
QTFP_plot <- plot_grid(QTFP_plot[[1]], QTFP_plot[[2]], QTFP_plot[[3]], QTFP_plot[[4]])
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/QTFP_plot.png", QTFP_plot, base_height=10, base_width=10)
###################################################################################################
###############TFP Over Time#######################################################
######################################################################################################
tau_t <- c(0.1, 0.25, 0.5, 0.9)
QLPestimatesTFP <- QLPestimates[tauvec%in%tau_t,]
All_ISIC_QLP <- data.frame(QLPestimatesTFP[(nrow(QLPestimatesTFP)-(length(tau_t)-1)):nrow(QLPestimatesTFP), c(2,4)])
All_ISIC_LP <- LPestimates[nrow(LPestimates), c(1,3)]
LP_TFP <- exp(COLdata$lnva-cbind(COLdata$lnk, COLdata$lnl)%*%as.numeric(All_ISIC_LP))
QLP_TFP <- data.frame(cbind(COLdata$id, COLdata$year, apply(All_ISIC_QLP, 1, function(x) exp(COLdata$lnva-cbind(COLdata$lnk, COLdata$lnl)%*%as.numeric(x))), LP_TFP))
colnames(QLP_TFP) <- c("id", "Year", paste("Q", tau_t, sep=""), "LP")
TFP_Data <- group_by(QLP_TFP, Year) %>% summarise_at(c(paste("Q", tau_t, sep=""), "LP"), mean, na.rm=TRUE) %>% mutate_at(vars(-Year), function(x) x/x[1L]*100)

TFP_Plot_Title <- ggdraw() + draw_label("Productivity Over Time", fontface="plain", size=22)
TFP <- melt(TFP_Data, "Year")
TFP_Plot <- ggplot(TFP, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("") + ggtitle("Productivity Over Time") + theme(plot.title=element_text(size=22, face="plain"))+ scale_colour_manual(name="", labels=c(tau_t, "LP"), values=c(pcolour, "red"))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/TFP_Plot.png", TFP_Plot, base_height=8, base_width=10)
################################################################################################
###################Coefficients over Time #######################################
################################################################################
T <- 5
dZ <- 2
time <- seq(min(COLdata$year), max(COLdata$year), by=T)
QLPT_Coef <- array(0, dim=c(length(tau_t), dZ, length(time)))
for (t in 1:length(time)){
  for (q in 1:length(tau_t)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/QLPT_COL_Q%s.RData", q))
    QLPT_Coef[,,t][q,] <- betahat[,t]
  }
}
KT <- data.frame(cbind(time, t(QLPT_Coef[,1,][,])))
colnames(KT) <- c("Year", paste("Q", tau_t, sep=" "))
KT <- melt(KT, "Year")
KTplot <- ggplot(KT, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Capital") + scale_colour_manual(name=expression(tau), labels=tau_t, values = pcolour)
LT <- data.frame(cbind(time, t(QLPT_Coef[,2,][,])))
colnames(LT) <- c("Year", paste("Q", tau_t, sep=""))
LT <- melt(LT, "Year")
LTplot <- ggplot(LT, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Labor") + scale_colour_manual(name=expression(tau), labels=tau_t, values = pcolour)
Plot_Title <- ggdraw() + draw_label("Output Elasticities Over Time", fontface="plain", size=22) 
Time_Plot <- plot_grid(Plot_Title, plot_grid(LTplot, KTplot), ncol=1, rel_heights = c(0.3, 1))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/Time_Plot.png", Time_Plot, base_height=6, base_width=10)
















