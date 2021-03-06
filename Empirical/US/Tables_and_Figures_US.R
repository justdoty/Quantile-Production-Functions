library(stringr)
library(dplyr)
library(xtable)
#Load US dataset
USdata <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/USdata.csv') %>% 
  select(id, year, Y, VA, K, L, M, naics3) %>% transmute(id=id, year=year, Y=log(Y), VA=log(VA), K=log(K), L=log(L), M=log(M), naics3=as.character(naics3), naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
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
NAICS_labels <- array(NA, 4*length(NAICS)); NAICS_labels[seq(1, 4*length(NAICS), by=4)] <- paste(NAICS, paste("(Total=", size$Total, ")", sep=""))
NAICS_labels[is.na(NAICS_labels)] <- ""
summary_table <- cbind(NAICS_labels, rep(c("Output", "Capital", "Labor", "Materials"), 4), sumstat)
colnames(summary_table) <- c("Industry (NAICS code)", " ", "1st Qu.", 'Median', "3rd Qu.", 'Mean', "sd")

summary_table <- xtable(summary_table, digits=c(2,2,0,4,4,2,2,2), type="latex")
align(summary_table) <- rep('c', 8)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline '
#For copy pasting into latex
print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#Saves to file
print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Estimates/US_Summary.tex")
############################################################################################################
#################################Load and prepare data frames for estimates#################################
############################################################################################################
alpha <- .1
#Vector of quantiles of firm-size
tauvec <- seq(5, 95, length.out=19)/100
#Vector of quantiles of TFP
tfptau <- seq(5, 95, length.out=19)/100
#Number of parameters
dZ <- 2
###############################################################################
#Store QLP Results
##############################################################################
#Store QLP Estimates and Standard Deviations
QLP_betahat <- array(0, c(length(tauvec), dZ, length(NAICS)))
QLP_betaSE <- array(0, c(length(tauvec), dZ, length(NAICS)))
#Store QR Estimates and Standard Deviations
QR_betahat <- array(0, c(length(tauvec), dZ, length(NAICS)))
QR_betaSE <- array(0, c(length(tauvec), dZ, length(NAICS)))
#Store QDIF Estimates and Standard Deviations
QDIF_hat <- array(0, c(length(tauvec), dZ, length(NAICS)))
QDIF_SE <- array(0, c(length(tauvec), dZ, length(NAICS)))
#Store QTFP Estimates and Standard Deviations
QTFP_hat <- array(0, c(length(tauvec), length(tfptau), length(NAICS)))
QTFP_SE <- array(0, c(length(tauvec), length(tfptau), length(NAICS)))
#Store QLP RTS Estimates and Standard Deviations
QLP_RTS <- array(0, c(length(tauvec), length(NAICS)))
QLP_RTS_SE <- array(0, c(length(tauvec), length(NAICS)))
#Store QLP Capital Intensity Estimates and Standard Deviations
QLP_IN <- array(0, c(length(tauvec), length(NAICS)))
QLP_IN_SE <- array(0, c(length(tauvec), length(NAICS)))
#############################################################################
#Store LP Results
#############################################################################
#Store LP Estimates and Standard Devitaions
LP_betahat <- array(0, c(length(NAICS), dZ))
LP_betaSE <- array(0, c(length(NAICS), dZ))
#Store LP RTS Standard Deviations
LP_RTS_SE <- array(0, c(length(NAICS), 1))
#Store LP Capital Intensity Standard Deviations
LP_IN_SE <- array(0, c(length(NAICS), 1))
##############################################################################
#Load LP and QLP Results
#############################################################################@
for (i in 1:length(NAICS)){
  for (j in 1:length(tauvec)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Environments/QLP_Boot_US_Q%s.RData", j))
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
    #Load the LP estimates from a single quantile environment: they should be the same across quantiles
    if (tauvec[j]==0.1){
      #LP Estimates and Standard Deviations
      LP_betahat[i,] <- LPhat[,i]
      LP_betaSE[i,] <- apply(LPboot[,,i], 2, sd)
      #Store LP RTS Standard Deviations
      LP_RTS_SE[i,] <- sd(apply(LPboot[,,i], 1, sum))
      #Store LP Capital Intensity Standard Deviations
      LP_IN_SE[i,] <- sd(apply(LPboot[,,i], 1, function(x) x[1]/x[2]))
    }
  }
}
#LP RTS Estimates
LP_RTS <- as.matrix(apply(LP_betahat, 1, sum))
#LP Capital Intensity Estimates
LP_IN <- as.matrix(apply(LP_betahat, 1, function(x) x[1]/x[2]))
#Make an estimates table for QLP Beta Estimates and Standard Deviations
QLP_betatable <- data.frame(cbind(rep(tauvec, length(NAICS)), cbind(do.call(rbind, lapply(seq(dim(QLP_betahat)[3]), function(x) QLP_betahat[ , , x])), do.call(rbind, lapply(seq(dim(QLP_betaSE)[3]), function(x) QLP_betaSE[ , , x])))[,c(rbind(c(1:dZ), dZ+(1:dZ)))]))
QLPestimates <- cbind(QLP_betatable, c(QLP_RTS), c(QLP_RTS_SE), c(QLP_IN), c(QLP_IN_SE))
colnames(QLPestimates) <- c('Tau','K',"se_K", 'L', "se_L", 'RTS', 'RTS_SE', 'In', 'In_SE')
#Make an estimates table for LP Beta Estimates and Standard Deviations
LP_betatable <- data.frame(cbind(LP_betahat, LP_betaSE)[,c(rbind(c(1:dZ), dZ+(1:dZ)))])
LPestimates <- cbind(LP_betatable, LP_RTS, LP_RTS_SE, LP_IN, LP_IN_SE)
colnames(LPestimates) <- c('K',"se_K", 'L', "se_L", 'RTS', 'RTS_SE', 'In', 'In_SE')
#Prepare estimates for table in paper/presentation
tau_table <- c(0.1, 0.25, 0.5, 0.9)

#Table Labels
NAICS_labels <- array(NA, length(tau_table)*length(NAICS)); NAICS_labels[seq(1, length(tau_table)*length(NAICS), by=length(tau_table))] <- NAICS
NAICS_labels[is.na(NAICS_labels)] <- ""

QLP_Table <- cbind(NAICS_labels, QLPestimates[rep(tauvec, length(NAICS))%in%tau_table, ])
colnames(QLP_Table) <- c("NAICS", "$\\tau$", rep(c("Coef.", "s.e"), dZ+2))
QLP_Table_X <- xtable(QLP_Table, digits=c(0,0,2,rep(c(3,4), dZ+2)), type="latex")
align(QLP_Table_X) <- rep('c', 7+dZ*2)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Returns to Scale} & \\multicolumn{2}{c}{Capital Intensity}\\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}'
#For copy pasting to latex
print(QLP_Table_X, hline.after=c(0,nrow(QLP_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#For saving to file
print(QLP_Table_X, hline.after=c(0,nrow(QLP_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Estimates/US_Beta_Estimates.tex")
##############################################################################################
############################Coefficicent Plots######################################
################################################################################################
require(ggplot2)
require(cowplot)
require(reshape2)
#Industry NAICS Code Plot Labels
QLP_Kplot <- list(); QLP_Lplot <- list(); QLP_RTSplot <- list()
QDIF_Kplot <- list(); QDIF_Lplot <- list(); QTFP_plot <- list()
pcolour <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
tau_t <- c(0.1, 0.25, 0.5, 0.9)
for (p in 1:length(NAICS)){
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
  QLP_Kplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=LB.K, ymax=UB.K), fill="grey70") + geom_line(aes(y=K)) + geom_hline(yintercept=LPplotdat$K, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.K, LPplotdat$UB.K), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(min(QLPplotdat$LB.K), LPplotdat$LB.K, LPplotdat$LB.L, min(QLPplotdat$LB.L)), max(max(QLPplotdat$UB.K), LPplotdat$UB.K, LPplotdat$UB.L, max(QLPplotdat$UB.L))))
  QLP_Lplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=LB.L, ymax=UB.L), fill="grey70") + geom_line(aes(y=L)) + geom_hline(yintercept=LPplotdat$L, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.L, LPplotdat$UB.L), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(min(QLPplotdat$LB.K), LPplotdat$LB.K, LPplotdat$LB.L, min(QLPplotdat$LB.L)), max(max(QLPplotdat$UB.K), LPplotdat$UB.K, LPplotdat$UB.L, max(QLPplotdat$UB.L))))
  QLP_RTSplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Returns to Scale") + geom_ribbon(aes(ymin=LB.RTS, ymax=UB.RTS), fill="grey70") + geom_line(aes(y=RTS)) + geom_hline(yintercept=LPplotdat$RTS, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.RTS, LPplotdat$UB.RTS), linetype='dashed', color='red')
  #QTFP Plots
  QTFP <- QTFP_hat[,,p][,tfptau%in%tau_t]
  taufac <- as.factor(rep(tau_t, each=length(tauvec)))
  QLP_QTFPdat <- data.frame(tau=rep(tauvec, length(tau_t)), tfptau=taufac, qtfp=c(QTFP))
  QTFP_plot[[p]] <- ggplot(QLP_QTFPdat, aes(x=tau, y=qtfp, group=tfptau)) + xlab(expression('percentile-'*tau)) + ylab("")+ ggtitle(paste("NAICS", NAICS[p], sep=" "))+ geom_line(aes(colour=tfptau)) + scale_colour_manual(name=expression(tau), labels=tau_t, values = pcolour)
  #Plotting data for QDIF Plots
  QDIF_dat <- data.frame(tau=tauvec, coef=QDIF_hat[,,p], LB=QDIF_hat[,,p]-QDIF_SE[,,p]*qnorm(1-.05/2), UB=QDIF_hat[,,p]+QDIF_SE[,,p]*qnorm(1-.05/2))
  #QDIF Plots
  QDIF_Kplot[[p]] <- ggplot(QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_point(aes(y=coef.1)) + geom_errorbar(aes(ymin=LB.1, ymax=UB.1)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(min(QDIF_dat$LB.1), min(QDIF_dat$LB.2), 0), max(max(QDIF_dat$UB.1), max(QDIF_dat$UB.2), 0)))
  QDIF_Lplot[[p]] <- ggplot(QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_point(aes(y=coef.2)) + geom_errorbar(aes(ymin=LB.2, ymax=UB.2)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(min(QDIF_dat$LB.1), min(QDIF_dat$LB.2), 0), max(max(QDIF_dat$UB.1), max(QDIF_dat$UB.2), 0)))
  #Combine Plots #######################################################
  #QLP Coefficient Plots
  NAICS_plots <- ggdraw() + draw_label(paste("NAICS", NAICS[p], sep=" "), fontface="plain", size=22) + theme(plot.title = element_text(hjust = 0.5))
  coef_row1 <- plot_grid(QLP_Kplot[[p]], QLP_Lplot[[p]])
  coef_row2 <- plot_grid(QDIF_Kplot[[p]], QDIF_Lplot[[p]])
  Coef_Plot <- plot_grid(NAICS_plots, coef_row1, coef_row2, ncol=1, align="h", rel_heights = c(0.3, 1, 1))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/Coef_Plot_NAICS_", NAICS[p], ".png", sep=""), Coef_Plot, base_height=8, base_width=7)
}
#Combine the QTFP plots over industries
QTFP_plot <- plot_grid(QTFP_plot[[1]], QTFP_plot[[2]], QTFP_plot[[3]], QTFP_plot[[4]])
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/QTFP_plot.png", QTFP_plot, base_height=10, base_width=10)
###################################################################################################
###############TFP Over Time#######################################################
######################################################################################################
QLPestimatesTFP <- QLPestimates[tauvec%in%tau_t,]
All_NAICS_QLP <- data.frame(QLPestimatesTFP[(nrow(QLPestimatesTFP)-(length(tau_t)-1)):nrow(QLPestimatesTFP), c(2,4)])
All_NAICS_LP <- LPestimates[nrow(LPestimates), c(1,3)]
LP_TFP <- exp(USdata$VA-cbind(USdata$K, USdata$L)%*%as.numeric(All_NAICS_LP))
QLP_TFP <- data.frame(cbind(USdata$id, USdata$year, apply(All_NAICS_QLP, 1, function(x) exp(USdata$VA-cbind(USdata$K, USdata$L)%*%as.numeric(x))), LP_TFP))
colnames(QLP_TFP) <- c("id", "Year", paste("Q", tau_t, sep=""), "LP")
TFP_Data <- group_by(QLP_TFP, Year) %>% summarise_at(c(paste("Q", tau_t, sep=""), "LP"), mean, na.rm=TRUE) %>% mutate_at(vars(-Year), function(x) x/x[1L]*100)

TFP_Plot_Title <- ggdraw() + draw_label("Productivity Over Time", fontface="plain", size=22)
TFP <- melt(TFP_Data, "Year")
TFP_Plot <- ggplot(TFP, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("") + ggtitle("Productivity Over Time") + theme(plot.title=element_text(size=22, face="plain"))+ scale_colour_manual(name="", labels=c(tau_t, "LP"), values=c(pcolour, "red"))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/TFP_Plot.png", TFP_Plot, base_height=8, base_width=10)
################################################################################################
###################Coefficients over Time #######################################
################################################################################
T <- 5
dZ <- 2
time <- seq(min(USdata$year), max(USdata$year), by=T)
QLPT_Coef <- array(0, dim=c(length(tau_t), dZ, length(time)))
LPT_Coef <- array(0, dim=c(length(time), dZ))
for (t in 1:length(time)){
  for (q in 1:length(tau_t)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Environments/QLPT_US_Q%s.RData", q))
    QLPT_Coef[,,t][q,] <- betahat[,t]
    if (tau_t[q]==0.5){
      LPT_Coef[t,] <- LPhat[,t]
    }
  }
}
KT <- data.frame(cbind(time, t(QLPT_Coef[,1,][,]), LPT_Coef[,1]))
colnames(KT) <- c("Year", paste("Q", tau_t, sep=" "), "LP")
KT <- melt(KT, "Year")
KTplot <- ggplot(KT, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Capital") + scale_colour_manual(name="", labels=c(tau_t, "LP"), values=c(pcolour, "red"))
LT <- data.frame(cbind(time, t(QLPT_Coef[,2,][,]), LPT_Coef[,2]))
colnames(LT) <- c("Year", paste("Q", tau_t, sep=""), "LP")
LT <- melt(LT, "Year")
LTplot <- ggplot(LT, aes(x=Year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Labor") + scale_colour_manual(name="", labels=c(tau_t, "LP"), values=c(pcolour, "red"))
Plot_Title <- ggdraw() + draw_label("Output Elasticities Over Time", fontface="plain", size=22) 
Time_Plot <- plot_grid(Plot_Title, plot_grid(LTplot, KTplot), ncol=1, rel_heights = c(0.3, 1))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/Time_Plot.png", Time_Plot, base_height=6, base_width=10)


###############################################################################
###############################################################################
###############################################################################
#APPENDIX PLOTS AND TABLES (OPTIONAL)#########################################
###############################################################################
###############################################################################
###############################################################################

###############################################################################
#Store QACF Results
##############################################################################
#Store QACF Estimates and Standard Deviations
QACF_betahat <- array(0, c(length(tauvec), dZ, length(NAICS)))
QACF_betaSE <- array(0, c(length(tauvec), dZ, length(NAICS)))
#Store QR Estimates and Standard Deviations
QR_betahat <- array(0, c(length(tauvec), dZ, length(NAICS)))
QR_betaSE <- array(0, c(length(tauvec), dZ, length(NAICS)))
#Store QDIF Estimates and Standard Deviations
QDIF_hat <- array(0, c(length(tauvec), dZ, length(NAICS)))
QDIF_SE <- array(0, c(length(tauvec), dZ, length(NAICS)))
#Store QTFP Estimates and Standard Deviations
QTFP_hat <- array(0, c(length(tauvec), length(tfptau), length(NAICS)))
QTFP_SE <- array(0, c(length(tauvec), length(tfptau), length(NAICS)))
#Store QACF RTS Estimates and Standard Deviations
QACF_RTS <- array(0, c(length(tauvec), length(NAICS)))
QACF_RTS_SE <- array(0, c(length(tauvec), length(NAICS)))
#Store QACF Capital Intensity Estimates and Standard Deviations
QACF_IN <- array(0, c(length(tauvec), length(NAICS)))
QACF_IN_SE <- array(0, c(length(tauvec), length(NAICS)))
#############################################################################
#Store ACF Results
#############################################################################
#Store ACF Estimates and Standard Devitaions
ACF_betahat <- array(0, c(length(NAICS), dZ))
ACF_betaSE <- array(0, c(length(NAICS), dZ))
#Store ACF RTS Standard Deviations
ACF_RTS_SE <- array(0, c(length(NAICS), 1))
#Store ACF Capital Intensity Standard Deviations
ACF_IN_SE <- array(0, c(length(NAICS), 1))
##############################################################################
#Load ACF and QACF Results
#############################################################################@
for (i in 1:length(NAICS)){
  for (j in 1:length(tauvec)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Environments/QACF/QACF_Boot_US_Q%s.RData", j))
    #QACF Estimates and Standard Deviations
    QACF_betahat[,,i][j,] <- betahat[,i]
    QACF_betaSE[,,i][j,] <- apply(betaboot[,,i], 2, sd)
    #QR Estimates and Standard Deviations
    QR_betahat[,,i][j,] <- qrhat[,i]
    QR_betaSE[,,i][j,] <- apply(qrboot[,,i], 2, sd)
    #QACF-QR Estimates and Standard Deviations
    QDIF_hat[,,i][j,] <- qdifhat[,i]
    QDIF_SE[,,i][j,] <- apply(qdifboot[,,i], 2, sd)
    #QACF QTFP Estimates and Standard Deviations
    QTFP_hat[,,i][j,] <- QTFPhat[,i]
    QTFP_SE[,,i][j,] <- apply(QTFPboot[,,i], 2, sd)
    #QACF RTS Estimates and Standard Deviations
    QACF_RTS[j,i] <- sum(betahat[,i])
    QACF_RTS_SE[j,i] <- sd(apply(betaboot[,,i], 1, sum))
    #QACF Capital Intensity Estimates and Standard Deviations
    QACF_IN[j,i] <- betahat[,i][1]/betahat[,i][2]
    QACF_IN_SE[j,i] <- sd(apply(betaboot[,,i], 1, function(x) x[1]/x[2]))
    #Load the ACF estimates from a single quantile environment: they should be the same across quantiles
    if (tauvec[j]==0.1){
      #ACF Estimates and Standard Deviations
      ACF_betahat[i,] <- ACFhat[,i]
      ACF_betaSE[i,] <- apply(ACFboot[,,i], 2, sd)
      #Store ACF RTS Standard Deviations
      ACF_RTS_SE[i,] <- sd(apply(ACFboot[,,i], 1, sum))
      #Store ACF Capital Intensity Standard Deviations
      ACF_IN_SE[i,] <- sd(apply(ACFboot[,,i], 1, function(x) x[1]/x[2]))
    }
  }
}
#ACF RTS Estimates
ACF_RTS <- as.matrix(apply(ACF_betahat, 1, sum))
#ACF Capital Intensity Estimates
ACF_IN <- as.matrix(apply(ACF_betahat, 1, function(x) x[1]/x[2]))
#Make an estimates table for QACF Beta Estimates and Standard Deviations
QACF_betatable <- data.frame(cbind(rep(tauvec, length(NAICS)), cbind(do.call(rbind, lapply(seq(dim(QACF_betahat)[3]), function(x) QACF_betahat[ , , x])), do.call(rbind, lapply(seq(dim(QACF_betaSE)[3]), function(x) QACF_betaSE[ , , x])))[,c(rbind(c(1:dZ), dZ+(1:dZ)))]))
QACFestimates <- cbind(QACF_betatable, c(QACF_RTS), c(QACF_RTS_SE), c(QACF_IN), c(QACF_IN_SE))
colnames(QACFestimates) <- c('Tau','K',"se_K", 'L', "se_L", 'RTS', 'RTS_SE', 'In', 'In_SE')
#Make an estimates table for ACF Beta Estimates and Standard Deviations
ACF_betatable <- data.frame(cbind(ACF_betahat, ACF_betaSE)[,c(rbind(c(1:dZ), dZ+(1:dZ)))])
ACFestimates <- cbind(ACF_betatable, ACF_RTS, ACF_RTS_SE, ACF_IN, ACF_IN_SE)
colnames(ACFestimates) <- c('K',"se_K", 'L', "se_L", 'RTS', 'RTS_SE', 'In', 'In_SE')
#Prepare estimates for table in paper/presentation
tau_table <- c(0.1, 0.25, 0.5, 0.9)

#Table Labels
NAICS_labels <- array(NA, length(tau_table)*length(NAICS)); NAICS_labels[seq(1, length(tau_table)*length(NAICS), by=length(tau_table))] <- NAICS
NAICS_labels[is.na(NAICS_labels)] <- ""

QACF_Table <- cbind(NAICS_labels, QACFestimates[rep(tauvec, length(NAICS))%in%tau_table, ])
colnames(QACF_Table) <- c("NAICS", "$\\tau$", rep(c("Coef.", "s.e"), dZ+2))
QACF_Table_X <- xtable(QACF_Table, digits=c(0,0,2,rep(c(3,4), dZ+2)), type="latex")
align(QACF_Table_X) <- rep('c', 7+dZ*2)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Returns to Scale} & \\multicolumn{2}{c}{Capital Intensity}\\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}'
#For copy pasting to latex
print(QACF_Table_X, hline.after=c(0,nrow(QACF_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#For saving to file
print(QACF_Table_X, hline.after=c(0,nrow(QACF_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Estimates/US_ACF_Estimates.tex")
