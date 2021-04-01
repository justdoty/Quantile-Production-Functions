library(stringr)
library(dplyr)
library(xtable)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
source('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Functions/QLP_aux.R')
#Load COL dataset
COLdata <- read.csv("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Colombia/COLdata.csv")
COLdata$year <- COLdata$year+1975
#Industries as listed in QLP_COL.R file
ISIC <- c("311", "322", "381", "All")
industries <- c("311", "322", "381", "^3")
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

summary_table <- xtable(summary_table, digits=c(2,2,0,4,4,2,2,2), type="latex")
align(summary_table) <- rep('c', 8)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline '
print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
print(summary_table, hline.after=c(0,nrow(summary_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Estimates/COL_Summary.tex")
############################################################################################################
#################################Load and prepare data frames for estimates#################################
############################################################################################################
alpha <- .1
#Vector of quantiles of firm-size
tauvec <- seq(5, 95, length.out=19)/100
#Prepare estimates for table in paper/presentation
tau_t <- c(0.1, 0.25, 0.5, 0.9)
#Number of parameters
dZ <- 2
###############################################################################
#Store QLP Results
##############################################################################
#Store QLP Estimates and Standard Deviations
QLP_betahat <- array(0, c(length(tauvec), dZ, length(ISIC)))
QLP_betaSE <- array(0, c(length(tauvec), dZ, length(ISIC)))
#Store QR Estimates and Standard Deviations
QLP_QR_betahat <- array(0, c(length(tauvec), dZ, length(ISIC)))
QLP_QR_betaSE <- array(0, c(length(tauvec), dZ, length(ISIC)))
#Store QDIF Estimates and Standard Deviations
QLP_QDIF_hat <- array(0, c(length(tauvec), dZ, length(ISIC)))
QLP_QDIF_SE <- array(0, c(length(tauvec), dZ, length(ISIC)))
#Store QTFP Estimates
QLP_QTFP_hat <- list()
#Store LP TFP Estimates
LP_TFP_hat <- list()
#Store Omega Estimates
LP_omega_hat <- list()
#Store Ex-Post Shock Estimates
LP_eta_hat <- list()
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
#Store LP RTS Standard Deviations
LP_RTS_SE <- array(0, c(length(ISIC), 1))
#Store LP Capital Intensity Standard Deviations
LP_IN_SE <- array(0, c(length(ISIC), 1))
##############################################################################
#Load LP and QLP Results
#############################################################################@
for (i in 1:length(ISIC)){
  QLP_QTFP_hat[[i]] <- list()
  for (j in 1:length(tauvec)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/QLP_Boot_COL_Q%s.RData", j))
    #QLP Estimates and Standard Deviations
    QLP_betahat[,,i][j,] <- QLPbetahat[,i]
    QLP_betaSE[,,i][j,] <- apply(QLPbetaboot[,,i], 2, sd)
    #QR Estimates and Standard Deviations
    QLP_QR_betahat[,,i][j,] <- QLPqrhat[,i]
    QLP_QR_betaSE[,,i][j,] <- apply(QLPqrboot[,,i], 2, sd)
    #QLP-QR Estimates and Standard Deviations
    QLP_QDIF_hat[,,i][j,] <- QLPqdifhat[,i]
    QLP_QDIF_SE[,,i][j,] <- apply(QLPqdifboot[,,i], 2, sd)
    #QLP RTS Estimates and Standard Deviations
    QLP_RTS[j,i] <- sum(QLPbetahat[,i])
    QLP_RTS_SE[j,i] <- sd(apply(QLPbetaboot[,,i], 1, sum))
    #QLP Capital Intensity Estimates and Standard Deviations
    QLP_IN[j,i] <- QLPbetahat[,i][1]/QLPbetahat[,i][2]
    QLP_IN_SE[j,i] <- sd(apply(QLPbetaboot[,,i], 1, function(x) x[1]/x[2]))
    #QTFP Estimates
    QLP_QTFP_hat[[i]][[j]] <- exp(QLPTFPhat[[i]])
    LP_eta_ecdf <- ecdf(LPexpost[[i]])
    LP_eta_hat[[i]] <- LP_eta_ecdf(LPexpost[[i]])
    #Load the LP estimates from a single quantile environment: they should be the same across quantiles
    if (tauvec[j]==0.1){
      #LP Estimates and Standard Deviations
      LP_betahat[i,] <- LPhat[,i]
      LP_betaSE[i,] <- apply(LPboot[,,i], 2, sd)
      LP_TFP_hat[[i]] <- exp(LPTFPhat[[i]])
      LP_omega_hat[[i]] <- exp(LPomegahat[[i]])
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
QLP_betatable <- data.frame(cbind(rep(tauvec, length(ISIC)), cbind(do.call(rbind, lapply(seq(dim(QLP_betahat)[3]), function(x) QLP_betahat[ , , x])), do.call(rbind, lapply(seq(dim(QLP_betaSE)[3]), function(x) QLP_betaSE[ , , x])))[,c(rbind(c(1:dZ), dZ+(1:dZ)))]))
QLPestimates <- cbind(QLP_betatable, c(QLP_RTS), c(QLP_RTS_SE), c(QLP_IN), c(QLP_IN_SE))
colnames(QLPestimates) <- c('Tau','K',"se_K", 'L', "se_L", 'RTS', 'RTS_SE', 'In', 'In_SE')
#Make an estimates table for LP Beta Estimates and Standard Deviations
LP_betatable <- data.frame(cbind(LP_betahat, LP_betaSE)[,c(rbind(c(1:dZ), dZ+(1:dZ)))])
LPestimates <- cbind(LP_betatable, LP_RTS, LP_RTS_SE, LP_IN, LP_IN_SE)
colnames(LPestimates) <- c('K',"se_K", 'L', "se_L", 'RTS', 'RTS_SE', 'In', 'In_SE')

#Table Labels
ISIC_labels <- array(NA, length(tau_t)*length(ISIC)); ISIC_labels[seq(1, length(tau_t)*length(ISIC), by=length(tau_t))] <- ISIC
ISIC_labels[is.na(ISIC_labels)] <- ""

QLP_Table <- cbind(ISIC_labels, QLPestimates[rep(tauvec, length(ISIC))%in%tau_t, ])
colnames(QLP_Table) <- c("ISIC", "$\\tau$", rep(c("Coef.", "s.e"), dZ+2))
QLP_Table_X <- xtable(QLP_Table, digits=c(0,0,2,rep(c(3,4), dZ+2)), type="latex")
align(QLP_Table_X) <- rep('c', 7+dZ*2)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Returns to Scale} & \\multicolumn{2}{c}{Capital Intensity}\\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}'
#For copy pasting to latex
print(QLP_Table_X, hline.after=c(0,nrow(QLP_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#For saving to file
print(QLP_Table_X, hline.after=c(0,nrow(QLP_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Estimates/COL_LP_Estimates.tex")
############Density Plots for TFP############################
LPTFPeta <- list()
for (p in 1:length(industries)){
  #Interpolation of TFP from RC Model
  LPparvec <- t(cbind(split(QLPestimates$K, ceiling(seq_along(QLPestimates$K)/length(tauvec)))[[p]], split(QLPestimates$L, ceiling(seq_along(QLPestimates$L)/length(tauvec)))[[p]]))
  COLeta <- subset(COLdata, str_detect(isic3, industries[p]))
  LPetafit <- rowSums(cbind(COLeta$lnk, COLeta$lnl)*lsplinemat(tau=tauvec, par=LPparvec, u=LP_eta_hat[[p]]))
  LPTFPeta[[p]] <- COLeta$lnva-LPetafit
  LPdensdata <- melt(data.frame(DS=LPTFPeta[[p]], LP=LPTFPhat[[p]]))
  LPdensplot <- ggplot(LPdensdata, aes(x=value, color=variable))+geom_density()+xlab("") + ylab("") + scale_colour_manual(name="", labels=c("DS", "LP"), values=c("black", "red"))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/TFP/LP/TFP_Dens_Plot_ISIC_", ISIC[p], ".png", sep=""), LPdensplot, base_height=8, base_width=7)
  #Plot TFP densities per quantile (0.25, 0.5, 0.75)
  LPTFPtau <- melt(data.frame(tau1=QLPTFPhat[[p]][[5]], tau2=QLPTFPhat[[p]][[10]], tau3=QLPTFPhat[[p]][[15]], LP=LPTFPhat[[p]]))
  LPdenstauplot <- ggplot(LPTFPtau, aes(x=value, color=variable))+geom_density()+xlab("") + ylab("")
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/TFP/LP/TFPtau_Dens_Plot_ISIC_", ISIC[p], ".png", sep=""), LPdenstauplot, base_height=8, base_width=7)
}
####################################################################################################
#TFP Data Formatting
####################################################################################################
#Combine quantiles of QTFP
QLP_QTFP <- lapply(QLP_QTFP_hat,  function(x) do.call(cbind, x))
#Combine with LP TFP and LP productivity
QLP_TFP_data <- lapply(1:length(ISIC), function(x) data.frame(cbind(subset(COLdata, str_detect(isic3, industries[x]))$id, subset(COLdata, str_detect(isic3, industries[x]))$year, QLP_QTFP[[x]], exp(LPTFPeta[[x]]), LP_TFP_hat[[x]], LP_omega_hat[[x]])))
QLP_TFP_data <- lapply(QLP_TFP_data, setNames, nm = c("id", "year", paste("Q", tauvec, sep=" "), "TFPeta", "TFP", "Omega"))
#Take un-weighted averages and set base year to 100
QLP_TFP_AVG <- lapply(QLP_TFP_data, function(x) group_by(x, year) %>% summarise_at(c(paste("Q", tauvec, sep=" "), "TFPeta", "TFP", "Omega"), mean, na.rm=TRUE) %>% mutate_at(vars(-year), function(z) z/z[1L]*100))
#Subset Quantiles of interest for plotting and group by quantile
QLP_TFP <- lapply(QLP_TFP_AVG, function(x) melt(x[,c("year", paste("Q", tau_t), "TFP")], "year"))
QLP_TFPeta <- lapply(QLP_TFP_AVG, function(x) melt(x[,c("year", "TFPeta", "TFP")], "year"))
#Plot entire industry sample
pcolour <- brewer.pal(n=length(tau_t), "Spectral")
LP_TFP_Plot <- ggplot(QLP_TFP[[length(ISIC)]], aes(x=year, y=value, group=variable, linetype=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("") + scale_colour_manual(name="", labels=c(paste("TFP" ,tau_t), "LP TFP"), values=c(pcolour, "black")) + theme(legend.text.align = 0) + scale_linetype_manual(name="", labels=c(paste("TFP" ,tau_t), "LP TFP"), values=c("TFP 0.1"="solid", "TFP 0.25"="solid", "TFP 0.5"="solid","TFP 0.9"="solid", "TFP"="longdash", "NA"))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/TFP/LP/LP_TFP_Plot.png", LP_TFP_Plot, base_height=8, base_width=10)
LP_TFPeta_Plot <- ggplot(QLP_TFPeta[[length(ISIC)]], aes(x=year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("") + scale_colour_manual(name="", labels=c("DS", "LP"), values=c("black", "red"))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/TFP/LP/LP_TFPeta_Plot.png", LP_TFPeta_Plot, base_height=8, base_width=10)
#Create table for productivity growth rates
TFP_growth <- lapply(QLP_TFP_AVG, function(x) cbind(x$year, apply(x[,-1], 2, function(y) ((lead(y)/y)-1)*100)))
#Create table for differences across columns
QLP_TFP_tau_t <- lapply(QLP_TFP_AVG, function(x) x[,c("year", paste("Q", rev(tau_t)), "TFP", "Omega")])
QLP_TFP_DIF <- lapply(QLP_TFP_tau_t, function(x) t(apply(x[,-1], 1, function(y) y-lead(y)))[,c(1:(length(tau_t)-1))])
#Coefficient Plots
#Industry ISIC Code Plot Labels
QLP_Kplot <- list(); QLP_Lplot <- list(); QLP_RTSplot <- list()
QLP_QDIF_Kplot <- list(); QLP_QDIF_Lplot <- list()
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
  QLP_Kplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=LB.K, ymax=UB.K), fill="grey70") + geom_line(aes(y=K)) + geom_hline(yintercept=LPplotdat$K, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.K, LPplotdat$UB.K), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(min(QLPplotdat$LB.K), LPplotdat$LB.K, LPplotdat$LB.L, min(QLPplotdat$LB.L)), max(max(QLPplotdat$UB.K), LPplotdat$UB.K, LPplotdat$UB.L, max(QLPplotdat$UB.L))))
  QLP_Lplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=LB.L, ymax=UB.L), fill="grey70") + geom_line(aes(y=L)) + geom_hline(yintercept=LPplotdat$L, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.L, LPplotdat$UB.L), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(min(QLPplotdat$LB.K), LPplotdat$LB.K, LPplotdat$LB.L, min(QLPplotdat$LB.L)), max(max(QLPplotdat$UB.K), LPplotdat$UB.K, LPplotdat$UB.L, max(QLPplotdat$UB.L))))
  QLP_RTSplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Returns to Scale") + geom_ribbon(aes(ymin=LB.RTS, ymax=UB.RTS), fill="grey70") + geom_line(aes(y=RTS)) + geom_hline(yintercept=LPplotdat$RTS, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.RTS, LPplotdat$UB.RTS), linetype='dashed', color='red')
  #Plotting data for QDIF Plots
  QLP_QDIF_dat <- data.frame(tau=tauvec, coef=QLP_QDIF_hat[,,p], LB=QLP_QDIF_hat[,,p]-QLP_QDIF_SE[,,p]*qnorm(1-alpha/2), UB=QLP_QDIF_hat[,,p]+QLP_QDIF_SE[,,p]*qnorm(1-alpha/2))
  #QDIF Plots
  QLP_QDIF_Kplot[[p]] <- ggplot(QLP_QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_point(aes(y=coef.1)) + geom_errorbar(aes(ymin=LB.1, ymax=UB.1)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(min(QLP_QDIF_dat$LB.1), min(QLP_QDIF_dat$LB.2), 0), max(max(QLP_QDIF_dat$UB.1), max(QLP_QDIF_dat$UB.2), 0)))
  QLP_QDIF_Lplot[[p]] <- ggplot(QLP_QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_point(aes(y=coef.2)) + geom_errorbar(aes(ymin=LB.2, ymax=UB.2)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(min(QLP_QDIF_dat$LB.1), min(QLP_QDIF_dat$LB.2), 0), max(max(QLP_QDIF_dat$UB.1), max(QLP_QDIF_dat$UB.2), 0)))
  #Combine Plots #######################################################
  #QLP Coefficient Plots
  QLP_coef_row1 <- plot_grid(QLP_Kplot[[p]], QLP_Lplot[[p]])
  QLP_coef_row2 <- plot_grid(QLP_QDIF_Kplot[[p]], QLP_QDIF_Lplot[[p]])
  QLP_Coef_Plot <- plot_grid(QLP_coef_row1, QLP_coef_row2, ncol=1, align="h", rel_heights = c(1, 1))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/Coefficients/LP/QLP_Coef_Plot_ISIC_", ISIC[p], ".png", sep=""), QLP_Coef_Plot, base_height=8, base_width=7)
}
################################################################################################
###################Coefficients over Time #######################################
################################################################################
T <- 2
dZ <- 2
time <- seq(min(COLdata$year), max(COLdata$year), by=T)
QLPT_Coef <- array(0, dim=c(length(tau_t), dZ, length(time)))
LPT_Coef <- array(0, dim=c(length(time), dZ))
for (t in 1:length(time)){
  for (q in 1:length(tau_t)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/QLPT_COL_Q%s.RData", q))
    QLPT_Coef[,,t][q,] <- betahat[,t]
    if (tau_t[q]==0.5){
      LPT_Coef[t,] <- LPhat[,t]
    }
  }
}
KT <- data.frame(cbind(time, t(QLPT_Coef[,1,][,]), LPT_Coef[,1]))
colnames(KT) <- c("Year", paste("Q", tau_t, sep=" "), "LP")
KT <- melt(KT, "Year")
KTplot <- ggplot(KT, aes(x=Year, y=value, group=variable, linetype=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Capital") + scale_colour_manual(name="", labels=c(tau_t, "LP"), values=c(pcolour, "black"))+theme(legend.text.align = 0) + scale_linetype_manual(name="", labels=c(tau_t, "LP"), values=c("0.1"="solid", "0.25"="solid", "0.5"="solid","0.9"="solid", "LP"="longdash", "NA"))
LT <- data.frame(cbind(time, t(QLPT_Coef[,2,][,]), LPT_Coef[,2]))
colnames(LT) <- c("Year", paste("Q", tau_t, sep=""), "LP")
LT <- melt(LT, "Year")
LTplot <- ggplot(LT, aes(x=Year, y=value, group=variable, linetype=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Labor") + scale_colour_manual(name="", labels=c(tau_t, "LP"), values=c(pcolour, "black"))+theme(legend.text.align = 0) + scale_linetype_manual(name="", labels=c(tau_t, "LP"), values=c("0.1"="solid", "0.25"="solid", "0.5"="solid","0.9"="solid", "LP"="longdash", "NA"))
Plot_Title <- ggdraw() + draw_label("Output Elasticities Over Time", fontface="plain", size=22) 
Time_Plot <- plot_grid(KTplot, LTplot, rel_heights = 0.7)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/Misc/QLP_Time_Plot.png", Time_Plot, base_height=6, base_width=10)

###############################################################################
###############################################################################
###############################################################################
#APPENDIX PLOTS AND TABLES (OPTIONAL)#########################################
###############################################################################
###############################################################################
###############################################################################
#Filter Data for T>=3
COLdata <- COLdata %>% group_by(id) %>% filter(n()>=3)
###############################################################################
#Store QACF Results
##############################################################################
#Store QACF Estimates and Standard Deviations
QACF_betahat <- array(0, c(length(tauvec), dZ, length(ISIC)))
QACF_betaSE <- array(0, c(length(tauvec), dZ, length(ISIC)))
#Store QR Estimates and Standard Deviations
QACF_QR_betahat <- array(0, c(length(tauvec), dZ, length(ISIC)))
QACF_QR_betaSE <- array(0, c(length(tauvec), dZ, length(ISIC)))
#Store QDIF Estimates and Standard Deviations
QACF_QDIF_hat <- array(0, c(length(tauvec), dZ, length(ISIC)))
QACF_QDIF_SE <- array(0, c(length(tauvec), dZ, length(ISIC)))
#Store QTFP Estimates
QACF_QTFP_hat <- list()
#Store ACF TFP Estimates
ACF_TFP_hat <- list()
#Store Omega Estimates
ACF_omega_hat <- list()
#Store Ex-Post Shock Estimates
ACF_eta_hat <- list()
#Store QACF RTS Estimates and Standard Deviations
QACF_RTS <- array(0, c(length(tauvec), length(ISIC)))
QACF_RTS_SE <- array(0, c(length(tauvec), length(ISIC)))
#Store QACF Capital Intensity Estimates and Standard Deviations
QACF_IN <- array(0, c(length(tauvec), length(ISIC)))
QACF_IN_SE <- array(0, c(length(tauvec), length(ISIC)))
#############################################################################
#Store ACF Results
#############################################################################
#Store ACF Estimates and Standard Devitaions
ACF_betahat <- array(0, c(length(ISIC), dZ))
ACF_betaSE <- array(0, c(length(ISIC), dZ))
#Store ACF RTS Standard Deviations
ACF_RTS_SE <- array(0, c(length(ISIC), 1))
#Store ACF Capital Intensity Standard Deviations
ACF_IN_SE <- array(0, c(length(ISIC), 1))
##############################################################################
#Load ACF and QACF Results
#############################################################################@
for (i in 1:length(ISIC)){
  QACF_QTFP_hat[[i]] <- list()
  for (j in 1:length(tauvec)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/QACF/QACF_Boot_COL_Q%s.RData", j))
    #QACF Estimates and Standard Deviations
    QACF_betahat[,,i][j,] <- QACFbetahat[,i]
    QACF_betaSE[,,i][j,] <- apply(QACFbetaboot[,,i], 2, sd)
    #QR Estimates and Standard Deviations
    QACF_QR_betahat[,,i][j,] <- QACFqrhat[,i]
    QACF_QR_betaSE[,,i][j,] <- apply(QACFqrboot[,,i], 2, sd)
    #QACF-QR Estimates and Standard Deviations
    QACF_QDIF_hat[,,i][j,] <- QACFqdifhat[,i]
    QACF_QDIF_SE[,,i][j,] <- apply(QACFqdifboot[,,i], 2, sd)
    #QACF RTS Estimates and Standard Deviations
    QACF_RTS[j,i] <- sum(QACFbetahat[,i])
    QACF_RTS_SE[j,i] <- sd(apply(QACFbetaboot[,,i], 1, sum))
    #QACF Capital Intensity Estimates and Standard Deviations
    QACF_IN[j,i] <- QACFbetahat[,i][1]/QACFbetahat[,i][2]
    QACF_IN_SE[j,i] <- sd(apply(QACFbetaboot[,,i], 1, function(x) x[1]/x[2]))
    #QTFP Estimates
    QACF_QTFP_hat[[i]][[j]] <- exp(QACFTFPhat[[i]])
    #Load the ACF estimates from a single quantile environment: they should be the same across quantiles
    if (tauvec[j]==0.1){
      #ACF Estimates and Standard Deviations
      ACF_betahat[i,] <- ACFhat[,i]
      ACF_betaSE[i,] <- apply(ACFboot[,,i], 2, sd)
      ACF_TFP_hat[[i]] <- exp(ACFTFPhat[[i]])
      ACF_omega_hat[[i]] <- exp(ACFomegahat[[i]])
      ACF_eta_ecdf <- ecdf(LPexpost[[i]])
      ACF_eta_hat[[i]] <- ACF_eta_ecdf(ACFexpost[[i]])
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
QACF_betatable <- data.frame(cbind(rep(tauvec, length(ISIC)), cbind(do.call(rbind, lapply(seq(dim(QACF_betahat)[3]), function(x) QACF_betahat[ , , x])), do.call(rbind, lapply(seq(dim(QACF_betaSE)[3]), function(x) QACF_betaSE[ , , x])))[,c(rbind(c(1:dZ), dZ+(1:dZ)))]))
QACFestimates <- cbind(QACF_betatable, c(QACF_RTS), c(QACF_RTS_SE), c(QACF_IN), c(QACF_IN_SE))
colnames(QACFestimates) <- c('Tau','K',"se_K", 'L', "se_L", 'RTS', 'RTS_SE', 'In', 'In_SE')
#Make an estimates table for ACF Beta Estimates and Standard Deviations
ACF_betatable <- data.frame(cbind(ACF_betahat, ACF_betaSE)[,c(rbind(c(1:dZ), dZ+(1:dZ)))])
ACFestimates <- cbind(ACF_betatable, ACF_RTS, ACF_RTS_SE, ACF_IN, ACF_IN_SE)
colnames(ACFestimates) <- c('K',"se_K", 'L', "se_L", 'RTS', 'RTS_SE', 'In', 'In_SE')

#Table Labels
ISIC_labels <- array(NA, length(tau_t)*length(ISIC)); ISIC_labels[seq(1, length(tau_t)*length(ISIC), by=length(tau_t))] <- ISIC
ISIC_labels[is.na(ISIC_labels)] <- ""

QACF_Table <- cbind(ISIC_labels, QACFestimates[rep(tauvec, length(ISIC))%in%tau_t, ])
colnames(QACF_Table) <- c("ISIC", "$\\tau$", rep(c("Coef.", "s.e"), dZ+2))
QACF_Table_X <- xtable(QACF_Table, digits=c(0,0,2,rep(c(3,4), dZ+2)), type="latex")
align(QACF_Table_X) <- rep('c', 7+dZ*2)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Returns to Scale} & \\multicolumn{2}{c}{Capital Intensity}\\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}'
#For copy pasting to latex
print(QACF_Table_X, hline.after=c(0,nrow(QACF_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#For saving to file
print(QACF_Table_X, hline.after=c(0,nrow(QACF_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H", file="/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Estimates/COL_ACF_Estimates.tex")
############Density Plots for TFP############################
ACFTFPeta <- list()
for (p in 1:length(industries)){
  #Interpolation of TFP from RC Model
  ACFparvec <- t(cbind(split(QACFestimates$K, ceiling(seq_along(QACFestimates$K)/length(tauvec)))[[p]], split(QACFestimates$L, ceiling(seq_along(QACFestimates$L)/length(tauvec)))[[p]]))
  COLeta <- subset(COLdata, str_detect(isic3, industries[p]))
  ACFetafit <- rowSums(cbind(COLeta$lnk, COLeta$lnl)*lsplinemat(tau=tauvec, par=ACFparvec, u=ACF_eta_hat[[p]]))
  ACFTFPeta[[p]] <- COLeta$lnva-ACFetafit
  ACFdensdata <- melt(data.frame(DS=ACFTFPeta[[p]], ACF=ACFTFPhat[[p]]))
  ACFdensplot <- ggplot(ACFdensdata, aes(x=value, color=variable))+geom_density()+xlab("") + ylab("") + scale_colour_manual(name="", labels=c("DS", "ACF"), values=c("black", "red"))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/TFP/ACF/TFP_Dens_Plot_ISIC_", ISIC[p], ".png", sep=""), ACFdensplot, base_height=8, base_width=7)
  #Plot TFP densities per quantile (0.25, 0.5, 0.75)
  ACFTFPtau <- melt(data.frame(tau1=QACFTFPhat[[p]][[5]], tau2=QACFTFPhat[[p]][[10]], tau3=QACFTFPhat[[p]][[15]], ACF=ACFTFPhat[[p]]))
  ACFdenstauplot <- ggplot(ACFTFPtau, aes(x=value, color=variable))+geom_density()+xlab("") + ylab("")
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/TFP/ACF/TFPtau_Dens_Plot_ISIC_", ISIC[p], ".png", sep=""), ACFdenstauplot, base_height=8, base_width=7)
}
####################################################################################################
#TFP Data Formatting
####################################################################################################
#Combine quantiles of QTFP
QACF_QTFP <- lapply(QACF_QTFP_hat,  function(x) do.call(cbind, x))
#Combine with ACF TFP and ACF productivity
QACF_TFP_data <- lapply(1:length(ISIC), function(x) data.frame(cbind(subset(COLdata, str_detect(isic3, industries[x]))$id, subset(COLdata, str_detect(isic3, industries[x]))$year, QACF_QTFP[[x]], exp(ACFTFPeta[[x]]), ACF_TFP_hat[[x]], ACF_omega_hat[[x]])))
QACF_TFP_data <- lapply(QACF_TFP_data, setNames, nm = c("id", "year", paste("Q", tauvec, sep=" "), "TFPeta", "TFP", "Omega"))
#Take un-weighted averages and set base year to 100
QACF_TFP_AVG <- lapply(QACF_TFP_data, function(x) group_by(x, year) %>% summarise_at(c(paste("Q", tauvec, sep=" "), "TFPeta", "TFP", "Omega"), mean, na.rm=TRUE) %>% mutate_at(vars(-year), function(z) z/z[1L]*100))
#Subset Quantiles of interest for plotting and group by quantile
QACF_TFP <- lapply(QACF_TFP_AVG, function(x) melt(x[,c("year", paste("Q", tau_t), "TFP")], "year"))
QACF_TFPeta <- lapply(QACF_TFP_AVG, function(x) melt(x[,c("year", "TFPeta", "TFP")], "year"))
#Plot entire industry sample
pcolour <- brewer.pal(n=length(tau_t), "Spectral")
ACF_TFP_Plot <- ggplot(QACF_TFP[[length(ISIC)]], aes(x=year, y=value, group=variable, linetype=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("") + scale_colour_manual(name="", labels=c(paste("TFP" ,tau_t), "ACF TFP"), values=c(pcolour, "black")) + theme(legend.text.align = 0) + scale_linetype_manual(name="", labels=c(paste("TFP" ,tau_t), "ACF TFP"), values=c("TFP 0.1"="solid", "TFP 0.25"="solid", "TFP 0.5"="solid","TFP 0.9"="solid", "TFP"="longdash", "NA"))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/TFP/ACF/ACF_TFP_Plot.png", ACF_TFP_Plot, base_height=8, base_width=10)
ACF_TFPeta_Plot <- ggplot(QACF_TFPeta[[length(ISIC)]], aes(x=year, y=value, group=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("") + scale_colour_manual(name="", labels=c("DS", "ACF"), values=c("black", "red"))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/TFP/ACF/ACF_TFPeta_Plot.png", ACF_TFPeta_Plot, base_height=8, base_width=10)
#Create table for productivity growth rates
TFP_growth <- lapply(QACF_TFP_AVG, function(x) cbind(x$year, apply(x[,-1], 2, function(y) ((lead(y)/y)-1)*100)))
#Create table for differences across columns
QACF_TFP_tau_t <- lapply(QACF_TFP_AVG, function(x) x[,c("year", paste("Q", rev(tau_t)), "TFP", "Omega")])
QACF_TFP_DIF <- lapply(QACF_TFP_tau_t, function(x) t(apply(x[,-1], 1, function(y) y-lead(y)))[,c(1:(length(tau_t)-1))])
#Coefficient Plots
#Industry ISIC Code Plot Labels
QACF_Kplot <- list(); QACF_Lplot <- list(); QACF_RTSplot <- list()
QACF_QDIF_Kplot <- list(); QACF_QDIF_Lplot <- list()
for (p in 1:length(ISIC)){
  #Plotting data for QACF
  QACFplotcoef <- apply(QACFestimates[c("K", "L", "RTS")], 2, function(x) split(x, ceiling(seq_along(x)/length(tauvec)))[[p]])
  QACFplotsd <- apply(QACFestimates[c("se_K", "se_L", "RTS_SE")], 2, function(x) split(x, ceiling(seq_along(x)/length(tauvec)))[[p]])
  QACFplotCI <- data.frame(LB=QACFplotcoef-QACFplotsd*qnorm(1-alpha/2), UB=QACFplotcoef+QACFplotsd*qnorm(1-alpha/2))
  QACFplotdat <- data.frame(tau=tauvec, QACFplotcoef, QACFplotsd, QACFplotCI)
  #Plotting data for ACF
  ACFplotcoef <- ACFestimates[c("K", "L", "RTS")][p,]
  ACFplotsd <- ACFestimates[c("se_K", "se_L", "RTS_SE")][p,]
  ACFplotCI <- data.frame(LB=ACFplotcoef-ACFplotsd*qnorm(1-alpha/2), UB=ACFplotcoef+ACFplotsd*qnorm(1-alpha/2))
  ACFplotdat <- data.frame(ACFplotcoef, ACFplotsd, ACFplotCI)
  #Plots
  QACF_Kplot[[p]] <- ggplot(QACFplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=LB.K, ymax=UB.K), fill="grey70") + geom_line(aes(y=K)) + geom_hline(yintercept=ACFplotdat$K, linetype='solid', color='red') + geom_hline(yintercept=c(ACFplotdat$LB.K, ACFplotdat$UB.K), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(min(QACFplotdat$LB.K), ACFplotdat$LB.K, ACFplotdat$LB.L, min(QACFplotdat$LB.L)), max(max(QACFplotdat$UB.K), ACFplotdat$UB.K, ACFplotdat$UB.L, max(QACFplotdat$UB.L))))
  QACF_Lplot[[p]] <- ggplot(QACFplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=LB.L, ymax=UB.L), fill="grey70") + geom_line(aes(y=L)) + geom_hline(yintercept=ACFplotdat$L, linetype='solid', color='red') + geom_hline(yintercept=c(ACFplotdat$LB.L, ACFplotdat$UB.L), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(min(QACFplotdat$LB.K), ACFplotdat$LB.K, ACFplotdat$LB.L, min(QACFplotdat$LB.L)), max(max(QACFplotdat$UB.K), ACFplotdat$UB.K, ACFplotdat$UB.L, max(QACFplotdat$UB.L))))
  QACF_RTSplot[[p]] <- ggplot(QACFplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Returns to Scale") + geom_ribbon(aes(ymin=LB.RTS, ymax=UB.RTS), fill="grey70") + geom_line(aes(y=RTS)) + geom_hline(yintercept=ACFplotdat$RTS, linetype='solid', color='red') + geom_hline(yintercept=c(ACFplotdat$LB.RTS, ACFplotdat$UB.RTS), linetype='dashed', color='red')
  #Plotting data for QDIF Plots
  QACF_QDIF_dat <- data.frame(tau=tauvec, coef=QACF_QDIF_hat[,,p], LB=QACF_QDIF_hat[,,p]-QACF_QDIF_SE[,,p]*qnorm(1-alpha/2), UB=QACF_QDIF_hat[,,p]+QACF_QDIF_SE[,,p]*qnorm(1-alpha/2))
  #QDIF Plots
  QACF_QDIF_Kplot[[p]] <- ggplot(QACF_QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_point(aes(y=coef.1)) + geom_errorbar(aes(ymin=LB.1, ymax=UB.1)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(min(QACF_QDIF_dat$LB.1), min(QACF_QDIF_dat$LB.2), 0), max(max(QACF_QDIF_dat$UB.1), max(QACF_QDIF_dat$UB.2), 0)))
  QACF_QDIF_Lplot[[p]] <- ggplot(QACF_QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_point(aes(y=coef.2)) + geom_errorbar(aes(ymin=LB.2, ymax=UB.2)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(min(QACF_QDIF_dat$LB.1), min(QACF_QDIF_dat$LB.2), 0), max(max(QACF_QDIF_dat$UB.1), max(QACF_QDIF_dat$UB.2), 0)))
  #Combine Plots #######################################################
  #QACF Coefficient Plots
  QACF_coef_row1 <- plot_grid(QACF_Kplot[[p]], QACF_Lplot[[p]])
  QACF_coef_row2 <- plot_grid(QACF_QDIF_Kplot[[p]], QACF_QDIF_Lplot[[p]])
  QACF_Coef_Plot <- plot_grid(QACF_coef_row1, QACF_coef_row2, ncol=1, align="h", rel_heights = c(1, 1))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/Coefficients/ACF/QACF_Coef_Plot_ISIC_", ISIC[p], ".png", sep=""), QACF_Coef_Plot, base_height=8, base_width=7)
}
################################################################################################
###################Coefficients over Time #######################################
################################################################################
T <- 2
dZ <- 2
time <- seq(min(COLdata$year), max(COLdata$year), by=T)
QACFT_Coef <- array(0, dim=c(length(tau_t), dZ, length(time)))
ACFT_Coef <- array(0, dim=c(length(time), dZ))
for (t in 1:length(time)){
  for (q in 1:length(tau_t)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Environments/QACF/QACFT_COL_Q%s.RData", q))
    QACFT_Coef[,,t][q,] <- QACFbetahat[,t]
    if (tau_t[q]==0.5){
      ACFT_Coef[t,] <- ACFhat[,t]
    }
  }
}
KT <- data.frame(cbind(time, t(QACFT_Coef[,1,][,]), ACFT_Coef[,1]))
colnames(KT) <- c("Year", paste("Q", tau_t, sep=" "), "ACF")
KT <- melt(KT, "Year")
KTplot <- ggplot(KT, aes(x=Year, y=value, group=variable, linetype=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Capital") + scale_colour_manual(name="", labels=c(tau_t, "ACF"), values=c(pcolour, "black"))+theme(legend.text.align = 0) + scale_linetype_manual(name="", labels=c(tau_t, "ACF"), values=c("0.1"="solid", "0.25"="solid", "0.5"="solid","0.9"="solid", "ACF"="longdash", "NA"))
LT <- data.frame(cbind(time, t(QACFT_Coef[,2,][,]), ACFT_Coef[,2]))
colnames(LT) <- c("Year", paste("Q", tau_t, sep=""), "ACF")
LT <- melt(LT, "Year")
LTplot <- ggplot(LT, aes(x=Year, y=value, group=variable, linetype=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("Labor") + scale_colour_manual(name="", labels=c(tau_t, "ACF"), values=c(pcolour, "black"))+theme(legend.text.align = 0) + scale_linetype_manual(name="", labels=c(tau_t, "ACF"), values=c("0.1"="solid", "0.25"="solid", "0.5"="solid","0.9"="solid", "ACF"="longdash", "NA"))
Plot_Title <- ggdraw() + draw_label("Output Elasticities Over Time", fontface="plain", size=22) 
Time_Plot <- plot_grid(KTplot, LTplot, rel_heights = 0.7)
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/Colombia/Plots/Misc/QACF_Time_Plot.png", Time_Plot, base_height=6, base_width=10)


















