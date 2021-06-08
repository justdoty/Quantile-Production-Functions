library(stringr)
library(dplyr)
library(xtable)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
source('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Functions/Aux_Fun.R')
#Load US dataset
USdata <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/USdata.csv') %>% 
  select(id, year, Y, VA, K, L, M, adv, rd, naics3) %>% transmute(id=id, year=year, Y=log(Y), VA=log(VA), K=log(K), L=log(L), M=log(M), adv=adv, rd=rd, naics3=as.character(naics3), naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
#Industries as listed in Estimation_US.R file
NAICS <- c("31", "32", "33", "All")
industries <- c("31", "32", "33", "^3")
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
tau_n <- length(tauvec)
#Prepare estimates for table in paper/presentation
tau_t <- c(0.1, 0.25, 0.5, 0.9)
#Number of parameters
dZLP <- 3
###############################################################################
#Store QLP Results
##############################################################################
#Store QLP Estimates and Standard Deviations
QLP_betahat <- array(0, c(length(tauvec), dZLP, length(NAICS)))
QLP_betaSE <- array(0, c(length(tauvec), dZLP, length(NAICS)))
#Store QR Estimates and Standard Deviations
QLP_QR_betahat <- array(0, c(length(tauvec), dZLP, length(NAICS)))
QLP_QR_betaSE <- array(0, c(length(tauvec), dZLP, length(NAICS)))
#Store QDIF Estimates and Standard Deviations
QLP_QDIF_hat <- array(0, c(length(tauvec), dZLP, length(NAICS)))
QLP_QDIF_SE <- array(0, c(length(tauvec), dZLP, length(NAICS)))
#Store QLP RTS Estimates and Standard Deviations
QLP_RTS <- array(0, c(length(tauvec), length(NAICS)))
QLP_RTS_SE <- array(0, c(length(tauvec), length(NAICS)))
#Store QLP Capital Intensity Estimates and Standard Deviations
QLP_IN <- array(0, c(length(tauvec), length(NAICS)))
QLP_IN_SE <- array(0, c(length(tauvec), length(NAICS)))
#Store QLP Productivity Estimates
QLP_DSPB <- array(0, c(length(tauvec), 2, length(NAICS)))
QLP_DSPB_SE <- array(0, c(length(tauvec), 2, length(NAICS)))
QLP_DSPC <- array(0, c(length(tauvec), 2, length(NAICS)))
QLP_DSPC_SE <- array(0, c(length(tauvec), 2, length(NAICS)))
#############################################################################
#Store LP Results
#############################################################################
#Store LP Estimates and Standard Devitaions
LP_betahat <- array(0, c(length(NAICS), dZLP))
LP_betaSE <- array(0, c(length(NAICS), dZLP))
#Store LP RTS Standard Deviations
LP_RTS_SE <- array(0, c(length(NAICS), 1))
#Store LP Capital Intensity Standard Deviations
LP_IN_SE <- array(0, c(length(NAICS), 1))
#Store LP Productivity Estimates
LP_PB <- array(0, c(length(NAICS), 2))
LP_PB_SE <- array(0, c(length(NAICS), 2))
LP_PC <- array(0, c(length(NAICS), 2))
LP_PC_SE <- array(0, c(length(NAICS), 2))
##############################################################################
#Load LP and QLP Results
#############################################################################@
for (i in 1:tau_n){
  for (j in 1:length(NAICS)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/QLP_Environments/QLP_Boot_US_NAICS%s.RData", j))
    #QLP Estimates and Standard Deviations
    QLP_betahat[,,j][i,] <- QLPbetahat[,i]
    QLP_betaSE[,,j][i,] <- apply(QLPbetaboot[,,i], 2, sd)
    #QR Estimates and Standard Deviations
    QLP_QR_betahat[,,j][i,] <- QLPqrhat[,i]
    QLP_QR_betaSE[,,j][i,] <- apply(QLPqrboot[,,i], 2, sd)
    #QLP-QR Estimates and Standard Deviations
    QLP_QDIF_hat[,,j][i,] <- QLPqdifhat[,i]
    QLP_QDIF_SE[,,j][i,] <- apply(QLPqdifboot[,,i], 2, sd)
    #QLP RTS Estimates and Standard Deviations
    QLP_RTS[i,j] <- sum(QLPbetahat[,i])
    QLP_RTS_SE[i,j] <- sd(apply(QLPbetaboot[,,i], 1, sum))
    #QLP Capital Intensity Estimates and Standard Deviations
    QLP_IN[i,j] <- QLPbetahat[,i][1]/QLPbetahat[,i][2]
    QLP_IN_SE[i,j] <- sd(apply(QLPbetaboot[,,i], 1, function(x) x[1]/x[2]))
    #QLP Productivity Estimates
    QLP_DSPB[,,j][i,] <- QLPDSPB[,i]
    QLP_DSPB_SE[,,j][i,] <- apply(QLPDSPB_boot[,,i], 2, sd)
    QLP_DSPC[,,j][i,] <- QLPDSPC[,i]
    QLP_DSPC_SE[,,j][i,] <- apply(QLPDSPC_boot[,,i], 2, sd)
    #Load the LP estimates
    #LP Estimates and Standard Deviations
    LP_betahat[j,] <- LPhat
    LP_betaSE[j,] <- apply(LPboot, 2, sd)
    #Store LP RTS Standard Deviations
    LP_RTS_SE[j,] <- sd(apply(LPboot, 1, sum))
    #Store LP Capital Intensity Standard Deviations
    LP_IN_SE[j,] <- sd(apply(LPboot, 1, function(x) x[1]/x[2]))
    #LP Productivity Estimates
    LP_PB[j,] <- LPPB
    LP_PB_SE[j,] <- apply(LPPB_boot, 2, sd)
    LP_PC[j,] <- LPPC
    LP_PC_SE[j,] <- apply(LPPC_boot, 2, sd)
  }
}
#LP RTS Estimates
LP_RTS <- as.matrix(apply(LP_betahat, 1, sum))
#LP Capital Intensity Estimates
LP_IN <- as.matrix(apply(LP_betahat, 1, function(x) x[1]/x[2]))
#Make an estimates table for QLP Beta Estimates and Standard Deviations
QLP_betatable <- data.frame(cbind(rep(tauvec, length(NAICS)), cbind(do.call(rbind, lapply(seq(dim(QLP_betahat)[3]), function(x) QLP_betahat[ , , x])), do.call(rbind, lapply(seq(dim(QLP_betaSE)[3]), function(x) QLP_betaSE[ , , x])))[,c(rbind(c(1:dZLP), dZLP+(1:dZLP)))]))
QLPestimates <- cbind(QLP_betatable, c(QLP_RTS), c(QLP_RTS_SE), c(QLP_IN), c(QLP_IN_SE))
colnames(QLPestimates) <- c('Tau','K',"se_K", 'L', "se_L", "M", "se_M", 'RTS', 'RTS_SE', 'In', 'In_SE')
#Make an estimates table for LP Beta Estimates and Standard Deviations
LP_betatable <- data.frame(cbind(LP_betahat, LP_betaSE)[,c(rbind(c(1:dZLP), dZLP+(1:dZLP)))])
LPestimates <- cbind(LP_betatable, LP_RTS, LP_RTS_SE, LP_IN, LP_IN_SE)
colnames(LPestimates) <- c('K',"se_K", 'L', "se_L", "M", "se_M", 'RTS', 'RTS_SE', 'In', 'In_SE')
#Make an estimates table for the productivity differentials (Continuous)
QLP_DSPC_estimates <- data.frame(cbind(rep(tauvec, length(NAICS)), cbind(do.call(rbind, lapply(seq(dim(QLP_DSPC)[3]), function(x) QLP_DSPC[ , , x])), do.call(rbind, lapply(seq(dim(QLP_DSPC_SE)[3]), function(x) QLP_DSPC_SE[ , , x])))[,c(rbind(c(1:2), 2+(1:2)))]))
colnames(QLP_DSPC_estimates) <- c("Tau", "R&D", "R&D_SE", "Adv", "Adv_SE")
LPPC_estimates <- data.frame(cbind(LP_PC, LP_PC_SE)[,c(rbind(c(1:2), 2+(1:2)))])
colnames(LPPC_estimates) <- c("R&D", "R&D_SE", "Adv", "Adv_SE")
#Table Labels
NAICS_labels <- array(NA, length(tau_t)*length(NAICS)); NAICS_labels[seq(1, length(tau_t)*length(NAICS), by=length(tau_t))] <- NAICS
NAICS_labels[is.na(NAICS_labels)] <- ""
#Tables for QLP Estimates
QLP_Table <- cbind(NAICS_labels, QLPestimates[rep(tauvec, length(NAICS))%in%tau_t, ])
colnames(QLP_Table) <- c("NAICS", "$\\tau$", rep(c("Coef.", "s.e"), dZLP+2))
QLP_Table_X <- xtable(QLP_Table, digits=c(0,0,2,rep(c(3,4), dZLP+2)), type="latex")
align(QLP_Table_X) <- rep('c', 7+dZLP*2)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Materials} & \\multicolumn{2}{c}{Returns to Scale} & \\multicolumn{2}{c}{Capital Intensity}\\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10} \\cmidrule(lr){11-12}'
#For copy pasting to latex
print(QLP_Table_X, hline.after=c(0,nrow(QLP_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#Tables for LP Estimates
LP_Table <- cbind(c("31", "32", "33", "All"), LPestimates)
colnames(LP_Table) <- c("NAICS", rep(c("Coef.", "s.e"), dZLP+2))
LP_Table_X <- xtable(LP_Table, digits=c(0,2,rep(c(3,4), dZLP+2)), type="latex")
align(LP_Table_X) <- rep('c', 6+dZLP*2)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & \\multicolumn{2}{c}{Capital} & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Materials} & \\multicolumn{2}{c}{Returns to Scale} & \\multicolumn{2}{c}{Capital Intensity}\\\\ \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \\cmidrule(lr){10-11}'
#For copy pasting to latex
print(LP_Table_X, hline.after=c(0,nrow(LP_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#Tables for Productivity Differentials
QLP_DSPC_table <- cbind(NAICS_labels, QLP_DSPC_estimates[rep(tauvec, length(NAICS))%in%tau_t, ])
colnames(QLP_DSPC_table) <- c("NAICS", "$\\tau$", rep(c("Coef.", "s.e"), 2))
QLP_DSPC_Table_X <- xtable(QLP_DSPC_table, digits=c(0,0,2,rep(c(3,4), 2)), type="latex")
align(QLP_DSPC_Table_X) <- rep('c', 7)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{R\\&D}  & \\multicolumn{2}{c}{Advertising} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6}'
print(QLP_DSPC_Table_X, hline.after=c(0,nrow(QLP_DSPC_Table_X)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#For LP
LPPC_table <- cbind(c("31", "32", "33", "All"), LPPC_estimates)
colnames(LPPC_table) <- c("NAICS", rep(c("Coef.", "s.e"), 2))
LPPC_Table_X <- xtable(LPPC_table, digits=c(0,2,rep(c(3,4), 2)), type="latex")
align(LPPC_Table_X) <- rep('c', 6)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & \\multicolumn{2}{c}{R\\&D}  & \\multicolumn{2}{c}{Advertising} \\\\ \\cmidrule(lr){2-3} \\cmidrule(lr){4-5}'
print(LPPC_Table_X, hline.after=c(0,nrow(LPPC_Table_X)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")


#Industry NAICS Code Plot Labels
QLP_Kplot <- list(); QLP_Lplot <- list(); QLP_Mplot <- list()
QLP_QDIF_Kplot <- list(); QLP_QDIF_Lplot <- list(); QLP_QDIF_Mplot <- list()
LPTFPM <- list(); LPTFPG <- list()
#Plotting Colours
pcolour <- brewer.pal(n=length(tau_t), "Spectral")
for (p in 1:length(NAICS)){
  #Plotting data for QLP
  QLPplotcoef <- apply(QLPestimates[c("K", "L", "M")], 2, function(x) split(x, ceiling(seq_along(x)/length(tauvec)))[[p]])
  QLPplotsd <- apply(QLPestimates[c("se_K", "se_L", "se_M")], 2, function(x) split(x, ceiling(seq_along(x)/length(tauvec)))[[p]])
  QLPplotCI <- data.frame(LB=QLPplotcoef-QLPplotsd*qnorm(1-alpha/2), UB=QLPplotcoef+QLPplotsd*qnorm(1-alpha/2))
  QLPplotdat <- data.frame(tau=tauvec, QLPplotcoef, QLPplotsd, QLPplotCI)
  #Plotting data for LP
  LPplotcoef <- LPestimates[c("K", "L", "M")][p,]
  LPplotsd <- LPestimates[c("se_K", "se_L", "se_M")][p,]
  LPplotCI <- data.frame(LB=LPplotcoef-LPplotsd*qnorm(1-alpha/2), UB=LPplotcoef+LPplotsd*qnorm(1-alpha/2))
  LPplotdat <- data.frame(LPplotcoef, LPplotsd, LPplotCI)
  #Coefficient Plots
  QLP_Kplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=LB.K, ymax=UB.K), fill="grey70") + geom_line(aes(y=K)) + geom_hline(yintercept=LPplotdat$K, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.K, LPplotdat$UB.K), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(min(QLPplotdat$LB.K), min(QLPplotdat$LB.L), min(QLPplotdat$LB.M), LPplotdat$LB.K, LPplotdat$LB.L, LPplotdat$LB.M), max(max(QLPplotdat$UB.K), max(QLPplotdat$UB.L), max(QLPplotdat$UB.M), LPplotdat$UB.K, LPplotdat$UB.L, LPplotdat$UB.M)))
  QLP_Lplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=LB.L, ymax=UB.L), fill="grey70") + geom_line(aes(y=L)) + geom_hline(yintercept=LPplotdat$L, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.L, LPplotdat$UB.L), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(min(QLPplotdat$LB.K), min(QLPplotdat$LB.L), min(QLPplotdat$LB.M), LPplotdat$LB.K, LPplotdat$LB.L, LPplotdat$LB.M), max(max(QLPplotdat$UB.K), max(QLPplotdat$UB.L), max(QLPplotdat$UB.M), LPplotdat$UB.K, LPplotdat$UB.L, LPplotdat$UB.M)))
  QLP_Mplot[[p]] <- ggplot(QLPplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Materials") + geom_ribbon(aes(ymin=LB.M, ymax=UB.M), fill="grey70") + geom_line(aes(y=M)) + geom_hline(yintercept=LPplotdat$M, linetype='solid', color='red') + geom_hline(yintercept=c(LPplotdat$LB.M, LPplotdat$UB.M), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(min(QLPplotdat$LB.K), min(QLPplotdat$LB.L), min(QLPplotdat$LB.M), LPplotdat$LB.K, LPplotdat$LB.L, LPplotdat$LB.M), max(max(QLPplotdat$UB.K), max(QLPplotdat$UB.L), max(QLPplotdat$UB.M), LPplotdat$UB.K, LPplotdat$UB.L, LPplotdat$UB.M)))
  #Plotting data for QDIF Plots (Difference between QR and DS)
  QLP_QDIF_dat <- data.frame(tau=tauvec, coef=QLP_QDIF_hat[,,p], LB=QLP_QDIF_hat[,,p]-QLP_QDIF_SE[,,p]*qnorm(1-alpha/2), UB=QLP_QDIF_hat[,,p]+QLP_QDIF_SE[,,p]*qnorm(1-alpha/2))
  #QDIF Plots
  QLP_QDIF_Kplot[[p]] <- ggplot(QLP_QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_point(aes(y=coef.1)) + geom_errorbar(aes(ymin=LB.1, ymax=UB.1)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(min(QLP_QDIF_dat$LB.1), min(QLP_QDIF_dat$LB.2), min(QLP_QDIF_dat$LB.3), 0), max(max(QLP_QDIF_dat$UB.1), max(QLP_QDIF_dat$UB.2), max(QLP_QDIF_dat$UB.3), 0)))
  QLP_QDIF_Lplot[[p]] <- ggplot(QLP_QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_point(aes(y=coef.2)) + geom_errorbar(aes(ymin=LB.2, ymax=UB.2)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(min(QLP_QDIF_dat$LB.1), min(QLP_QDIF_dat$LB.2), min(QLP_QDIF_dat$LB.3), 0), max(max(QLP_QDIF_dat$UB.1), max(QLP_QDIF_dat$UB.2), max(QLP_QDIF_dat$UB.3), 0)))
  QLP_QDIF_Mplot[[p]] <- ggplot(QLP_QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Materials") + geom_point(aes(y=coef.3)) + geom_errorbar(aes(ymin=LB.3, ymax=UB.3)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(min(QLP_QDIF_dat$LB.1), min(QLP_QDIF_dat$LB.2), min(QLP_QDIF_dat$LB.3), 0), max(max(QLP_QDIF_dat$UB.1), max(QLP_QDIF_dat$UB.2), max(QLP_QDIF_dat$UB.3), 0)))
  #Combine Plots #######################################################
  #QLP Coefficient Plots
  QLP_coef_row1 <- plot_grid(QLP_Kplot[[p]], QLP_Lplot[[p]], QLP_Mplot[[p]], nrow=1)
  QLP_coef_row2 <- plot_grid(QLP_QDIF_Kplot[[p]], QLP_QDIF_Lplot[[p]], QLP_QDIF_Mplot[[p]], nrow=1)
  QLP_Coef_Plot <- plot_grid(QLP_coef_row1, QLP_coef_row2, ncol=1, nrow=2, align="h", rel_heights = c(1, 1))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/Coefficients/LP/QLP_Coef_Plot_NAICS_", NAICS[p], ".png", sep=""), QLP_Coef_Plot, base_height=8, base_width=10)
  # Plot TFP densities at median and mean
  # TFP data for selected quantiles 
  US <- filter(USdata, str_detect(naics2, industries[p]))
  QLPTFP <- apply(QLPplotcoef[match(tau_t, tauvec), ], 1, function(x) US$Y-cbind(US$K, US$L, US$M)%*%x)
  LPTFP <- US$Y-cbind(US$K, US$L, US$M)%*%as.matrix(as.numeric(LPplotcoef))
  LPMdat <- melt(data.frame('TFP=0.1'=QLPTFP[,match(0.1, tau_t)], 'TFP=0.5'=QLPTFP[,match(0.5, tau_t)], 'TFP=0.9'=QLPTFP[,match(0.9, tau_t)] , LP=LPTFP))
  LPTFPM[[p]] <- ggplot(LPMdat, aes(x=value, fill=variable))+geom_density(alpha=0.5) + xlab("TFP") + ylab("")+ guides(fill=guide_legend(title="Estimator"))+ggtitle(paste("NAICS", NAICS[p]))+ theme(plot.title=element_text(face='plain'))
  #Plot TFP growth over time
  LPGdat <- data.frame(US$id, US$year, QLPTFP, LPTFP)
  colnames(LPGdat) <- c("id", "year", paste("QLP", tau_t, sep=""), "LPTFP")
  LPGdat <- LPGdat %>% group_by(year) %>% summarise_at(c(paste("QLP", tau_t, sep=""), "LPTFP"), mean, na.rm=TRUE) %>% mutate_at(vars(-year), function(z) z/z[1L]*100)
  LPGdat <- melt(LPGdat[,c("year", paste("QLP", tau_t, sep=""), "LPTFP")], "year")
  LPTFPG[[p]] <- ggplot(LPGdat, aes(x=year, y=value, group=variable, linetype=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("") + scale_colour_manual(name="", labels=c(paste("TFP" ,tau_t), "LP TFP"), values=c(pcolour, "black")) + theme(legend.text.align = 0) + scale_linetype_manual(name="", labels=c(paste("TFP" ,tau_t), "LP TFP"), values=c("TFP 0.1"="solid", "TFP 0.25"="solid", "TFP 0.5"="solid","TFP 0.9"="solid", "TFP"="longdash", "NA"))+ggtitle(paste("NAICS", NAICS[p]))+ theme(plot.title=element_text(face='plain'))
}
#Combine TFP median vs mean plots over industries and save
LPTFProw1 <- plot_grid(LPTFPM[[1]], LPTFPM[[2]])
LPTFProw2 <- plot_grid(LPTFPM[[3]], LPTFPM[[4]])
LPTFPplot <- plot_grid(LPTFProw1, LPTFProw2, ncol=1, align="h", rel_heights = c(1, 1))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/TFP/QLP_TFP_Plot.png", LPTFPplot, base_height=8, base_width=10)
#Combine TFP growth plots over industries and save
LPTFPgrowthrow1 <- plot_grid(LPTFPG[[1]], LPTFPG[[2]])
LPTFPgrowthrow2 <- plot_grid(LPTFPG[[3]], LPTFPG[[4]])
QLPTFPgrowthplot <- plot_grid(LPTFPgrowthrow1, LPTFPgrowthrow2, ncol=1, align="h", rel_heights = c(1, 1))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/TFP/QLP_TFPgrowth_Plot.png", QLPTFPgrowthplot, base_height=8, base_width=10)
###############################################################################
###############################################################################
###############################################################################
#ACF#########################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#Store QACF Results
##############################################################################
#Number of parameters
dZ <- 2
###############################################################################
#Store QACF Results
##############################################################################
#Store QACF Estimates and Standard Deviations
QACF_betahat <- array(0, c(length(tauvec), dZ, length(NAICS)))
QACF_betaSE <- array(0, c(length(tauvec), dZ, length(NAICS)))
#Store QR Estimates and Standard Deviations
QACF_QR_betahat <- array(0, c(length(tauvec), dZ, length(NAICS)))
QACF_QR_betaSE <- array(0, c(length(tauvec), dZ, length(NAICS)))
#Store QDIF Estimates and Standard Deviations
QACF_QDIF_hat <- array(0, c(length(tauvec), dZ, length(NAICS)))
QACF_QDIF_SE <- array(0, c(length(tauvec), dZ, length(NAICS)))
#Store QACF RTS Estimates and Standard Deviations
QACF_RTS <- array(0, c(length(tauvec), length(NAICS)))
QACF_RTS_SE <- array(0, c(length(tauvec), length(NAICS)))
#Store QACF Capital Intensity Estimates and Standard Deviations
QACF_IN <- array(0, c(length(tauvec), length(NAICS)))
QACF_IN_SE <- array(0, c(length(tauvec), length(NAICS)))
#Store QACF Productivity Estimates
QACF_DSPB <- array(0, c(length(tauvec), 2, length(NAICS)))
QACF_DSPB_SE <- array(0, c(length(tauvec), 2, length(NAICS)))
QACF_DSPC <- array(0, c(length(tauvec), 2, length(NAICS)))
QACF_DSPC_SE <- array(0, c(length(tauvec), 2, length(NAICS)))
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
#Store ACF Productivity Estimates
ACF_PB <- array(0, c(length(NAICS), 2))
ACF_PB_SE <- array(0, c(length(NAICS), 2))
ACF_PC <- array(0, c(length(NAICS), 2))
ACF_PC_SE <- array(0, c(length(NAICS), 2))
##############################################################################
#Load ACF and QACF Results
#############################################################################@
for (i in 1:tau_n){
  for (j in 1:length(NAICS)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/QACF_Environments/QACF_Boot_US_NAICS%s.RData", j))
    #QACF Estimates and Standard Deviations
    QACF_betahat[,,j][i,] <- QACFbetahat[,i]
    QACF_betaSE[,,j][i,] <- apply(QACFbetaboot[,,i], 2, sd)
    #QR Estimates and Standard Deviations
    QACF_QR_betahat[,,j][i,] <- QACFqrhat[,i]
    QACF_QR_betaSE[,,j][i,] <- apply(QACFqrboot[,,i], 2, sd)
    #QACF-QR Estimates and Standard Deviations
    QACF_QDIF_hat[,,j][i,] <- QACFqdifhat[,i]
    QACF_QDIF_SE[,,j][i,] <- apply(QACFqdifboot[,,i], 2, sd)
    #QACF RTS Estimates and Standard Deviations
    QACF_RTS[i,j] <- sum(QACFbetahat[,i])
    QACF_RTS_SE[i,j] <- sd(apply(QACFbetaboot[,,i], 1, sum))
    #QACF Capital Intensity Estimates and Standard Deviations
    QACF_IN[i,j] <- QACFbetahat[,i][1]/QACFbetahat[,i][2]
    QACF_IN_SE[i,j] <- sd(apply(QACFbetaboot[,,i], 1, function(x) x[1]/x[2]))
    #QACF Productivity Estimates
    QACF_DSPB[,,j][i,] <- QACFDSPB[,i]
    QACF_DSPB_SE[,,j][i,] <- apply(QACFDSPB_boot[,,i], 2, sd)
    QACF_DSPC[,,j][i,] <- QACFDSPC[,i]
    QACF_DSPC_SE[,,j][i,] <- apply(QACFDSPC_boot[,,i], 2, sd)
    #Load the ACF estimates
    #ACF Estimates and Standard Deviations
    ACF_betahat[j,] <- ACFhat
    ACF_betaSE[j,] <- apply(ACFboot, 2, sd)
    #Store ACF RTS Standard Deviations
    ACF_RTS_SE[j,] <- sd(apply(ACFboot, 1, sum))
    #Store ACF Capital Intensity Standard Deviations
    ACF_IN_SE[j,] <- sd(apply(ACFboot, 1, function(x) x[1]/x[2]))
    #ACF Productivity Estimates
    ACF_PB[j,] <- ACFPB
    ACF_PB_SE[j,] <- apply(ACFPB_boot, 2, sd)
    ACF_PC[j,] <- ACFPC
    ACF_PC_SE[j,] <- apply(ACFPC_boot, 2, sd)
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
#Make an estimates table for the productivity differentials (Continuous)
QACF_DSPC_estimates <- data.frame(cbind(rep(tauvec, length(NAICS)), cbind(do.call(rbind, lapply(seq(dim(QACF_DSPC)[3]), function(x) QACF_DSPC[ , , x])), do.call(rbind, lapply(seq(dim(QACF_DSPC_SE)[3]), function(x) QACF_DSPC_SE[ , , x])))[,c(rbind(c(1:2), 2+(1:2)))]))
colnames(QACF_DSPC_estimates) <- c("Tau", "R&D", "R&D_SE", "Adv", "Adv_SE")
ACFPC_estimates <- data.frame(cbind(ACF_PC, ACF_PC_SE)[,c(rbind(c(1:2), 2+(1:2)))])
colnames(ACFPC_estimates) <- c("R&D", "R&D_SE", "Adv", "Adv_SE")
#Table Labels
NAICS_labels <- array(NA, length(tau_t)*length(NAICS)); NAICS_labels[seq(1, length(tau_t)*length(NAICS), by=length(tau_t))] <- NAICS
NAICS_labels[is.na(NAICS_labels)] <- ""

QACF_Table <- cbind(NAICS_labels, QACFestimates[rep(tauvec, length(NAICS))%in%tau_t, ])
colnames(QACF_Table) <- c("NAICS", "$\\tau$", rep(c("Coef.", "s.e"), dZ+2))
QACF_Table_X <- xtable(QACF_Table, digits=c(0,0,2,rep(c(3,4), dZ+2)), type="latex")
align(QACF_Table_X) <- rep('c', 7+dZ*2)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Returns to Scale} & \\multicolumn{2}{c}{Capital Intensity}\\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}'
#For copy pasting to latex
print(QACF_Table_X, hline.after=c(0,nrow(QACF_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#Tables for Productivity Differentials
QACF_DSPC_table <- cbind(NAICS_labels, QACF_DSPC_estimates[rep(tauvec, length(NAICS))%in%tau_t, ])
colnames(QACF_DSPC_table) <- c("NAICS", "$\\tau$", rep(c("Coef.", "s.e"), 2))
QACF_DSPC_Table_X <- xtable(QACF_DSPC_table, digits=c(0,0,2,rep(c(3,4), 2)), type="latex")
align(QACF_DSPC_Table_X) <- rep('c', 7)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{R\\&D}  & \\multicolumn{2}{c}{Advertising} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6}'
print(QACF_DSPC_Table_X, hline.after=c(0,nrow(QACF_DSPC_Table_X)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#For ACF
ACFPC_table <- cbind(c("31", "32", "33", "All"), ACFPC_estimates)
colnames(ACFPC_table) <- c("NAICS", rep(c("Coef.", "s.e"), 2))
ACFPC_Table_X <- xtable(ACFPC_table, digits=c(0,2,rep(c(3,4), 2)), type="latex")
align(ACFPC_Table_X) <- rep('c', 6)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & \\multicolumn{2}{c}{R\\&D}  & \\multicolumn{2}{c}{Advertising} \\\\ \\cmidrule(lr){2-3} \\cmidrule(lr){4-5}'
print(ACFPC_Table_X, hline.after=c(0,nrow(ACFPC_Table_X)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")


#Tables for ACF Estimates
ACF_Table <- cbind(c("31", "32", "33", "All"), ACFestimates)
colnames(ACF_Table) <- c("NAICS", rep(c("Coef.", "s.e"), dZ+2))
ACF_Table_X <- xtable(ACF_Table, digits=c(0,2,rep(c(3,4), dZ+2)), type="latex")
align(ACF_Table_X) <- rep('c', 6+dZ*2)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & \\multicolumn{2}{c}{Capital} & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Returns to Scale} & \\multicolumn{2}{c}{Capital Intensity}\\\\ \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}'
#For copy pasting to latex
print(ACF_Table_X, hline.after=c(0,nrow(ACF_Table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x, table.placement="H")
#Industry NAICS Code Plot Labels
QACF_Kplot <- list(); QACF_Lplot <- list(); 
QACF_QDIF_Kplot <- list(); QACF_QDIF_Lplot <- list(); 
ACFTFPM <- list(); ACFTFPG <- list()
#Plotting Colours
pcolour <- brewer.pal(n=length(tau_t), "Spectral")
for (p in 1:length(NAICS)){
  #Plotting data for QACF
  QACFplotcoef <- apply(QACFestimates[c("K", "L")], 2, function(x) split(x, ceiling(seq_along(x)/length(tauvec)))[[p]])
  QACFplotsd <- apply(QACFestimates[c("se_K", "se_L")], 2, function(x) split(x, ceiling(seq_along(x)/length(tauvec)))[[p]])
  QACFplotCI <- data.frame(LB=QACFplotcoef-QACFplotsd*qnorm(1-alpha/2), UB=QACFplotcoef+QACFplotsd*qnorm(1-alpha/2))
  QACFplotdat <- data.frame(tau=tauvec, QACFplotcoef, QACFplotsd, QACFplotCI)
  #Plotting data for ACF
  ACFplotcoef <- ACFestimates[c("K", "L")][p,]
  ACFplotsd <- ACFestimates[c("se_K", "se_L")][p,]
  ACFplotCI <- data.frame(LB=ACFplotcoef-ACFplotsd*qnorm(1-alpha/2), UB=ACFplotcoef+ACFplotsd*qnorm(1-alpha/2))
  ACFplotdat <- data.frame(ACFplotcoef, ACFplotsd, ACFplotCI)
  #Coefficient Plots
  QACF_Kplot[[p]] <- ggplot(QACFplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=LB.K, ymax=UB.K), fill="grey70") + geom_line(aes(y=K)) + geom_hline(yintercept=ACFplotdat$K, linetype='solid', color='red') + geom_hline(yintercept=c(ACFplotdat$LB.K, ACFplotdat$UB.K), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(min(QACFplotdat$LB.K), min(QACFplotdat$LB.L), ACFplotdat$LB.K, ACFplotdat$LB.L), max(max(QACFplotdat$UB.K), max(QACFplotdat$UB.L), ACFplotdat$UB.K, ACFplotdat$UB.L)))
  QACF_Lplot[[p]] <- ggplot(QACFplotdat, aes(x=tau)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=LB.L, ymax=UB.L), fill="grey70") + geom_line(aes(y=L)) + geom_hline(yintercept=ACFplotdat$L, linetype='solid', color='red') + geom_hline(yintercept=c(ACFplotdat$LB.L, ACFplotdat$UB.L), linetype='dashed', color='red') + coord_cartesian(ylim=c(min(min(QACFplotdat$LB.K), min(QACFplotdat$LB.L), ACFplotdat$LB.K, ACFplotdat$LB.L), max(max(QACFplotdat$UB.K), max(QACFplotdat$UB.L), ACFplotdat$UB.K, ACFplotdat$UB.L)))
  #Plotting data for QDIF Plots (Difference between QR and DS)
  QACF_QDIF_dat <- data.frame(tau=tauvec, coef=QACF_QDIF_hat[,,p], LB=QACF_QDIF_hat[,,p]-QACF_QDIF_SE[,,p]*qnorm(1-alpha/2), UB=QACF_QDIF_hat[,,p]+QACF_QDIF_SE[,,p]*qnorm(1-alpha/2))
  #QDIF Plots
  QACF_QDIF_Kplot[[p]] <- ggplot(QACF_QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_point(aes(y=coef.1)) + geom_errorbar(aes(ymin=LB.1, ymax=UB.1)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(min(QACF_QDIF_dat$LB.1), min(QACF_QDIF_dat$LB.2), 0), max(max(QACF_QDIF_dat$UB.1), max(QACF_QDIF_dat$UB.2), 0)))
  QACF_QDIF_Lplot[[p]] <- ggplot(QACF_QDIF_dat, aes(x=tau))+ xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_point(aes(y=coef.2)) + geom_errorbar(aes(ymin=LB.2, ymax=UB.2)) + geom_hline(yintercept=0, linetype='dashed', color='red') + coord_cartesian(ylim =c(min(min(QACF_QDIF_dat$LB.1), min(QACF_QDIF_dat$LB.2), 0), max(max(QACF_QDIF_dat$UB.1), max(QACF_QDIF_dat$UB.2), 0)))
  #Combine Plots #######################################################
  #QACF Coefficient Plots
  QACF_coef_row1 <- plot_grid(QACF_Kplot[[p]], QACF_Lplot[[p]])
  QACF_coef_row2 <- plot_grid(QACF_QDIF_Kplot[[p]], QACF_QDIF_Lplot[[p]])
  QACF_Coef_Plot <- plot_grid(QACF_coef_row1, QACF_coef_row2, ncol=1, align="h", rel_heights = c(1, 1))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/Coefficients/ACF/QACF_Coef_Plot_NAICS_", NAICS[p], ".png", sep=""), QACF_Coef_Plot, base_height=8, base_width=10)
  # Plot TFP densities at median and mean
  # TFP data for selected quantiles 
  US <- filter(USdata, str_detect(naics2, industries[p]))
  QACFTFP <- apply(QACFplotcoef[match(tau_t, tauvec), ], 1, function(x) US$VA-cbind(US$K, US$L)%*%x)
  ACFTFP <- US$VA-cbind(US$K, US$L)%*%as.matrix(as.numeric(ACFplotcoef))
  ACFMdat <- melt(data.frame('TFP=0.1'=QACFTFP[,match(0.1, tau_t)], 'TFP=0.5'=QACFTFP[,match(0.5, tau_t)], 'TFP=0.9'=QACFTFP[,match(0.9, tau_t)] , ACF=ACFTFP))
  ACFTFPM[[p]] <- ggplot(ACFMdat, aes(x=value, fill=variable))+geom_density(alpha=0.5) + xlab("TFP") + ylab("")+ guides(fill=guide_legend(title="Estimator"))+ggtitle(paste("NAICS", NAICS[p]))+ theme(plot.title=element_text(face='plain'))
  #Plot TFP growth over time
  ACFGdat <- data.frame(US$id, US$year, QACFTFP, ACFTFP)
  colnames(ACFGdat) <- c("id", "year", paste("QACF", tau_t, sep=""), "ACFTFP")
  ACFGdat <- ACFGdat %>% group_by(year) %>% summarise_at(c(paste("QACF", tau_t, sep=""), "ACFTFP"), mean, na.rm=TRUE) %>% mutate_at(vars(-year), function(z) z/z[1L]*100)
  ACFGdat <- melt(ACFGdat[,c("year", paste("QACF", tau_t, sep=""), "ACFTFP")], "year")
  ACFTFPG[[p]] <- ggplot(ACFGdat, aes(x=year, y=value, group=variable, linetype=variable)) + geom_line(aes(colour=variable)) + xlab("Year") + ylab("") + scale_colour_manual(name="", labels=c(paste("TFP" ,tau_t), "ACF TFP"), values=c(pcolour, "black")) + theme(legend.text.align = 0) + scale_linetype_manual(name="", labels=c(paste("TFP" ,tau_t), "ACF TFP"), values=c("TFP 0.1"="solid", "TFP 0.25"="solid", "TFP 0.5"="solid","TFP 0.9"="solid", "TFP"="longdash", "NA"))+ggtitle(paste("NAICS", NAICS[p]))+ theme(plot.title=element_text(face='plain'))
}
#Combine TFP median vs mean plots over industries and save
ACFTFProw1 <- plot_grid(ACFTFPM[[1]], ACFTFPM[[2]])
ACFTFProw2 <- plot_grid(ACFTFPM[[3]], ACFTFPM[[4]])
ACFTFPplot <- plot_grid(ACFTFProw1, ACFTFProw2, ncol=1, align="h", rel_heights = c(1, 1))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/TFP/QACF_TFP_Plot.png", ACFTFPplot, base_height=8, base_width=10)
#Combine TFP growth plots over industries and save
ACFTFPgrowthrow1 <- plot_grid(ACFTFPG[[1]], ACFTFPG[[2]])
ACFTFPgrowthrow2 <- plot_grid(ACFTFPG[[3]], ACFTFPG[[4]])
QACFTFPgrowthplot <- plot_grid(ACFTFPgrowthrow1, ACFTFPgrowthrow2, ncol=1, align="h", rel_heights = c(1, 1))
save_plot("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/TFP/QACF_TFPgrowth_Plot.png", QACFTFPgrowthplot, base_height=8, base_width=10)
###########################################################################################






