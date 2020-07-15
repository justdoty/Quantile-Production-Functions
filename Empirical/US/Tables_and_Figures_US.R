library(qwraps2)
library(stringr)
library(dplyr)
#Load US dataset
USdata <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/USdata.csv') %>% 
  select(id, year, Y, K, L, M, naics3) %>% transmute(id=id, year=year, Y=Y/1e6, K=K/1e6, L=L, M=M/1e6, naics3=as.character(naics3), naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
#Industries as listed in Estimation_US.R file
All <- "^3"
NAICS <- c("^31", "^32", "^33", All)
########################################################################################################
##########################################Summary Statistics############################################
########################################################################################################
#Formatting for summary statistics using the package qwraps2
format <- 
  list("Output" =
       list("mean" = ~ round(mean(.data$Y), 2),
            "median" = ~ round(median(.data$Y), 2),
            "sd" = ~ round(sd(.data$Y), 2)),
       "Capital" =
       list("mean" = ~ round(mean(.data$K), 2),
            "median" = ~ round(median(.data$K), 2),
            "sd" = ~ round(sd(.data$K), 2)),
       "Labor" =
       list("mean" = ~ round(mean(.data$L)),
            "median" = ~ round(median(.data$L)),
            "sd" = ~ round(sd(.data$L), 2)),
       "Materials" =
       list("mean" = ~ round(mean(.data$M), 2),
            "median" = ~ round(median(.data$M), 2),
            "sd" = ~ round(sd(.data$M), 2)),
       "Size"=
       list("Firms"= ~length(unique(.data$id)),
            "Total"= ~n())
       ) 
sumind <- NAICS
summary <- do.call(cbind, lapply(sumind, function(x) summary_table(dplyr::filter(USdata, str_detect(naics3, x)), format)))
sumnames <- sumind %>% str_replace_all("[|]", ",") %>% str_replace_all("[\\.^]", "")
sumnames[length(sumnames)] <- "All"
print(summary, cnames=sumnames, rtitle="NAICS", align=c("l", rep("c", length(sumind))))
#Industry size breakdown
naicssize <- group_by(USdata, naics2) %>% summarise(Firms=length(unique(id)), Total=n())
print(data.frame(naicssize))
############################################################################################################
#################################Load and prepare data frames for estimates#################################
############################################################################################################
#Vector of quantiles
tau <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.75, 0.8, 0.9)
#Number of parameters
dZ <- 2
R <- 500
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
QLP_Coef <- array(0, dim=c(2, length(NAICS), length(tau)))
# QLP_Coef <- lapply(QLP_true, function(x) t(x))
QLP_SE <- array(0, dim=c(2, length(NAICS), length(tau)))
# QLP_RtS <- lapply(QLP_results, colSums)
QLP_RtS <- array(0, dim=c(length(tau), length(NAICS)))
QLP_RtS_SE <- array(0, dim=c(length(tau), length(NAICS)))
QLP_Lower <- array(0, dim=c(2, length(NAICS), length(tau)))
QLP_Upper <- array(0, dim=c(2, length(NAICS), length(tau)))
for (i in 1:length(tau)){
	for (j in 1:length(NAICS)){
		QLP_Coef[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, mean)
		QLP_SE[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, sd)
    QLP_RtS[i,j] <- mean(apply(QLP_results[[i]][,,j], 1, sum))
    QLP_RtS_SE[i,j] <- sd(apply(QLP_results[[i]][,,j], 1, sum))
		QLP_Lower[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, function(x) quantile(x, 0.05))
		QLP_Upper[,,i][,j] <- apply(QLP_results[[i]][,,j], 2, function(x) quantile(x, 0.95))
	}
}
#A little bit of reformating to obtain estimates over industries instead of quantiles
#Bootstrapped QLP Coefficient Estimates
QLP_Coef <- aperm(QLP_Coef, c(3, 1, 2))
# QLP_Coef <- array(as.numeric(unlist(QLP_Coef)), dim=c(3, 2, 10))
QLP_SE <- aperm(QLP_SE, c(3, 1, 2))
QLP_RtS <- c(t(array(as.numeric(unlist(QLP_RtS)), dim=c(length(NAICS), length(tau)))))
QLP_RtS_SE <- c(QLP_RtS_SE)
QLP_Lower <- aperm(QLP_Lower, c(3, 1, 2))
QLP_Upper <- aperm(QLP_Upper, c(3, 1, 2))
#Prepare estimates for LP
#Bootstrapped LP Coefficient Estimates
LP_Coef <- array(0, dim=c(2, length(NAICS)))
# LP_Coef <- array(as.numeric(unlist(LP_true)), dim=c(2,3))
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
QLP_Coef <- lapply(seq(dim(QLP_Coef)[3]), function(x) QLP_Coef[ , , x])
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
estimates <- data.frame(cbind(rep(tau, length(NAICS)), cbind(do.call(rbind, QLP_Coef), do.call(rbind, QLP_SE))[,c(rbind(c(1:2), 2+(1:2)))]))
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
require(xtable)
tau_table <- c(0.1, 0.25, 0.5, 0.75)

#Table Labels
NAICS[length(NAICS)] <- "All"
NAICS_relabel <- NAICS %>% str_replace_all("[|]", ",") %>% str_replace_all("[\\.^]", "")
NAICS_labels <- array(NA, length(tau_table)*length(NAICS)); NAICS_labels[seq(1, length(tau_table)*length(NAICS), by=length(tau_table))] <- NAICS_relabel
NAICS_labels[is.na(NAICS_labels)] <- ""

estimates_table <- cbind(NAICS_labels, estimates[rep(tau, length(NAICS))%in%tau_table, ])
colnames(estimates_table) <- c("Industry (NAICS code)", "$\\tau$", "Coef.", 's.e.', "Coef.", 's.e.', "Coef.", "s.e")

estimates_table <- xtable(estimates_table, digits=c(0,0,2,3,4,3,4,3,4))
align(estimates_table) <- rep('c', 9)
addtorow <- list()
addtorow$pos <- list(-1)
addtorow$command <- '\\hline\\hline & & \\multicolumn{2}{c}{Capital}  & \\multicolumn{2}{c}{Labor} & \\multicolumn{2}{c}{Returns to Scale} \\\\ \\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}'
print(estimates_table, hline.after=c(0,nrow(estimates_table)), add.to.row=addtorow, auto=FALSE, include.rownames=FALSE, sanitize.text.function=function(x) x)

############################Coefficicent Plots######################################
require(ggplot2)
require(cowplot)
#Industry NAICS Code Plot Labels
QLP_K_plot <- list(); QLP_L_plot <- list()
QR_K_plot <- list(); QR_L_plot <- list()
for (p in 1:length(NAICS)){
  ###################################Plots for QLP and LP ############################
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
  QLP_K_plot[[p]] <- ggplot(QLPK_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QLPK_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QLPK_data$lower_LP, QLPK_data$upper_LP), linetype='dashed', color='red')
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
  QLP_L_plot[[p]] <- ggplot(QLPL_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QLPL_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QLPL_data$lower_LP, QLPL_data$upper_LP), linetype='dashed', color='red')
  ###################################Plots for QR and OLS ############################
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
  QR_K_plot[[p]] <- ggplot(QRK_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QRK_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QRK_data$lower_OLS, QRK_data$upper_OLS), linetype='dashed', color='red')
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
  QR_L_plot[[p]] <- ggplot(QRL_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_hline(yintercept=QRL_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(QRL_data$lower_OLS, QRL_data$upper_OLS), linetype='dashed', color='red')
  ###############################Combine Plots ##############################################
  ############################################################################################
  NAICS_plots <- ggdraw() + draw_label(paste("NAICS", NAICS_relabel[p], sep=" "), fontface="plain") + theme(plot.title = element_text(hjust = 0.5))
  Lrow <- plot_grid(QLP_L_plot[[p]], QR_L_plot[[p]])
  Krow <- plot_grid(QLP_K_plot[[p]], QR_K_plot[[p]])
  Coef_Plot <- plot_grid(NAICS_plots, Lrow, Krow, ncol=1, align="h", rel_heights = c(0.3, 1, 1))
  save_plot(paste("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Plots/Coef_Plot_NAICS_", NAICS_relabel[p], ".png", sep=""), Coef_Plot, base_height=8, base_width=7)
}
##################################Prepare Plots over Time##############################
tau_t <- c(0.1, 0.3, 0.5, 0.7, 0.9)
T <- 5
time <- seq(min(USdata$year), max(USdata$year), by=T)
QLPT_Coef <- array(0, dim=c(length(tau), dZ, length(time)))
QLPT_True <- array(0, dim=c(length(tau), dZ, length(time)))
for (t in 1:length(time)){
  for (q in 1:length(tau)){
    load(sprintf("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US/Environments/QLP_US_Q%s.RData", q))
    QLPT_Coef[,,t][q,] <- colMeans(results_T[,,t])
    QLPT_True[,,t][q,] <- true.beta_T[,t]
  }
}



















