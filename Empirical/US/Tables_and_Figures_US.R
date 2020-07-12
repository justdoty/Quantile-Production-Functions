setwd('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Code/Empirical/US')
library(qwraps2)
library(stringr)
library(dplyr)
#Load US dataset
USdata <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/USdata.csv') %>% 
  select(id, year, Y, K, L, M, naics3) %>% transmute(id=id, year=year, Y=Y/1e6, K=K/1e6, L=L, M=M/1e6, naics3=as.character(naics3), naics2=as.numeric(str_extract(as.character(naics3), "^.{2}")))
#Industries as listed in Estimation_US.R file
All <- "^3"
industries <- c("311|312", "313|314|315|316", "^31", "321", "322|323", "324|325", "326|327", "^32", "331", "332", "333", "334", "335", "336", "337|339", "^33", All)
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
sumind <- c("^31", "^32", "^33", All)
summary <- do.call(cbind, lapply(sumind, function(x) summary_table(dplyr::filter(USdata, str_detect(naics3, x)), format)))
sumnames <- sumind %>% str_replace_all("[|]", ",") %>% str_replace_all("[\\.^]", "")
sumnames[length(sumnames)] <- "All"
print(summary, cnames=sumnames, rtitle="NAICS", align=c("l", rep("c", length(sumind))))
#Industry size breakdown
naicssize <- group_by(USdata, naics2) %>% summarise(Firms=length(unique(id)), Total=n())
print(data.frame(naicssize))
#Choose which industry to select
#Vector of quantiles
tau <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.75, 0.8, 0.9)
#Number of parameters
dZ <- 2
require(abind)
#Load the results from QLP, output is a list of estimates over quantiles for each industry
QLP_results <- list()
QLP_true <- list()
for (i in 1:length(tau)){
	load(sprintf("QLP_Estimation_US%s.RData", i))
	QLP_results[[i]] <- results
  QLP_true[[i]] <- true.beta

}
#Load the results from LP, output is a list of estimates over industries
LP_results <- list()
LP_true <- list()
for (i in 1:length(industries)){
	load(sprintf("US_LP_Estimation%s.RData", i))
	LP_results[[i]] <- results
  LP_true[[i]] <- true.beta.LP
}
#Load the results from OLS and QR, output is a list of estimates over industries
QR_OLS_results <- load('US_QR_OLS.Rdata')
#Prepare estimates for QLP: Outputs a list of quantile estimates over industry (columns)
#and estimates (rows)
#Bootstrapped QLP Coefficient Estimates
QLP_Coef <- array(0, dim=c(2, length(industries), length(tau)))
# QLP_Coef <- lapply(QLP_true, function(x) t(x))
QLP_SE <- array(0, dim=c(2, length(industries), length(tau)))
# QLP_RtS <- lapply(QLP_results, colSums)
QLP_RtS <- array(0, dim=c(length(tau), length(industries)))
QLP_RtS_SE <- array(0, dim=c(length(tau), length(industries)))
QLP_Lower <- array(0, dim=c(2, length(industries), length(tau)))
QLP_Upper <- array(0, dim=c(2, length(industries), length(tau)))
for (i in 1:length(tau)){
	for (j in 1:length(industries)){
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
QLP_RtS <- c(t(array(as.numeric(unlist(QLP_RtS)), dim=c(3, 10))))
QLP_RtS_SE <- c(QLP_RtS_SE)
QLP_Lower <- aperm(QLP_Lower, c(3, 1, 2))
QLP_Upper <- aperm(QLP_Upper, c(3, 1, 2))
#Prepare estimates for LP
#Bootstrapped LP Coefficient Estimates
LP_Coef <- array(0, dim=c(2, length(industries)))
# LP_Coef <- array(as.numeric(unlist(LP_true)), dim=c(2,3))
LP_SE <- array(0, dim=c(2, length(industries)))
LP_Lower <- array(0, dim=c(2, length(industries)))
LP_Upper <- array(0, dim=c(2, length(industries)))
for (i in 1:length(industries)){
	LP_Coef[,i] <- colMeans(LP_results[[i]]) 
	LP_SE[,i] <- apply(LP_results[[i]], 2, sd) 
	LP_Lower[,i] <- apply(LP_results[[i]], 2, function(x) quantile(x, 0.05))
	LP_Upper[,i] <- apply(LP_results[[i]], 2, function(x) quantile(x, 0.95)) 
}
#Prepare estimates for OLS
OLS_Coef <- lapply(lm.soln, function(x) coefficients(x))
OLS_Coef <- lapply(OLS_Coef, function(x) as.numeric(x[2:4]))
OLS_CI <- lapply(lm.soln, function(x) confint(x, c('US$lnk', 'US$lnl'), level=0.9))
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
#For Quantile Regression
QR_Coef <- lapply(seq(dim(qr.soln)[3]), function(x) qr.soln[ , , x])


#Make an estimates table for Quantile GMM
estimates <- data.frame(cbind(rep(tau, length(industries)), cbind(do.call(rbind, QLP_Coef), do.call(rbind, QLP_SE))[,c(rbind(c(1:2), 2+(1:2)))]))
estimates <- cbind(estimates, QLP_RtS, QLP_RtS_SE)
colnames(estimates) <- c('Tau','K',"se_K", 'L', "se_L", 'RtS', 'RtS_SE')
#Make a Confidence Interval Table for Quantile GMM
QLP_CI <- data.frame(cbind(rep(tau, length(industries)), cbind(do.call(rbind, QLP_Lower), do.call(rbind, QLP_Upper))[,c(rbind(c(1:2), 2+(1:2)))]))
colnames(QLP_CI) <- c('Tau', 'Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
#Make an estimates table for LP
estimates_LP <- data.frame(cbind(do.call(rbind, LP_Coef), do.call(rbind, LP_SE))[,c(rbind(c(1:2), 2+(1:2)))])
colnames(estimates_LP) <- c('K', 'se_K', 'L', 'se_L')
#Make a Confidence Interval Table for LP
LP_CI <- data.frame(cbind(do.call(rbind, LP_Lower), do.call(rbind, LP_Upper))[,c(rbind(c(1:2), 2+(1:2)))])
colnames(LP_CI) <- c('Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
#Make an estimates table for OLS
estimates_OLS <- data.frame(do.call(rbind, OLS_Coef))
colnames(estimates_OLS) <- c('K', 'L')
#Make a Confidence Interval Table for OLS
OLS_CI <- data.frame(cbind(do.call(rbind, OLS_Lower), do.call(rbind, OLS_Upper))[,c(rbind(c(1:2), 2+(1:2)))])
colnames(OLS_CI) <- c('Lower_K', 'Upper_K', 'Lower_L', 'Upper_L')
#Make an estimates table for QR
estimates_QR <- data.frame(cbind(rep(tau, length(industries)), do.call(rbind, QR_Coef)))
colnames(estimates_QR) <- c("Tau", "K", "L")
#Prepare estimates for table in paper/presentation
require(xtable)
tau_table <- c(0.1, 0.25, 0.5, 0.75)

#Table Labels
NAICS_labels <- array(NA, length(tau_table)*length(industries)); NAICS_labels[seq(1, length(tau_table)*length(industries), by=length(tau_table))] <- industries
NAICS_labels[is.na(NAICS_labels)] <- ""

estimates_table <- cbind(NAICS_labels, estimates[rep(tau, length(industries))%in%tau_table, ])
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
K_plot <- list(); L_plot <- list()
for (p in 1:length(industries)){
  ky <- split(estimates$K, ceiling(seq_along(estimates$K)/length(tau)))[[p]]
  klow <- split(QLP_CI$Lower_K, ceiling(seq_along(QLP_CI$Lower_K)/length(tau)))[[p]]
  klow_LP <- LP_CI$Lower_K[p]
  # klow_LP <- kcoef[p]-qnorm(0.95)*kse[p]
  kup <- split(QLP_CI$Upper_K, ceiling(seq_along(QLP_CI$Upper_K)/length(tau)))[[p]]
  kup_LP <- LP_CI$Upper_K[p]
  # kup_LP <- kcoef[p]+qnorm(0.95)*kse[p]
  kqr <- split(estimates_QR$K, ceiling(seq_along(estimates_QR$K)/length(tau)))[[p]]
  K_data <- data.frame(x=tau, y=ky, z=estimates_LP$K[p], qr=kqr, lower=klow, upper=kup, lower_LP=klow_LP, upper_LP=kup_LP)
  K_plot[[p]] <- ggplot(K_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Capital") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) + geom_line(aes(y=kqr), linetype="twodash") +  geom_hline(yintercept=K_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(K_data$lower_LP, K_data$upper_LP), linetype='dashed', color='red')
  #Labor
  ly <- split(estimates$L, ceiling(seq_along(estimates$L)/length(tau)))[[p]]
  llow <- split(QLP_CI$Lower_L, ceiling(seq_along(QLP_CI$Lower_L)/length(tau)))[[p]]
  llow_LP <- LP_CI$Lower_L[p]
  # lwlow_LP <- wcoef[p]-qnorm(0.95)*wse[p]
  lup <- split(QLP_CI$Upper_L, ceiling(seq_along(QLP_CI$Upper_L)/length(tau)))[[p]]
  lup_LP <- LP_CI$Upper_L[p]
  # lwup_LP <- wcoef[p]+qnorm(0.95)*wse[p]
  lqr <- split(estimates_QR$L, ceiling(seq_along(estimates_QR$L)/length(tau)))[[p]]
  L_data <- data.frame(x=tau, y=ly, z=estimates_LP$L[p], qr=lqr, lower=llow, upper=lup, lower_LP=llow_LP, upper_LP=lup_LP)
  L_plot[[p]] <- ggplot(L_data, aes(x=x)) + xlab(expression('percentile-'*tau)) + ylab("Labor") + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey70") + geom_line(aes(y=y)) +  geom_line(aes(y=lqr), linetype="twodash") +  geom_hline(yintercept=L_data$z, linetype='solid', color='red') + geom_hline(yintercept=c(L_data$lower_LP, L_data$upper_LP), linetype='dashed', color='red')

}
#Combine plots across industries
#Industry ISIC Code Plot Labels
Industry_1 <- ggdraw() + draw_label("NAICS 31", fontface='plain')
Industry_2 <- ggdraw() + draw_label("NAICS 32", fontface='plain')
Industry_3 <- ggdraw() + draw_label("NAICS 33", fontface='plain')


#Coefficient Plot for 1st industry
Coefficient_Plot_1 <- plot_grid(Industry_1, L_plot[[1]], K_plot[[1]], ncol=1)
save_plot("Coefficient_Plot_1.png", Coefficient_Plot_1, base_height=10, base_width=5)
#Coefficient Plot for 2nd industry
Coefficient_Plot_2 <- plot_grid(Industry_2, L_plot[[2]], K_plot[[2]], ncol=1)
save_plot("Coefficient_Plot_2.png", Coefficient_Plot_2, base_height=10, base_width=5)

#Coefficient Plot for 3rd industry
Coefficient_Plot_3 <- plot_grid(Industry_3, L_plot[[3]], K_plot[[3]], ncol=1)
save_plot("Coefficient_Plot_3.png", Coefficient_Plot_3, base_height=10, base_width=5)


















