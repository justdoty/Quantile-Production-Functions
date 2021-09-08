library(dplyr)
library(stringr)
library(reshape2)
library(purrr)
library(zoo)
#Data for other deflators from NBER
nbernaics <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/naicsdef.csv', header=TRUE)
macro <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Nonlinear_Production_Function_QR/Data/macro_vars.csv')[,-1]
########################################################################################################################
#Clean the NBER data:
#########################################################################################################################
#I aggregate price deflators up to the 3-digit. I also aggregate up to the 2-digit levels for firms in the Compustat sample that only have the first two-digits of their NAICS code listed
naics3def <- nbernaics %>% select(naics, year, emp, pay, cap, invest, piinv, piship, pimat, vship, matcost, invest) %>% na.omit() %>%  mutate(Y=vship/piship, M=matcost/pimat, I=invest/piinv) %>% 
group_by(naics) %>% mutate(dep=ifelse(year==first(year), 0, (cap-(lag(invest/piinv)))/lag(cap)), avgpay=(pay*1e6)/(emp*1e3), naics3=str_extract(as.character(naics), "^.{3}")) %>% ungroup() %>%
group_by(naics3, year) %>% summarise(lprice=mean(avgpay), drate=mean(dep), yprice=sum(vship)/sum(Y), mprice=sum(matcost)/sum(M), iprice=sum(invest)/sum(I))
naics3def <- data.frame(naics3def)
names(naics3def) <- c("naics3", "year", "lprice3", "drate3", "yprice3", "mprice3", "iprice3")
#For the two-digit level:
naics2def <- nbernaics %>% select(naics, year, emp, pay, cap, invest, piinv, piship, pimat, vship, matcost, invest) %>% na.omit() %>%  mutate(Y=vship/piship, M=matcost/pimat, I=invest/piinv) %>% 
group_by(naics) %>% mutate(dep=ifelse(year==first(year), 0, (cap-(lag(invest/piinv)))/lag(cap)), avgpay=(pay*1e6)/(emp*1e3), naics2=str_extract(as.character(naics), "^.{2}")) %>% ungroup() %>%
group_by(naics2, year) %>% summarise(lprice=mean(avgpay), drate=mean(dep), yprice=sum(vship)/sum(Y), mprice=sum(matcost)/sum(M), iprice=sum(invest)/sum(I))
naics2def <- data.frame(naics2def)
names(naics2def) <- c("naics2", "year", "lprice2", "drate2", "yprice2", "mprice2", "iprice2")
#For deflating capital stock, care needs to be taken in deflating to the year investment was made which may not line up with the current year
#Similar to other authors, I delfate capital by the average age calculated as the ratio dpact/dp which is smoothed by a 3-year rolling average
#I create a separate data file for this
kdef3 <- naics3def %>% select(naics3, year, iprice3) %>% rename(kyear=year, kprice3=iprice3)
kdef2 <- naics2def %>% select(naics2, year, iprice2) %>% rename(kyear=year, kprice2=iprice2)
#Cleaning and Merging with the Compustat Data (renaming employ not to be confused with emp from the price deflator data)
wind <- 0.01
compstat <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/compustat.csv', header=TRUE) %>% rename(id=gvkey, year=fyear, employ=emp, rd=xrd, adv=xad) %>% 
	#Standard screening conditions, only firms incorporated in the US
	filter(indfmt=="INDL", consol=="C", datafmt=="STD", fic=="USA", curcd=="USD", str_detect(naics, "^3"), fic=="USA", year>=1958, year<=2016) %>%
	#Selecting only key variables and constructing age of the firm
	select(id, year, sale, employ, ppegt, ppent, cogs, xsga, capx, naics, oibdp, dp, dpact, rd, adv) %>% group_by(id) %>% mutate(age=year-first(year)+1) %>% ungroup() %>%
	#Truncating 3 and 2 digit NAICS codes
	mutate(naics3=ifelse(nchar(as.character(naics))>=3, substr(as.character(naics), 1, 3), NA), naics2=ifelse(nchar(as.character(naics))==2, substr(as.character(naics), 1, 2), NA)) %>%
	#Merge the deflator data
	left_join(naics3def, by=c("year", "naics3")) %>% left_join(naics2def, by=c("year", "naics2")) %>%
	#Coalesce the deflator data
	mutate(yprice=coalesce(yprice3, yprice2), lprice=coalesce(lprice3, lprice2), mprice=coalesce(mprice3, mprice2), iprice=coalesce(iprice3, iprice2)) %>%
	#Selecting only key variables and filtering depreciation expenses
	select(id, year, sale, employ, ppegt, ppent, cogs, xsga, capx, oibdp, dp, dpact, age, yprice, lprice, mprice, iprice, naics, naics3, naics2, rd, adv) %>% filter(!is.na(dpact), !is.na(dp)) %>%
	#Calculate Average Age of Capital
	mutate(kage=dpact/dp) %>% group_by(id) %>% mutate(ksmooth=if(n()>=3) rollmean(kage, 3, fill=first(kage), align="right") else kage) %>% ungroup() %>% mutate(kyear=year-round(ksmooth)) %>%
	#Merge in capital deflators by average age of capital
	left_join(kdef3, by=c("kyear", "naics3")) %>% left_join(kdef2, by=c("kyear", "naics2")) %>%
	#Coalesce
	mutate(kprice=coalesce(kprice3, kprice2)) %>% mutate(naics=coalesce(naics3, naics2)) %>%
	#Selecting only key variables
	select(id, year, sale, employ, ppegt, ppent, cogs, xsga, capx, oibdp, age, yprice, lprice, mprice, iprice, kprice, naics, rd, adv) %>%
	#Put variable into comparabe units (in thousands)
	mutate(sale=sale*1e3, ppegt=ppegt*1e3, ppent=ppent*1e3, capx=capx*1e3, employ=employ, lexp=lprice, oibdp=oibdp*1e3, rd=rd*1e3, adv=adv*1e3) %>%
	#Filter employees less than 10
	filter(employ>0.01) %>%
	#Omitting NA data and deflating output, capital, and investment and constructing labor expense: COGS and XSGA are filtered in the construction of intermediate inputs
	mutate(Y=sale/yprice, K1=ppegt/kprice, K2=ppent/kprice, L=employ, I=capx/iprice, lexp=employ*lprice, oibdp=oibdp) %>% filter(cogs>0, xsga>0) %>%	
	#Construct interemdiate input expense, deflate with material price deflator and construct deflated value added. Construct 2-digit NAICS codes
	mutate(mexp=(sale-oibdp)-lexp) %>% mutate(M=mexp/mprice) %>% mutate(VA=(oibdp+lexp)/yprice, naics2=str_extract(as.character(naics), "^.{2}")) %>% 
	#Create indicator for RnD and Adv status
	mutate(rdB=ifelse(rd==0, 0, 1), advB=ifelse(adv==0, 0, 1)) %>%
	#Filtering only non-negative observations
	group_by(id) %>% filter(Y>0, !is.na(Y), K1>0, !is.na(K1), K2>0, !is.na(K2), L>0, !is.na(L), M>0, !is.na(M), I>0, !is.na(I), VA>0, !is.na(VA)) %>% 
	select(id, year, Y, VA, K1, K2, L, M, I, age, naics2, rd, adv, rdB, advB)
#The main sample I consider contains firms with more than 3 firm-year observations
compstat <- compstat %>% group_by(id) %>% filter(n()>=2)
USdata <- compstat %>% transmute(id=id, year=year, lny=log(Y), lnva=log(VA), lnk1=log(K1), lnk2=log(K2), lnl=log(L), lnm=log(M), lni=log(I), age=age, rd=rd, rdB=rdB, adv=adv, advB=advB, naics2=naics2)
####################################################################################################
#Summary statistics for the cleaned data set
print(summary(USdata))
#Average number of firms per year
avg_firms <- round(as.numeric(colMeans(group_by(USdata, year) %>% summarise(firms=n()))[2]))
print(avg_firms)
#Total number of firms in the dataset
unique_firms <- length(unique(USdata$id))
print(unique_firms)
#Panel time length
panelT <- max(USdata$year)-min(USdata$year)+1
print(panelT)
#Average number of years a firm in the data
avg_time <- summary(group_by(USdata, id) %>% summarise(time=n()))
print(avg_time)
#Total observations
print(nrow(USdata))
#Industry Specific Sample Sizes
naicsfirms <- group_by(USdata, naics2) %>% summarise(firms=n())
print(naicsfirms)
#Prodest Test
USind <- USdata %>% filter(naics2=="32")
require(prodest)
VA <- prodestLP(Y=USind$lnva, fX=USind$lnl, sX=USind$lnk1, pX=USind$lnm, idvar=USind$id, timevar=USind$year)
print(VA)
# TFP <- USind$lnva-cbind(USind$lnl, USind$lnk1)%*%as.numeric(VA@Estimates$pars)
# print(summary(TFP))
#Save clean data
path_out <- '/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/'
fileName <- paste(path_out, 'USdata.csv',sep = '')
write.csv(USdata, fileName, row.names=FALSE)
