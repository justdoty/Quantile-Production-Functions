library(dplyr)
library(stringr)
library(reshape2)
library(purrr)
#US GDP Deflator and Capital User Cost
macro <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/macro_vars.csv')[,-1]
#Data for other deflators from NBER
nber <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/naicsdef.csv', header=TRUE)
#Clean the NBER data a bit
#I keep emp and pay to construct industry level average wage bill since XLR (staff expense) is often missing in compustat
#EMP is measured in thousdands
#PAY is measured in millions
#I keep cap, invest, and piinv to construct industry level average depreciation rates on capital
nberdef <- nber %>% select(naics, year, emp, pay, cap, invest, piinv) %>% na.omit() %>% group_by(naics) %>% 
mutate(dep=ifelse(year==first(year), 0, (cap-(lag(invest/piinv)))/lag(cap)), avgpay=(pay)/(emp), naics3=str_extract(as.character(naics), "^.{3}")) %>%
group_by(naics3, year) %>% summarise(lprice=mean(avgpay), drate=mean(dep))
#lprice is measured in thousands of dollars per worker (in a year)
#Cleaning and merging computstat data with price deflator data
#This is the main dataset, it includes all variables needed for production function estimation using OP, LP, ACF, or GNR
compstat <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/compustat.csv', header=TRUE) %>% rename(year=fyear, employ=emp) %>% 
	select(gvkey, year, sale, employ, ppegt, ppent, cogs, xsga, capx, naics, fic, oibdp, dp, xad, xrd) %>%
	filter(str_detect(naics, "^31|^32|^33"), fic=="USA", year>=1961, year<=2018) %>%
	transmute(id=gvkey, year=year, sale=sale*1e3, oibdp=oibdp*1e3, cogs=cogs*1e3, xsga=xsga*1e3, ppegt=ppegt*1e3, ppent=ppent*1e3, employ=employ*1e3, capx=capx*1e3, dp=dp*1e3, adv=xad*1e3, rd=xrd*1e3, naics3=str_extract(as.character(naics), "^.{3}")) %>%
	#Merge with US GDP deflator and NBER Data
	inner_join(macro, "year") %>% inner_join(nberdef, c("naics3", "year")) %>% mutate(lexp=employ*lprice) %>% mutate(mexp=cogs+xsga-lexp) %>%
	mutate(id=id, year=year, Y=(sale/USGDP)*100, K=(ppegt/USGDP)*100, K2=(ppent/USGDP)*100, M=(mexp/USGDP)*100, I=(capx/USGDP)*100, L=employ, dp=dp, adv=adv, rd=rd) %>% 
	mutate(VA=((sale-mexp)/USGDP)*100, S=mexp/sale) %>%
	group_by(id) %>% na.omit() %>% filter(Y>0, K>0, K2>0, M>0, I>0, L>0, VA>0, rd>=0, adv>0) %>% group_by(year)
compstat <- compstat %>% select(id, year, Y, K, K2, M, I, L, VA, S, adv, rd, advB, rdB, naics3) %>% group_by(id) %>% filter(n()>=2)
####################################################################################################
#Summary statistics for the cleaned data set
print(summary(compstat))
#Average number of firms per year
avg_firms <- round(as.numeric(colMeans(group_by(compstat, year) %>% summarise(firms=n()))[2]))
print(avg_firms)
#Total number of firms in the dataset
unique_firms <- length(unique(compstat$id))
print(unique_firms)
#Panel time length
panelT <- max(compstat$year)-min(compstat$year)+1
print(panelT)
#Average number of years a firm in the data
avg_time <- summary(group_by(compstat, id) %>% summarise(time=n()))
print(avg_time)
#Total observations
print(nrow(compstat))
#Save clean data
path_out <- '/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/'
fileName <- paste(path_out, 'USdata.csv',sep = '')
write.csv(compstat,fileName, row.names=FALSE)
