library(dplyr)
library(stringr)
library(reshape2)
library(purrr)
#Data for output, materials, and labor deflators from NBER Manufacturing Productivity Database
nberprod <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/naicsdef.csv', header=TRUE)
#Data for capital deflators from the BLS: start date 1987 and 2010 base year. 
# blsprod <- read.csv('cap_details.csv', header=TRUE, check.names=FALSE)
#See which industries have NA's in the NBER database
na.naics <- unique(nberprod[is.na(nberprod$emp),]$naics)
print(na.naics)
######################################################################################################
#Clean the deflator data from NBER and aggregate to 3-digit NAICS
nberdeflators <- select(nberprod, naics, year, emp, pay, invest, cap, piship, pimat, piinv) %>% na.omit() %>%
	#Back out depreciation rates using real capital and investment (Becker et. al. (2016))
	group_by(naics) %>% mutate(dep=ifelse(year==first(year), 0, (cap-(lag(invest/piinv)))/lag(cap))) %>%
	#Calculate average industry wage, real investment flows, and extract first 3 digits of NAICS code
	mutate(avgpay=(pay*1e6)/(emp*1e3), rinvest=invest/piinv, naics=str_extract(as.character(naics), "^.{3}")) %>% group_by(naics, year) %>% 
	#Average deflators according to the 3 digit NAICS code for each year
	summarise(yprice=mean(piship), mprice=mean(pimat), lprice=mean(avgpay), iprice=mean(piinv), dep=mean(dep))

########################################################################################################
#This is earlier code used to prepare a merge with the capital deflators from the BLS
#This was abandoned due to lack of clarity of the deflators and differing base years
#Clean the data with naics code beginning with 31
# naics31 <- select(blsprod, NAICS, Measure, "Asset Category", "Duration Title", paste0(1987:2011)) %>%
# 	rename(naics=NAICS, measure=Measure, category="Asset Category", duration="Duration Title") %>%
# 	filter(str_detect(naics, "^31"), str_detect(measure, "^Investment"), str_detect(category, "^All"), str_detect(duration, "^Indexes")) %>%
# 	mutate(naics=as.character(naics))
# unlist31 <- cbind(paste0(311:316), naics31[rep(1:nrow(naics31), each=2),-1]) 
# colnames(unlist31)[1] <- "naics"
# #Clean the deflator data from BLS and aggregate to 3-digit NAICS
# blsdeflators <- select(blsprod, NAICS, Measure, "Asset Category", "Duration Title", paste0(1987:2011)) %>%
# 	rename(naics=NAICS, measure=Measure, category="Asset Category", duration="Duration Title") %>%
# 	filter(str_detect(naics, "^32|^33"), str_detect(measure, "^Investment"), str_detect(category, "^All"), str_detect(duration, "^Indexes")) %>% 
# 	add_row(unlist31) %>% select(naics, paste0(1987:2011)) %>% melt("naics") %>% arrange(naics) %>% 
# 	mutate(variable=as.numeric(as.character(variable))) %>% rename(year=variable, blskprice=value)

####################################################################################################

#Cleaning and merging computstat data with price deflator data
compstat <- read.csv('/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/Compustat.csv', header=TRUE) %>% rename(year=fyear, employ=emp) %>% 
	select(gvkey, year, sale, employ, ppegt, ppent, cogs, dpact, dp, xsga, capx, naics, fic) %>% group_by(gvkey) %>%
	#remove firms with nonpositive values
	filter(!any(sale<=0), !any(employ<0.01), !any(ppegt<=0), !any(ppent<=0), !any(cogs<=0), !any(xsga<=0), !any(dp<=0), !any(capx<=0)) %>%
	#Only manufacturing firms
	ungroup() %>% filter(str_detect(naics, "^31|^32|^33"), fic=="USA", year>=1961, year<=2010) %>%
	#Merge with NBER price deflator data
	mutate(naics=str_extract(as.character(naics), "^.{3}")) %>% inner_join(nberdeflators, c("naics", "year")) %>%
	#Use PPI method to calculate capital stocks
	group_by(gvkey) %>% mutate(realcap=ifelse(year==first(year), ppent/iprice, 0)) %>% 
	mutate(realcap=ifelse(year==first(year), first(realcap), lag(dep)*lag(realcap)+lag(capx/iprice))) %>% ungroup() %>%
	#Unit changes
	transmute(id=gvkey, year=year, Y=(sale/yprice)*1e6, VA=(sale*1e6-(cogs*1e6+xsga*1e6-dp*1e6-employ*lprice*1e3))/yprice, K=realcap*1e6, L=employ*1e3, M=(cogs*1e6+xsga*1e6-dp*1e6-employ*lprice*1e3)/mprice, I=capx/iprice, dep=dep, naics3=naics) %>% group_by(id) %>%
	#Year-to-year changes for output and inputs, firms with extreme changes are dropped
	mutate(Yratio=ifelse(year==first(year), 0, abs((Y-lag(Y))/lag(Y))), Kratio=ifelse(year==first(year), 0, abs((K-lag(Y))/lag(K))), Lratio=ifelse(year==first(year), 0, abs((L-lag(L))/lag(L))), Mratio=ifelse(year==first(year), 0, abs((M-lag(M))/lag(M))), Iratio=ifelse(year==first(year), 0, abs((I-lag(I))/lag(I)))) %>%
	#Year-to-year changes for total input/output ratio, firms with values extremely different from one are dropped
	mutate(IOratio=ifelse(year==first(year), 0, ((K+L+M-lag(K+L+M))/lag(K+L+M))/((Y-lag(Y))/lag(Y)))) %>%
	filter(!any(M<=0), !any(VA<=0), !any(abs(IOratio)>500),!any(Yratio>500), !any(Kratio>500), !any(Lratio>500), !any(Mratio>500)) %>% ungroup() %>%
	select(id, year, Y, VA, K, L, M, I, naics3)
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
panelT <- max(compstat$year)-min(compstat$year)
print(panelT)
#Total observations
print(nrow(compstat))
#Save clean data
path_out <- '/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/US/'
fileName <- paste(path_out, 'USdata.csv',sep = '')
write.csv(compstat,fileName, row.names=FALSE)

