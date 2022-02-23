library(readstata13)
library(dplyr)
#####################For the data from Greenstreet########################################
#Note that dropping the NA's in Exports reduced sample size greatly
CHL <- read.dta13("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Chile/Raw_Data/chile_data.dta") %>%
select(plantID, year, salesTotal, salesDeflTotal, totalEmp, totlab, payrollTotal, matDeflTotal, capDeflTotal, matTotal, matPayments, salesTotal, exports, rawmatsi, adverts, ind3Digit) %>%
mutate(id=plantID, year=year, Y=salesDeflTotal, sale=salesTotal, vsales=salesTotal-matPayments, VA=salesDeflTotal-matDeflTotal, K=capDeflTotal, L=totlab, lexp=payrollTotal, M=matDeflTotal, mexp=matPayments) %>%
	mutate(ExportB=ifelse(exports!=0, 1, 0), ImportB=ifelse(rawmatsi!=0, 1, 0), AdvB=ifelse(adverts!=0, 1, 0), isic3=ind3Digit) %>%
	filter(Y>0, VA>0, vsales>0, K>0, M>0, mexp>0, lexp>0) %>% mutate(VA=Y-M) %>% select(id, year, Y, sale, VA, vsales, K, L, lexp, M, mexp, ExportB, ImportB, AdvB, exports, rawmatsi, adverts, isic3) %>% arrange(id, year)
CHL <- CHL %>% group_by(id) %>% filter(n()>2)
####################################################################################################
#Summary statistics for the cleaned data set
print(summary(CHL))
#Average number of firms per year
avg_firms <- round(as.numeric(colMeans(group_by(CHL, year) %>% summarise(firms=n()))[2]))
print(avg_firms)
#Total number of firms in the dataset
unique_firms <- length(unique(CHL$id))
print(unique_firms)
#Panel time length
panelT <- max(CHL$year)-min(CHL$year)+1
print(panelT)
#Average number of years a firm in the data
avg_time <- summary(group_by(CHL, id) %>% summarise(time=n()))
print(avg_time)
#Total observations
print(nrow(CHL))

write.csv(CHL, "/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Chile/CHLdata.csv", row.names=FALSE)

