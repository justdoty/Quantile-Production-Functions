library(readstata13)
library(dplyr)
####################For merging data from GNR (2020)#################################################
#RGO: Real Gross Output
#RVA: Real Value-Added
#si: log material share of revenue
#K: Capital
#L: Labor
#RI: Real Intermediate Inputs
#Import/Export: Whether a firm exports output or inputs intermediate inputs
#adv: Wheter a firm advertises or not
#hiwag: Whether a firm pays higher than median wages
CHLlist <- list()
datind <- c("311", "321", "322", "331", "381")
for (i in 1:5){
	CHLGNR <- read.dta13(sprintf("/Users/justindoty/Documents/Research/Dissertation/Data/GNR/Chile/data_chi_%s.dta", datind[i])) %>% 
	select(id, year, RGO, RVA, si, K, L, RI, import, export, adv, hiwag) %>% group_by(id) %>%
	filter(RGO!=-1000, RVA!=-1000, si!=-1000, K!=-1000, L!=-1000, RI!=-1000, import!=-1000, export!=-1000, hiwag!=-1000, adv!=-1000) %>%
	transmute(id=id, year=year, Y=RGO, VA=RVA, share=si, K=K, L=L, M=RI, import=import, export=export, adv=adv, hiwag=hiwag, isic3=datind[i]) %>% arrange(id, year)
	CHLlist[[i]] <- data.frame(CHLGNR)
}
CHLGNRdata <- do.call(rbind, CHLlist)
write.csv(CHLdata, "/Users/justindoty/Documents/Research/Dissertation/Data/GNR/Chile/GNRCHLdata.csv", row.names=FALSE)
#####################For the data from Greenstreet########################################
#Note that dropping the NA's in Exports reduced sample size greatly
CHL <- read.dta13("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Chile/Raw_Data/chile_data.dta") %>%
select(plantID, year, salesTotal, salesDeflTotal, totalEmp, totlab, matDeflTotal, capDeflTotal, matTotal, salesTotal, exports, rawmatsi, adverts, ind3Digit) %>%
mutate(id=plantID, year=year, Y=salesDeflTotal, VA=salesDeflTotal-matDeflTotal, K=capDeflTotal, L=totlab, M=matDeflTotal) %>%
	mutate(ExportB=ifelse(exports!=0, 1, 0), ImportB=ifelse(rawmatsi!=0, 1, 0), AdvB=ifelse(adverts!=0, 1, 0), isic3=ind3Digit) %>%
	filter(Y>0, VA>0, K>0, M>0) %>% mutate(VA=Y-M) %>% select(id, year, Y, VA, K, L, M, ExportB, ImportB, AdvB, exports, rawmatsi, adverts, isic3) %>% arrange(id, year)
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

