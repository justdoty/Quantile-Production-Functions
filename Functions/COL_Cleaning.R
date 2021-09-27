library(readstata13)
library(dplyr)
#####################For the data from Greenstreet########################################
#Note that dropping the NA's in Exports reduced sample size greatly
COL <- read.dta13("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Colombia/Raw Data/data/colombia_panel.dta") %>%
select(id, year, salesDeflTotal, capDeflTotal, totalEmp, matDeflTotal, matTotal, salesTotal, s4, s7, c13, isic3) %>% 
mutate(id=id, year=year, Y=salesDeflTotal, K=capDeflTotal, L=totalEmp, M=matDeflTotal) %>% group_by(id) %>%
filter(!is.na(Y), Y>0, Y>M, !is.na(K), K>0, !is.na(L), L>=1, !is.na(M), M>0) %>% mutate(ExC=s4, ExB=ifelse(s4!=0, 1, 0), ImC=s7, ImB=ifelse(s7!=0, 1, 0), AdvC=c13, AdvB=ifelse(c13!=0, 1, 0), isic3=isic3) %>%
mutate(VA=Y-M) %>% select(id, year, Y, VA, K, L, M, ExC, ImC, AdvC, ExB, ImB, AdvB, isic3) %>% arrange(id, year) %>% filter(n()>=3)
# ####################################################################################################
#Summary statistics for the cleaned data set
print(summary(COL))
#Average number of firms per year
avg_firms <- round(as.numeric(colMeans(group_by(COL, year) %>% summarise(firms=n()))[2]))
print(avg_firms)
#Total number of firms in the dataset
unique_firms <- length(unique(COL$id))
print(unique_firms)
#Panel time length
panelT <- max(COL$year)-min(COL$year)+1
print(panelT)
#Average number of years a firm in the data
avg_time <- summary(group_by(COL, id) %>% summarise(time=n()))
print(avg_time)
#Total observations
print(nrow(COL))

write.csv(COL, "/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Colombia/COLdata.csv", row.names=FALSE)












