library(readstata13)
library(dplyr)
chile_data <- read.dta13("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Chile/Raw Data/chile_data.dta") %>% select(plantID, year, isic3, salesDeflTotal, whitecnt, bluecnt, totalEmp, matDeflTotal, capDeflTotal) %>%
	transmute(id=plantID, year=year, isic3=isic3, Y=salesDeflTotal, VA=salesDeflTotal-matDeflTotal, Lw=whitecnt, Lb=bluecnt, L=totalEmp, M=matDeflTotal, K=capDeflTotal) %>%
	group_by(id) %>% filter(!any(VA<=0), !any(L<10), !any(M<=0), !any(K<=0))

#Average number of firms per year
avg_firms <- round(as.numeric(colMeans(group_by(chile_data, year) %>% summarise(firms=n()))[2]))
print(avg_firms)
#Total number of firms in the dataset
unique_firms <- length(unique(chile_data$id))
print(unique_firms)
#Panel time length
panelT <- max(chile_data$year)-min(chile_data$year)
print(panelT)

write.csv(chile_data, "/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Chile/CHLdata.csv", row.names=FALSE)


























