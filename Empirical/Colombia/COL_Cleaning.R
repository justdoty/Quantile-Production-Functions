library(dplyr)
COL_data <- read.csv("/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Colombia/Raw Data/data.csv") %>%
	select(firm, time, ind, Y, K, L, M, ind) %>% transmute(id=firm, year=time, lny=Y, lnva=log(exp(Y)-exp(M)), lnk=K, lnl=L, lnm=M, isic3=ind) %>%
	group_by(id) %>% filter(!is.na(lnva))

#Average number of firms per year
avg_firms <- round(as.numeric(colMeans(group_by(COL_data, year) %>% summarise(firms=n()))[2]))
print(avg_firms)
#Total number of firms in the dataset
unique_firms <- length(unique(COL_data$id))
print(unique_firms)
#Panel time length
panelT <- max(COL_data$year)-min(COL_data$year)
print(panelT)
print(nrow(COL_data))


path_out <- '/Users/justindoty/Documents/Research/Dissertation/Production_QR_Proxy/Data/Colombia/'
fileName <- paste(path_out, 'COLdata.csv',sep = '')
write.csv(COL_data, fileName, row.names=FALSE)




























