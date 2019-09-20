library(readxl)
library(reshape2)

setwd('~/overflow_dropbox/mh-execute/')

# get data
cardata <- read.xlsx('190918_data_from_RTS.xlsx',sheetIndex=2,rowIndex = 3:43)
# tidy data
colnames(cardata) <- c('region','road',1:6,'total')
cardata <- cardata[,c(1,2,9)]

# reformat data
cities <- unique(cardata$region)
roads <- unique(cardata$road)
carmatrix <- matrix(0,nrow=length(cities),ncol=length(roads),dimnames=list(cities,roads))
for(i in 1:length(cities))
  for(j in 1:length(roads))
    carmatrix[i,j] <- cardata$total[cardata$region==cities[i]&cardata$road==roads[j]]

# save
saveRDS(carmatrix,'car_million_km_2010_to_2015.Rds')
write.csv(carmatrix,'car_million_km_2010_to_2015.csv')

