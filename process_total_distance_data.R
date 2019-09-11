library(readxl)
library(reshape2)

setwd('~/overflow_dropbox/mh-execute/')

# get data
cardata <- read_xlsx('190905_data_from_RTS.xlsx',sheet=2)
# tidy data
cardata <- cardata[-c(1:2),]
colnames(cardata) <- c('region','road',1:6)
cardata$total <- rowSums(apply(cardata[,-c(1:2)],2,as.numeric))
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

