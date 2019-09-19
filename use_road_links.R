#setwd('~/Desktop/')
## AADFs
raw_aadf <- read.csv('dft_traffic_counts_aadf.csv',stringsAsFactors = F)
dim(raw_aadf)
colnames(raw_aadf)
unique(raw_aadf$year)
raw_aadf$road_letter <- sapply(raw_aadf$road_name,function(x)strsplit(x,'')[[1]][1])
raw_aadf$road_letter[raw_aadf$road_name=='A1'] <- 'M'
raw_aadf$link_length_km <- as.numeric(raw_aadf$link_length_km)
raw_aadf$distance <- raw_aadf$link_length_km*raw_aadf$cars_and_taxis


## city regions definitions
la_table <- read.csv('inputs/mh_regions_lad_lookup.csv')
regions <- unique(la_table$cityregion)
regions <- regions[regions!='']

## get sum of travel for A and M for 2010-2015
tab <- sapply(c('A','M'),function(y)sapply(regions,function(x)
  sum(subset(raw_aadf,year%in%2010:2015&road_letter==y&local_authority_code%in%subset(la_table,cityregion==x)$lad11cd)$distance)))*365/1000
rownames(tab) <- regions


## compare to RTS
rts_estimates <- read.xlsx('190918_data_from_RTS.xlsx',sheetIndex=2,rowIndex = 3:43)
rts_estimates$NA. <- tolower(rts_estimates$NA.)
rts_estimates$NA.[rts_estimates$NA.=='greater manchester combined authority'] <- 'greatermanchester'
rts_estimates$NA.[rts_estimates$NA.=='liverpool city region combined authority'] <- 'liverpool'
rts_estimates$NA.[rts_estimates$NA.=='north east combined authority'] <- 'northeast'
rts_estimates$NA.[rts_estimates$NA.=='sheffield city region combined authority'] <- 'sheffield'
rts_estimates$NA.[rts_estimates$NA.=='west midlands combined authority'] <- 'westmidlands'
rts_tab <- sapply(list(c('Urban A','Rural A'),'Motorway') ,function(y) sapply(regions,function(x) sum(subset(rts_estimates,NA.==x&Road.Type%in%y)$sum)))
colnames(rts_tab) <- c('A','M')
rownames(rts_tab) <- regions
cols <- rainbow(9)
{
x11(height=10,width=5); 
par(mfrow=c(2,1)); 
plot(log(rts_tab[,1]),log(tab[,1]),main='A',col=cols,pch=16,ylim=c(15,18.2),xlim=c(15,18.2),xlab='RTS',ylab='Links'); 
plot(log(rts_tab[,2]),log(tab[,2]),main='M',col=cols,pch=16,ylim=c(11,17),xlim=c(11,17),xlab='RTS',ylab='Links'); 
legend(x=11,y=15,col=cols,legend=regions,pch=16)
}

raw_2010 <- subset(raw_aadf,year==2010&local_authority_code%in%c('E06000022','E06000023','E06000024','E06000025')&road_type=='major')
sapply(c('A','M'),function(x)sum(subset(raw_2010,road_letter==x)$distance))


require(rgdal)
shape <- readOGR(dsn = ".", layer = "2018-MRDB-minimal")

#require(sf)
#shape <- read_sf(dsn = ".", layer = "2018-MRDB-minimal")