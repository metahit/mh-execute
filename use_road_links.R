
library(dplyr)
library(xlsx)
setwd('~/overflow_dropbox/mh-execute')
## AADFs
raw_aadf <- read.csv('inputs/dft_traffic_counts_aadf.csv',stringsAsFactors = F)
dim(raw_aadf)
colnames(raw_aadf)
unique(raw_aadf$year)
raw_aadf$road_letter <- sapply(raw_aadf$road_category,function(x)strsplit(x,'')[[1]][2])
raw_aadf$link_length_km <- as.numeric(raw_aadf$link_length_km)

## city regions definitions
la_table <- read.csv('inputs/mh_regions_lad_lookup.csv')
regions <- unique(la_table$cityregion)
regions <- regions[regions!='']

## la distances
la_dist <- read.xlsx('inputs/VehicleType_LALevel.xlsx',sheetIndex = 1,rowIndex = 6:1667)
la_dist$LA_Name <- as.character(la_dist$LA_Name)
la_dist$LA_Name[la_dist$LA_Name=='Bristol'] <- 'Bristol, City of'## compare to RTS

## names
aadf_names <- c("pedal_cycles","two_wheeled_motor_vehicles","cars_and_taxis","buses_and_coaches","lgvs","all_hgvs")
la_names <- c("Pedal.Cycles","Two.Wheeled.Motor.Vehicles", "Car","Bus","LGV","HGV")
mh_names <- c('bicycle','motorcycle','car','bus','lgv','hgv')
rts_indices <- c(3,5,6)

## get most recent RTS values
for(i in 1:length(rts_indices)){
  rts_estimates <- read.xlsx('inputs/190918_data_from_RTS.xlsx',sheetIndex=i+1,rowIndex = 3:48)
  #rownames(road_dist) <- sapply(rownames(road_dist),function(x)tolower(gsub(' ','',x)))
  #citymap <- list(bristol='bristol',
  #                nottingham='',
  #                liverpool='liverpoolcityregioncombinedauthority',
  #                northeast='northeastcombinedauthority',
  #                greatermanchester='greatermanchestercombinedauthority',
  #                sheffield='sheffieldcityregioncombinedauthority',
  #                westmidlands='westmidlandscombinedauthority',
  #                leeds='westyorkshirecombinedauthority',
  #                london='london')
  if(i==1){
    rts_estimates$NA. <- tolower(rts_estimates$NA.)
    rts_estimates$Road.Type <- as.character(rts_estimates$Road.Type)
    rts_estimates$NA.[rts_estimates$NA.=='greater manchester combined authority'] <- 'greatermanchester'
    rts_estimates$NA.[rts_estimates$NA.=='liverpool city region combined authority'] <- 'liverpool'
    rts_estimates$NA.[rts_estimates$NA.=='north east combined authority'] <- 'northeast'
    rts_estimates$NA.[rts_estimates$NA.=='sheffield city region combined authority'] <- 'sheffield'
    rts_estimates$NA.[rts_estimates$NA.=='west midlands combined authority'] <- 'westmidlands'
    rts_estimates$NA.[rts_estimates$NA.=='west yorkshire combined authority'] <- 'leeds'
    rts_estimates$Road.Type[rts_estimates$Road.Type=='Rural B,C or Unclassified'] <- 'Rural minor'
    rts_estimates$Road.Type[rts_estimates$Road.Type=='Urban B,C or Unclassified'] <- 'Urban minor'
    rts_est <- rts_estimates[,c(1,2,9)]
    colnames(rts_est) <- c('city','road',mh_names[rts_indices[i]])
  }else {
    rts_est[[mh_names[rts_indices[i]]]] <- rts_estimates[,9]
  }
}

#######################################################
## urban fraction of A roads

buff <- 0
if(file.exists(paste0('inputs/urban_road_fraction_',buff,'.Rds'))&file.exists(paste0('inputs/urban_road_points_',buff,'.Rds'))){
  road_df <- readRDS(paste0('inputs/urban_road_fraction_',buff,'.Rds'))
  point_df <- readRDS(paste0('inputs/urban_road_points_',buff,'.Rds'))
}else{
  library(rgdal)
  library(raster)
  library(rgeos)
  library(spatialEco)
  road_shape <- readOGR(dsn = "shapefiles", layer = "2018-MRDB-minimal")
  urban_shape <- readOGR(dsn = "shapefiles", layer = "Builtup_Areas_December_2011_Boundaries_V2")
  
  crs_string <- "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs"
  
  urban_shp <- spTransform(urban_shape,CRS(crs_string))
  road_shp <- spTransform(road_shape,CRS(crs_string))
  urban_shape_urban <- urban_shp[urban_shp$urban_bua=='Yes',]
  urban_shp_buffered <- buffer(urban_shape_urban,buff)
  
  minor_road_coords <- raw_aadf[,c(1,16,17)]
  coordinates(minor_road_coords) <- c('longitude','latitude')
  proj4string(minor_road_coords) <- CRS("+proj=longlat")
  minor_road_coords <- spTransform(minor_road_coords,CRS(crs_string))
  point_shape <- point.in.poly(minor_road_coords,urban_shp_buffered)
  point_df <- point_shape@data
  colnames(point_df)[2] <- 'urban_point'
  point_df$urban_point[is.na(point_df$urban_point)] <- 0
  saveRDS(point_df,paste0('inputs/urban_road_points_',buff,'.Rds'))
  
  pdf('buffered_urban_area.pdf'); par(mar=c(1,1,1,1))
  plot(urban_shp_buffered,xlim=c(520000 , 550000),ylim=c( 150000,  240000))
  lines(urban_shape_urban,xlim=c(520000 , 550000),ylim=c( 150000,  240000),col='red',lty=2)
  #points(minor_road_coords[!is.na(point_shape@data$poly.ids),],cex=0.5,pch=16,col='grey')
  dev.off()
  
  urban_road <- raster::intersect(road_shp,urban_shp_buffered)
  
  urban_df <- as.data.frame(urban_road)
  urban_df$urban_length <- gLength(urban_road,byid=T)
  road_df <- as.data.frame(road_shp)
  road_df$length <- gLength(road_shp,byid=T)
  road_df <- left_join(road_df,urban_df,by=c('CP_Number','RoadNumber'))
  road_df$urban_length[is.na(road_df$urban_length)] <- 0
  #road_df <- road_df[,c(1,2,3,12,13)]
  road_df$rural_length <- road_df$length - road_df$urban_length
  road_df$urban_fraction <- road_df$urban_length/road_df$length
  colnames(road_df)[1] <- 'count_point_id'
  saveRDS(road_df,paste0('inputs/urban_road_fraction_',buff,'.Rds'))
}
raw_aadf <- left_join(raw_aadf,road_df,by='count_point_id')
raw_aadf <- left_join(raw_aadf,point_df,by='count_point_id')

##########################################################
## compute for modes

tabs_list <- list()
for(mode_number in c(rts_indices,c(1:length(mh_names))[-rts_indices])){
  mh_name <- mh_names[mode_number]
  la_name <- la_names[mode_number]
  aadf_name <- aadf_names[mode_number]
  
  raw_aadf$distance <- raw_aadf$link_length_km*raw_aadf[[aadf_name]]
  raw_aadf$urban_distance <- raw_aadf$distance * raw_aadf$urban_fraction
  raw_aadf$rural_distance <- raw_aadf$distance * (1-raw_aadf$urban_fraction)
  
  
  ## get sum of travel for A and M for 2010-2015
  tab <- t(sapply(regions,function(x) #sapply(c('A','M'),function(y)
  {
    subtab <- subset(raw_aadf,year%in%2010:2015&(local_authority_code%in%subset(la_table,cityregion==x)$lad14cd|
                                                   local_authority_code%in%subset(la_table,cityregion==x)$lad11cd))
    m_dist <- sum(subset(subtab,road_letter=='M')$distance,na.rm=T)
    r_dist <- sum(subset(subtab,road_letter=='A')$rural_distance,na.rm=T)
    u_dist <- sum(subset(subtab,road_letter=='A')$distance,na.rm=T) - r_dist
    c(m_dist,
      u_dist,
      r_dist)
  }
  ))*365/1000
  rownames(tab) <- regions
  colnames(tab) <- c('Motorway','Urban A','Rural A')
  
  
  #########################################################
  ## get minor distances
  if(!mode_number%in%rts_indices){
    la_totals <- list()
    for(city in regions){
      las <- subset(la_table,cityregion==city)$lad11nm
      missing <- las[!las%in%la_dist$LA_Name]
      #if(length(missing)>0){
      #  cat(paste0('Unmatched LAs, distances not extracted for ',city,':\n'))
      #  cat(paste0(missing,'\n'))
      #}else{
      la_sum <- sum(subset(la_dist,Year>2009&LA_Name%in%las)[[la_name]])
      la_totals[[city]] <- la_sum*1.6/1000
      #}
    }
    
    remaining <- unlist(la_totals)-rowSums(tab)
    ## use car change in ratio for bike, bus, motorcycle
    A_car <- tabs_list$car[,2]/(tabs_list$car[,2]+tabs_list$car[,3])
    minor_car <- tabs_list$car[,4]/(tabs_list$car[,4]+tabs_list$car[,5])
    extra <- (minor_car-A_car)/(1-A_car)
    A_ratio <- tab[,2]/(tab[,2]+tab[,3])
    ratio <- A_ratio + extra * (1-A_ratio)
    urban_m <- ratio*remaining
    rural_m <- (1-ratio)*remaining
  }else{
    urban_m <- rts_est[sapply(regions,function(x)which(rts_est$city==x&rts_est$road=='Urban minor')),which(colnames(rts_est)==mh_name)]
    rural_m <- rts_est[sapply(regions,function(x)which(rts_est$city==x&rts_est$road=='Rural minor')),which(colnames(rts_est)==mh_name)]
  }
  tab <- cbind(tab,urban_m,rural_m)
  colnames(tab)[4:5] <- c('Urban minor','Rural minor')
  
  saveRDS(tab,paste0('outputs/',mh_name,'dist2010to2015.Rds'))
  tabs_list[[mh_name]] <- tab
}
write.csv(do.call(rbind,lapply(1:length(tabs_list),function(x)cbind(names(tabs_list)[x],tabs_list[[x]]))),'outputs/mode_road_city.csv')



##########################################################

## replace minor roads


## get minor link lengths
library(readODS)
link_lengths_sheet <- readODS::read.ods('inputs/rdl0202.ods',sheet=13)
link_lengths <- link_lengths_sheet[-c(1:5),-c(4:14,21:26)]
colnames(link_lengths) <- c('lad11cd','region','lad11nm','rural_b','urban_b','rural_c','urban_c','rural_u','urban_u')
link_lengths <- link_lengths[-c(1:2),]
link_lengths <- left_join(link_lengths,la_table,by='lad11cd')#    lad11cd                      lad11nm   lad14cd cityregion                   gordet
link_lengths <- subset(link_lengths,!is.na(cityregion)&cityregion!='')
link_lengths <- link_lengths[,-c(10,13)]
urban_cols <- paste0('urban_',c('b','c','u'))
rural_cols <- paste0('rural_',c('b','c','u'))
for(column in c(urban_cols,rural_cols)) link_lengths[[column]] <- as.numeric(link_lengths[[column]])
#link_lengths$urban_minor <- rowSums(link_lengths[,colnames(link_lengths)%in%urban_cols])
#link_lengths$rural_minor <- rowSums(link_lengths[,colnames(link_lengths)%in%rural_cols])
colnames(link_lengths)[10] <- 'local_authority_code'
link_lengths <- link_lengths[,c(1,3:11)]

link_aadf <- left_join(raw_aadf,link_lengths[,c(3:10)],by='local_authority_code')

minor_aadf <- subset(link_aadf,road_type=='minor'&!is.na(cityregion))
minor_aadf <- minor_aadf[,colnames(minor_aadf)%in%c('count_point_id','year','local_authority_name','local_authority_code','road_name','cars_and_taxis','road_letter',
                                                    'RoadNumber','length','urban_length','rural_length','urban_fraction','urban_point','rural_b','urban_b',
                                                    'rural_c','urban_c','rural_u','urban_u','cityregion')]#c(1,2,7,8,9,24,35:43)]

minor_aadf$road_letter2 <- sapply(minor_aadf$road_name,function(x)strsplit(x,'')[[1]][1])
#x11(); par(mfrow=c(2,2)); for(x in c('U','C','B')) hist(subset(minor_aadf,road_letter2==x)$cars_and_taxis,main=x)
sapply(c('U','C','B'),function(x)nrow(subset(minor_aadf,road_letter2==x)))
sapply(c('U','C','B'),function(x)mean(subset(minor_aadf,road_letter2==x)$cars_and_taxis))

cities <- unique(minor_aadf$cityregion)
for(city in cities){
  print(city)
  lminor <- subset(minor_aadf,cityregion==city)
  sapply(unique(lminor$local_authority_name),function(x)c(nrow(subset(lminor,local_authority_name==x&urban_point==0)),nrow(subset(lminor,local_authority_name==x))))
  print(lapply(c(0,1),function(x) (
    sapply(2010:2015,function(y)
      sum(sapply(unique(lminor$local_authority_name),function(z){
        tab <- subset(lminor,urban_point==x&year==y&local_authority_name==z)
        sapply(1:3,function(w){
          subtab <- subset(tab,road_letter2==c('B','C','U')[w])
          lab <- list(rural_cols,urban_cols)[[x+1]]
          mean(subtab$cars_and_taxis*subtab[[lab[w]]],na.rm=T)*365/1000000
        })
      }),na.rm=T)
    ))))
}

raw_counts <- list()
for(city in cities){
  raw_counts[[city]] <- list()
  lminor <- subset(minor_aadf,cityregion==city)
  for(x in 0:1){
    raw_counts[[city]][[c('rural','urban')[x+1]]]
    for(y in 2010:2015){
      raw_counts[[city]][[c('rural','urban')[x+1]]][[as.character(y)]] <- 
        sapply(unique(lminor$local_authority_name),function(z){
          tab <- subset(lminor,urban_point==x&year==y&local_authority_name==z)
          sapply(1:3,function(w)
            mean(subset(tab,road_letter2==c('B','C','U')[w])$cars_and_taxis)
          )
        })
    }
  }
}

minor_df <- as.data.frame(do.call(rbind,lapply(as.character(cities),function(city){
  lminor <- subset(minor_aadf,cityregion==city)
  do.call(rbind,lapply(0:1,function(x){
    do.call(rbind,lapply(2000:2018,function(y){
      do.call(rbind,lapply(unique(lminor$local_authority_name),function(z){
        tab <- subset(lminor,urban_point==x&year==y&local_authority_name==z)
        t(sapply(1:3,function(w)
          c(city,z,y,x,c('B','C','U')[w],mean(subset(tab,road_letter2==c('B','C','U')[w])$cars_and_taxis))
        ))
      }))
    }))
  }))
})),stringsAsFactors=F)

colnames(minor_df) <- c('city','local_authority_name','year','urban','road','count')
sapply(minor_df,class)
minor_df$year <- as.numeric(minor_df$year)
minor_df$count <- as.numeric(minor_df$count)
minor_df$logcount <- log(minor_df$count)
pred_model <- glm(logcount ~ city+local_authority_name*road*urban+year,data=minor_df)
plot(minor_df$logcount[!is.na(minor_df$logcount)],pred_model$fitted.values)
minor_df$predlogcount <- stats::predict(pred_model,newdata=minor_df)
minor_df$predcount <- exp(minor_df$predlogcount)
plot(minor_df$count,minor_df$predcount)
minor_df <- minor_df[,-c(6:8)]

minor_df_count <- cbind(minor_df[minor_df$urban==0,c(1:3,5)],rural=minor_df$predcount[minor_df$urban==0],urban=minor_df$predcount[minor_df$urban==0])
urban_count <- sapply(c('B','C','U'),function(x) minor_df_count$urban[minor_df_count$road==x])
colnames(urban_count) <- c('urban_b_count','urban_c_count','urban_u_count')
minor_df_count <- cbind(minor_df_count[minor_df_count$road=='B',c(1:3)],sapply(c('B','C','U'),function(x) minor_df_count$rural[minor_df_count$road==x]))
colnames(minor_df_count)[4:6] <- c('rural_b_count','rural_c_count','rural_u_count')
minor_df_count <- cbind(minor_df_count,urban_count)




minor_aadf$local_authority_name[!minor_aadf$local_authority_name%in%minor_df_count$local_authority_name]
link_lengths$local_authority_code[!link_lengths$local_authority_code%in%minor_aadf$local_authority_code]
la_code_name <- unique(minor_aadf[,colnames(minor_aadf)%in%c('local_authority_code','local_authority_name')])


minor_df_count <- left_join(minor_df_count,la_code_name,by='local_authority_name')
minor_df_count <- left_join(minor_df_count,link_lengths[,colnames(link_lengths)%in%c('local_authority_code',"rural_b","urban_b","rural_c","urban_c","rural_u","urban_u")],by='local_authority_code')

for(lab in c('rural','urban')) for(roadletter in c('b','c','u')){
  length_col <- paste0(lab,'_',roadletter)
  count_col <- paste0(length_col,'_count')
  dis_col <- paste0(length_col,'_dist')
  minor_df_count[[dis_col]] <- minor_df_count[[length_col]]*minor_df_count[[count_col]]
}


smooth_counts <- list()
for(cit in cities){
  smooth_counts[[cit]] <- list()
  lminor <- subset(minor_df_count,city==cit)
  for(x in 0:1){
    lab <- c('rural','urban')[x+1]
    smooth_counts[[cit]][[lab]]
    for(y in 2010:2015){
      smooth_counts[[cit]][[lab]][[as.character(y)]] <- 
        sum(sapply(unique(lminor$local_authority_name),function(z){
          tab <- subset(lminor,year==y&local_authority_name==z)
          dis_col <- paste0(lab,'_',c('b','c','u'),'_dist')
          sum(tab[,colnames(tab)%in%dis_col])*365/1e6
        }))
    }
  }
}



###########################################################

## diagnostic plots

for(j in rts_indices){
  rts_tab <- sapply(colnames(tab) ,function(y) sapply(regions,function(x) sum(subset(rts_est,city==x&road%in%y)[[mh_names[j]]])))
  rownames(rts_tab) <- regions
  tab <- tabs_list[[mh_names[j]]]
  cbind(tab,rts_tab)*1e-6
  cols <- rainbow(9)
  {
    pdf(paste0('outputs/Road_vs_link_',mh_names[j],'.pdf'),height=3,width=15); 
    par(mfrow=c(1,5)); 
    for(i in 1:5){
      limits <- range(c(log(rts_tab[,i]),log(tab[,i])))
      plot(log(rts_tab[,i]),log(tab[,i]),main=colnames(tab)[i],col=cols,pch=16,ylim=limits,xlim=limits,xlab='RTS',ylab='Links',cex=2,cex.axis=1.5,cex.lab=1.5); 
    }
    legend(x=mean(limits),y=mean(limits)*1.05,col=cols,legend=regions,pch=16,bty='n')
    dev.off()
  }
  major <- tab[,2]/tab[,3]
  minor <- tab[,4]/tab[,5]
  print(rbind(major,minor,major/minor))
}

