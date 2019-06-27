rm(list=ls())
library(ithimr)
library(splines)
#setwd('~/overflow_dropbox/mh-scenarios/')
## overwrite some functions for METAHIT's pp_summary use (instead of TIGTHAT's tripset use)
## in general, the overwriting functions are from ithimr's uncertain_travel branch
## in general, ithimr functions are written ithimr::function()
source('metahit_functions.R')

#####################################################################
## 1 ##
## set variables, which are TIGTHAT studies' input parameters.

## general settings
CITY <- 'bristol'
setup_call_summary_filename <- 'setup_call_summary.txt'
AGE_RANGE <- c(0,150)
REFERENCE_SCENARIO <- 'Baseline'
PATH_TO_LOCAL_DATA <- 'data/'

## placeholders for uncertain parameters
NSAMPLES <- 1
MMET_CYCLING <- 4.63
MMET_WALKING <- 2.53
#PM_CONC_BASE <- 10 
PM_TRANS_SHARE <- 0.225
PA_DOSE_RESPONSE_QUANTILE <- F
AP_DOSE_RESPONSE_QUANTILE <- F
BACKGROUND_PA_SCALAR <- 1
BACKGROUND_PA_CONFIDENCE <- 1
INJURY_REPORTING_RATE <- 1
CHRONIC_DISEASE_SCALAR <- 1
INJURY_LINEARITY <- 1
CASUALTY_EXPONENT_FRACTION <- 0.5
EMISSION_INVENTORY_CONFIDENCE <- 1
DISTANCE_SCALAR_CAR_TAXI <- 1
DISTANCE_SCALAR_WALKING <- 1
DISTANCE_SCALAR_PT <- 1
DISTANCE_SCALAR_CYCLING <- 1
DISTANCE_SCALAR_MOTORCYCLE <- 1

## things we need for ithim-r to work
ADD_TRUCK_DRIVERS <<- F
ADD_BUS_DRIVERS <<- F

## SUMMARY OF INPUTS
# seed = double. sets seed to allow some reproducibility.
# CITY = string. used to identify input files.

# speeds = named list of doubles. average mode speeds.
# emission_inventory = named list of doubles. vehicle emission factors.
# setup_call_summary_filename = string. Where to write input call summary.
# DIST_CAT = vector of strings. defines distance categories for scenario generation (5 accra scenarios)

# AGE_RANGE = vector of length 2, specifying the minimum and maximum ages to be used in the model. Note that the actual 
# maximum and minimum will coincide with boundaries in the population and GBD files.

# ADD_WALK_TO_BUS_TRIPS = logic. T: adds walk trips to all bus trips whose duration exceeds BUS_WALK_TIME. F: no trips added
# ADD_BUS_DRIVERS = logic. T: adds `ghost trips', i.e. trips not taken by any participant. F: no trips added
# ADD_TRUCK_DRIVERS = logic. T: adds `ghost trips', i.e. trips not taken by any participant. F: no trips added

# TEST_WALK_SCENARIO = logic. T: run `scenario 0', one simple scenario where everyone takes one (extra) ten-minute walk trip. F: 5 Accra scenarios.
# TEST_CYCLE_SCENARIO = logic. F: 5 Accra scenarios.
# MAX_MODE_SHARE_SCENARIO = logic. T: run scenarios where we take the maximum mode share across cities and distance categories. F: 5 Accra scenarios.

# REFERENCE_SCENARIO = string: at present, one of 'Baseline' or 'Scenario N' where N is an integer
# PATH_TO_LOCAL_DATA = string: path to input files, if not one of the default case studies 

# NSAMPLES = integer: number of samples to take for each parameter to be sampled

# BUS_WALK_TIME = parameter. double: time taken to walk to bus. vector: samples from distribution.
# MMET_CYCLING = parameter. double: sets cycling (M)METs. vector: samples from distribution.
# MMET_WALKING = parameter. double: sets walking (M)METs. vector: samples from distribution.
# PM_CONC_BASE = parameter. double: sets background PM. vector: samples from distribution.
# PM_TRANS_SHARE = parameter. double: sets PM proportion that comes from transport. vector: samples from distribution.

# PA_DOSE_RESPONSE_QUANTILE = logic. T: PA dose--response relationship is sampled. F: relationship is fixed.
# AP_DOSE_RESPONSE_QUANTILE = logic. T: AP dose--response relationship is sampled. F: relationship is fixed.
# CHRONIC_DISEASE_SCALAR = parameter. double: sets scalar for chronic disease background burden. vector: samples from distribution.

# BACKGROUND_PA_SCALAR = parameter. double: sets scalar for background PA. vector: samples from distribution.
# BACKGROUND_PA_CONFIDENCE = parameter. double between 0 and 1. 1 = use PA data as they are.
# INJURY_REPORTING_RATE = parameter. double: sets scalar for injury counts (inverse). vector: samples from distribution.
# INJURY_LINEARITY = parameter. double: sets scalar. vector: samples from distribution.
# CASUALTY_EXPONENT_FRACTION = parameter. double: sets scalar. vector: samples from distribution.

# DAY_TO_WEEK_TRAVEL_SCALAR = parameter. double: sets scalar for extrapolation from day to week. vector: samples from distribution.
# MOTORCYCLE_TO_CAR_RATIO = parameter. double: sets motorcycle distance relative to car. vector: samples from distribution.
# BUS_TO_PASSENGER_RATIO = parameter. double: sets bus distance relative to bus passenger distance. vector: samples from distribution.
# TRUCK_TO_CAR_RATIO = parameter. double: sets truck distance relative to car. vector: samples from distribution.
# EMISSION_INVENTORY_CONFIDENCE = parameter. double between 0 and 1. 1 = use emission data as they are.
# DISTANCE_SCALAR_CAR_TAXI = double: sets scalar. vector: samples from distribution.
# DISTANCE_SCALAR_WALKING = double: sets scalar. vector: samples from distribution.
# DISTANCE_SCALAR_PT = double: sets scalar. vector: samples from distribution.
# DISTANCE_SCALAR_CYCLING = double: sets scalar. vector: samples from distribution.
# DISTANCE_SCALAR_MOTORCYCLE = double: sets scalar. vector: samples from distribution.

#####################################################################
## 2 ##
## setting all the global variables at the beginning to minimise ITHIM computation
## copied from ithimr::run_ithim_setup

## SET GLOBAL VALUES
## PROGRAMMING VARIABLES

## MODEL VARIABLES
if(is.null(PATH_TO_LOCAL_DATA)){
  PATH_TO_LOCAL_DATA <<- file.path(find.package('ithimr',lib.loc=.libPaths()), 'extdata/local/',CITY,'/') 
}else{
  PATH_TO_LOCAL_DATA <<- PATH_TO_LOCAL_DATA
}

## fixed parameters for AP inhalation
BASE_LEVEL_INHALATION_RATE <<- 1
CLOSED_WINDOW_PM_RATIO <<- 0.5
CLOSED_WINDOW_RATIO <<- 0.9
ROAD_RATIO_MAX <<- 3.216
ROAD_RATIO_SLOPE <<- 0.379
SUBWAY_PM_RATIO <<- 0.8

## default speeds that can be edited by input. 
default_speeds <- list(
  bus=15,
  bus_driver=15,
  car=21,
  taxi=21,
  walking=4.8,
  bicycle=14.5,
  motorcycle=25,
  truck=21,
  van=15,
  subway=28,
  rail=35,
  shared_taxi=21
)
TRAVEL_MODES <<- tolower(names(default_speeds))
MODE_SPEEDS <<- data.frame(stage_mode = TRAVEL_MODES, speed = unlist(default_speeds), stringsAsFactors = F)

## default emission contributions that can be edited by input. 
default_emission_inventory <- list(
  bus=0,
  bus_driver=0.82,
  car=0.228,
  taxi=0.011,
  walking=0,
  bicycle=0,
  motorcycle=0.011,
  truck=0.859,
  big_truck=0.711,
  other=0.082
)
#names(default_emission_inventory) <- tolower(names(default_emission_inventory))
EMISSION_INVENTORY <<- default_emission_inventory

#####################################################################
## LOAD DATA
## copied from ithimr ithim_load_data
global_path <- file.path(find.package('ithimr',lib.loc=.libPaths()), 'extdata/global/')
## for windows??
global_path <- paste0(global_path, "/")

## DATA FILES FOR MODEL  
DISEASE_INVENTORY <<- read.csv(paste0(global_path,"dose_response/disease_outcomes_lookup.csv"))
DR_AP <<- read.csv(paste0(global_path,"dose_response/drap/dose_response.csv"))
# root of list_of_files matches DISEASE_INVENTORY$pa_acronym
list_of_files <- list.files(path = paste0(global_path,"dose_response/drpa/extdata/"), recursive = TRUE, pattern = "\\.csv$", full.names = TRUE)
for (i in 1:length(list_of_files)){
  assign(stringr::str_sub(basename(list_of_files[[i]]), end = -5),
         readr::read_csv(list_of_files[[i]],col_types = cols()),
         pos = 1)
}

background_pollution <- read.csv('../mh-air-pollution/02_DataCreated/1_apmeans.csv')
PM_CONC_BASE <- background_pollution$apmean_bpm25[grepl(CITY,tolower(background_pollution$apgroup_name))]

#####################################################################
## these datasets are all local, saved in local folder.
## there will be one folder per city. this block will have to loop over CITY.
## OR we have one file with, e.g., all the GBD data in.
local_path <- PATH_TO_LOCAL_DATA


demography <- readxl::read_xlsx('190330_sp_ind_codebook.xlsx',sheet=2,col_names=F)
demogindex_to_numerical <- unlist(demography[,3])
demography[,3] <- 1:nrow(demography)
demo_indices <- unlist(demography[,3])
age_table <- readxl::read_xlsx('190330_sp_ind_codebook.xlsx',sheet=1,col_names=F)
age_category <- unlist(age_table[,1])
age_lower_bounds <- as.numeric(sapply(age_category,function(x)strsplit(x,' to ')[[1]][1]))


## DATA FILES FOR CITY
filename <- paste0(local_path,CITY,"/synth_pop.Rds")
synth_pop <- readRDS(filename)

# GBD file needs to have the following columns: 
# age (=label, e.g. 15-49)
# sex (=male or female)
# measure
# cause (GBD_DATA$cause matches DISEASE_INVENTORY$GBD_name)
# metric
# burden
filename <- paste0(local_path,CITY,"/gbd.csv")
GBD_DATA <- read_csv(filename,col_types = cols())
filename <- paste0(local_path,CITY,"/population.csv")
demographic <- read_csv(filename,col_types = cols())
demographic$dem_index <- 1:nrow(demographic)

## find min and max age from AGE_RANGE, trips, and demographic.
age_category <- demographic$age
max_age <- max(as.numeric(sapply(age_category,function(x)strsplit(x,'-')[[1]][2])))
max_age <- min(max_age,AGE_RANGE[2])
min_age <- min(as.numeric(sapply(age_category,function(x)strsplit(x,'-')[[1]][1])))
min_age <- max(min_age,AGE_RANGE[1])
demographic <- demographic[as.numeric(sapply(age_category,function(x)strsplit(x,'-')[[1]][1]))<=max_age&
                             as.numeric(sapply(age_category,function(x)strsplit(x,'-')[[1]][2]))>=min_age,]
POPULATION <<- demographic
demographic <- demographic[,names(demographic)!='population']
names(demographic)[which(names(demographic)=='age')] <- 'age_cat'
DEMOGRAPHIC <<- demographic

# get age-category details from (modified) population data
AGE_CATEGORY <<- unique(POPULATION$age)
AGE_LOWER_BOUNDS <<- as.numeric(sapply(AGE_CATEGORY,function(x)strsplit(x,'-')[[1]][1]))
MAX_AGE <<- max(as.numeric(sapply(AGE_CATEGORY,function(x)strsplit(x,'-')[[1]][2])))

## now process GBD_DATA
# keep named subset of diseases
disease_names <- c(as.character(DISEASE_INVENTORY$GBD_name),'Road injuries')
GBD_DATA <- subset(GBD_DATA,cause_name%in%disease_names)
# keep entries in correct age range
GBD_DATA$min_age <- as.numeric(sapply(GBD_DATA$age_name,function(x)str_split(x,' to ')[[1]][1]))
GBD_DATA$max_age <- as.numeric(sapply(GBD_DATA$age_name,function(x)str_split(x,' to ')[[1]][2]))
GBD_DATA <- subset(GBD_DATA,max_age>=AGE_LOWER_BOUNDS[1])
GBD_DATA <- subset(GBD_DATA,min_age<=MAX_AGE)
##!! hard-coded rename...
names(GBD_DATA)[c(1,3,4,5)] <- c('measure','sex','age','cause')
# ensure lower case
GBD_DATA$sex <- tolower(GBD_DATA$sex)

## get burden of disease for each city by scaling according to population
burden_of_disease <- expand.grid(measure=unique(GBD_DATA$measure),sex=unique(POPULATION$sex),age=unique(POPULATION$age),
                                 cause=disease_names,stringsAsFactors = F)
burden_of_disease <- left_join(burden_of_disease,POPULATION,by=c('age','sex'))
burden_of_disease$min_age <- as.numeric(sapply(burden_of_disease$age,function(x)str_split(x,'-')[[1]][1]))
burden_of_disease$max_age <- as.numeric(sapply(burden_of_disease$age,function(x)str_split(x,'-')[[1]][2]))
## when we sum ages, we assume that all age boundaries used coincide with the GBD age boundaries.
##!! this isn't the case for metahit: age category 15-19 vs 16-19. therefore, have added '-1' for now.
burden_of_disease$rate <- apply(burden_of_disease,1,
                                function(x){
                                  subtab <- subset(GBD_DATA,measure==as.character(x[1])&sex==as.character(x[2])&cause==as.character(x[4])&
                                                     min_age>=as.numeric(x[7])-1&max_age<=as.numeric(x[8])); 
                                  sum(subtab$val)/sum(subtab$population)
                                }
)
burden_of_disease$burden <- burden_of_disease$population*burden_of_disease$rate
##!! if an entry is missing in GBD, we set it to zero. we should also issue a warning.
burden_of_disease$burden[is.na(burden_of_disease$burden)] <- 0
DISEASE_BURDEN <<- burden_of_disease

## for tigthat, use GBD to scale from fatalities to YLL. calculate this ratio here.
gbd_injuries <- DISEASE_BURDEN[which(DISEASE_BURDEN$cause == "Road injuries"),]
gbd_injuries$sex_age <- paste0(gbd_injuries$sex,"_",gbd_injuries$age)
## calculating the ratio of YLL to deaths for each age and sex group
gbd_injuries <- arrange(gbd_injuries, measure)
gbd_inj_yll <- gbd_injuries[which(gbd_injuries$measure == "YLLs (Years of Life Lost)"),]
gbd_inj_dth <- gbd_injuries[which(gbd_injuries$measure == "Deaths"),]
gbd_inj_yll$yll_dth_ratio <- gbd_inj_yll$burden/gbd_inj_dth$burden 
GBD_INJ_YLL <<- gbd_inj_yll


#####################################################################
## 3 ##
## SET PARAMETERS
parameters <- ithimr::ithim_setup_parameters(NSAMPLES=NSAMPLES,
                                                  MMET_CYCLING=MMET_CYCLING,
                                                  MMET_WALKING=MMET_WALKING,
                                                  PM_CONC_BASE=PM_CONC_BASE,  
                                                  PM_TRANS_SHARE=PM_TRANS_SHARE,
                                                  PA_DOSE_RESPONSE_QUANTILE=PA_DOSE_RESPONSE_QUANTILE,
                                                  AP_DOSE_RESPONSE_QUANTILE=AP_DOSE_RESPONSE_QUANTILE,
                                                  BACKGROUND_PA_SCALAR=BACKGROUND_PA_SCALAR,
                                                  BACKGROUND_PA_CONFIDENCE=BACKGROUND_PA_CONFIDENCE,
                                                  INJURY_REPORTING_RATE=INJURY_REPORTING_RATE,
                                                  CHRONIC_DISEASE_SCALAR=CHRONIC_DISEASE_SCALAR,
                                                  INJURY_LINEARITY=INJURY_LINEARITY,
                                                  CASUALTY_EXPONENT_FRACTION=CASUALTY_EXPONENT_FRACTION,
                                                  EMISSION_INVENTORY_CONFIDENCE=EMISSION_INVENTORY_CONFIDENCE,
                                                  DISTANCE_SCALAR_CAR_TAXI=DISTANCE_SCALAR_CAR_TAXI,
                                                  DISTANCE_SCALAR_WALKING=DISTANCE_SCALAR_WALKING,
                                                  DISTANCE_SCALAR_PT=DISTANCE_SCALAR_PT,
                                                  DISTANCE_SCALAR_CYCLING=DISTANCE_SCALAR_CYCLING,
                                                  DISTANCE_SCALAR_MOTORCYCLE=DISTANCE_SCALAR_MOTORCYCLE)

ithimr::set_vehicle_inventory() # sets vehicle inventory

#####################################################################
## start metahit
#####################################################################
## 4 ##

## convert synth pop to ithim-r style
names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('wkhr','dur',x))
names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('wkkm','dist',x))
names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('walk_','walking_',x))
names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('_cycle','_bicycle',x))
names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('cardrive','car',x))
names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('mbikedrive','motorcycle',x))
names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('mbike','motorcycle',x))
names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('tube','subway',x))
names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('train','rail',x))
synth_pop$participant_id <- 1:nrow(synth_pop)
synthetic_pop <- synth_pop[,names(synth_pop)%in%c('participant_id','dem_index')]
demographic <- DEMOGRAPHIC
demographic$age <- sapply(demographic$age_cat,function(x)strsplit(x,'-')[[1]][1])
##!! not sure we need this as a separate object but, for now...
SYNTHETIC_POPULATION <<- left_join(synthetic_pop,demographic[,names(demographic)%in%c('dem_index','age')],by='dem_index')

#####################################################################
## set scenario variables. these can (should) be determined from input data rather than hard coded.
NSCEN <<- 1
SCEN_SHORT_NAME <<- c('base','scen')
SCEN <<- c('Baseline','Scenario 1')

## we effectively have a "SYNTHETIC_POPULATION" per scenario.
pp_summary <- list()
for(scenario in SCEN_SHORT_NAME){
  scenario_name_flag <- sapply(names(synth_pop),function(x)grepl(paste0(scenario,'_'),x))
  scenario_names <- names(synth_pop)[scenario_name_flag]
  # choose subset for each scenario per person summary
  pp_summary[[scenario]] <- synth_pop[,names(synth_pop)%in%c('participant_id','dem_index',scenario_names)]
  ##!! assuming mmetwk is total work and leisure (= nontravel) mmets -- but is it different in different scenarios?
  names(pp_summary[[scenario]])[names(pp_summary[[scenario]])==paste0(scenario,'_mmetwk')] <- 'work_ltpa_marg_met'
  names(pp_summary[[scenario]]) <- sapply(names(pp_summary[[scenario]]),function(x)gsub(paste0(scenario,'_'),'',x))
}

# Generate distance and duration matrices
dist_and_dir <- dist_dur_tbls(pp_summary)
dist <- dist_and_dir$dist
dur <- dist_and_dir$dur

# set as data.table for speed
for(scenario in SCEN_SHORT_NAME) pp_summary[[scenario]] <- setDT(pp_summary[[scenario]])

#####################################################################
## 5 ## 
## ITHIM health calculation

## (1) AP PATHWAY
# Calculate PM2.5 concentrations
system.time(pm_conc <- scenario_pm_calculations(dist,pp_summary))
scenario_pm <- pm_conc$scenario_pm
pm_conc_pp <- pm_conc$pm_conc_pp
# Air pollution DR calculation
system.time(RR_AP_calculations <- ithimr::gen_ap_rr(pm_conc_pp))


#####################################################################
## (2) PA PATHWAY
# Calculate total mMETs
## pp_summary and SYNTHETIC_POPULATION are basically the same thing.
# Only difference is pp_summary is a list for scenarios. This could be more efficient.
# this function differs from ithim-r because mmets differ in baseline and scenario
system.time(mmets_pp <- total_mmet(pp_summary))
# Physical activity calculation
system.time(RR_PA_calculations <- ithimr::gen_pa_rr(mmets_pp))


#####################################################################
## (3) COMBINE (1) AND (2)
# Physical activity and air pollution combined
system.time(RR_PA_AP_calculations <- combined_rr_ap_pa(RR_PA_calculations,RR_AP_calculations))


#####################################################################
## (4) INJURIES
# get data and model
path_to_injury_model_and_data <- 'injury_processing/rds_storage/'
injury_table <- readRDS(paste0(path_to_injury_model_and_data,'processed_injuries_9.Rds'))
baseline_injury_model <- list()
for(i in 1:2){
  baseline_injury_model[[i]] <- list()
  for(j in 1:2){
    baseline_injury_model[[i]][[j]] <- readRDS(paste0(path_to_injury_model_and_data,'city_region',i,j,'.Rds'))
  }
}

# get city data
city_table <- injury_table
for(i in 1:2)
  for(j in 1:2)
    city_table[[i]][[j]] <- injury_table[[i]][[j]][injury_table[[i]][[j]]$region==CITY,]
## for each scenario, add/subtract distance

##!! copied from add_distance_to_injury
# get indices for fast matching data
roads <- unique(injury_table[[1]][[1]]$road)
model_modes <- c('pedestrian','cyclist','motorcycle','car/taxi')
mode_proportions <- matrix(0,nrow=length(roads),ncol=length(model_modes))
colnames(mode_proportions) <- model_modes
rownames(mode_proportions) <- roads
mode_proportions[,colnames(mode_proportions)%in%'pedestrian'] <- c(0,  0.05,0.1,0.1,0.75)
mode_proportions[,colnames(mode_proportions)%in%'cyclist'] <-    c(0,  0.05,0.1,0.1,0.75)
mode_proportions[,colnames(mode_proportions)%in%'motorcycle'] <- c(0.2,0.1, 0.1,0.3,0.3)
mode_proportions[,colnames(mode_proportions)%in%'car/taxi'] <-   c(0.2,0.1, 0.1,0.3,0.3)

# get scenario distance differences
modes <- c('walking','bicycle','motorcycle','car')
mode_labels <- paste0(modes,'_dist')
distance_columns <- colnames(pp_summary[[1]])%in%c(mode_labels)

injury_deaths <- secondary_deaths <- list()
# get prediction for baseline (using smoothed data, not raw data)
for(i in 1:2)
  for(j in 1:2)
    city_table[[i]][[j]]$pred <- predict(baseline_injury_model[[i]][[j]],newdata=city_table[[i]][[j]],type='response')
injury_predictions <- predict_injuries(city_table)
injury_deaths[[1]] <- injury_predictions[[1]] 
secondary_deaths[[1]] <- injury_predictions[[2]] 
# store baseline data
baseline_city_table <- city_table
# for each scenario, add/subtract observed change in travel to/from smoothed baseline data
for(scen in 1:NSCEN+1){
  city_table <- baseline_city_table
  # get difference between scenario and baseline
  scen_diff <- pp_summary[[scen]][,distance_columns,with=F] - pp_summary[[1]][,distance_columns,with=F]
  scen_diff$dem_index <- pp_summary[[1]]$dem_index
  # sum over modes and groups
  scen_diff_dem <- scen_diff[,.(pedestrian=sum(walking_dist),cyclist=sum(bicycle_dist),'car/taxi'=sum(car_dist),motorcycle=sum(motorcycle_dist)),by='dem_index']
  
  # casualty distances
  for(j in 1:2){
    # get indices to match
    road_index <- match(city_table[[1]][[j]]$road,roads)
    mode_index <- match(city_table[[1]][[j]]$cas_mode,model_modes)
    # add/subtract mode-group distances according to road proportions
    ##!! test this
    new_cas_dist <- city_table[[1]][[j]]$cas_distance + as.data.frame(scen_diff_dem)[cbind(city_table[[1]][[j]]$cas_index,mode_index+1)] * mode_proportions[cbind(road_index,mode_index)]
    ##!! some negative distances.
    new_cas_dist[new_cas_dist<0] <- 0
    # edit dataset with new distances
    city_table[[1]][[j]]$cas_distance <- new_cas_dist
  }
  
  # striker distances
  for(i in 1:2){
    # get indices to match
    road_index <- match(city_table[[i]][[1]]$road,roads)
    mode_index <- match(city_table[[i]][[1]]$strike_mode,model_modes)
    # add/subtract mode-group distances according to road proportions
    ##!! test this
    new_str_dist <- city_table[[i]][[1]]$strike_distance + as.data.frame(scen_diff_dem)[cbind(city_table[[i]][[1]]$strike_index,mode_index+1)] * mode_proportions[cbind(road_index,mode_index)]
    ##!! some negative distances.
    new_str_dist[new_str_dist<0] <- 0
    # edit dataset with new distances
    city_table[[i]][[1]]$strike_distance <- new_str_dist
  }
  # get prediction for scenario using modified smoothed data, not raw data
  for(i in 1:2)
    for(j in 1:2)
      city_table[[i]][[j]]$pred <- predict(baseline_injury_model[[i]][[j]],newdata=city_table[[i]][[j]],type='response')
  # summarise predicted fatalities
  injury_predictions <- predict_injuries(city_table)
  # store results
  injury_deaths[[scen]] <- injury_predictions[[1]] 
  secondary_deaths[[scen]] <- injury_predictions[[2]] 
}
# convert to ithimr format
injuries <- cbind(do.call(rbind,injury_deaths),rep(SCEN,each=nrow(injury_deaths[[1]])))
names(injuries) <- c('dem_index','Deaths','scenario')
# compute ylls from deaths
(deaths_yll_injuries <- injury_death_to_yll(injuries))
# store reference number of deaths and ylls
ref_injuries <- deaths_yll_injuries$ref_injuries
##TODO report by mode. covert to burden. then sum.


#####################################################################
## (5) COMBINE (3) AND (4)
# Combine health burden from disease and injury
(hb <- health_burden(RR_PA_AP_calculations,deaths_yll_injuries$deaths_yll_injuries))
pathway_hb <- NULL
constant_mode <- T
if(constant_mode) pathway_hb <- health_burden(RR_PA_AP_calculations,deaths_yll_injuries$deaths_yll_injuries,combined_AP_PA=F)


#####################################################################
## plot
{x11(width=9,height=6); par(mar=c(6,5,1,1),mfrow=c(1,2))
  for(type in c('deaths','ylls')){
    plot_cols <- sapply(names(hb[[type]]),function(x)grepl('scen',x))
    outcomes <- colSums(hb[[type]][,plot_cols])
    names(outcomes) <- sapply(names(outcomes),function(x)last(strsplit(x,'_')[[1]]))
    barplot(outcomes,las=2,cex.axis=1.5,cex.lab=1.5,ylab=paste0('Number of ',type,' averted in Scenario'),xlab='',cex.names=1.5)
  }
}

