{
rm(list=ls())
library(ithimr)
library(splines)
require(tidyverse)
require(knitr)
require(kableExtra)
require(citr)
library(compiler)
library(earth)
library(doParallel)
registerDoParallel(4)
#setwd('~/overflow_dropbox/mh-execute/')
## overwrite some functions for METAHIT's pp_summary use (instead of TIGTHAT's tripset use)
## in general, the overwriting functions are from ithimr's uncertain_travel branch
## in general, ithimr functions are written ithimr::function()
source('metahit_functions.R')
source('mslt_functions.R')
enableJIT(3)

`%myinfix%` <- ifelse(Sys.info()[['sysname']] == "Windows", `%do%`, `%dopar%`)


## 1 SET GLOBAL VARIABLES ##########################################

## set variables, which are TIGTHAT studies' input parameters.

## general settings
setup_call_summary_filename <- 'setup_call_summary.txt'
AGE_RANGE <- c(0,150)
REFERENCE_SCENARIO <- 'Baseline'

## placeholders for uncertain parameters
NSAMPLES <- 1024
MMET_CYCLING <- c(log(4.63),log(1.2)) # 4.63 # 
MMET_WALKING <- c(log(2.53),log(1.1)) # 2.53 # 
PM_CONC_BASE_QUANTILE <- T
PM_TRANS_SHARE_QUANTILE <- T#F 
PA_DOSE_RESPONSE_QUANTILE <- T#F
AP_DOSE_RESPONSE_QUANTILE <- T#F
BACKGROUND_PA_SCALAR <- c(log(1),log(1.1)) # 1 
BACKGROUND_PA_CONFIDENCE <- 1
INJURY_REPORTING_RATE <- c(40,5) # 1
CHRONIC_DISEASE_SCALAR <- c(log(1),log(1.1)) #1
SIN_EXPONENT_SUM <- c(log(1.9),log(1.03)) #2
CASUALTY_EXPONENT_FRACTION <- c(20,20) # 0.5 # 
EMISSION_INVENTORY_CONFIDENCE <- 0.9
DISTANCE_SCALAR_CAR_TAXI <- c(log(1),log(1.1)) # 1
DISTANCE_SCALAR_WALKING <- c(log(1),log(1.1)) # 1
DISTANCE_SCALAR_PT <- c(log(1),log(1.1)) # 1
DISTANCE_SCALAR_CYCLING <- c(log(1),log(1.1)) # 1
DISTANCE_SCALAR_MOTORCYCLE <- c(log(1),log(1.1)) # 1

DIABETES_IHD_RR_F <<- 2.82 ## 2.35
DIABETES_STROKE_RR_F <<- 2.28 ## 1.93
DIABETES_IHD_RR_M <<- 2.16 ## 2.16
DIABETES_STROKE_RR_M <<- 1.83 ## 1.6

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

# REFERENCE_SCENARIO = string: at present, one of 'Baseline' or 'Scenario N' where N is an integer

# NSAMPLES = integer: number of samples to take for each parameter to be sampled

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

# EMISSION_INVENTORY_CONFIDENCE = parameter. double between 0 and 1. 1 = use emission data as they are.
# DISTANCE_SCALAR_CAR_TAXI = double: sets scalar. vector: samples from distribution.
# DISTANCE_SCALAR_WALKING = double: sets scalar. vector: samples from distribution.
# DISTANCE_SCALAR_PT = double: sets scalar. vector: samples from distribution.
# DISTANCE_SCALAR_CYCLING = double: sets scalar. vector: samples from distribution.
# DISTANCE_SCALAR_MOTORCYCLE = double: sets scalar. vector: samples from distribution.


## setting all the global variables at the beginning to minimise ITHIM computation
## copied from ithimr::run_ithim_setup

## SET GLOBAL VALUES
## PROGRAMMING VARIABLES

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
#EMISSION_INVENTORY <<- default_emission_inventory
## 2 GET GLOBAL DATA ##################################################

## copied from ithimr ithim_load_data
global_path <- file.path(find.package('ithimr',lib.loc=.libPaths()), 'extdata/global/')
## for windows??
global_path <- paste0(global_path, "/")

## DATA FILES FOR MODEL  
DISEASE_INVENTORY <<- read.csv(paste0("inputs/dose_response/disease_outcomes_lookup.csv"))
DR_AP <<- read.csv(paste0(global_path,"dose_response/drap/dose_response.csv"))
# root of list_of_files matches DISEASE_INVENTORY$pa_acronym
list_of_files <- list.files(path = paste0(global_path,"dose_response/drpa/extdata/"), recursive = TRUE, pattern = "\\.csv$", full.names = TRUE)
for (i in 1:length(list_of_files)){
  assign(stringr::str_sub(basename(list_of_files[[i]]), end = -5),
         readr::read_csv(list_of_files[[i]],col_types = cols()),
         pos = 1)
}

BACKGROUND_POLLUION_TABLE <<- read.csv('inputs/background-air-pollution/1_apmeans.csv')

disease_short_names <- read.csv("inputs/mslt/disease_names.csv")
DISEASE_SHORT_NAMES <<- disease_short_names


demography <- readxl::read_xlsx('inputs/scenarios/190330_sp_ind_codebook.xlsx',sheet=2,col_names=F)
demogindex_to_numerical <- unlist(demography[,3])
demography[,3] <- 1:nrow(demography)
demo_indices <- unlist(demography[,3])
age_table <- readxl::read_xlsx('inputs/scenarios/190330_sp_ind_codebook.xlsx',sheet=1,col_names=F)
age_category <- unlist(age_table[,1])
age_lower_bounds <- as.numeric(sapply(age_category,function(x)strsplit(x,' to ')[[1]][1]))

## 3 GET MULTI-CITY DATA #################################################
## set scenario variables. these can (should) be determined from input data rather than hard coded.
NSCEN <<- 1
SCEN_SHORT_NAME <<- c('base','scen')
SCEN <<- c('Baseline','Scenario 1')
all_distances <- list()
for(i in 1:length(SCEN)){
  scen_name <- SCEN_SHORT_NAME[i]
  all_distances[[scen_name]] <- list()
  for(file_name in c('emissions_distances','pa_distances'))
    all_distances[[scen_name]][[file_name]] <- readRDS(paste0('inputs/distances/',scen_name,'_',file_name,'.Rds'))
  # if(i==1&&as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE)) > 1e7){
  #   ##!! hack for Rob - laptop can't compute london inh
  #   all_distances[[scen_name]]$inh_distances$london <- readRDS(paste0('inputs/distances/',scen_name,'_london_inh_distances.Rds'))
  #   INCLUDE_LONDON <- T
  #   if(i==1) cat('Including London.\n')
  # }else if(i==1){
     INCLUDE_LONDON <- F
  #   cat('Excluding London.\n')
  # }
}

## get city distances for e.g. bus mode
city_total_distances <- read.csv('inputs/distances/mode_road_city.csv',stringsAsFactors = F)
for(i in 3:ncol(city_total_distances)) city_total_distances[,i] <- as.numeric(city_total_distances[,i])

## injury model / preprocessed data
# get data and model
path_to_injury_model_and_data <- 'inputs/injury/'
injury_table <- readRDS(paste0(path_to_injury_model_and_data,'processed_injuries_9.Rds'))
baseline_injury_model <- list()
for(i in 1:2){
  baseline_injury_model[[i]] <- list()
  for(j in 1:2){
    baseline_injury_model[[i]][[j]] <- readRDS(paste0(path_to_injury_model_and_data,'city_region',i,j,'.Rds'))
    if(INCLUDE_LONDON==F) injury_table[[i]][[j]] <- subset(injury_table[[i]][[j]],region!='london')
  }
}

## use bristol data to define demography etc

filename <- 'inputs/populations/bristol.csv'
demographic <- read_csv(filename,col_types = cols())
demographic$dem_index <- 1:nrow(demographic)

## find min and max age from AGE_RANGE, trips, and demographic.
##!! a lot of this is actually global, but it's coded within cities. It can be brought outside the city loop (to join demography code) by re-writing.
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
demographic$age <- sapply(demographic$age_cat,function(x)strsplit(x,'-')[[1]][1])
DEMOGRAPHIC <<- demographic

# get age-category details from (modified) population data
AGE_CATEGORY <<- unique(POPULATION$age)
AGE_LOWER_BOUNDS <<- as.numeric(sapply(AGE_CATEGORY,function(x)strsplit(x,'-')[[1]][1]))
MAX_AGE <<- max(as.numeric(sapply(AGE_CATEGORY,function(x)strsplit(x,'-')[[1]][2])))

## 4 PREPARE LOCAL (CITY) DATA ##########################################

city_regions_table <- read.csv('inputs/mh_regions_lad_lookup.csv',stringsAsFactors = F)
city_regions <- unique(city_regions_table$cityregion)
city_regions <- city_regions[city_regions!='']
city_las <- city_regions_table$lad11cd[city_regions_table$cityregion%in%city_regions]
la_city_indices <- sapply(city_las,function(x) which(city_regions==city_regions_table$cityregion[city_regions_table$lad11cd==x]))
city_regions_dt <- setDT(city_regions_table[city_regions_table$cityregion%in%city_regions,1:4])
city_regions_dt$la <- 1:nrow(city_regions_dt)
city_regions_dt$city_index <- la_city_indices
if(INCLUDE_LONDON==F) city_regions <- city_regions[city_regions!='london']


## DATA FILES FOR CITY
##!! what are we doing with modes tube, train?
synth_pop_path <- 'inputs/scenarios/'
synth_pop_files <- list.files(synth_pop_path)
synth_pop_files <- synth_pop_files[sapply(synth_pop_files,function(x)grepl('SPind_E[[:digit:]]+.Rds',x))]
la_names <- sapply(synth_pop_files,function(x)gsub('SPind_','',x))
la_names <- sapply(la_names,function(x)gsub('.Rds','',x))
synth_pop_list_in_la_order <- match(la_names,city_regions_dt$lad14cd)
##!! check they're in the right order
print(synth_pop_list_in_la_order)




inventory <- read.csv('inputs/background-air-pollution/emission_inventory.csv')
emission_inventories <- list()
for(city in city_regions){
  row_index <- grepl(city,tolower(BACKGROUND_POLLUION_TABLE$apgroup_name))
  col_indices <- (which(colnames(inventory)=='apgroup_name')+1):ncol(inventory)
  emission_inventories[[city]] <- list()
  for(i in col_indices)
    emission_inventories[[city]][[colnames(inventory)[i]]] <- inventory[row_index,i] 
}

EMISSION_INVENTORIES <<- emission_inventories


## 5 SET PARAMETERS ################################################

parameters <- ithim_setup_parameters(NSAMPLES=NSAMPLES,
                                     MMET_CYCLING=MMET_CYCLING,
                                     MMET_WALKING=MMET_WALKING,
                                     PM_CONC_BASE_QUANTILE=PM_CONC_BASE_QUANTILE,  
                                     PM_TRANS_SHARE_QUANTILE=PM_TRANS_SHARE_QUANTILE,
                                     PA_DOSE_RESPONSE_QUANTILE=PA_DOSE_RESPONSE_QUANTILE,
                                     AP_DOSE_RESPONSE_QUANTILE=AP_DOSE_RESPONSE_QUANTILE,
                                     BACKGROUND_PA_SCALAR=BACKGROUND_PA_SCALAR,
                                     BACKGROUND_PA_CONFIDENCE=BACKGROUND_PA_CONFIDENCE,
                                     INJURY_REPORTING_RATE=INJURY_REPORTING_RATE,
                                     CHRONIC_DISEASE_SCALAR=CHRONIC_DISEASE_SCALAR,
                                     SIN_EXPONENT_SUM=SIN_EXPONENT_SUM,
                                     CASUALTY_EXPONENT_FRACTION=CASUALTY_EXPONENT_FRACTION,
                                     EMISSION_INVENTORY_CONFIDENCE=EMISSION_INVENTORY_CONFIDENCE,
                                     DISTANCE_SCALAR_CAR_TAXI=DISTANCE_SCALAR_CAR_TAXI,
                                     DISTANCE_SCALAR_WALKING=DISTANCE_SCALAR_WALKING,
                                     DISTANCE_SCALAR_PT=DISTANCE_SCALAR_PT,
                                     DISTANCE_SCALAR_CYCLING=DISTANCE_SCALAR_CYCLING,
                                     DISTANCE_SCALAR_MOTORCYCLE=DISTANCE_SCALAR_MOTORCYCLE)


## start metahit
## 6 START LOOP OVER CITIES #################################################

city_results <- list()
}
for(city_ind in 1:length(city_regions)){
  
  ## 7 GET LOCAL (city) DATA ###############################################
  
  #city_ind <- 3
  CITY <<- city_regions[city_ind]
  city_results[[CITY]] <- list()
  
  ## these datasets are all local, saved in local folder.
  ## there will be one folder per city. this block will have to loop over CITY.
  ## OR we have one file with, e.g., all the GBD data in.
  
  # GBD file needs to have the following columns: 
  # age (=label, e.g. 15-49)
  # sex (=male or female)
  # measure
  # cause (GBD_DATA$cause matches DISEASE_INVENTORY$GBD_name)
  # metric
  # burden
  
  ## now process GBD_DATA
  filename <- paste0('inputs/gbd/',CITY,".csv")
  GBD_DATA <- read_csv(filename,col_types = cols())
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
  
  ## get mslt tables
  mslt_df <- read.csv(paste0('inputs/mslt/',CITY, "_mslt.csv"))
  MSLT_DF <<- mslt_df
  
  ## get inh distance
  file_name <- 'inh_distances'
  for(i in 1:length(SCEN)){
    scen_name <- SCEN_SHORT_NAME[i]
    all_distances[[scen_name]][[file_name]] <- readRDS(paste0('inputs/distances/',scen_name,'_',CITY,'_',file_name,'.Rds'))
  }
  
  ## 8 GET/SET CITY SYNTH POP #########################################
  
  # select city LAs
  la_indices <- synth_pop_list_in_la_order[city_regions_dt$city_index==city_ind]
  # remove na (nottinghamshire)
  la_indices <- la_indices[!is.na(la_indices)]
  # set to data table
  synth_pops <- list()
  for(i in 1:length(la_indices)) synth_pops[[i]] <- setDT(readRDS(paste0(synth_pop_path,synth_pop_files[la_indices[i]])))
  # take subset of columns
  for(i in 1:length(synth_pops)) synth_pops[[i]] <- 
    synth_pops[[i]][,sapply(colnames(synth_pops[[i]]),
                            function(x)x%in%c('census_id','demogindex','sport_wkmmets')
    ),with=F]
  # rename
  names(synth_pops) <- la_names[la_indices]
  number_city_las <- length(synth_pops)
  synth_pop <- do.call(rbind,synth_pops)
  synth_pops <- NULL
  
  ## convert synth pop to ithim-r style
  synth_pop$participant_id <- 1:nrow(synth_pop)
  demog_to_dem <- data.table(demogindex=demogindex_to_numerical,dem_index=1:length(demogindex_to_numerical))
  synth_pop <- synth_pop[demog_to_dem,on='demogindex']
  synthetic_pop <- synth_pop[,names(synth_pop)%in%c('participant_id','dem_index'),with=F]
  ##!! not sure we need this as a separate object but, for now...
  SYNTHETIC_POPULATION <<- left_join(synthetic_pop,DEMOGRAPHIC[,names(DEMOGRAPHIC)%in%c('dem_index','age')],by='dem_index')
  synthetic_pop <- NULL
  
  ## we effectively have a "SYNTHETIC_POPULATION" per scenario.
  pp_summary <- list()
  dist_mode_names <- c('walk','cycle','mbikedrive','cardrive','vandrive','subway','bus')
  function_mode_names <- c('walking','bicycle','motorcycle','car','van','subway','bus')
  for(scenario in SCEN_SHORT_NAME){
    #scenario_name_flag <- sapply(names(synth_pop),function(x)grepl(paste0(scenario,'_'),x))
    #scenario_names <- names(synth_pop)[scenario_name_flag]
    # choose subset for each scenario per person summary
    pp_summary[[scenario]] <- synth_pop[,names(synth_pop)%in%c('participant_id','dem_index','census_id','sport_wkmmets'),with=F]
    ## pa
    pp_summary[[scenario]][all_distances[[scenario]]$pa_distances,on='census_id',bicycle_dur_pa:=i.cycle_dur_pa]
    pp_summary[[scenario]][all_distances[[scenario]]$pa_distances,on='census_id',walking_dur_pa:=i.walking_dur_pa]
    ## inh
    for(modenumber in 1:length(dist_mode_names)){
      cols <- sapply(colnames(all_distances[[scenario]]$inh_distances),function(x)grepl(dist_mode_names[modenumber],x))
      pp_summary[[scenario]][,c(paste0(function_mode_names[modenumber],'_dur')):=0]
      pp_summary[[scenario]][match(all_distances[[scenario]]$inh_distances$census_id,pp_summary[[scenario]]$census_id),paste0(function_mode_names[modenumber],'_dur'):=rowSums(all_distances[[scenario]]$inh_distances[,cols,with=F])]
      }
    names(pp_summary[[scenario]])[names(pp_summary[[scenario]])=='sport_wkmmets'] <- 'work_ltpa_marg_met'
  }
  true_pops <- pp_summary[[1]][,.N,by='dem_index']
  POPULATION$population <- true_pops$N[match(POPULATION$dem_index,true_pops$dem_index)]
  synth_pop <- NULL
  INH_NAMES <<- colnames(pp_summary[[1]])%in%paste0(function_mode_names,'_dur')
  PA_NAMES <<- colnames(pp_summary[[1]])%in%c('bicycle_dur_pa','walking_dur_pa')
  
  ##!! use all_distances$...$inh_distances and use all_distances$...$emissions_distances$distance_for_emission. Sum over LA and road type. Join to pp_summary.
  ##!! at present the duration from pa is used, for cycle and walk only.
  ##!! hard coded to maintain naming conventions etc
  DIST <- matrix(0,nrow=3,ncol=NSCEN+1)
  rownames(DIST) <- c('car','motorcycle','bus')
  colnames(DIST) <- SCEN
  n_roads <- ncol(all_distances[[SCEN_SHORT_NAME[1]]]$emissions_distances$distance_for_emission)-2
  for(scen in 1:(NSCEN+1))
    DIST[,scen] <- c(sum(all_distances[[SCEN_SHORT_NAME[scen]]]$emissions_distances$distance_for_emission[mode_name=='cardrive'&la%in%la_names[la_indices],2:(n_roads+1)]),
                     sum(all_distances[[SCEN_SHORT_NAME[scen]]]$emissions_distances$distance_for_emission[mode_name=='mbikedrive'&la%in%la_names[la_indices],2:(n_roads+1)]),
                     ##!! assume total bus travel doesn't change in scenario
                     sum(city_total_distances[city_total_distances[,1]==CITY&city_total_distances[,2]=='bus',3:ncol(city_total_distances)]))
  
  ## 9 ITHIM ########################################
  
  ##!! start loop over parameters here
  
  ## set city-specific parameters
  # background pm2.5
  pm_conc_base <- BACKGROUND_POLLUION_TABLE$apmean_bpm25[grepl(CITY,tolower(BACKGROUND_POLLUION_TABLE$apgroup_name))]
  if(PM_CONC_BASE_QUANTILE==F){
    PM_CONC_BASE <<- pm_conc_base
  }else{
    pm_sd <- BACKGROUND_POLLUION_TABLE$apsd_bpm25[grepl(CITY,tolower(BACKGROUND_POLLUION_TABLE$apgroup_name))]
    lnorm_params <- get_lnorm_params(pm_conc_base,pm_sd)
    parameters$PM_CONC_BASE <- qlnorm(parameters$PM_CONC_BASE_QUANTILE,lnorm_params[1],lnorm_params[2])
  }
  
  # transport portion of pm2.5
  if(PM_TRANS_SHARE_QUANTILE==F){
    pm_transport_share <- BACKGROUND_POLLUION_TABLE$transport_fraction[grepl(CITY,tolower(BACKGROUND_POLLUION_TABLE$apgroup_name))]
    PM_TRANS_SHARE <<- pm_transport_share
  }else{
    pm_share_alpha <- BACKGROUND_POLLUION_TABLE$alpha[grepl(CITY,tolower(BACKGROUND_POLLUION_TABLE$apgroup_name))]
    pm_share_beta <- BACKGROUND_POLLUION_TABLE$beta[grepl(CITY,tolower(BACKGROUND_POLLUION_TABLE$apgroup_name))]
    parameters$PM_TRANS_SHARE <- qbeta(parameters$PM_TRANS_SHARE_QUANTILE,pm_share_alpha,pm_share_beta)
  }
  
  if(EMISSION_INVENTORY_CONFIDENCE<1){
    total <- sum(unlist(EMISSION_INVENTORIES[[CITY]]))
    parameters$EMISSION_INVENTORY <- list()
    for(n in 1:NSAMPLES){
      quantiles <- parameters$EMISSION_INVENTORY_QUANTILE[[n]]
      samples <- lapply(names(quantiles),function(x) qgamma(quantiles[[x]],shape=EMISSION_INVENTORIES[[CITY]][[x]]/total*dirichlet_pointiness(EMISSION_INVENTORY_CONFIDENCE),scale=1))
      names(samples) <- names(quantiles)
      new_total <- sum(unlist(samples))
      parameters$EMISSION_INVENTORY[[n]] <- lapply(samples,function(x)x/new_total)
    }
  }else{
    EMISSION_INVENTORY <<- emission_inventories[[CITY]]
  }
  
  # other parameters to set by city:
  #DISTANCE_SCALAR_CAR_TAXI 
  #DISTANCE_SCALAR_WALKING 
  #DISTANCE_SCALAR_PT
  #DISTANCE_SCALAR_CYCLING 
  #DISTANCE_SCALAR_MOTORCYCLE 
  
  parameters <<- parameters
  DIST <<- DIST
  pp_summary <<- pp_summary
  injury_table <<- injury_table
  baseline_injury_model <<- baseline_injury_model
  
  city_results[[CITY]] <- foreach(sampl = 1:NSAMPLES) %myinfix% {
    for(i in 1:length(parameters))
      assign(names(parameters)[i],parameters[[i]][[sampl]],pos=1)
    CAS_EXPONENT <<- CASUALTY_EXPONENT_FRACTION * SIN_EXPONENT_SUM
    STR_EXPONENT <<- SIN_EXPONENT_SUM - CAS_EXPONENT
    
    ## instead of ithimr::set_vehicle_inventory() # sets vehicle inventory
    vehicle_inventory <- MODE_SPEEDS
    vehicle_inventory$emission_inventory <- 0
    for(m in names(EMISSION_INVENTORY))
      vehicle_inventory$emission_inventory[vehicle_inventory$stage_mode%in%m] <- EMISSION_INVENTORY[[m]]
    VEHICLE_INVENTORY <<- vehicle_inventory
    
    ## (1) AP PATHWAY ######################################
    # Calculate PM2.5 concentrations
    ##!! using pa durations for now, which don't differentiate between road types and las.
    ##!! we don't have durations by road type and la. We could map from distances.
    pm_conc <- scenario_pm_calculations(DIST,pp_summary)
    ## change inh column names
    for(i in 1:length(pp_summary)) colnames(pp_summary[[i]])[INH_NAMES] <- paste0(colnames(pp_summary[[i]])[INH_NAMES],'_inh')
    scenario_pm <- pm_conc$scenario_pm
    pm_conc_pp <- pm_conc$pm_conc_pp
    pm_conc <- NULL
    # Air pollution DR calculation
    RR_AP_calculations <- ithimr::gen_ap_rr(pm_conc_pp)
    pm_conc_pp <- NULL
    
    ## (2) PA PATHWAY ##############################################
    
    # Calculate total mMETs
    ## pp_summary and SYNTHETIC_POPULATION are basically the same thing.
    # Only difference is pp_summary is a list for scenarios. This could be more efficient.
    # this function differs from ithim-r because mmets differ in baseline and scenario
    ##!! check these look sensible
    ## rename pa columns
    for(i in 1:length(pp_summary)) colnames(pp_summary[[i]]) <- sapply(colnames(pp_summary[[i]]),function(x) gsub('_pa','',x))
    mmets_pp <- total_mmet(pp_summary)
    ## change names back
    ##!! alternatively, re-write ITHIM-R functions within metahit_functions.R so that scenario_pm_calculations and total_mmet look for different columns, e.g. _dur_inh and _dur_pa.
    for(i in 1:length(pp_summary)) colnames(pp_summary[[i]])[PA_NAMES] <- paste0(colnames(pp_summary[[i]])[PA_NAMES],'_pa')
    for(i in 1:length(pp_summary)) colnames(pp_summary[[i]]) <- sapply(colnames(pp_summary[[i]]),function(x) gsub('_inh','',x))
    
    # Physical activity calculation
    RR_PA_calculations <- ithimr::gen_pa_rr(mmets_pp)
    mmets_pp <- NULL
    
    ## (3) COMBINE (1) AND (2) #################################################
    
    # Physical activity and air pollution combined
    RR_PA_AP_calculations <- combined_rr_ap_pa(RR_PA_calculations,RR_AP_calculations)
    RR_PA_calculations <- RR_AP_calculations <- NULL
    
    ## (4) INJURIES ##############################################
    
    
    # get city data
    city_table <- injury_table
    for(i in 1:2)
      for(j in 1:2)
        city_table[[i]][[j]] <- injury_table[[i]][[j]][injury_table[[i]][[j]]$region==CITY,]
    ## for each scenario, add/subtract distance
    
    # get indices for fast matching data
    roads <- unique(injury_table[[1]][[1]]$road)
    
    model_modes <- c('pedestrian','cyclist','motorcycle','car/taxi')
    
    injury_deaths <- secondary_deaths <- list()
    # get prediction for baseline (using smoothed data, not raw data)
    for(i in 1:2)
      for(j in 1:2){
        city_table[[i]][[j]]$cas_distance <- city_table[[i]][[j]]$base_cas_distance
        city_table[[i]][[j]]$strike_distance <- city_table[[i]][[j]]$base_strike_distance
        city_table[[i]][[j]]$cas_distance_sum <- city_table[[i]][[j]]$base_cas_distance_sum
        city_table[[i]][[j]]$strike_distance_sum <- city_table[[i]][[j]]$base_strike_distance_sum
        city_table[[i]][[j]]$pred <- predict(baseline_injury_model[[i]][[j]],newdata=city_table[[i]][[j]],type='response')
      }
    injury_predictions <- predict_injuries(city_table)
    injury_deaths[[1]] <- injury_predictions[[1]] 
    secondary_deaths[[1]] <- injury_predictions[[2]] 
    injury_predictions_for_bz_baseline <- predict_injuries_for_bz(city_table)
    # store baseline data
    baseline_city_table <- city_table
    injury_ratios_for_bz <- list()
    injury_ratios_for_bz[[1]] <- injury_predictions_for_bz_baseline
    injury_ratios_for_bz[[1]][,c(1:ncol(injury_ratios_for_bz[[1]]))[-1]] <- injury_ratios_for_bz[[1]][,-1]/injury_predictions_for_bz_baseline[,-1]
    
    ##!!
    for(scen in 1:NSCEN+1){
      scen_name <- SCEN_SHORT_NAME[scen]
      
      city_table <- baseline_city_table
      
      # casualty distances
      for(j in 1:2){
        # edit dataset with new distances
        city_table[[1]][[j]]$cas_distance <- city_table[[1]][[j]][[paste0(scen_name,'_cas_distance')]]
        city_table[[1]][[j]]$cas_distance_sum <- city_table[[1]][[j]][[paste0(scen_name,'_cas_distance_sum')]]
      }
      
      # striker distances
      for(i in 1:2){
        # edit dataset with new distances
        city_table[[i]][[1]]$strike_distance <- city_table[[i]][[1]][[paste0(scen_name,'_strike_distance')]]
        city_table[[i]][[1]]$strike_distance_sum <- city_table[[i]][[1]][[paste0(scen_name,'_strike_distance_sum')]]
      }
      # get prediction for scenario using modified smoothed data, not raw data
      for(i in 1:2)
        for(j in 1:2)
          city_table[[i]][[j]]$pred <- predict(baseline_injury_model[[i]][[j]],newdata=city_table[[i]][[j]],type='response')
      # summarise predicted fatalities
      injury_predictions <- predict_injuries(city_table)
      injury_ratios_for_bz[[scen]] <- predict_injuries_for_bz(city_table)
      # store results
      injury_deaths[[scen]] <- injury_predictions[[1]] 
      secondary_deaths[[scen]] <- injury_predictions[[2]] 
    }
    city_table <- baseline_city_table <- scen_diff <- NULL
    # convert to ithimr format
    injuries <- cbind(do.call(rbind,injury_deaths),rep(SCEN,each=nrow(injury_deaths[[1]])))
    names(injuries) <- c('dem_index','Deaths','scenario')
    # compute ylls from deaths
    (deaths_yll_injuries <- injury_death_to_yll(injuries))
    # store reference number of deaths and ylls
    ref_injuries <- deaths_yll_injuries$ref_injuries
    ##TODO report by mode. covert to burden. then sum.
    
    
    ## (5) COMBINE (3) AND (4)###########################################
    
    # Combine health burden from disease and injury
    (hb <- health_burden(RR_PA_AP_calculations,deaths_yll_injuries$deaths_yll_injuries))
    (pif_table <- health_burden_2(RR_PA_AP_calculations))
    for(scen in 1:NSCEN+1) {
      for(i in 2:ncol(injury_ratios_for_bz[[scen]])) {
        injury_col_name <- colnames(injury_ratios_for_bz[[scen]])[i]
        pif_table[[paste0(SCEN_SHORT_NAME[scen],'_',injury_col_name)]] <- injury_ratios_for_bz[[scen]][[i]]/injury_ratios_for_bz[[1]][[i]]
      }
    }
    
    ## add in population column
    for(i in 1:length(hb))
      hb[[i]] <- left_join(hb[[i]],POPULATION[,c(colnames(POPULATION)%in%c('population','dem_index'))],by='dem_index')
    
    
    pathway_hb <- NULL
    constant_mode <- F
    if(constant_mode) {
      pathway_hb <- health_burden(RR_PA_AP_calculations,deaths_yll_injuries$deaths_yll_injuries,combined_AP_PA=F)
      pathway_pif_table <- health_burden_2(RR_PA_AP_calculations,combined_AP_PA=F)
      x11(); plot(pif_table$scen_pif_pa_ap_noise_no2_ihd,1-(1-pathway_pif_table$scen_pif_pa_ihd)*(1-pathway_pif_table$scen_pif_ap_ihd))
      lines(c(0,1),c(0,1))
    }
    
    RR_PA_AP_calculations <- NULL
    
    #profvis(hb_2 <- belens_function(pif_table) )
    #sort(sapply(ls(),function(x)object.size(get(x))))
    
    ## Rob, added this line to save to my repo, but not sure if you have it too, so I commented it out. 
    # write_csv(hb_2, '../mh-mslt/data/pif.csv')
    
    hb
    
  }
  ## clear memory
  SYNTHETIC_POPULATION <<- NULL
  pp_summary <- NULL
  # reduce size of injury table
  for(i in 1:2)
    for(j in 1:2)
      injury_table[[i]][[j]] <- injury_table[[i]][[j]][injury_table[[i]][[j]]$region!=CITY,]
  
  
  saveRDS(city_results[[CITY]],paste0('outputs/files/',CITY,'_results.Rds'))
  city_results[[CITY]] <- c()
}

for(city_ind in 1:length(city_regions)){
  CITY <<- city_regions[city_ind]
  city_results[[CITY]] <- readRDS(paste0('outputs/files/',CITY,'_results.Rds'))
}
saveRDS(city_results,'outputs/files/city_results.Rds')

## 10 EXTRACT RESULTS AND PLOT ############################################################

outcomes <- list()
plot_cols <- sapply(names(city_results[[1]][[1]][[1]]),function(x)grepl('scen',x))
col_names <- sapply(names(city_results[[1]][[1]][[1]])[plot_cols],function(x)last(strsplit(x,'_')[[1]]))
for(type in c('deaths','ylls')){
  outcomes[[type]] <- list()
  outcomes[[type]]$lower <- matrix(0,nrow=length(city_regions),ncol=length(col_names))
  colnames(outcomes[[type]]$lower) <- col_names
  rownames(outcomes[[type]]$lower) <- city_regions
  outcomes[[type]]$upper <- outcomes[[type]]$median <- outcomes[[type]]$lower
  for(i in 1:length(city_regions)){
    CITY <- city_regions[i]
    sum_results <- sapply(city_results[[CITY]],function(x)colSums(x[[type]][,plot_cols]))
    outcomes[[type]]$median[i,] <- apply(sum_results,1,function(x)median(x))/sum(city_results[[CITY]][[1]][[type]]$population)*1e3
    outcomes[[type]]$lower[i,] <- apply(sum_results,1,quantile,0.05)/sum(city_results[[CITY]][[1]][[type]]$population)*1e3
    outcomes[[type]]$upper[i,] <- apply(sum_results,1,quantile,0.95)/sum(city_results[[CITY]][[1]][[type]]$population)*1e3
  }
}
cols <- rainbow(length(city_regions))
for(type in c('deaths','ylls')){
  pdf(paste0('outputs/figures/',type,'.pdf'),width=9,height=6); 
  par(mar=c(6,5,1,1))
  x<-barplot(outcomes[[type]]$median,las=2,cex.axis=1.5,cex.lab=1.5,ylab=paste0('Thousand ',type,' pp averted in Scenario'),xlab='',cex.names=1.5,beside=T,col=cols)
  legend(fill=cols,bty='n',legend=city_regions,x=prod(dim(outcomes[[type]][[1]])-1),y=max(outcomes[[type]]$median))
  dev.off()
}

for(type in c('deaths','ylls')){
  pdf(paste0('outputs/figures/',type,'.pdf'),width=9,height=6); 
  par(mar=c(6,5,1,1))
  plot(x,outcomes[[type]]$median,las=2,cex.axis=1.5,cex.lab=1.5,ylab=paste0('Thousand ',type,' pp averted in Scenario'),xlab='',xaxt='n',
       cex=1.5,col=cols,pch=15,frame=F,ylim=c(min(outcomes[[type]]$lower),max(outcomes[[type]]$upper)))
  abline(h=0)
  #legend(fill=cols,bty='n',legend=city_regions,x=prod(dim(outcomes[[type]][[1]])-1),y=max(outcomes[[type]]$upper))
  legend(fill=cols,bty='n',legend=city_regions,x=prod(dim(outcomes[[type]][[1]])-2),y=min(outcomes[[type]]$median))
  for(i in 1:nrow(x)) for(j in 1:ncol(x)) 
    lines(c(x[i,j],x[i,j]),c(outcomes[[type]]$lower[i,j],outcomes[[type]]$upper[i,j]),col=cols[i],lwd=2)
  axis(1,at=x[5,],labels=col_names,las=2)
  dev.off()
}


## 11 VOI ############################################################

#if('EMISSION_INVENTORY'%in%names(parameters)){
#  for(i in 1:length(parameters$EMISSION_INVENTORY[[1]])){
#    extract_vals <- sapply(parameters$EMISSION_INVENTORY,function(x)x[[i]])
#    if(sum(extract_vals)!=0)
#      parameters[[paste0('EMISSION_INVENTORY_',names(parameters$EMISSION_INVENTORY[[1]])[i])]] <- extract_vals
#  }
#}

parameter_store <- parameters
for(list_names in c('DR_AP_LIST','PM_CONC_BASE_QUANTILE','PM_TRANS_SHARE_QUANTILE','EMISSION_INVENTORY','EMISSION_INVENTORY_QUANTILES'))
  parameters[[list_names]] <- NULL
parameter_samples <- do.call(cbind,parameters)
saveRDS(parameter_samples,'outputs/files/parameter_samples.Rds')
parameter_samples <- readRDS('outputs/files/parameter_samples.Rds')
#parameter_samples <- parameter_samples[,!colnames(parameter_samples)%in%c('DR_AP_LIST','PM_CONC_BASE_QUANTILE','PM_TRANS_SHARE_QUANTILE','EMISSION_INVENTORY','EMISSION_INVENTORY_QUANTILES')]

plot_cols <- sapply(names(city_results[[1]][[1]][[1]]),function(x)grepl('scen',x)&!(grepl('ac',x)|grepl('neo',x)))
col_names <- sapply(names(city_results[[1]][[1]][[1]])[plot_cols],function(x)last(strsplit(x,'_')[[1]]))

outcome <- list()
type <- 'ylls'
for(i in 1:length(city_regions)){
  CITY <- city_regions[i]
  outcome[[CITY]] <- t(sapply(city_results[[CITY]],function(x)colSums(x[[type]][,plot_cols])))
}

## get basic evppi matrix
numcores <- 4
evppi <- mclapply(1:ncol(parameter_samples), 
                  FUN = compute_evppi,
                  as.data.frame(parameter_samples),
                  outcome, 
                  nscen=NSCEN,
                  all=T,
                  multi_city_outcome=F,
                  mc.cores = ifelse(Sys.info()[['sysname']] == "Windows",  1,  numcores))
evppi <- do.call(rbind,evppi)
colnames(evppi) <- apply(expand.grid(SCEN_SHORT_NAME[2:length(SCEN_SHORT_NAME)],names(outcome)),1,function(x)paste0(x,collapse='_'))
rownames(evppi) <- colnames(parameter_samples)

## replace some rows of evppi if some parameters should be combined
## add four-dimensional EVPPI if AP_DOSE_RESPONSE is uncertain.
numcores <- 1
if("AP_DOSE_RESPONSE_QUANTILE_ALPHA_lri"%in%names(parameters)&&NSAMPLES>=1024){
  AP_names <- sapply(names(parameters),function(x)length(strsplit(x,'AP_DOSE_RESPONSE_QUANTILE_ALPHA')[[1]])>1)
  diseases <- sapply(names(parameters)[AP_names],function(x)strsplit(x,'AP_DOSE_RESPONSE_QUANTILE_ALPHA_')[[1]][2])
  sources <- list()
  for(di in diseases){
    col_names <- sapply(colnames(parameter_samples),function(x)grepl('AP_DOSE_RESPONSE_QUANTILE',x)&grepl(di,x))
    sources[[di]] <- parameter_samples[,col_names]
  }
  evppi_for_AP <- mclapply(1:length(sources), 
                           FUN = compute_evppi,
                           sources,
                           outcome, 
                           all=T,
                           multi_city_outcome=F,
                           mc.cores = ifelse(Sys.info()[['sysname']] == "Windows",  1,  numcores))
  names(evppi_for_AP) <- paste0('AP_DOSE_RESPONSE_QUANTILE_',diseases)
  evppi <- rbind(evppi,do.call(rbind,evppi_for_AP))
  ## get rows to remove
  keep_names <- sapply(rownames(evppi),function(x)!any(c('ALPHA','BETA','GAMMA','TMREL')%in%strsplit(x,'_')[[1]]))
  evppi <- evppi[keep_names,]
}

if("EMISSION_INVENTORY_QUANTILES"%in%names(parameter_store)&&NSAMPLES>=1024){
  sources <- list()
  for(ci in 1:length(city_regions)){
    city <- city_regions[ci]
    sources[[ci]] <- matrix(0,nrow=NSAMPLES,ncol=length(parameter_store$EMISSION_INVENTORY_QUANTILES[[1]]))
    total <- sum(unlist(EMISSION_INVENTORIES[[city]]))
    parameter_store$EMISSION_INVENTORY <- list()
    for(n in 1:NSAMPLES){
      quantiles <- parameter_store$EMISSION_INVENTORY_QUANTILE[[n]]
      samples <- sapply(names(quantiles),function(x) qgamma(quantiles[[x]],shape=EMISSION_INVENTORIES[[city]][[x]]/total*dirichlet_pointiness(EMISSION_INVENTORY_CONFIDENCE),scale=1))
      new_total <- sum(unlist(samples))
      sources[[ci]][n,] <- samples/new_total
    }
  }
  evppi_for_emissions <- mclapply(1:length(sources),
                                  FUN = compute_evppi,
                                  sources,
                                  outcome,
                                  all=F,
                                  multi_city_outcome=F,
                                  mc.cores = ifelse(Sys.info()[['sysname']] == "Windows",  1,  numcores))
  
  #names(evppi_for_emissions) <- paste0('EMISSION_INVENTORY_',city_regions)
  #sapply(evppi_for_emissions,function(x)x[x>0])
  ## get rows to remove
  keep_names <- sapply(rownames(evppi),function(x)!grepl('EMISSION_INVENTORY_',x))
  evppi <- evppi[keep_names,]
  
  evppi <- rbind(evppi,sapply(evppi_for_emissions,function(x)x[x>0]))
  rownames(evppi)[nrow(evppi)] <- 'EMISSION_INVENTORY'
}
print(evppi)

## PA
if(sum(c("BACKGROUND_PA_SCALAR","BACKGROUND_PA_ZEROS")%in%names(parameters))==2&&NSAMPLES>=1024){
  sources <- list()
  for(ci in 1:length(city_regions)){
    city <- city_regions[ci]
    pa_names <- sapply(colnames(parameter_samples),function(x)(grepl('BACKGROUND_PA_SCALAR',x)||grepl('BACKGROUND_PA_ZEROS',x)))
    sources[[ci]] <- parameter_samples[,pa_names]
  }
  evppi_for_pa <- mclapply(1:length(sources), 
                           FUN = compute_evppi,
                           sources, 
                           outcome, 
                           all=F,
                           multi_city_outcome=F,
                           mc.cores = ifelse(Sys.info()[['sysname']] == "Windows",  1,  numcores))
  
  #names(evppi_for_pa) <- paste0('BACKGROUND_PA_',city_regions)
  ## get rows to remove
  keep_names <- sapply(rownames(evppi),function(x)!grepl('BACKGROUND_PA_',x))
  evppi <- evppi[keep_names,]
  evppi <- rbind(evppi,sapply(evppi_for_pa,function(x)x[x>0]))
  #evppi <- rbind(evppi,do.call(rbind,evppi_for_pa))
  rownames(evppi)[nrow(evppi)] <- 'BACKGROUND_PA'
}

## plot evppi
library(RColorBrewer)
library(plotrix)

evppi <- apply(evppi,2,function(x){x[is.na(x)]<-0;x})
{pdf('outputs/figures/evppi.pdf',height=15,width=8);
  par(mar=c(6,20,3.5,5.5))
  labs <- rownames(evppi)
  get.pal=colorRampPalette(brewer.pal(9,"Reds"))
  redCol=rev(get.pal(12))
  bkT <- seq(max(evppi)+1e-10, 0,length=13)
  cex.lab <- 1.5
  maxval <- round(bkT[1],digits=1)
  col.labels<- c(0,maxval/2,maxval)
  cellcolors <- vector()
  for(ii in 1:length(unlist(evppi)))
    cellcolors[ii] <- redCol[tail(which(unlist(evppi[ii])<bkT),n=1)]
  color2D.matplot(evppi,cellcolors=cellcolors,main="",xlab="",ylab="",cex.lab=2,axes=F,border='white')
  fullaxis(side=1,las=2,at=NSCEN*0:(length(outcome)-1)+NSCEN/2,labels=names(outcome),line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=1)
  fullaxis(side=2,las=1,at=(length(labs)-1):0+0.5,labels=labs,line=NA,pos=NA,outer=FALSE,font=NA,lwd=0,cex.axis=0.8)
  mtext(3,text='By how much (%) could we reduce uncertainty in\n the outcome if we knew this parameter perfectly?',line=1)
  color.legend(NSCEN*length(outcome)+0.5,0,NSCEN*length(outcome)+0.8,length(labs),col.labels,rev(redCol),gradient="y",cex=1,align="rb")
  for(i in seq(0,NSCEN*length(outcome),by=NSCEN)) abline(v=i)
  for(i in seq(0,length(labs),by=NSCEN)) abline(h=i)
  dev.off()
  }
