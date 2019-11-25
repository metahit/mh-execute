{
rm(list=ls())
library(ithimr)
library(splines)
require(tidyverse)
require(knitr)
require(kableExtra)
require(citr)
library(compiler)
#setwd('~/overflow_dropbox/mh-execute/')
## overwrite some functions for METAHIT's pp_summary use (instead of TIGTHAT's tripset use)
## in general, the overwriting functions are from ithimr's uncertain_travel branch
## in general, ithimr functions are written ithimr::function()
source('metahit_functions.R')
source('mslt_functions.R')
enableJIT(3)


## 1 SET GLOBAL VARIABLES ##########################################

## set variables, which are TIGTHAT studies' input parameters.

## general settings
setup_call_summary_filename <- 'setup_call_summary.txt'
AGE_RANGE <- c(0,150)
REFERENCE_SCENARIO <- 'Baseline'

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
SIN_EXPONENT_SUM <- 2
CASUALTY_EXPONENT_FRACTION <- 0.5
EMISSION_INVENTORY_CONFIDENCE <- 1
BUS_TO_PASSENGER_RATIO <- 0.22 # estimated from SP vs RTS totals, focussing on urban roads
DISTANCE_SCALAR_CAR_TAXI <- 1
DISTANCE_SCALAR_WALKING <- 1
DISTANCE_SCALAR_PT <- 1
DISTANCE_SCALAR_CYCLING <- 1
DISTANCE_SCALAR_MOTORCYCLE <- 1

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

# ADD_WALK_TO_BUS_TRIPS = logic. T: adds walk trips to all bus trips whose duration exceeds BUS_WALK_TIME. F: no trips added
# ADD_BUS_DRIVERS = logic. T: adds `ghost trips', i.e. trips not taken by any participant. F: no trips added
# ADD_TRUCK_DRIVERS = logic. T: adds `ghost trips', i.e. trips not taken by any participant. F: no trips added

# TEST_WALK_SCENARIO = logic. T: run `scenario 0', one simple scenario where everyone takes one (extra) ten-minute walk trip. F: 5 Accra scenarios.
# TEST_CYCLE_SCENARIO = logic. F: 5 Accra scenarios.
# MAX_MODE_SHARE_SCENARIO = logic. T: run scenarios where we take the maximum mode share across cities and distance categories. F: 5 Accra scenarios.

# REFERENCE_SCENARIO = string: at present, one of 'Baseline' or 'Scenario N' where N is an integer

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

# MOTORCYCLE_TO_CAR_RATIO = parameter. double: sets motorcycle distance relative to car. vector: samples from distribution.
# BUS_TO_PASSENGER_RATIO = parameter. double: sets bus distance relative to bus passenger distance. vector: samples from distribution.
# TRUCK_TO_CAR_RATIO = parameter. double: sets truck distance relative to car. vector: samples from distribution.
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
EMISSION_INVENTORY <<- default_emission_inventory

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
  for(file_name in c('emissions_distances','inh_distances','injury_distances','pa_distances'))
    all_distances[[scen_name]][[file_name]] <- readRDS(paste0('inputs/distances/',scen_name,'_',file_name,'.Rds'))
  if(i==1&&as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE)) > 1e7){
    ##!! hack for Rob - laptop can't compute london inh
    all_distances[[scen_name]]$inh_distances$london <- readRDS(paste0('inputs/distances/',scen_name,'_london_inh_distances.Rds'))
    INCLUDE_LONDON <- T
    if(i==1) cat('Including London.\n')
  }else if(i==1){
    INCLUDE_LONDON <- F
    cat('Excluding London.\n')
  }
}
##!! inh distance is by distance, but should be duration. 
##!! however, we don't have duration by road type.
##!! if we get inh duration, we can delete pa and sum inh to get pa.

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

## start metahit
## 5 START LOOP OVER CITIES #################################################

city_results <- list()
}
for(city_ind in 1:length(city_regions)){
  
  ## 6 GET LOCAL (city) DATA ###############################################
  
  #city_ind <- 3
  CITY <<- city_regions[city_ind]
  
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
  
  PM_CONC_BASE <<- BACKGROUND_POLLUION_TABLE$apmean_bpm25[grepl(CITY,tolower(BACKGROUND_POLLUION_TABLE$apgroup_name))]
  
  mslt_df <- read.csv(paste0('inputs/mslt/',CITY, "_mslt.csv"))
  MSLT_DF <<- mslt_df
  
  
  ## 7 SET PARAMETERS ################################################
  
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
                                               SIN_EXPONENT_SUM=SIN_EXPONENT_SUM,
                                               CASUALTY_EXPONENT_FRACTION=CASUALTY_EXPONENT_FRACTION,
                                               EMISSION_INVENTORY_CONFIDENCE=EMISSION_INVENTORY_CONFIDENCE,
                                               DISTANCE_SCALAR_CAR_TAXI=DISTANCE_SCALAR_CAR_TAXI,
                                               DISTANCE_SCALAR_WALKING=DISTANCE_SCALAR_WALKING,
                                               DISTANCE_SCALAR_PT=DISTANCE_SCALAR_PT,
                                               DISTANCE_SCALAR_CYCLING=DISTANCE_SCALAR_CYCLING,
                                               DISTANCE_SCALAR_MOTORCYCLE=DISTANCE_SCALAR_MOTORCYCLE)
  
  ithimr::set_vehicle_inventory() # sets vehicle inventory
  
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
  ##!! this can be done in mh-distance/process_distances_for_execute.R
  #names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('wkhr','dur',x))
  #names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('wkkm','dist',x))
  #names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('walk_','walking_',x))
  #names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('_cycle','_bicycle',x))
  #names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('cardrive','car',x))
  #names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('mbikedrive','motorcycle',x))
  #names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('mbike','motorcycle',x))
  #names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('tube','subway',x))
  #names(synth_pop) <- sapply(names(synth_pop),function(x)gsub('train','rail',x))
  synth_pop$participant_id <- 1:nrow(synth_pop)
  demog_to_dem <- data.table(demogindex=demogindex_to_numerical,dem_index=1:length(demogindex_to_numerical))
  synth_pop <- synth_pop[demog_to_dem,on='demogindex']
  synthetic_pop <- synth_pop[,names(synth_pop)%in%c('participant_id','dem_index'),with=F]
  ##!! not sure we need this as a separate object but, for now...
  SYNTHETIC_POPULATION <<- left_join(synthetic_pop,DEMOGRAPHIC[,names(DEMOGRAPHIC)%in%c('dem_index','age')],by='dem_index')
  synthetic_pop <- NULL
  
  ## we effectively have a "SYNTHETIC_POPULATION" per scenario.
  pp_summary <- list()
  for(scenario in SCEN_SHORT_NAME){
    #scenario_name_flag <- sapply(names(synth_pop),function(x)grepl(paste0(scenario,'_'),x))
    #scenario_names <- names(synth_pop)[scenario_name_flag]
    # choose subset for each scenario per person summary
    pp_summary[[scenario]] <- synth_pop[,names(synth_pop)%in%c('participant_id','dem_index','census_id','sport_wkmmets'),with=F]
    pp_summary[[scenario]] <- pp_summary[[scenario]][all_distances[[scenario]]$pa_distances[all_distances[[scenario]]$pa_distances$census_id%in%pp_summary[[scenario]]$census_id],on='census_id']
    names(pp_summary[[scenario]])[names(pp_summary[[scenario]])=='sport_wkmmets'] <- 'work_ltpa_marg_met'
    names(pp_summary[[scenario]]) <- sapply(names(pp_summary[[scenario]]),function(x)gsub('wkhr','dur',x))
    names(pp_summary[[scenario]]) <- sapply(names(pp_summary[[scenario]]),function(x)gsub('walk','walking',x))
    names(pp_summary[[scenario]]) <- sapply(names(pp_summary[[scenario]]),function(x)gsub('cycle','bicycle',x))
  }
  synth_pop <- NULL
  
  # Generate distance and duration matrices
  #dist_and_dir <- dist_dur_tbls(pp_summary)
  #dist <- dist_and_dir$dist
  #dur <- dist_and_dir$dur
  
  ##!! use all_distances$...$inh_distances and use all_distances$...$emissions_distances$distance_for_emission. Sum over LA and road type. Join to pp_summary.
  ##!! at present the duration from pa is used, for cycle and walk only.
  ##!! hard coded to maintain naming conventions etc
  dist <- matrix(0,nrow=3,ncol=2)
  rownames(dist) <- c('car','motorcycle','bus')
  colnames(dist) <- c('Baseline','Scenario 1')
  dist[,1] <- c(sum(all_distances$base$emissions_distances$distance_for_emission[mode_name=='cardrive',2:7]),
                sum(all_distances$base$emissions_distances$distance_for_emission[mode_name=='mbikedrive',2:7]),
                BUS_TO_PASSENGER_RATIO*sum(all_distances$base$emissions_distances$distance_for_emission[mode_name=='bus',2:7]))
  dist[,2] <- c(sum(all_distances$scen$emissions_distances$distance_for_emission[mode_name=='cardrive',2:7]),
                sum(all_distances$scen$emissions_distances$distance_for_emission[mode_name=='mbikedrive',2:7]),
                BUS_TO_PASSENGER_RATIO*sum(all_distances$scen$emissions_distances$distance_for_emission[mode_name=='bus',2:7]))
  
  
  ## 9 ITHIM ########################################
  
  ##!! start loop over parameters here
  
  
  ## (1) AP PATHWAY ######################################
  # Calculate PM2.5 concentrations
  ##!! using pa durations for now, which don't differentiate between road types and las.
  ##!! we don't have durations by road type and la. We could map from distances.
  system.time(pm_conc <- scenario_pm_calculations(dist,pp_summary))
  SYNTHETIC_POPULATION <<- NULL
  scenario_pm <- pm_conc$scenario_pm
  pm_conc_pp <- pm_conc$pm_conc_pp
  pm_conc <- NULL
  # Air pollution DR calculation
  system.time(RR_AP_calculations <- ithimr::gen_ap_rr(pm_conc_pp))
  pm_conc_pp <- NULL
  
  ## (2) PA PATHWAY ##############################################
  
  # Calculate total mMETs
  ## pp_summary and SYNTHETIC_POPULATION are basically the same thing.
  # Only difference is pp_summary is a list for scenarios. This could be more efficient.
  # this function differs from ithim-r because mmets differ in baseline and scenario
  ##!! check these look sensible
  system.time(mmets_pp <- total_mmet(pp_summary))
  # Physical activity calculation
  system.time(RR_PA_calculations <- ithimr::gen_pa_rr(mmets_pp))
  mmets_pp <- NULL
  
  ## (3) COMBINE (1) AND (2) #################################################
  
  # Physical activity and air pollution combined
  system.time(RR_PA_AP_calculations <- combined_rr_ap_pa(RR_PA_calculations,RR_AP_calculations))
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
  # reduce size of injury table
  for(i in 1:2)
    for(j in 1:2)
      injury_table[[i]][[j]] <- injury_table[[i]][[j]][injury_table[[i]][[j]]$region!=CITY,]
  
  model_modes <- c('pedestrian','cyclist','motorcycle','car/taxi')
  
  injury_deaths <- secondary_deaths <- list()
  # get prediction for baseline (using smoothed data, not raw data)
  for(i in 1:2)
    for(j in 1:2){
      city_table[[i]][[j]]$cas_distance <- city_table[[i]][[j]]$base_cas_distance
      city_table[[i]][[j]]$strike_distance <- city_table[[i]][[j]]$base_strike_distance
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
  scenarios <- c('base_','scen_')
  for(scen in 1:NSCEN+1){
    scen_name <- scenarios[scen]
    
    city_table <- baseline_city_table
    
    # casualty distances
    for(j in 1:2){
      # edit dataset with new distances
      city_table[[1]][[j]]$cas_distance <- city_table[[1]][[j]][[paste0(scen_name,'cas_distance')]]
    }
    
    # striker distances
    for(i in 1:2){
      # edit dataset with new distances
      city_table[[i]][[1]]$strike_distance <- city_table[[i]][[1]][[paste0(scen_name,'strike_distance')]]
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
  city_table <- baseline_city_table <- scen_diff <- pp_summary <- NULL
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
  for(scen in 1:NSCEN+1) 
    for(i in 2:ncol(injury_ratios_for_bz[[scen]])) {
      injury_col_name <- colnames(injury_ratios_for_bz[[scen]])[i]
      pif_table[[paste0(SCEN_SHORT_NAME[scen],'_',injury_col_name)]] <- injury_ratios_for_bz[[scen]][[i]]/injury_ratios_for_bz[[1]][[i]]
    }
  
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
  
  city_results[[CITY]] <- hb
  
}

saveRDS(city_results,'outputs/files/city_results.Rds')

## plot ############################################################

outcomes <- list()
plot_cols <- sapply(names(city_results[[1]][[1]]),function(x)grepl('scen',x))
col_names <- sapply(names(city_results[[1]][[1]])[plot_cols],function(x)last(strsplit(x,'_')[[1]]))
for(type in c('deaths','ylls')){
  outcomes[[type]] <- matrix(0,nrow=length(city_regions),ncol=length(col_names))
  colnames(outcomes[[type]]) <- col_names
  rownames(outcomes[[type]]) <- city_regions
  for(i in 1:length(city_regions)){
    CITY <- city_regions[i]
    outcomes[[type]][i,] <- colSums(city_results[[CITY]][[type]][,plot_cols])
  }
}
cols <- rainbow(length(city_regions))
for(type in c('deaths','ylls')){
  pdf(paste0('outputs/figures/',type,'.pdf'),width=9,height=6); par(mar=c(6,5,1,1))
  barplot(outcomes[[type]],las=2,cex.axis=1.5,cex.lab=1.5,ylab=paste0('Number of ',type,' averted in Scenario'),xlab='',cex.names=1.5,beside=T,col=cols)
  legend(fill=cols,bty='n',legend=city_regions,x=prod(dim(outcomes[[type]])-1),y=max(outcomes[[type]]))
  dev.off()
}
