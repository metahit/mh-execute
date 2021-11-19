get_lnorm_params <- function(mean_val,sd_val){
  mu <- -log(((sd_val/mean_val)^2+1)/(mean_val^2))/2
  sig2 <- 2*(log(mean_val)-mu)
  c(mu,sqrt(sig2))
}

ithim_setup_parameters <- function(NSAMPLES = 1,
         MMET_CYCLING = 4.63,
         MMET_WALKING = 2.53,
         PM_CONC_BASE_QUANTILE = F,  
         PM_TRANS_SHARE_QUANTILE = F,
         PA_DOSE_RESPONSE_QUANTILE = F,
         AP_DOSE_RESPONSE_QUANTILE = F,
         BACKGROUND_PA_SCALAR = 1,
         BACKGROUND_PA_CONFIDENCE = 1,
         INJURY_REPORTING_RATE = 1,
         CHRONIC_DISEASE_SCALAR = 1,
         SIN_EXPONENT_SUM = 2,
         CASUALTY_EXPONENT_FRACTION = 0.5,
         PM_EMISSION_INVENTORY_CONFIDENCE = 1,
         DISTANCE_SCALAR_CAR_TAXI = 1,
         DISTANCE_SCALAR_WALKING = 1,
         DISTANCE_SCALAR_PT = 1,
         DISTANCE_SCALAR_CYCLING = 1,
         DISTANCE_SCALAR_MOTORCYCLE = 1){
  
  ## PARAMETERS
  ##RJ parameters are assigned to the environment and so are set for every function. They are over-written when sample_parameters is called.
  MMET_CYCLING <<- MMET_CYCLING
  MMET_WALKING <<- MMET_WALKING
  PM_CONC_BASE_QUANTILE <<- PM_CONC_BASE_QUANTILE
  PM_TRANS_SHARE_QUANTILE <<- PM_TRANS_SHARE_QUANTILE
  PA_DOSE_RESPONSE_QUANTILE <<- PA_DOSE_RESPONSE_QUANTILE
  BACKGROUND_PA_SCALAR <<- BACKGROUND_PA_SCALAR
  BACKGROUND_PA_CONFIDENCE <<- BACKGROUND_PA_CONFIDENCE
  INJURY_REPORTING_RATE <<- INJURY_REPORTING_RATE
  CHRONIC_DISEASE_SCALAR <<- CHRONIC_DISEASE_SCALAR
  SIN_EXPONENT_SUM <<- SIN_EXPONENT_SUM
  CASUALTY_EXPONENT_FRACTION <<- CASUALTY_EXPONENT_FRACTION
  DISTANCE_SCALAR_CAR_TAXI <<- DISTANCE_SCALAR_CAR_TAXI
  DISTANCE_SCALAR_WALKING <<- DISTANCE_SCALAR_WALKING
  DISTANCE_SCALAR_PT <<- DISTANCE_SCALAR_PT
  DISTANCE_SCALAR_CYCLING <<-  DISTANCE_SCALAR_CYCLING
  DISTANCE_SCALAR_MOTORCYCLE <<- DISTANCE_SCALAR_MOTORCYCLE
  parameters <- list()
  
  ##Variables with normal distribution
  normVariables <- c("MMET_CYCLING",
                     "MMET_WALKING",
                     "BACKGROUND_PA_SCALAR",
                     "CHRONIC_DISEASE_SCALAR",
                     "SIN_EXPONENT_SUM",
                     "DISTANCE_SCALAR_CAR_TAXI",
                     "DISTANCE_SCALAR_WALKING",
                     "DISTANCE_SCALAR_PT",
                     "DISTANCE_SCALAR_CYCLING",
                     "DISTANCE_SCALAR_MOTORCYCLE")
  for (i in 1:length(normVariables)) {
    name <- normVariables[i]
    val <- get(normVariables[i])
    if (length(val) == 1) {
      assign(name, val, envir = .GlobalEnv)
    } else {
      parameters[[name]] <-
        rlnorm(NSAMPLES, val[1], val[2])
    }
  }
  
  ##Variables with beta distribution
  betaVariables <- c("INJURY_REPORTING_RATE",
                     "CASUALTY_EXPONENT_FRACTION")
  for (i in 1:length(betaVariables)) {
    name <- betaVariables[i]
    val <- get(betaVariables[i])
    if (length(val) == 1) {
      assign(name, val, envir = .GlobalEnv)
    } else {
      parameters[[name]] <-
        rbeta(NSAMPLES, val[1], val[2])
    }
  }
  
  ##Variables with uniform distribution
  unifVariables <- c("PM_CONC_BASE_QUANTILE",
                     "PM_TRANS_SHARE_QUANTILE")
  for (i in 1:length(unifVariables)) {
    name <- unifVariables[i]
    val <- get(unifVariables[i])
    if (val == F) {
      assign(name, val, envir = .GlobalEnv)
    } else {
      parameters[[name]] <- runif(NSAMPLES,0,1)
    }
  }
  
  if(BACKGROUND_PA_CONFIDENCE<1){
    parameters$BACKGROUND_PA_ZEROS <- runif(NSAMPLES,0,1)
  }
  
  if(PM_EMISSION_INVENTORY_CONFIDENCE<1){
    parameters$PM_EMISSION_INVENTORY_QUANTILES <- list()
    for(n in 1:NSAMPLES){
      parameters$PM_EMISSION_INVENTORY_QUANTILES[[n]] <- lapply(PM_EMISSION_INVENTORIES[[1]],function(x) runif(1))
    }
  }
  
  ## PA DOSE RESPONSE
  if(PA_DOSE_RESPONSE_QUANTILE == T ) {
    pa_diseases <- subset(DISEASE_INVENTORY,physical_activity==1)
    dr_pa_list <- list()
    for(disease in pa_diseases$pa_acronym)
      parameters[[paste0('PA_DOSE_RESPONSE_QUANTILE_',disease)]] <- runif(NSAMPLES,0,1)
  }
  
  #### AP DOSE RESPONSE
  AP_DOSE_RESPONSE_QUANTILE <<- AP_DOSE_RESPONSE_QUANTILE
  ## shortcut: use saved median values
  if(!AP_DOSE_RESPONSE_QUANTILE){
    global_path <- file.path(find.package('ithimr',lib.loc=.libPaths()), 'extdata/global/')
    global_path <- paste0(global_path, "/")
    DR_AP_LIST <<- readRDS(paste0(global_path,"dose_response/drap/dr_ap_list.Rds"))
  }else{
    dr_ap_list <- list()
    ap_diseases <- subset(DISEASE_INVENTORY,air_pollution==1)
    ap_parameters <- list()
    for(disease in ap_diseases$ap_acronym){ 
      for(letter in c('ALPHA_','BETA_','GAMMA_','TMREL_')){
        if(AP_DOSE_RESPONSE_QUANTILE){
          ap_parameters[[paste0('AP_DOSE_RESPONSE_QUANTILE_',letter,disease)]] <- runif(NSAMPLES,0,1)
          parameters[[paste0('AP_DOSE_RESPONSE_QUANTILE_',letter,disease)]] <- ap_parameters[[paste0('AP_DOSE_RESPONSE_QUANTILE_',letter,disease)]]
        } else {
          ap_parameters[[paste0('AP_DOSE_RESPONSE_QUANTILE_',letter,disease)]] <- 0.5
        }
      }
      dr_ap <- subset(DR_AP,cause_code==disease)
      dr_ap_list[[disease]] <- list()
      quant1 <- ap_parameters[[paste0('AP_DOSE_RESPONSE_QUANTILE_GAMMA_',disease)]]
      quant2 <- ap_parameters[[paste0('AP_DOSE_RESPONSE_QUANTILE_BETA_',disease)]]
      quant3 <- ap_parameters[[paste0('AP_DOSE_RESPONSE_QUANTILE_ALPHA_',disease)]]
      quant4 <- ap_parameters[[paste0('AP_DOSE_RESPONSE_QUANTILE_TMREL_',disease)]]
      for(age in unique(dr_ap$age_code)){
        dr_ap_age <- subset(dr_ap,age_code==age)
        #######################################
        lbeta <- log(dr_ap_age$beta)
        lgamma <- log(dr_ap_age$gamma)
        gamma_val <- quantile(density(lgamma),quant1)
        beta_val <- c()
        for(i in 1:ifelse(AP_DOSE_RESPONSE_QUANTILE,NSAMPLES,1)){
          den <- kde2d(lgamma,lbeta,n=c(1,100),h=0.2,lims=c(gamma_val[i],gamma_val[i],min(lbeta)-1,max(lbeta)+1))
          beta_val[i] <- approx(x=cumsum(den$z)/sum(den$z),y=den$y,xout=quant2[i])$y
        }
        mod <- gam(log(alpha)~te(log(gamma),log(beta)),data=dr_ap_age)
        pred_val <- predict(mod, newdata=data.frame(beta=exp(beta_val),gamma=exp(gamma_val)),se.fit=T)
        alpha_val <- qnorm(quant3,pred_val$fit,sqrt(mod$sig2))
        # generate a value for tmrel given alpha, beta and gamma
        mod <- gam(log(tmrel)~ns(log(gamma),df=8)+ns(log(beta),df=8)+ns(log(alpha),df=8),data=dr_ap_age)
        pred_val <- predict(mod, newdata=data.frame(alpha=exp(alpha_val),beta=exp(beta_val),gamma=exp(gamma_val)),se.fit=T)
        tmrel_val <- qnorm(quant4,pred_val$fit,sqrt(mod$sig2))
        dr_ap_list[[disease]][[as.character(age)]] <- data.frame(alpha=exp(alpha_val),beta=exp(beta_val),gamma=exp(gamma_val),tmrel=exp(tmrel_val))
      }
      if(AP_DOSE_RESPONSE_QUANTILE){ 
        # turn list inside out, so it's indexed first by sample
        parameters$DR_AP_LIST <- lapply(1:NSAMPLES,function(x)lapply(dr_ap_list,function(y) lapply(y,function(z)z[x,])))
      }else{
        DR_AP_LIST <<- dr_ap_list
      }
    }
    
  }
  parameters
}

#' @export
scenario_pm_calculations <- function(dist,pp_summary){
  
  # concentration contributed by non-transport share (remains constant across the scenarios)
  non_transport_pm_conc <- PM_CONC_BASE*(1 - PM_TRANS_SHARE)  
  
  ## adding in travel not covered in the synthetic trip set, based on distances travelled relative to car, set in VEHICLE_INVENTORY
  emission_dist <- dist
  
  ## get emission factor by dividing inventory by baseline distance. (We don't need to scale to a whole year, as we are just scaling the background concentration.)
  ordered_efs <- VEHICLE_INVENTORY$pm_emission_inventory[match(rownames(emission_dist),VEHICLE_INVENTORY$stage_mode)]/emission_dist[,'Baseline']
  ## get new emission by multiplying emission factor by scenario distance.
  trans_emissions <- emission_dist*t(repmat(ordered_efs,NSCEN+1,1))
  ## augment with travel emission contributions that aren't included in distance calculation
  for(mode_type in which(!VEHICLE_INVENTORY$stage_mode%in%rownames(emission_dist))){
    em <- VEHICLE_INVENTORY$pm_emission_inventory[mode_type]
    if(em>0){
      trans_emissions <- rbind(trans_emissions,rep(em,ncol(trans_emissions)))
      rownames(trans_emissions)[nrow(trans_emissions)] <- VEHICLE_INVENTORY$stage_mode[mode_type]
    }
  }
  
  ## scenario travel pm2.5 calculated as relative to the baseline
  ##!! as we divide scenario by baseline, we do not need to multiply through by distance scalars
  baseline_sum <- sum(trans_emissions[,SCEN[1]])
  conc_pm <- c()
  ## in this sum, the non-transport pm is constant; the transport emissions scale the transport contribution (PM_TRANS_SHARE) to the base level (PM_CONC_BASE)
  for(i in 1:length(SCEN_SHORT_NAME))
    conc_pm[i] <- non_transport_pm_conc + PM_TRANS_SHARE*PM_CONC_BASE*sum(trans_emissions[,SCEN[i]])/baseline_sum
  
  ##RJ rewriting ventilation as a function of MMET_CYCLING and MMET_WALKING, loosely following de Sa's SP model.
  vent_rates <- data.frame(stage_mode=VEHICLE_INVENTORY$stage_mode,stringsAsFactors = F) 
  vent_rates$vent_rate <- BASE_LEVEL_INHALATION_RATE  # L / min
  vent_rates$vent_rate[vent_rates$stage_mode=='cycle'] <- BASE_LEVEL_INHALATION_RATE + MMET_CYCLING
  vent_rates$vent_rate[vent_rates$stage_mode%in%c('pedestrian','walk_to_bus')] <- BASE_LEVEL_INHALATION_RATE + MMET_WALKING
  
  ##RJ rewriting exposure ratio as function of ambient PM2.5, as in Goel et al 2015
  ##!! five fixed parameters: BASE_LEVEL_INHALATION_RATE (10), CLOSED_WINDOW_PM_RATIO (0.5), CLOSED_WINDOW_RATIO (0.5), ROAD_RATIO_MAX (3.216), ROAD_RATIO_SLOPE (0.379)
  ##RJ question for RG: should this function account for PM_TRANS_SHARE?
  on_road_off_road_ratio <- ROAD_RATIO_MAX - ROAD_RATIO_SLOPE*log(conc_pm)
  ##RJ question for RG: why is 'in car' twice better than 'away from road'?
  # averaging over windows open and windows closed
  in_vehicle_ratio <- (1-CLOSED_WINDOW_RATIO)*on_road_off_road_ratio + CLOSED_WINDOW_RATIO*CLOSED_WINDOW_PM_RATIO 
  # subway ratio is a constant
  subway_ratio <- rep(SUBWAY_PM_RATIO,length(conc_pm))
  # open vehicles experience the ``on_road_off_road_ratio'', and closed vehicles experience the ``in_vehicle_ratio''
  ratio_by_mode <- rbind(on_road_off_road_ratio,in_vehicle_ratio,subway_ratio)
  # assign rates according to the order of the ratio_by_mode array: 1 is open vehicle, 2 is closed vehicle, 3 is subway
  open_vehicles <- c('pedestrian','walk_to_bus','cycle','motorcycle','auto_rickshaw','shared_auto','cycle_rickshaw')
  rail_vehicles <- c('subway','rail') 
  vent_rates$vehicle_ratio_index <- sapply(vent_rates$stage_mode,function(x) ifelse(x%in%rail_vehicles,3,ifelse(x%in%open_vehicles,1,2)))
  
  pp_summary2 <- pp_summary#lapply(pp_summary,function(y)y[,sapply(colnames(y),function(x)!grepl('_dist',x)),with=F])
  for(i in 1:length(pp_summary2)) colnames(pp_summary2[[i]]) <- sapply(colnames(pp_summary2[[i]]),function(x)gsub('_dur','',x))
  ## multiply through by distance scalars
  for(i in 1:length(pp_summary2)){
    pp_summary2[[i]][,pedestrian := pedestrian * DISTANCE_SCALAR_WALKING]
    pp_summary2[[i]][,cycle := cycle * DISTANCE_SCALAR_CYCLING]
    pp_summary2[[i]][,motorcycle := motorcycle * DISTANCE_SCALAR_MOTORCYCLE]
    pp_summary2[[i]][,car := car * DISTANCE_SCALAR_CAR_TAXI]
    pp_summary2[[i]][,bus := bus * DISTANCE_SCALAR_PT]
    pp_summary2[[i]][,subway := subway * DISTANCE_SCALAR_PT]
  }
  travel_indices <- which(colnames(pp_summary2[[1]])%in%vent_rates$stage_mode)
  travel_modes <- colnames(pp_summary2[[1]])[travel_indices]
  vent_modes <- match(travel_modes,vent_rates$stage_mode)
  
  # prepare individual-level dataset
  pm_conc_pp <- SYNTHETIC_POPULATION
  vent_multiplier <- repmat(vent_rates$vent_rate[vent_modes],nrow(pm_conc_pp),1)
  vent_and_ratio_multiplier <- vent_multiplier*repmat(ratio_by_mode[vent_rates$vehicle_ratio_index[vent_modes],1],nrow(pm_conc_pp),1)
  # compute individual-level pm scenario by scenario
  for (i in 1:length(SCEN)){
    
    scen_travel <- pp_summary2[[i]]
    # duration is per week
    scen_travel[, on_road_dur := Reduce(`+`, .SD), .SDcols=travel_indices]
    #vent_travel <- scen_travel[,travel_indices,with=F] * vent_and_ratio_multiplier
    scen_travel[, on_road_pm := Reduce(`+`, lapply(seq_along(.SD),function(x)(.SD[[x]]*vent_and_ratio_multiplier[,x]))), .SDcols=names(scen_travel)[travel_indices]]
    #vent_travel[, on_road_pm := Reduce(`+`, .SD), .SDcols=names(vent_travel)]
    
    ## PM2.5 inhalation = total mg inhaled / total volume inhaled
    # calculate non-travel air inhalation
    non_transport_air_inhaled <- (24*7-scen_travel$on_road_dur)*BASE_LEVEL_INHALATION_RATE
    # concentration of pm inhaled = total pm inhaled / total air inhaled
    pm_conc <- ((non_transport_air_inhaled * as.numeric(conc_pm[i])) + scen_travel$on_road_pm)#/(non_transport_air_inhaled+individual_data$air_inhaled)
    # match individual ids to set per person pm exposure
    pm_conc_pp[[paste0('pm_conc_',SCEN_SHORT_NAME[i])]] <- pm_conc/24/7 #* conc_pm[i]
  }
  
  #####PM normalise
  ## Rahul made changes here/./-- no normalisation  
  ## calculating means of individual-level concentrations
  #mean_conc <- mean(pm_conc_pp[[paste0("pm_conc_", SCEN_SHORT_NAME[1])]])
  
  #normalise <- as.numeric(conc_pm[1])/as.numeric(mean_conc)
  #for (i in 1: length(SCEN_SHORT_NAME))
  #pm_conc_pp[[paste0("pm_conc_", SCEN_SHORT_NAME[i])]] <- normalise*pm_conc_pp[[paste0("pm_conc_", SCEN_SHORT_NAME[i])]]
  
  pm_conc_pp$participant_id <- as.integer(pm_conc_pp$participant_id)
  
  list(scenario_pm=conc_pm, pm_conc_pp=pm_conc_pp)
  
}

#' @export
total_mmet <- function(pp_summary){
  
  
  ##!! maybe we don't need individual distance and can remove it from pp_summary?
  pp_summary2 <- pp_summary#lapply(pp_summary,function(y)y[,sapply(colnames(y),function(x)!grepl('_dist',x)),with=F])
  for(i in 1:length(pp_summary2)) colnames(pp_summary2[[i]]) <- sapply(colnames(pp_summary2[[i]]),function(x)gsub('_dur','',x))
  # Get total individual level walking and cycling and sport mmets 
  synth_pop_return <- pp_summary2[[1]]
  for (i in 1:length(SCEN)){
    synth_pop_temp <- pp_summary2[[i]]
    synth_pop_return[[paste0(SCEN_SHORT_NAME[i],'_mmet')]] <- ifelse(any(names(synth_pop_temp) == 'work_ltpa_marg_met'), synth_pop_temp$work_ltpa_marg_met * BACKGROUND_PA_SCALAR, 0)
    
    scen_travel <- subset(pp_summary2[[i]],participant_id%in%synth_pop_return$participant_id)
    ##!! check units: duration is in hours per week, and mmets multiply hours?
    scen_travel$cycling_mmet <- scen_travel$cycle * MMET_CYCLING * DISTANCE_SCALAR_CYCLING
    if('walk_to_bus'%in%names(scen_travel)) scen_travel$pedestrian <- scen_travel$pedestrian+scen_travel$walk_to_bus
    scen_travel$walking_mmet <- scen_travel$pedestrian * MMET_WALKING * DISTANCE_SCALAR_WALKING
    
    individual_data <- scen_travel
    
    part_id <- match(individual_data$participant_id,synth_pop_return$participant_id)
    synth_pop_return[[paste0(SCEN_SHORT_NAME[i],'_mmet')]][part_id] <- 
      synth_pop_return[[paste0(SCEN_SHORT_NAME[i],'_mmet')]][part_id] + individual_data$cycling_mmet + individual_data$walking_mmet
  }
  
  mmets <- dplyr::select(synth_pop_return, any_of(c('participant_id', 'sex', 'age', 'dem_index', paste0(SCEN_SHORT_NAME,'_mmet'))))
  mmets
  
}

#' @export
combined_rr_ap_pa <- function(ind_pa,ind_ap){
  
  # Replace NaNs with 1
  ind_ap[is.na(ind_ap)] <- 1
  
  # Replace Na with 1
  ind_pa[is.na(ind_pa)] <- 1
  
  # join pa and ap datasets
  ind_ap_pa <- left_join(ind_pa, ind_ap, by = c('participant_id','dem_index'))
  
  ### iterating over all all disease outcomes
  for ( j in c(1:nrow(DISEASE_INVENTORY))[DISEASE_INVENTORY$physical_activity == 1 & DISEASE_INVENTORY$air_pollution == 1]){
    ac <- as.character(DISEASE_INVENTORY$acronym[j])
    for (scen in SCEN_SHORT_NAME){
      ind_ap_pa[[paste('RR_pa_ap', scen, ac, sep = '_')]] <- ind_ap_pa[[paste('RR_pa', scen, ac, sep = '_')]] * ind_ap_pa[[paste('RR_ap', scen, ac, sep = '_')]]
    }
  }
  
  ind_ap_pa
}

#' @export
predict_without_model <- function(x,newdata,type='response'){
  x <- newdata$base_pred*
    (newdata$cas_distance_sum/newdata$base_cas_distance_sum)^(CAS_EXPONENT-1)*
    (newdata$strike_distance_sum/newdata$base_strike_distance_sum)^(STR_EXPONENT-1)*
    newdata$cas_distance/newdata$base_cas_distance*
    newdata$strike_distance/newdata$base_strike_distance
  x[is.na(x)] <- 0
  x
}

#' @export
summarise_injuries <- function(city_table){
  fatal_data <- list()
  for(i in 1:2){
    fatal_data[[i]] <- list()
    for(j in 1:2){
      fatal_data[[i]][[j]] <- city_table[[i]][[j]][city_table[[i]][[j]]$cas_severity=='Fatal',]
    }
  }
  #cas_modes <- unique(fatal_data[[1]][[1]]$cas_mode)
  cas_fatal1 <- setDT(fatal_data[[1]][[1]])[,.(Deaths=sum(pred)),by=c('cas_index')]
  cas_fatal2 <- setDT(fatal_data[[1]][[2]])[,.(Deaths=sum(pred)),by=c('cas_index')]
  cas_fatal <- cas_fatal1
  ##!! check this line
  cas_fatal$Deaths <- cas_fatal$Deaths + cas_fatal2$Deaths[match(cas_fatal1$cas_index,cas_fatal2$cas_index)]
  
  nonspecific_fatalities <- sum(fatal_data[[2]][[2]]$pred) + sum(fatal_data[[2]][[1]]$pred)
  
  return(list(cas_fatal,nonspecific_fatalities))
  
}

#' @export
summarise_injuries_for_bz <- function(city_table){
  fatal_data <- list()
  for(i in 1:2){
    fatal_data[[i]] <- list()
    for(j in 1:2){
      fatal_data[[i]][[j]] <- city_table[[i]][[j]][city_table[[i]][[j]]$cas_severity%in%c('Serious','Fatal'),]
    }
  }
  #cas_modes <- unique(fatal_data[[1]][[1]]$cas_mode)
  cas_fatal1 <- setDT(fatal_data[[1]][[1]])[,.(Deaths=sum(pred)),by=c('cas_index','cas_mode','cas_severity')]
  cas_fatal2 <- setDT(fatal_data[[1]][[2]])[,.(Deaths=sum(pred)),by=c('cas_index','cas_mode','cas_severity')]
  cas_fatal <- cas_fatal1
  ##!! check this line
  cas_fatal[cas_fatal2,Deaths2:=i.Deaths,on=c('cas_index','cas_mode','cas_severity')]
  cas_fatal[,burden:=Deaths+Deaths2]
  cas_fatal$burden[is.na(cas_fatal$burden)] <- 0
  injury_by_mode_and_demo <- dcast(cas_fatal,cas_index ~ cas_mode+cas_severity,value.var = 'burden')
  #for (i in names(injury_by_mode_and_demo))
  #  injury_by_mode_and_demo[is.na(get(i)), (i):=0]
  for (j in seq_len(ncol(injury_by_mode_and_demo)))
    set(injury_by_mode_and_demo,which(is.na(injury_by_mode_and_demo[[j]])),j,0)
  #print(injury_by_mode_and_demo)
  
  nonspecific_fatalities <- sum(fatal_data[[2]][[2]]$pred) + sum(fatal_data[[2]][[1]]$pred)
  
  return(injury_by_mode_and_demo)
}

#' @export
injury_death_to_yll <- function(injuries){
  joined_injury <- left_join(injuries, GBD_INJ_YLL[,c('dem_index','yll_dth_ratio')], by="dem_index")
  
  joined_injury$YLL <- joined_injury$Deaths*joined_injury$yll_dth_ratio
  death_and_yll <- dplyr::select(joined_injury, c('dem_index','scenario','Deaths','YLL'))
  
  x_deaths <- dplyr::select(death_and_yll, -YLL)
  x_deaths <- spread(x_deaths,scenario, Deaths) %>% as.data.frame()
  x_yll <- dplyr::select(death_and_yll, -Deaths)
  x_yll <- spread(x_yll,scenario, YLL) %>% as.data.frame()
  
  ref_scen <- REFERENCE_SCENARIO
  ref_scen_index <- which(SCEN==ref_scen)
  calc_scen <- SCEN[SCEN!=ref_scen]
  calc_scen_index <- which(colnames(x_deaths)%in%calc_scen)
  
  ref_injuries <- as.data.frame(cbind(dem_index=x_deaths$dem_index,deaths=x_deaths[[ref_scen]],ylls=x_yll[[ref_scen]]))
  deaths <- t(repmat(unlist(ref_injuries$deaths),NSCEN,1)) - x_deaths[,calc_scen_index]
  ylls <- t(repmat(unlist(ref_injuries$ylls),NSCEN,1)) - x_yll[,calc_scen_index]
  deaths_yll_injuries <- as.data.frame(cbind(dem_index=x_deaths$dem_index,deaths, ylls))
  
  metric <- c("deaths", "yll")
  k <- 1
  for  (i in 1: 2)
    for (j in c(1:(NSCEN+1))[-ref_scen_index]){
      names(deaths_yll_injuries)[1+k] <- paste0(SCEN_SHORT_NAME[j],"_",metric[i],"_inj")
      k<-k+1
    }
  
  list(deaths_yll_injuries=deaths_yll_injuries,ref_injuries=ref_injuries)
}

#' @export
health_burden <- function(ind_ap_pa,inj,combined_AP_PA=T){
  
  # subset gbd data for outcome types
  gbd_data_scaled <- DISEASE_BURDEN
  #gbd_data_scaled$burden[gbd_data_scaled$cause%in%c("Neoplasms","Ischemic heart disease","Tracheal, bronchus, and lung cancer","Breast cancer","Colon and rectum cancer","Uterine cancer")] <- 
  #  gbd_data_scaled$burden[gbd_data_scaled$cause%in%c("Neoplasms","Ischemic heart disease","Tracheal, bronchus, and lung cancer","Breast cancer","Colon and rectum cancer","Uterine cancer")]*CHRONIC_DISEASE_SCALAR
  ## chronic disease scalar scales all diseases
  gbd_data_scaled$burden <- gbd_data_scaled$burden*CHRONIC_DISEASE_SCALAR
  gbd_deaths <- subset(gbd_data_scaled,measure=='Deaths')
  gbd_ylls <- subset(gbd_data_scaled,measure=='YLLs (Years of Life Lost)')
  pop_details <- DEMOGRAPHIC
  deaths <- ylls <- pop_details
  # set up reference (scen1)
  reference_scenario <- SCEN_SHORT_NAME[which(SCEN==REFERENCE_SCENARIO)]
  scen_names <- SCEN_SHORT_NAME[SCEN_SHORT_NAME!=reference_scenario]
  ### iterating over all all disease outcomes
  for ( j in 1:nrow(DISEASE_INVENTORY)){
    # Disease acronym and full name
    ac <- as.character(DISEASE_INVENTORY$acronym[j])
    gbd_dn <- as.character(DISEASE_INVENTORY$GBD_name[j])
    # calculating health outcome, or independent pathways?
    pathways_to_calculate <- ifelse(combined_AP_PA,1,DISEASE_INVENTORY$physical_activity[j]+DISEASE_INVENTORY$air_pollution[j])
    for(path in 1:pathways_to_calculate){
      # set up column names
      if(combined_AP_PA){
        middle_bit <-
          paste0(
            ifelse(DISEASE_INVENTORY$physical_activity[j] == 1, 'pa_', ''),
            ifelse(DISEASE_INVENTORY$air_pollution[j] == 1, 'ap_', '')
          )
      }else{
        # if independent, choose which one
        middle_bit <- c('pa_','ap_')[which(c(DISEASE_INVENTORY$physical_activity[j],DISEASE_INVENTORY$air_pollution[j])==1)[path]]
      }
      base_var <- paste0('RR_', middle_bit, reference_scenario, '_', ac)
      scen_vars <- paste0('RR_', middle_bit, scen_names, '_', ac)
      # subset gbd data
      gbd_deaths_disease <- subset(gbd_deaths,cause==gbd_dn)
      gbd_ylls_disease <- subset(gbd_ylls,cause==gbd_dn)
      # set up pif tables
      pif_table <- setDT(ind_ap_pa[,colnames(ind_ap_pa)%in%c(base_var,'dem_index')])
      setnames(pif_table,base_var,'outcome')
      pif_ref <- pif_table[,.(sum(outcome)),by='dem_index']
      ## sort pif_ref
      setorder(pif_ref,dem_index)
      for (index in 1:length(scen_vars)){
        # set up naming conventions
        scen <- scen_names[index]
        scen_var <- scen_vars[index]
        yll_name <- paste0(scen, '_ylls_',middle_bit,ac)
        deaths_name <- paste0(scen, '_deaths_',middle_bit,ac)
        # Calculate PIFs for selected scenario
        pif_table <- setDT(ind_ap_pa[,colnames(ind_ap_pa)%in%c(scen_var,'dem_index')])
        setnames(pif_table,scen_var,'outcome')
        pif_temp <- pif_table[,.(sum(outcome)),by='dem_index']
        ## sort pif_temp
        setorder(pif_temp,dem_index)
        pif_scen <- (pif_ref[,2] - pif_temp[,2]) / pif_ref[,2]
        # Calculate ylls 
        yll_dfs <- combine_health_and_pif(pif_values=pif_scen, hc = gbd_ylls_disease)
        ylls[[yll_name]] <- yll_dfs[,V1]
        # Calculate deaths 
        death_dfs <- combine_health_and_pif(pif_values=pif_scen,hc=gbd_deaths_disease)
        deaths[[deaths_name]] <- death_dfs[,V1]
      }
    }
  }
  # Select deaths columns
  inj_deaths <- dplyr::select(inj, c(dem_index, contains("deaths")))
  # Select yll columns
  inj_ylls <- dplyr::select(inj, c(dem_index, contains("yll")))
  # Join injuries data to global datasets
  deaths <- left_join(deaths, inj_deaths, by = c("dem_index"))
  ylls <- left_join(ylls, inj_ylls, by = c("dem_index"))
  list(deaths=deaths,ylls=ylls)
}

#' @export
health_burden_2 <- function(ind_ap_pa,combined_AP_PA=T){
  pop_details <- DEMOGRAPHIC
  pif_scen <- pop_details
  # set up reference (scen1)
  reference_scenario <- SCEN_SHORT_NAME[which(SCEN==REFERENCE_SCENARIO)]
  scen_names <- SCEN_SHORT_NAME[SCEN_SHORT_NAME!=reference_scenario]
  ### iterating over all all disease outcomes
  for ( j in 1:nrow(DISEASE_INVENTORY)){
    # Disease acronym and full name
    ac <- as.character(DISEASE_INVENTORY$acronym[j])
    gbd_dn <- as.character(DISEASE_INVENTORY$GBD_name[j])
    # calculating health outcome, or independent pathways?
    pathways_to_calculate <- ifelse(combined_AP_PA,1,DISEASE_INVENTORY$physical_activity[j]+DISEASE_INVENTORY$air_pollution[j])
    for(path in 1:pathways_to_calculate){
      # set up column names
      if(combined_AP_PA){
        middle_bit <-
          paste0(
            ifelse(DISEASE_INVENTORY$physical_activity[j] == 1, 'pa_', ''),
            ifelse(DISEASE_INVENTORY$air_pollution[j] == 1, 'ap_', '')
          )
        middle_bit_plus <-
          paste0(
            ifelse(DISEASE_INVENTORY$physical_activity[j] == 1, 'pa_', ''),
            ifelse(DISEASE_INVENTORY$air_pollution[j] == 1, 'ap_', ''),
            ifelse(DISEASE_INVENTORY$noise[j] == 1, 'noise_', ''),
            ifelse(DISEASE_INVENTORY$nitrogen_dioxide[j] == 1, 'no2_', '')
          )
      }else{
        # if independent, choose which one
        middle_bit <- middle_bit_plus <- c('pa_','ap_')[which(c(DISEASE_INVENTORY$physical_activity[j],DISEASE_INVENTORY$air_pollution[j])==1)[path]]
      }
      base_var <- paste0('RR_', middle_bit, reference_scenario, '_', ac)
      scen_vars <- paste0('RR_', middle_bit, scen_names, '_', ac)
      # set up pif tables
      pif_table <- setDT(ind_ap_pa[,colnames(ind_ap_pa)%in%c(base_var,'dem_index')])
      setnames(pif_table,base_var,'outcome')
      pif_ref <- pif_table[,.(sum(outcome)),by='dem_index']
      ## sort pif_ref
      setorder(pif_ref,dem_index)
      for (index in 1:length(scen_vars)){
        # set up naming conventions
        scen <- scen_names[index]
        scen_var <- scen_vars[index]
        pif_name <- paste0(scen, '_pif_',middle_bit_plus,ac)
        # Calculate PIFs for selected scenario
        pif_table <- setDT(ind_ap_pa[,colnames(ind_ap_pa)%in%c(scen_var,'dem_index')])
        setnames(pif_table,scen_var,'outcome')
        pif_temp <- pif_table[,.(sum(outcome)),by='dem_index']
        ## sort pif_temp
        setorder(pif_temp,dem_index)
        pif_scen[[pif_name]] <- (pif_ref[,V1] - pif_temp[,V1]) / pif_ref[,V1]
      }
    }
  }
  return(pif_scen)
}

#' @export
# ---- RunMSLT --- ## used to be called belens_function
#' @export
#' 
RunMSLT <- function(mslt_df, i_sex, i_age_cohort, disease_names, pif) {
  
  source("~/mh-execute/mslt_functions.R")
  
  # mslt_df=read_csv("~/mh-execute/inputs/mslt/bristol_mslt.csv")
  # disease_names=readRDS("~/mh-mslt/output/parameters/DISEASE_SHORT_NAMES.rds")
  # i_sex=c("male", "female")
  # i_age_cohort=seq(from=17, to=97, by =5)
  # pif=read_csv("~/mh-mslt/input/pif_place_holder.csv") ## PLACE HOLDER (from Aus)
  
  # 
  # ### Relative risks diabetes (BZ (19/11/21) move to parameters
  DIABETES_IHD_RR_F <- 2.82 ## c(2.82, CI (2.35, 3.38) get SD from CI
  DIABETES_STROKE_RR_F <- 2.28 ## c(2.28) CI (1.93, 2.69) get SD from CI
  DIABETES_IHD_RR_M <- 2.16 ## c(2.16, CI (1.82, 2.56) get SD from CI
  DIABETES_STROKE_RR_M <- 1.83 ## c(1.83) CI (1.60, 2.08) get SD from CI
  
  DISEASE_SHORT_NAMES <- disease_names
  
  
  ### USE Australian pifs for diseases as not ready for city regions
  ### Get pif, generated in mh-execute and saved in input mh-mslt
  pif_expanded <- pif %>%
    mutate(pif_cyclist_deaths=1.1,
           pif_pedestrian_deaths=1.1,
           pif_cyclist_ylds=1.1,
           pif_pedestrian_ylds=1.1,
           pif_motor_deaths=1.1,
           pif_motorcyclist_deaths=1.1,
           pif_motor_ylds=1.1,
           pif_motorcyclist_ylds=1.1,
           pif_road_deaths=1.1,
           pif_road_ylds=1.1,
           pif_lri_ylds=0.03,
           pif_lri_deaths=0.03,
           pif_copd=0.03) %>%
    dplyr::slice(rep(1:dplyr::n(), each = 5)) %>% ### expand to 1-yr group (repeat values)
    mutate(age=rep(seq(16,100,1), times = 2))
  
  
  # dataframe of the age and sex cohorts (crossing just does a cross product) for loop below
  age_sex_cohorts <- crossing(data.frame(age=i_age_cohort),
                              data.frame(sex=c('male', 'female'))) %>%
    dplyr::mutate(cohort=paste0(age,"_",sex))
  
  ### Run baseline general dataframes
  general_life_table_list_bl <- list()
  
  
  for (i in 1:nrow(age_sex_cohorts)){
    suppressWarnings(
      general_life_table_list_bl[[i]] <- RunLifeTable(
        in_idata    = mslt_df,
        in_sex      = age_sex_cohorts$sex[i],
        in_mid_age  = age_sex_cohorts$age[i],
        death_rates = NA ## mortality trends data if available, not for now. We would need trends for each city region.
      ))
    names(general_life_table_list_bl)[i] <- age_sex_cohorts$cohort[i]
  }
  
  # convert the list of dataframes to single dataframes
  general_life_table_bl <- bind_rows(general_life_table_list_bl, .id = "age_group") %>%
    mutate(age_group = as.numeric(gsub("_.*","",age_group)))
  
  ### Run disease life tables
  
  disease_cohorts <- DISEASE_SHORT_NAMES %>%
    # Exclude non-diseases, road injuries, and diseases with no pif
    dplyr::filter(is_not_dis == 0 & acronym != 'no_pif' & acronym != 'other' ) %>%
    dplyr::select(sname,acronym,males,females)
  
  # adding the age and sex cohorts:
  age_sex_disease_cohorts <- crossing(age_sex_cohorts,disease_cohorts) %>%
    mutate(cohort=paste0(age,'_',sex,'_',sname)) %>%
    # Exclude non-male diseases (and non-female if there were any)
    filter( (sex=='male' & males==1) | (sex=='female' & females==1)) %>%
    dplyr::select(age,sex,sname,acronym,cohort) %>%
    # ishd and strk have the prerequisite disease dmt2
    mutate(prerequsite=ifelse(sname %in% c("ishd","strk"),paste0(age,"_",sex,"_dmt2"),0)) %>%
    # ensuring prequisites are calculated first
    arrange(age,sex,prerequsite,sname)
  
  
  disease_life_table_list_bl <- list()
  
  for (i in 1:nrow(age_sex_disease_cohorts)){
    disease_life_table_list_bl[[i]] <- RunDisease(
      in_idata         = mslt_df,
      in_mid_age       = age_sex_disease_cohorts$age[i],
      in_sex           = age_sex_disease_cohorts$sex[i],
      in_disease       = age_sex_disease_cohorts$sname[i],
      incidence_trends = NA,
      mortality_trends = NA
    )
    names(disease_life_table_list_bl)[i] <- age_sex_disease_cohorts$cohort[i]
  }
  
  ### Run non diseases (injuries and lower respiratory disease)
  
  non_disease_cohorts <- DISEASE_SHORT_NAMES %>%
    # Exclude non-diseases, road injuries, and diseases with no pif
    dplyr::filter(is_not_dis == 1 & acronym != 'no_pif' & acronym != 'other' ) %>%
    dplyr::select(sname,acronym,males,females)
  
  # adding the age and sex cohorts:
  age_sex_non_disease_cohorts <- crossing(age_sex_cohorts,non_disease_cohorts) %>%
    mutate(cohort=paste0(age,'_',sex,'_',sname)) %>%
    dplyr::select(age,sex,sname,acronym,cohort)
  
  non_disease_life_table_list_bl <- list()
  
  for (i in 1:nrow(age_sex_non_disease_cohorts )){ 
    
    non_disease_life_table_list_bl[[i]] <- RunNonDisease(
      in_idata         = mslt_df,
      in_sex           = age_sex_non_disease_cohorts$sex[i],
      in_mid_age       = age_sex_non_disease_cohorts$age[i],
      in_non_disease   = age_sex_non_disease_cohorts$sname[i]
    )
    names(non_disease_life_table_list_bl)[i] <- age_sex_non_disease_cohorts$cohort[i]
  }  
  
  ### Run scenario disease life tables
  
  disease_relative_risks <- tribble(
    ~sex    , ~prerequsite, ~disease , ~relative_risk       ,
    "male"  ,  "dmt2"     ,  "ishd"  ,  DIABETES_IHD_RR_M   ,
    "female",  "dmt2"     ,  "ishd"  ,  DIABETES_IHD_RR_F   ,
    "male"  ,  "dmt2"     ,  "strk"  ,  DIABETES_STROKE_RR_M,
    "female",  "dmt2"     ,  "strk"  ,  DIABETES_STROKE_RR_F
  )
  
  disease_life_table_list_sc <- list()
  
  for (i in 1:nrow(age_sex_disease_cohorts)){
    # i=39
    td1_age_sex <- mslt_df %>% ### new mslt dataframe with modified incidence rates
      filter(age >= age_sex_disease_cohorts$age[i] & sex == age_sex_disease_cohorts$sex[i])
    
    pif_colname <- paste0('pif_',age_sex_disease_cohorts$acronym[i])
    
    pif_disease <- pif_expanded %>%
      filter(age >= age_sex_disease_cohorts$age[i] & sex == age_sex_disease_cohorts$sex[i]) %>%
      dplyr::select(age,sex,pif_colname)
    
    # adjustment for diabetes effect on ihd and stroke
    if(age_sex_disease_cohorts$prerequsite[i] != 0){
      # get name for pif column
      target_disease <- paste0("pif_",age_sex_disease_cohorts$acronym[i])
      # get prerequisite disease cohort name (i.e., age_sex_dmt2 for diabetes)
      dia_col <- age_sex_disease_cohorts$prerequsite[i]
      # select relative risk of disease given diabetes (depends on sex, not age)
      relative_risk <- disease_relative_risks %>%
        filter(sex == age_sex_disease_cohorts$sex[i] &
                 disease == age_sex_disease_cohorts$sname[i]) %>%
        pull(relative_risk)
      # (store old pif)
      # old_pif <- pif_disease[[target_disease]]
      # diabetes pif = - { scenario prevalence - baseline prevalence } * (RR - 1)  / { baseline prevalence * (RR - 1) + 1 }
      scenario_prevalence <- disease_life_table_list_sc[[dia_col]]$px
      baseline_prevalence <- disease_life_table_list_bl[[dia_col]]$px
      pif_dia <- -(scenario_prevalence - baseline_prevalence)*(relative_risk-1)/
        (baseline_prevalence * (relative_risk-1) + 1)
      # modify pif for target disease: new pif =  (1 - old pif) * (1 - diabetes pif)
      pif_disease[[target_disease]] <- 1- (1-pif_disease[[target_disease]]) * (1-pif_dia)
      # print(sum(old_pif-pif_disease[[target_disease]]))
    }
    
    incidence_colname <- paste0('incidence_', age_sex_disease_cohorts$sname[i])
    new_col <- td1_age_sex%>%pull(incidence_colname) * (1 - (pif_disease%>%pull(pif_colname)))
    new_col[is.na(new_col)] <- 0
    td1_age_sex[[incidence_colname]] <- new_col
    
    ## Instead of idata, feed td to run scenarios. Now all diseases are run again, with the effect of diabetes
    ## on cardiovascular diseases taken into account. 
    
    disease_life_table_list_sc[[i]] <- RunDisease(
      in_idata         = td1_age_sex,
      in_sex           = age_sex_disease_cohorts$sex[i],
      in_mid_age       = age_sex_disease_cohorts$age[i],
      in_disease       = age_sex_disease_cohorts$sname[i],
      incidence_trends = NA,
      mortality_trends = NA
    )
    names(disease_life_table_list_sc)[i] <- age_sex_disease_cohorts$cohort[i]
  }
  
  
  
  for (cohort in age_sex_disease_cohorts$cohort) {
    disease_life_table_list_sc[[cohort]]$diff_inc_disease <-
      disease_life_table_list_sc[[cohort]]$incidence_disease - disease_life_table_list_bl[[cohort]]$incidence_disease
    
    disease_life_table_list_sc[[cohort]]$diff_prev_disease <-
      disease_life_table_list_sc[[cohort]]$px - disease_life_table_list_bl[[cohort]]$px
    
    disease_life_table_list_sc[[cohort]]$diff_mort_disease <-
      disease_life_table_list_sc[[cohort]]$mx - disease_life_table_list_bl[[cohort]]$mx
    
    disease_life_table_list_sc[[cohort]]$diff_pylds_disease <-
      (disease_life_table_list_sc[[cohort]]$px - disease_life_table_list_bl[[cohort]]$px) * 
      (disease_life_table_list_bl[[cohort]]$dw_disease)
  }
  
  
  # convert the list of dataframes to single dataframes
  disease_life_table_bl <- bind_rows(disease_life_table_list_bl, .id = "age_sex_disease_cohort") %>%
    mutate(age_sex_disease_cohort = as.numeric(gsub("_.*","",age_sex_disease_cohort))) %>%
    dplyr::rename(age_group=age_sex_disease_cohort,
                  cause=disease)
  
  disease_life_table_sc <- bind_rows(disease_life_table_list_sc, .id = "age_sex_disease_cohort") %>%
    mutate(age_sex_disease_cohort = as.numeric(gsub("_.*","",age_sex_disease_cohort))) %>%
    dplyr::rename(age_group=age_sex_disease_cohort,
                  cause=disease)
  
  
  ### Run scenario non diseases
  
  non_disease_life_table_list_sc <- list()
  
  for (i in 1:nrow(age_sex_non_disease_cohorts)){
    # i=6
    td1_age_sex <- mslt_df %>% ### new mslt dataframe with modified injuries and lri deaths and ylds rates
      filter(age >= age_sex_non_disease_cohorts$age[i] & sex == age_sex_non_disease_cohorts$sex[i])
    
    pif_colname_deaths <- paste0('pif_',age_sex_non_disease_cohorts$acronym[i], '_deaths')
    pif_colname_ylds <- paste0('pif_',age_sex_non_disease_cohorts$acronym[i], '_ylds')
    
    pif_non_disease <- pif_expanded %>%
      filter(age >= age_sex_non_disease_cohorts$age[i] & sex == age_sex_non_disease_cohorts$sex[i]) %>%
      dplyr::select(age,sex,pif_colname_deaths,pif_colname_ylds)
    
    
    death_colname <- paste0('deaths_rate_', age_sex_non_disease_cohorts$sname[i])
    ylds_colname <- paste0('ylds_rate_', age_sex_non_disease_cohorts$sname[i])
    
    
    new_deaths <- td1_age_sex%>%pull(death_colname) * (1 - (pif_non_disease%>%pull(pif_colname_deaths)))
    new_deaths[is.na(new_deaths)] <- 0
    
    new_ylds <- td1_age_sex%>%pull(ylds_colname) * (1 - (pif_non_disease%>%pull(pif_colname_ylds)))
    new_ylds[is.na(new_ylds)] <- 0
    
    td1_age_sex[[death_colname]] <- new_deaths
    td1_age_sex[[ylds_colname]] <- new_ylds
    
    ## Instead of idata, feed td to run scenarios. Now all non_diseases are run again
    
    
    non_disease_life_table_list_sc[[i]] <- RunNonDisease(
      in_idata         = td1_age_sex,
      in_sex           = age_sex_non_disease_cohorts$sex[i],
      in_mid_age       = age_sex_non_disease_cohorts$age[i],
      in_non_disease   = age_sex_non_disease_cohorts$sname[i]
    )
    names(non_disease_life_table_list_sc)[i] <- age_sex_non_disease_cohorts$cohort[i]
    
  }
  
  ### Difference rates
  
  for (cohort in age_sex_non_disease_cohorts$cohort) {
    non_disease_life_table_list_sc[[cohort]]$diff_mort <-
      non_disease_life_table_list_sc[[cohort]]$deaths_rate -  non_disease_life_table_list_bl[[cohort]]$deaths_rate
    
    non_disease_life_table_list_sc[[cohort]]$diff_pylds <-
      non_disease_life_table_list_sc[[cohort]]$ylds_rate -  non_disease_life_table_list_bl[[cohort]]$ylds_rate }
  
  
  
  # convert the list of dataframes to single dataframes
  non_disease_life_table_bl <- bind_rows(non_disease_life_table_list_bl, .id = "age_sex_non_disease_cohort") %>%
    mutate(age_sex_non_disease_cohort = as.numeric(gsub("_.*","",age_sex_non_disease_cohort))) %>%
    dplyr::rename(age_group=age_sex_non_disease_cohort,
                  cause=non_disease,
                  mx=deaths_rate)
  
  non_disease_life_table_sc <- bind_rows(non_disease_life_table_list_sc, .id = "age_sex_non_disease_cohort") %>%
    mutate(age_sex_non_disease_cohort = as.numeric(gsub("_.*","",age_sex_non_disease_cohort))) %>%
    dplyr::rename(age_group=age_sex_non_disease_cohort,
                  cause=non_disease,
                  mx=deaths_rate)
  
  ### Add up mortality and pylds changes
  
  ### Diseases: Sum mortality rate and pylds change scenarios
  mx_pylds_sc_total_disease_df <- disease_life_table_sc %>%
    group_by(age_group,sex,age) %>%
    dplyr::summarise(mortality_sum=sum(diff_mort_disease,na.rm=T),
                     pylds_sum=sum(diff_pylds_disease,na.rm=T)) %>%
    ungroup() %>%
    mutate(age_sex_cohort=paste0(age_group,'_',sex))
  
  ### Non-diseases
  mx_pylds_sc_total_non_disease_df <- non_disease_life_table_sc %>%
    group_by(age_group,sex,age) %>%
    dplyr::summarise(mortality_sum=sum(diff_mort,na.rm=T),
                     pylds_sum=sum(diff_pylds,na.rm=T)) %>%
    ungroup() %>%
    mutate(age_sex_cohort=paste0(age_group,'_',sex))
  
  
  ### Run scenario life tables
  
  general_life_table_list_sc <- list()
  
  for (i in 1:nrow(age_sex_cohorts)){
    # modify idata's mortality and pyld total for the said scenario
    mx_pylds_sc_total_disease_df_cohort <- mx_pylds_sc_total_disease_df %>%
      filter(age_sex_cohort==age_sex_cohorts$cohort[i]) %>%
      dplyr::select(age,mortality_sum,pylds_sum)
    
    mx_pylds_sc_total_non_disease_df_cohort <- mx_pylds_sc_total_non_disease_df %>%
      filter(age_sex_cohort==age_sex_cohorts$cohort[i]) %>%
      dplyr::select(age,mortality_sum,pylds_sum)
    
    ### Modify rates in static MSLT  (pylds are always static, mx can include future trends)
    #### With diseases changes in mortality and ylds
    td2 <- mslt_df %>%
      filter(sex==age_sex_cohorts$sex[i]) %>%
      left_join(mx_pylds_sc_total_disease_df_cohort,by="age") %>%
      mutate(mx=mx+replace_na(mortality_sum,0),
             pyld_rate=pyld_rate+replace_na(pylds_sum,0)) %>%
      dplyr::select(-mortality_sum,-pylds_sum)
    
    #### With diseases changes in mortality and ylds
    td3 <-  td2 %>%
      filter(sex==age_sex_cohorts$sex[i]) %>%
      left_join(mx_pylds_sc_total_non_disease_df_cohort,by="age") %>%
      mutate(mx=mx+replace_na(mortality_sum,0),
             pyld_rate=pyld_rate+replace_na(pylds_sum,0)) %>%
      dplyr::select(-mortality_sum,-pylds_sum)
    
    # ### Modify death rates with future trends NOT AVAILABLE FOR METAHIT
    # td3 <- death_projections %>%
    #   mutate(cohort=paste(age_cohort, sex, sep = "_")) %>% # variable to match change in mortality rates df
    #   filter(cohort==age_sex_cohorts$cohort[i]) %>%
    #   left_join(mx_pylds_sc_total_disease_df_cohort) %>%
    #   mutate(rate=rate+replace_na(mortality_sum,0))%>%
    #   dplyr::select(-mortality_sum,-pylds_sum)   
    
    
    
    suppressWarnings(
      general_life_table_list_sc[[i]] <- RunLifeTable(
        in_idata    = td3,
        in_sex      = age_sex_cohorts$sex[i],
        in_mid_age  = age_sex_cohorts$age[i],
        death_rates = NA
      ))
    names(general_life_table_list_sc)[i] <- age_sex_cohorts$cohort[i]
  }
  
  # convert the list of dataframes to single dataframes
  general_life_table_sc <- bind_rows(general_life_table_list_sc, .id = "age_group") %>%
    mutate(age_group = as.numeric(gsub("_.*","",age_group)))
  
  
  ## In the following list 'output_life_table', 34 data frames are nested per age and sex cohort
  ## Outputs are generated following the index order of disease life tables baseline and scenarios where ##diabetes is     firstcalculated as it impacts on cardiovascular diseases. 
  
  ### Generate outputs dataframe
  #### Diseases life tables
  dia_index <- which(DISEASE_SHORT_NAMES$sname=='dmt2')
  dia_order <- c(dia_index,c(1:nrow(DISEASE_SHORT_NAMES))[-dia_index])
  
  ### Combine diseases and general life tables for scenarios
  ### Step needed to calculate numbers (rates*people cohort)
  scenario_d <- inner_join(disease_life_table_sc %>%
                             dplyr::select(age_group,sex,age,cause,incidence_disease,mx,px),
                           general_life_table_sc %>%
                             dplyr::select(age_group,sex,age,Lx,ex,Lwx,ewx),
                           by=c("age","sex","age_group")) %>%
    mutate(intervention="sc")
  
  baseline_d <- inner_join(disease_life_table_bl %>%
                             dplyr::select(age_group,sex,age,cause,incidence_disease,mx,px),
                           general_life_table_bl %>%
                             dplyr::select(age_group,sex,age,Lx,ex,Lwx,ewx),
                           by=c("age","sex","age_group")) %>%
    mutate(intervention="bl")
  
  disease_combined <- bind_rows(scenario_d,baseline_d) %>%
    pivot_wider(names_from  = intervention,
                values_from = c(incidence_disease,mx,px,Lx,ex,Lwx,ewx))  %>%
    mutate(inc_num_bl   = incidence_disease_bl*(1-px_bl)*Lx_bl,
           inc_num_sc   = incidence_disease_sc*(1-px_sc)*Lx_sc,
           inc_num_diff = inc_num_sc-inc_num_bl,
           mx_num_bl    = mx_bl*Lx_bl,
           mx_num_sc    = mx_sc*Lx_sc,
           mx_num_diff  = mx_num_sc-mx_num_bl) %>%
    dplyr::select(c(-mx_sc, -mx_bl, -px_sc, -px_bl, -Lx_sc, -Lx_bl, -ex_sc, -ex_bl, -Lwx_sc, -Lwx_bl, -ewx_sc, -ewx_bl)) %>%    dplyr::select(-starts_with(c( "incidence_"))) %>%
    pivot_wider(names_from  = cause,
                values_from = inc_num_bl:mx_num_diff)
  
  
  ### Non_diseases life tables
  
  scenario_nd <- inner_join(non_disease_life_table_sc %>%
                              dplyr::select(age_group,sex,age,cause,mx,ylds_rate, diff_mort, diff_pylds),
                            general_life_table_sc %>%
                              dplyr::select(age_group,sex,age,Lx),
                            by=c("age","sex","age_group")) %>%
    mutate(intervention="sc")
  
  baseline_nd <- inner_join(non_disease_life_table_sc %>%
                              dplyr::select(age_group,sex,age,cause,mx,ylds_rate, diff_mort, diff_pylds),
                            general_life_table_bl %>%
                              dplyr::select(age_group,sex,age,Lx),
                            by=c("age","sex","age_group")) %>%
    mutate(intervention="bl")
  
  non_disease_combined <- bind_rows(scenario_nd,baseline_nd) %>%
    pivot_wider(names_from  = intervention,
                values_from = c(mx, ylds_rate, Lx,
                                values_fill=0)) %>%
    mutate(ylds_num_bl   = ylds_rate_bl*Lx_bl,
           ylds_num_sc   = ylds_rate_sc*Lx_sc,
           ylds_num_diff = ylds_num_sc-ylds_num_bl,
           mx_num_bl    = mx_bl*Lx_bl,
           mx_num_sc    = mx_sc*Lx_sc,
           mx_num_diff  = mx_num_sc-mx_num_bl) %>%
    dplyr::select(c(cause, sex, age, age_group, ylds_num_bl, ylds_num_sc, ylds_num_diff, mx_num_bl, mx_num_sc, mx_num_diff)) %>%  
    pivot_wider(names_from  = cause,
                values_from = ylds_num_bl:mx_num_diff)
  
  ### General life tables
  
  general_lf <- bind_rows(
    general_life_table_sc %>%
      dplyr::select(age_group,sex,age,Lx,ex,Lwx,ewx) %>%
      mutate(intervention="sc"),
    general_life_table_bl %>%
      dplyr::select(age_group,sex,age,Lx,ex,Lwx,ewx) %>%
      mutate(intervention="bl")) %>%
    pivot_wider(names_from  = intervention,
                values_from = c(Lx,ex,Lwx,ewx), 
                values_fill=0) %>%
    mutate(Lx_diff  = Lx_sc-Lx_bl,
           Lwx_diff = Lwx_sc-Lwx_bl,
           ex_diff  = ex_sc-ex_bl,
           ewx_diff = ewx_sc-ewx_bl)
  
  
  ######## Dataframe with all outputs by age and sex cohort over the simulation years (years of the cohort)
  output_df <- inner_join(disease_combined,
                          general_lf,
                          by=c("age","sex","age_group")) %>%
    inner_join(non_disease_combined, 
               by=c("age","sex","age_group"))
  
  ### Generte outputs
  
  #### Summary data frame by age and sex and total 
  
  ######## Dataframe with all outputs aggregated by year of simlation by sex
  output_df_agg_sex  <- output_df   %>% ### Create a simulation year columns
    group_by(age_group, sex, .add=TRUE) %>%
    dplyr::mutate(year = 1:dplyr::n()) %>%
    ungroup() %>%
    dplyr::select(sex, year, Lx_bl, Lx_sc, Lx_diff, Lwx_bl, Lwx_sc, Lwx_diff, contains("num")) %>%
    group_by(year, sex, .add=TRUE) %>%
    dplyr::summarise_if(is.numeric, funs(sum)) %>%
    ungroup()
  # 
  # ######## Dataframe with all outputs aggregated by year of simulation all
  output_df_agg_all  <- output_df   %>% ### Create a simulation year columns
    group_by(age_group, sex, .add=TRUE) %>%
    dplyr::mutate(year = 1:dplyr::n()) %>%
    ungroup() %>%
    dplyr::select(sex, year, Lx_bl, Lx_sc, Lx_diff, Lwx_bl, Lwx_sc, Lwx_diff, contains("num")) %>%
    group_by(year, .add=TRUE) %>%
    dplyr::summarise_if(is.numeric, funs(sum)) %>%
    ungroup()
  
  ### Create age groups variable, easier to read
  
  output_df <- output_df %>%
    mutate(age_group_2 = case_when(
      age_group == 17 ~ "16-19",
      age_group == 22 ~ "20-24",
      age_group == 27 ~ "25-29",
      age_group == 32 ~ "30-34",
      age_group == 37 ~ "35-39",
      age_group == 42 ~ "40-44",
      age_group == 47 ~ "45-49",
      age_group == 52 ~ "50-54",
      age_group == 57 ~ "55-59",
      age_group == 62 ~ "60-64",
      age_group == 67 ~ "65-69",
      age_group == 72 ~ "70-74",
      age_group == 77 ~ "75-79",
      age_group == 82 ~ "80-84",
      age_group == 87 ~ "85-89",
      age_group == 92 ~ "90-94",
      age_group == 97 ~ "95 plus")) %>%
    mutate(cohort=paste(sex, age_group, sep = "_")) 
  
  
  population <- mslt_df %>% dplyr::select(age_cat, sex, population_number) %>%
    dplyr::filter(population_number!=0) %>%
    dplyr::rename(age_group=age_cat)
  
  output_df <- output_df %>% left_join(population)%>%
    dplyr::rename(`Age group` = age_group_2)
  
  
  ##################### Below outcomes for presentation ####################################################
  
  # Table: Life expectancy and health adjusted life expectancy 
  output_life_expectancy_change <- output_df[!duplicated(output_df$cohort), c("Age group", "cohort", "sex", "ex_bl", "ex_sc", "ewx_bl", "ewx_sc", "ex_diff", "ewx_diff", "population_number")] %>%
    dplyr::rename(`Life expectancy at baseline` = ex_bl, 
                  `Life expectancy scenario` = ex_sc, 
                  `Health adjusted life expectancy baseline` = ewx_bl, 
                  `Health adjusted life expectancy scenario` = ewx_sc) %>%
    dplyr::mutate(`Difference in life expectancy in days` = ex_diff * 365, 
                  `Difference in health adjusted life expectancy in days` = ewx_diff* 365) %>% 
    mutate_if(is.numeric, round, digits = 3) %>%
    dplyr::select(-c(ex_diff, ewx_diff, cohort)) %>%
    relocate(population_number, .after = sex)%>%
    dplyr::rename('Population cohort'=population_number)
  
  output_life_expectancy_change <- output_life_expectancy_change[order(output_life_expectancy_change$sex),] 
  
  # Table: Life years and health adjusted life years ----
  
  output_life_years_change <- output_df %>% 
    group_by(sex, `Age group`, cohort, .add=TRUE) %>%
    dplyr::summarise_if(is.numeric, funs(sum)) %>%
    ungroup() %>%
    dplyr::select(`Age group`,cohort, sex,Lx_diff, Lwx_diff, population_number) %>%
    dplyr::rename(`Life years` = Lx_diff, 
                  `Health adjusted life years` = Lwx_diff)  %>% 
    mutate_if(is.numeric, round) %>%
    relocate(population_number, .after = sex)%>%
    dplyr::rename('Population cohort'=population_number) %>%
    dplyr::select(-cohort)
  
  # Table: Diseases deaths, incidence and ylds ----
  
  output_diseases_change <- output_df %>% 
    group_by(sex, `Age group`, cohort, .add=TRUE) %>%
    dplyr::summarise_if(is.numeric, funs(sum)) %>% 
    mutate_if(is.numeric, round) %>%
    relocate(population_number, .after = sex)%>%
    dplyr::rename('Population cohort'=population_number) %>%
    dplyr::select(-cohort) %>%
    dplyr::select(`Age group`, sex, contains("diff")) %>%
    dplyr::select(-contains(c("Lx", "Lwx", "ex")))
  
  return(list(output_df=output_df,LifeYears=output_life_years_change, LifeExpectancy=output_life_expectancy_change, Cause=output_diseases_change, 
              output_df_agg_sex=output_df_agg_sex, output_df_agg_all=output_df_agg_all))
  
}

#' @export
combine_health_and_pif <- function(pif_values, hc=DISEASE_BURDEN){
  # pif_values are already ordered as in pop; reorder hc values to match.
  setorder(hc,dem_index)
  hm_cn_values <- hc$burden
  return_values <- hm_cn_values * pif_values
  round(as.vector(return_values),5)
}

gen_ap_rr <- function(pm_conc_pp){
  
  pm_rr_pp <- pm_conc_pp 
  
  ## assigning air pollution age band to the individual_level data
  min_ages <- c(seq(24,94,by=5),200)
  pm_rr_pp$age <- as.numeric(pm_rr_pp$age)
  pm_rr_pp$ap_age <- 0
  for(i in 1:length(min_ages)) pm_rr_pp$ap_age[pm_rr_pp$age>min_ages[i]] <- min_ages[i]+1
  
  pm_indices <- sapply(SCEN_SHORT_NAME,function(x)which(colnames(pm_rr_pp)==paste0("pm_conc_",x)))
  ### iterating over all all disease outcomes
  for ( j in c(1:nrow(DISEASE_INVENTORY))[DISEASE_INVENTORY$air_pollution == 1]){
    # initialise lists
    for (x in 1:length(SCEN_SHORT_NAME))
      pm_rr_pp[[paste0("RR_ap_", SCEN_SHORT_NAME[x])]] <- 1
    cause <- as.character(DISEASE_INVENTORY$ap_acronym[j])
    dr_ap_disease <- DR_AP[DR_AP$cause_code == cause,]
    # apply by age groups
    ages <- unique(dr_ap_disease$age_code)
    for(age in ages){
      dr_ap_sub <- dr_ap_disease[dr_ap_disease$age_code == age,]
      if(age==99){
        i <-1:nrow(pm_rr_pp)
      }else{
        i <- which(pm_rr_pp$ap_age==age)
      }
      
      age <- as.character(age)
      # print(paste(cause, age, sep = " - "))
      
      
      # get parameters
      alpha <- DR_AP_LIST[[cause]][[age]]$alpha
      beta <- DR_AP_LIST[[cause]][[age]]$beta
      gamma <- DR_AP_LIST[[cause]][[age]]$gamma
      tmrel <- DR_AP_LIST[[cause]][[age]]$tmrel
      # calculate AP and apply to all in age group
      for(x in 1: length(SCEN_SHORT_NAME)) 
        pm_rr_pp[[paste0("RR_ap_", SCEN_SHORT_NAME[x])]][i] <- ap_dose_response_curve(pm_rr_pp[[pm_indices[x]]][i],alpha,beta,gamma,tmrel)
    }
    ## change the names of the columns as per the disease
    for (n in 1: length(SCEN_SHORT_NAME)){
      col <- which(names(pm_rr_pp)== paste0("RR_ap_",SCEN_SHORT_NAME[n]))
      names(pm_rr_pp)[col]<- paste0("RR_ap_",SCEN_SHORT_NAME[n],"_",DISEASE_INVENTORY$acronym[j])
    }
  }
  pm_rr_pp
}
