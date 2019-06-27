#' @export
dist_dur_tbls <- function(pp_summary){
  
  durations <- list()
  
  ## calculate all distances & durations for each scenario
  for (i in 1:length(SCEN)){
    
    scen_travel <- pp_summary[[SCEN_SHORT_NAME[i]]][,colnames(pp_summary[[SCEN_SHORT_NAME[i]]])%in%paste0(VEHICLE_INVENTORY$stage_mode,'_dur')]
    total_travel <- colSums(scen_travel)
    
    # add Short Walking, if Short Walking has been added
    if("walk_to_bus"%in%names(total_travel)){
      total_travel[['walking_dur']] <-  total_travel[['walking_dur']] +  total_travel[['walk_to_bus_dur']]
    }
    
    durations[[i]] <- total_travel
  }
  dur <- do.call('cbind',durations)
  colnames(dur) <- SCEN
  rownames(dur) <- sapply(rownames(dur),function(x)strsplit(x,'_dur')[[1]][1])
  dur <- dur[rownames(dur)!='walk_to_bus',]
  
  mode_indices <- match(rownames(dur),VEHICLE_INVENTORY$stage_mode)
  mode_speeds <- VEHICLE_INVENTORY$speed[mode_indices]
  mode_speeds[is.na(mode_speeds)] <- 0
  dist <- dur * matrix(rep(mode_speeds,NSCEN+1),ncol=NSCEN+1) / 60
  
  ## bus travel is linear in bus passenger travel
  bus_passenger_row <- which(rownames(dur)=='bus')
  dist <- rbind(dist,dist[bus_passenger_row,] * BUS_TO_PASSENGER_RATIO)
  dur <- rbind(dur,dur[bus_passenger_row,] * BUS_TO_PASSENGER_RATIO) 
  rownames(dist)[nrow(dist)] <- 'bus_driver'
  rownames(dur)[nrow(dur)] <- 'bus_driver'
  ## truck travel is linear in car travel
  if(ADD_TRUCK_DRIVERS){
    car_row <- which(rownames(dur)=='car')
    dist <- rbind(dist,dist[car_row,] * TRUCK_TO_CAR_RATIO)
    dur <- rbind(dur,dur[car_row,] * TRUCK_TO_CAR_RATIO) 
    rownames(dist)[nrow(dist)] <- 'truck'
    rownames(dur)[nrow(dur)] <- 'truck'
  }
  
  return(list(dist=dist,dur=dur))
}

#' @export
scenario_pm_calculations <- function(dist,pp_summary){
  
  # concentration contributed by non-transport share (remains constant across the scenarios)
  non_transport_pm_conc <- PM_CONC_BASE*(1 - PM_TRANS_SHARE)  
  
  ## adding in travel not covered in the synthetic trip set, based on distances travelled relative to car, set in VEHICLE_INVENTORY
  emission_dist <- dist
  
  ## get emission factor by dividing inventory by baseline distance. (We don't need to scale to a whole year, as we are just scaling the background concentration.)
  ordered_efs <- VEHICLE_INVENTORY$emission_inventory[match(rownames(emission_dist),VEHICLE_INVENTORY$stage_mode)]/emission_dist[,'Baseline']
  ## get new emission by multiplying emission factor by scenario distance.
  trans_emissions <- emission_dist*t(repmat(ordered_efs,NSCEN+1,1))
  ## augment with travel emission contributions that aren't included in distance calculation
  ##!! need to sort out bus driver. Should be sorted in dist/dur? Or, if not sorted in dist/dur, sort here.
  for(mode_type in which(!VEHICLE_INVENTORY$stage_mode%in%rownames(emission_dist))){
    em <- VEHICLE_INVENTORY$emission_inventory[mode_type]
    if(em>0){
      trans_emissions <- rbind(trans_emissions,rep(em,ncol(trans_emissions)))
      rownames(trans_emissions)[nrow(trans_emissions)] <- VEHICLE_INVENTORY$stage_mode[mode_type]
    }
  }
  
  ## scenario travel pm2.5 calculated as relative to the baseline
  baseline_sum <- sum(trans_emissions[,SCEN[1]])
  conc_pm <- c()
  ## in this sum, the non-transport pm is constant; the transport emissions scale the transport contribution (PM_TRANS_SHARE) to the base level (PM_CONC_BASE)
  for(i in 1:length(SCEN_SHORT_NAME))
    conc_pm[i] <- non_transport_pm_conc + PM_TRANS_SHARE*PM_CONC_BASE*sum(trans_emissions[,SCEN[i]])/baseline_sum
  
  ##RJ rewriting ventilation as a function of MMET_CYCLING and MMET_WALKING, loosely following de Sa's SP model.
  vent_rates <- data.frame(stage_mode=VEHICLE_INVENTORY$stage_mode,stringsAsFactors = F) 
  vent_rates$vent_rate <- BASE_LEVEL_INHALATION_RATE  # L / min
  vent_rates$vent_rate[vent_rates$stage_mode=='bicycle'] <- BASE_LEVEL_INHALATION_RATE + MMET_CYCLING
  vent_rates$vent_rate[vent_rates$stage_mode%in%c('walking','walk_to_bus')] <- BASE_LEVEL_INHALATION_RATE + MMET_WALKING
  
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
  open_vehicles <- c('walking','walk_to_bus','bicycle','motorcycle','auto_rickshaw','shared_auto','cycle_rickshaw')
  rail_vehicles <- c('subway','rail')
  vent_rates$vehicle_ratio_index <- sapply(vent_rates$stage_mode,function(x) ifelse(x%in%rail_vehicles,3,ifelse(x%in%open_vehicles,1,2)))
  
  pp_summary2 <- lapply(pp_summary,function(y)y[,sapply(colnames(y),function(x)!grepl('_dist',x)),with=F])
  for(i in 1:length(pp_summary2)) colnames(pp_summary2[[i]]) <- sapply(colnames(pp_summary2[[i]]),function(x)gsub('_dur','',x))
  travel_indices <- which(colnames(pp_summary2[[1]])%in%vent_rates$stage_mode)
  travel_modes <- colnames(pp_summary2[[1]])[travel_indices]
  vent_modes <- match(travel_modes,vent_rates$stage_mode)
  
  # prepare individual-level dataset
  pm_conc_pp <- SYNTHETIC_POPULATION
  vent_multiplier <- repmat(vent_rates$vent_rate[vent_modes],nrow(pm_conc_pp),1)/60
  vent_and_ratio_multiplier <- vent_multiplier*repmat(ratio_by_mode[vent_rates$vehicle_ratio_index[vent_modes],1],nrow(pm_conc_pp),1)
  # compute individual-level pm scenario by scenario
  for (i in 1:length(SCEN)){
    
    scen_travel <- pp_summary[[i]]
    scen_travel[, on_road_dur := Reduce(`+`, .SD), .SDcols=travel_indices]
    #vent_travel <- scen_travel[,travel_indices,with=F] * vent_and_ratio_multiplier
    scen_travel[, on_road_pm := Reduce(`+`, lapply(seq_along(.SD)[travel_indices],function(x)(.SD[[x]]*vent_and_ratio_multiplier[,x]))), .SDcols=names(scen_travel)[travel_indices]]
    #vent_travel[, on_road_pm := Reduce(`+`, .SD), .SDcols=names(vent_travel)]
    
    ## PM2.5 inhalation = total mg inhaled / total volume inhaled
    # calculate non-travel air inhalation
    non_transport_air_inhaled <- (24-scen_travel$on_road_dur/60)*BASE_LEVEL_INHALATION_RATE
    # concentration of pm inhaled = total pm inhaled / total air inhaled
    pm_conc <- ((non_transport_air_inhaled * as.numeric(conc_pm[i])) + scen_travel$on_road_pm)#/(non_transport_air_inhaled+individual_data$air_inhaled)
    # match individual ids to set per person pm exposure
    pm_conc_pp[[paste0('pm_conc_',SCEN_SHORT_NAME[i])]] <- pm_conc/24 * conc_pm[i]
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
  pp_summary2 <- lapply(pp_summary,function(y)y[,sapply(colnames(y),function(x)!grepl('_dist',x)),with=F])
  for(i in 1:length(pp_summary2)) colnames(pp_summary2[[i]]) <- sapply(colnames(pp_summary2[[i]]),function(x)gsub('_dur','',x))
  # Get total individual level walking and cycling and sport mmets 
  synth_pop_return <- pp_summary2[[1]]
  for (i in 1:length(SCEN)){
    synth_pop_temp <- pp_summary2[[i]]
    synth_pop_return[[paste0(SCEN_SHORT_NAME[i],'_mmet')]] <- synth_pop_temp$work_ltpa_marg_met * BACKGROUND_PA_SCALAR
    
    scen_travel <- subset(pp_summary2[[i]],participant_id%in%synth_pop_return$participant_id)
    scen_travel$cycling_mmet <- scen_travel$bicycle/60*DAY_TO_WEEK_TRAVEL_SCALAR * MMET_CYCLING
    if('walk_to_bus'%in%names(scen_travel)) scen_travel$walking <- scen_travel$walking+scen_travel$walk_to_bus
    scen_travel$walking_mmet <- scen_travel$walking/60*DAY_TO_WEEK_TRAVEL_SCALAR * MMET_WALKING
    
    individual_data <- scen_travel
    
    part_id <- match(individual_data$participant_id,synth_pop_return$participant_id)
    synth_pop_return[[paste0(SCEN_SHORT_NAME[i],'_mmet')]][part_id] <- 
      synth_pop_return[[paste0(SCEN_SHORT_NAME[i],'_mmet')]][part_id] + individual_data$cycling_mmet + individual_data$walking_mmet
  }
  
  name_indices <- which(colnames(synth_pop_return)%in%c('participant_id', 'sex', 'age', 'dem_index', paste0(SCEN_SHORT_NAME,'_mmet')))
  mmets <- tbl_df(synth_pop_return)[,name_indices]
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
predict_injuries <- function(city_table){
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
injury_death_to_yll <- function(injuries){
  
  joined_injury <- left_join(injuries, GBD_INJ_YLL[,c('dem_index','yll_dth_ratio')], by="dem_index")
  
  joined_injury$YLL <- joined_injury$Deaths*joined_injury$yll_dth_ratio
  death_and_yll <- dplyr::select(joined_injury, c('dem_index','scenario','Deaths','YLL'))
  
  x_deaths <- dplyr::select(death_and_yll, -YLL)
  x_deaths <- spread(x_deaths,scenario, Deaths)
  x_yll <- dplyr::select(death_and_yll, -Deaths)
  x_yll <- spread(x_yll,scenario, YLL)
  
  ref_scen <- REFERENCE_SCENARIO
  ref_scen_index <- which(SCEN==ref_scen)
  calc_scen <- SCEN[SCEN!=ref_scen]
  calc_scen_index <- which(colnames(x_deaths)%in%calc_scen)
  
  ref_injuries <- as.data.frame(cbind(dem_index=x_deaths$dem_index,deaths=x_deaths[[ref_scen]],ylls=x_yll[[ref_scen]]))
  deaths <- t(repmat(unlist(ref_injuries$deaths),NSCEN,1)) - x_deaths[,calc_scen_index,drop=F]
  ylls <- t(repmat(unlist(ref_injuries$ylls),NSCEN,1)) - x_yll[,calc_scen_index,drop=F]
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
combine_health_and_pif <- function(pif_values, hc=DISEASE_BURDEN){
  # pif_values are already ordered as in pop; reorder hc values to match.
  setorder(hc,dem_index)
  hm_cn_values <- hc$burden
  return_values <- hm_cn_values * pif_values
  round(as.vector(return_values),5)
}
