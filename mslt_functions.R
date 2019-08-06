# ---- Packages for functions ----

#require(dplyr)
#require(tidyverse)
#require(knitr)
#require(kableExtra)
#require(citr)

# ---- Explanation method ---- 

# Briefly, the proportional multi-state multi cohort life table consists of a 
# general life table and a life table for each of the modelled diseases.
# The diseases are those associated to the studied risk factor/s. 
# The link between the general life table and disease life tables is via the 
# potential impact fraction (pif), also called paf (population attributable fraction)
# The pif combines exposure to the risk factor and relative risks. The pif is 
# appleid to modify incidence in the individual disease life tables, which in turn
# modify prevalence and mortality. Changes in mortality and prevalence rates
# feed bacak into the general life table to modify total mortality and disability. 
# Changes in total mortality impact on life years and changes in disability impact 
# the disability adjusment of life years. 

# Method reference: 1.	Barendregt JJ, Oortmarssen vGJ, Murray CJ, Vos T. A generic model for the assessment of disease epidemiology: the computational basis of DisMod II. Popul Health Metr. 2003;1(1):4-.
# Naming convention for functions: Function.Name

# ---- Functions ----

# SortGbdInput, RunLocDf, RunLifeTable, RunDisease, RunPif, RunOutput run_life_table, run_disease, run_pif (temp), run_output 

# ---- Sort.Gbd.Input ----
#' @export
SortGbdInput <- function(in_data, in_year, in_locality) {
  data <- in_data[which(in_data$year== in_year & in_data$location == in_locality),]

}

## Selects year and localities from GBD data frame dowloaded from: http://ghdx.healthdata.org/gbd-results-tool

# --- RunLocDf ----
#' @export
RunLocDf <- function(i_data) {
  
  gbd_df <- NULL 
  
  
  for (ag in 1:length(unique(i_data$age))){
    for (gender in c("Male", "Female")){
      age_sex_df <- NULL
      for (dm in 1:length(disease_measures_list)){
        for (d in 1:nrow(disease_short_names)){
          dn <- disease_short_names$disease[d]
          dmeasure <- disease_measures_list[dm] %>% as.character()
          # gender <- "Male"
          
          agroup <- unique(i_data$age)[ag]
          
          idf <- filter(i_data, sex == gender & age == agroup & measure == dmeasure & cause == dn) 
        
          
          
          if (nrow(idf) > 0){
            
            population_numbers <- filter(idf, metric == "Number") %>% dplyr::select("val")
            
            idf_rate <- filter(idf, metric == "Rate") %>% dplyr::select("val") 
            
            current_idf_rate <- idf_rate
            
            current_population_numbers <- population_numbers
            
            
            idf$population_number <- 0
            
            if (idf_rate$val != 0 && population_numbers$val != 0)
              idf$population_number <- (100000 * population_numbers$val) / idf_rate$val
            
            else{
              
              current_idf_rate <- idf_rate
              
              current_population_numbers <- population_numbers
              
              idf <- filter(i_data, sex == gender & age == agroup & measure == dmeasure & val > 0) 
              
              idf <- filter(idf, cause == unique(idf$cause)[1])
              
              idf$cause <- dn
              
              population_numbers <- filter(idf, metric == "Number") %>% dplyr::select("val")
              
              idf_rate <- filter(idf, metric == "Rate") %>% dplyr::select("val") 
              
              # browser()
              
              idf$population_number <- 0
              
              if (idf_rate$val != 0 && population_numbers$val != 0)
                idf$population_number <- (100000 * population_numbers$val) / idf_rate$val
              
            }
            
            # if (is.nan(idf$population_number)){
            #   browser()
            # }
            
            idf$rate_per_1 <- round(current_idf_rate$val / 100000, 6)
            
            
            idf[[tolower(paste(dmeasure, "rate", disease_short_names$sname[d], sep = "_"))]] <- idf$rate_per_1
            
            idf[[tolower(paste(dmeasure, "number", disease_short_names$sname[d], sep = "_"))]] <- current_population_numbers$val
            
            #print(unique(i_data$age)[ag])
            #print(gender)
            
            # if (ag == "Under 5" && gender == "Female"){
            #   browser()
            # }
            #idf$rate_per_1 <- NULL
            
            idf <- filter(idf, metric == "Number")
            
            if (is.null(age_sex_df)){
              #browser()
              # print("if")
              # print(names(idf)[ncol(idf) - 1])
              # print(names(idf)[ncol(idf)])
              age_sex_df <- dplyr::select(idf, age, sex, population_number, location, names(idf)[ncol(idf) - 1] , names(idf)[ncol(idf)])
            }
            else{
              #browser()
              # print("else")
              # print(names(idf)[ncol(idf) - 1])
              # print(names(idf)[ncol(idf)])
              
              age_sex_df <- cbind(age_sex_df, dplyr::select(idf, names(idf)[ncol(idf) - 1] , names(idf)[ncol(idf)]))
            }
            
          }
          
          # age_range <- years %>% str_match_all("[0-9]+") %>% unlist %>% as.numeric
          
        }
      }
      
      # browser()
      
      if (is.null(gbd_df)){
        # browser()
        gbd_df <- age_sex_df
      }
      else{
        # browser()
        age_sex_df[setdiff(names(gbd_df), names(age_sex_df))] <- 0
        gbd_df[setdiff(names(age_sex_df), names(gbd_df))] <- 0
        
        gbd_df <- rbind(gbd_df, age_sex_df)
        
        
      }
    }
  }
  return(gbd_df)
}


# --- InterFunc --- (Belen: within data preparation, if not in data preparation interpolation code, it does not work)

### Embeded in loop

# InterFunc <- stats::splinefun(x, y, method = "monoH.FC", ties = mean)

# --- IsNanDataFrame ---
#' @export
IsNanDataFrame <- function(x)
  do.call(cbind, lapply(x, is.nan))

# --- IsInfDataFrame ---
#' @export
IsInfDataFrame <- function(x)
  do.call(cbind, lapply(x, is.infinite))



## Sorts out each locality data frame into a list with column names for age, sex, and each cause and disease combination and calculates population numbers.

# ---- RunLifeTable ----

#####Function to generate age and sex life tables. The function is then use in the model.R script
#####to calculate variables for the baseline and scenario life tables. 
#' @export
RunLifeTable <- function(in_idata, in_sex, in_mid_age)
{
  
  # Create a life table data frame
  
  # Filter in_idata by age and dplyr::select columns for lifetable calculations
  #lf_df <- filter(in_idata, age >= in_mid_age & sex == in_sex) %>% dplyr::select(sex, age, pyld_rate, mx)
  lf_df <- in_idata[in_idata$age >= in_mid_age & in_idata$sex == in_sex,] 
  lf_df <- lf_df[,colnames(lf_df)%in%c('sex', 'age', 'pyld_rate', 'mx')]
  #col_names <- c('sex','age','pyld_rate','mx')
  #lf_df <- lf_df[,colnames(lf_df)%in%col_names]
  
  # Create list of required columns (variables)
  
  # probability of dying
  #lf_df$qx <-  ifelse(lf_df$age < 100, 1 - exp(-1 * lf_df$mx), 1)
  qx <-  ifelse(lf_df$age < 100, 1 - exp(-1 * lf_df$mx), 1)
  
  # number of survivors
  #lf_df$lx <- 0
  num_row <- nrow(lf_df)
  lx <- rep(0,num_row)
  # Create it for males population
  #lf_df$lx[1] <- as.numeric(in_idata$population_number[in_idata$age == in_mid_age & in_idata$sex == in_sex]) # filter(in_idata, age == in_mid_age & sex == in_sex) %>% dplyr::select(population_number)
  #lf_df$lx <- as.numeric(lf_df$lx)
  lx[1] <- as.numeric(in_idata$population_number[in_idata$age == in_mid_age & in_idata$sex == in_sex]) # filter(in_idata, age == in_mid_age & sex == in_sex) %>% dplyr::select(population_number)
  
  # number died
  #lf_df$dx <- 0
  dx <- rep(0,num_row)
  
  # Create it for males population
  #lf_df$dx[1] <- lf_df$lx[1] * lf_df$qx[1]
  dx[1] <- lx[1] * qx[1]
  
  for (i in 2:num_row){
    #lf_df$lx[i] <- lf_df$lx[i - 1] - lf_df$dx[i - 1]
    #lf_df$dx[i] <- lf_df$lx[i] * lf_df$qx[i]
    lx[i] <- lx[i - 1] - dx[i - 1]
    dx[i] <- lx[i] * qx[i]
  }
  
  # number of persons lived by cohort to age x + 1/2 (average people)
  #lf_df$Lx <- 0
  Lx <- rep(0,num_row)
  
  for (i in 1:(num_row-1))
    Lx[i] <- (lx[i] + lx[i + 1]) / 2
    #lf_df$Lx[i] <- (lf_df$lx[i] + lf_df$lx[i + 1]) / 2
  Lx[num_row] <- lx[num_row] / lf_df$mx[num_row]
  #lf_df$Lx[num_row] <- lf_df$lx[num_row] / lf_df$mx[num_row]
  
  
  # create life expectancy variable
  ex <- rep(0,num_row)
  for (i in 1:num_row){
    #lf_df$ex[i] <- sum(lf_df$Lx[i:nrow(lf_df)]) / lf_df$lx[i]
    ex[i] <- sum(Lx[i:num_row]) / lx[i]
  }
  
  # create health adjusted life years variable 
  
  #lf_df$Lwx <- lf_df$Lx * (1 - lf_df$pyld_rate)
  Lwx <- Lx * (1 - lf_df$pyld_rate)
  
  # create health adjusted life expectancy variable
  ewx <- rep(0,num_row)
  for (i in 1:num_row){
    #lf_df$ewx[i] <- sum(lf_df$Lwx[i:nrow(lf_df)]) / lf_df$lx[i]
    ewx[i] <- sum(Lwx[i:num_row]) / lx[i]
  }
  
  lf_df$qx <- qx
  lf_df$lx <- lx
  lf_df$dx <- dx
  lf_df$Lx <- Lx
  lf_df$ex <- ex
  lf_df$Lwx <- Lwx
  lf_df$ewx <- ewx
  lf_df
}

# ---- RunDisease ----(ADD REMISSION FOR CANCERS?)

RunDisease <- function(in_idata, in_mid_age, in_sex, in_disease) 
  
{
  
  # Uncomment the variables below to debug your code  
  # in_idata = sub_idata
  # in_sex = "males"
  # in_mid_age = 22
  # in_disease = "ihd"
  
  # create disease variable for the disease life table function 
  dw_disease <- paste("dw_adj", in_disease, sep = "_")
  incidence_disease <- paste("incidence", in_disease, sep = "_")
  case_fatality_disease <- paste("case_fatality", in_disease, sep = "_")
  
  ## add generic variable names to the source data frame (in_idata)
  in_idata$dw_disease <- in_idata[[dw_disease]]
  in_idata$incidence_disease <- in_idata[[incidence_disease]]
  in_idata$case_fatality_disease <- in_idata[[case_fatality_disease]]
  
  # filter in_idata by age and select columns for lifetable calculations 
  ##!!RJ filtered before passed to function
  #dlt_df <- in_idata[in_idata$age >= in_mid_age & in_idata$sex == in_sex,] #%>% 
  #print(c(nrow(dlt_df),nrow(in_idata)))
  dlt_df <- in_idata[,colnames(in_idata)%in%c('sex', 'age', 'dw_disease', 'incidence_disease', 'case_fatality_disease')] # dplyr::select(sex, age, dw_disease, incidence_disease, case_fatality_disease)
  
  dlt_df$disease <- in_disease
  
  # create list of required columns
  ## intermediate variables lx, qx, wx and vx
  ###lx
  
  #browser()
  lx <- dlt_df$incidence_disease + dlt_df$case_fatality_disease
  ###qx
  qx <-  sqrt((dlt_df$incidence_disease - dlt_df$case_fatality_disease) * (dlt_df$incidence_disease - dlt_df$case_fatality_disease))
  ### wx
  wx <- exp(-1*(lx+qx)/2)
  ### vx
  vx <- exp(-1*(lx-qx)/2)
  
  ## Healthy (Sx), Disease (Cx) and Death (Dx), total (Tx) (control check, has to be 1000), total alive (Ax)
  ## persons years live at risk (PYx), prevalence rate (px), mortality rate (mx)
  ### Remission and mortality from other causes were replaced by zero in the formulas (as we assume no remission and independence of disease mortality with total mortlaity). 
  
  #### first create empty variables
  
  number_of_ages <- nrow(dlt_df)
  Sx <- Cx <- Dx <- Tx  <- Ax <- PYx <- px <- mx <- rep(0,number_of_ages)
  cfds <- dlt_df$case_fatality_disease
  ages <- dlt_df$age
  ## set initial conditions
  # Dx, Cx, PYx, px, mx stay zero
  Sx[1] <- Ax[1] <- 1000
  ##### start with variables without calculation exceptions
  
  ##### variables with exceptions  
  for (i in 2:(number_of_ages-1)){ ##!! this can go to "number_of_ages" now (?)
    if(qx[i-1] > 0){
      vxmwx <- vx[i-1] - wx[i-1]
      SxpCx <- Sx[i-1]+Cx[i-1]
      dqx <- 2 * qx[i-1]
      qxmlx <- qx[i-1] - lx[i-1]
      qxplx <- qx[i-1] + lx[i-1]
      Sx[i] <- Sx[i-1] * (2*vxmwx * cfds[i-1]  + (vx[i-1] * qxmlx + wx[i-1] * qxplx)) / dqx
      Cx[i] <- -1*(vxmwx*(2*(cfds[i-1]  * SxpCx - lx[i-1] * Sx[i-1]) - Cx[i-1] * lx[i-1]) - Cx[i-1] * qx[i-1] * (vx[i-1]+wx[i-1])) / dqx
      Dx[i] <- (vxmwx * (2 * cfds[i-1] * Cx[i-1] - lx[i-1]*SxpCx)- qx[i-1] * SxpCx*(vx[i-1]+wx[i-1]) + dqx * (SxpCx+Dx[i-1]) ) / dqx
    }else{
      Sx[i] <- Sx[i - 1] 
      Cx[i] <- Cx[i - 1]
      Dx[i] <- Dx[i - 1]
    }
  }
  Tx   <- Sx + Cx + Dx 
  Ax <- Sx + Cx
  first_indices <- 1:(number_of_ages-1)
  last_indices <- 2:number_of_ages
  PYx <- (Ax[first_indices] + Ax[last_indices])/2
  mx[first_indices] <- (Dx[last_indices] - Dx[first_indices])/PYx[first_indices]
  mx[mx<0] <- 0
  px[first_indices] <- (Cx[last_indices] + Cx[first_indices])/2/PYx[first_indices]
  #for (i in 1:(number_of_ages-1)){
    #if ((Dx[i+1] - Dx[i]) < 0){
    #  mx[i] <- 0
    #}else{
    #  mx[i] <- (Dx[i+1] - Dx[i])/PYx[i]
    #}
  #  Csum <- Cx[i]+Cx[i+1]
  #  px[i] <- Csum/2/ PYx[i]   
  #}
  
  dlt_df$Tx <- Tx
  dlt_df$mx <- mx
  dlt_df$px <- px
  dlt_df
}


# Run non_diseases
#' @export
RunNonDisease <- function(in_idata, in_sex, in_mid_age, in_non_disease)

  {

  # deaths_rate <- paste("deaths_rate", in_non_disease, sep = "_")
  # pyld_rate <- paste("ylds (years lived with disability)_rate", in_non_disease, sep = "_")
  # 
  # 
  # # df$deaths_rate <- df[[deaths_rate]]
  # # df$pyld_rate <- df[[pyld_rate]]

  ##!!RJ filtered before passing
  #df <- filter(in_idata, age >= in_mid_age & sex == in_sex) #%>%
  df <- in_idata[,colnames(in_idata)%in%c('sex', 'age',  paste0("deaths_rate_", in_non_disease), paste0("ylds_rate_", in_non_disease))]
  
  
  
  return(df)
}


# GetPif (for Metahit)
#' @export
GetPif <- function(in_pif, in_age, in_sex, pif_name){
  p <- as.data.frame(in_pif)
  ##!!RJ filtered before passed to function
  #p <- df[df$age >= in_age & df$sex == in_sex,]# %>% 
  #print(c(dim(df),dim(p)))
  p <- p[,colnames(p)%in%c('age', pif_name)] # dplyr::select(age, pif_name)
  
  ## Expand to repeat values between age groups, for example, same value from 17 to 21
  
  outage <- min(p$age):100
  
  ind <- findInterval(outage, p$age)
  p <- p[ind,]
  p$age <- outage
  
  return(p)
}

# 
# test_pif <- GetPif(pif, 27, "male", "scen_pif_pa_ac")


# RunPif (temp) ----

# The code for PIFs will depend on the data sources. 

#' @export
RunPif <- function(in_idata, i_irr, i_exposure, in_mid_age, in_sex, in_disease, in_met_sc) 
  # 
{
  
  ## uncomment to debug function
  
  # in_idata = idata
  # i_irr = irr
  # i_exposure = edata
  # in_sex = "females"
  # in_mid_age = 22
  # in_disease = "breast_cancer"
  # in_met_sc <- effect
  
  ### filter data to use in pif calculations (age and sex). Add rrs, ee and calculations
  
  pif_df <- filter(in_idata, age >= in_mid_age & sex == in_sex) %>%
    dplyr::select(sex, age)
  
  ### add varaibles to data.frame (different age category for breast cancer)
  
  pif_df$disease <- in_disease
  
  if(in_disease == "breast_cancer") {
    pif_df$age_cat [pif_df$age <=30] <- 30
    pif_df$age_cat [pif_df$age >30 & pif_df$age <=45 ] <- 45
    pif_df$age_cat [pif_df$age >45 & pif_df$age <=70 ] <- 70
    pif_df$age_cat [pif_df$age >70 & pif_df$age <=100 ] <- 80
  }
  else {
    pif_df$age_cat [pif_df$age <=30] <- 30
    pif_df$age_cat [pif_df$age >30 & pif_df$age <=70 ] <- 70
    pif_df$age_cat [pif_df$age >70 & pif_df$age <=100 ] <- 80
  }
  
  pif_df$sex_cat <- ifelse(in_disease == "breast_cancer", "female", "female_male")
  
  # create concatenated variables to match pif_df with i_irr
  pif_df$sex_age_dis_cat <- paste(pif_df$disease,pif_df$age_cat, pif_df$sex_cat, sep = "_"  )
  i_irr$sex_age_dis_cat <- paste(i_irr$disease,i_irr$age, i_irr$sex, sep = "_"  )
  
  # remove sex, age and disease variables from i_irr df, as they are not needed
  i_irr <- dplyr::select(i_irr, -one_of('sex','age', 'disease'))
  
  # the code below is working but copies age, sex and disease for x and y, how can this be avoided?
  pif_df <-  inner_join(pif_df, i_irr, by = c("sex_age_dis_cat" = "sex_age_dis_cat") , copy = FALSE)
  
  # creation of splineFun which uses baseline's RR and EE to use to estimate intervention RRs
  for (i in 1:nrow(pif_df)){
    sp_obj <-  splinefun(y = c(pif_df$rr_inactive[i], 
                               pif_df$rr_insufficiently_active[i], 
                               pif_df$rr_recommended_level_active[i], 
                               pif_df$rr_highly_active[i]), 
                         x = c(pif_df$ee_inactive[i], 
                               pif_df$ee_insufficiently_active[i], 
                               pif_df$ee_recommended_level_active[i], 
                               pif_df$ee_highly_active[i]), 
                         method = "hyman")
    
    # use created spline function above to estimate intervention RRs
    
    pif_df$sc_rr_inactive[i] <- sp_obj(x = pif_df$ee_inactive[i] + in_met_sc)
    pif_df$sc_rr_insufficiently_active[i] <-  sp_obj(x = pif_df$ee_insufficiently_active[i] + in_met_sc)
    pif_df$sc_rr_recommended_level_active[i] <-  sp_obj(x = pif_df$ee_recommended_level_active[i] + in_met_sc)
    pif_df$sc_rr_highly_active[i] <-  sp_obj(x = pif_df$ee_highly_active[i] + in_met_sc)
    
    # plot(sp_obj, xlab = "RR", ylab = "EE", main = paste("Spline ", i))
  }
  
  
  ## round sc_rr_highly_active column - it should be 1
  pif_df$sc_rr_highly_active <- round(pif_df$sc_rr_highly_active)
  
  ##Calculate PIFs. I already process the data to generate categories in stata.
  ##First add PA categories to pif_df
  
  pif_df$sex_age_cat <- paste(pif_df$sex, pif_df$age, sep = "_"  )
  i_exposure$sex_age_cat <- paste(i_exposure$sex, i_exposure$age, sep = "_"  )
  
  # remove sex, age and disease variables from i_irr df, as they are not needed
  i_exposure <- dplyr::select(i_exposure, -one_of('sex','age'))
  
  # join edata (PA prevalence to pif_df)
  
  pif_df <-  inner_join(pif_df, i_exposure, by = c("sex_age_cat" = "sex_age_cat") , copy = FALSE)
  
  
  # we need to adapt to ITHIMR developments. REPLACE DATA FRAME FROM WHICH PREVALENCE OF PA IS TAKEN
  
  pif_df$pif <- 1-(pif_df$sc_rr_inactive *pif_df$inactive +
                     pif_df$sc_rr_insufficiently_active*pif_df$insufficiently_active +
                     pif_df$sc_rr_recommended_level_active*pif_df$recommended_level_active +
                     pif_df$sc_rr_highly_active *pif_df$highly_active)/
    (pif_df$rr_inactive *pif_df$inactive  +
       pif_df$rr_insufficiently_active *pif_df$insufficiently_active +
       pif_df$rr_recommended_level_active *pif_df$recommended_level_active +
       pif_df$rr_highly_active *pif_df$highly_active)
  
  
  pif_df
  
}


# ---- PlotOutput ----

# Function to generate graphs by age and sex, per outcome of interest. 

#' @export
PlotOutput <- function(in_data, in_age, in_population, in_outcomes, in_legend = "", in_disease = ""){
  
  # in_data <- output_df
  # in_population <- "male"
  # in_age <- 22
  # in_outcomes <- c('age', 'inc_num_bl_ihd', 'inc_num_sc_ihd')
  # in_legend <- "none"
  # in_cols <- c('alpha', 'beta')
  
  data <- in_data
  
  if (in_population != "total")
    data <- filter(data, sex == in_population)
  if (length(in_age) > 0)
    data <- filter(data, age_cohort == in_age)
  if (length(in_outcomes) > 0)
    data <- dplyr::select(data, in_outcomes)
  
  td <- data
  p <- ggplot(data = td, aes (x = td[[in_outcomes[[1]]]]))
  
  # loop
  for (i in 2:length(in_outcomes)) {
    # use aes_string with names of the data.frame
    p <- p + geom_line(aes_string(y = td[[in_outcomes[i]]], color = as.factor(in_outcomes[i])), size = 0.8) +
      
      theme_classic() 
    
    
  }
  
  p <- p + scale_color_discrete(name = paste(in_legend), labels = c("Baseline", "Difference", "Scenario")) +
    theme(legend.title = element_text(size = 9))
  
  p <- p + xlab ('Age') + ylab ('Cases') + labs (title = ifelse(length(in_disease) > 0, 
                                                                in_disease, paste('Cohort', in_age, "years old", in_population, sep = " "))) +
    theme(plot.title = element_text(hjust = 0.5, size = 9)) +
    theme(legend.text = element_text(size = 9)) +
    # theme(axis.title.x = element_text(size = 7)) +
    xlim(in_age, 100) +
    geom_hline(yintercept=0, linetype="dashed", color = "black")
  
  
  return(p)
  
  last_plot()
  
  
}


# ---- GenAggregate ----
# Function to aggreate outcomes by age an sex

#' @export
GenAggregate <- function(in_data, in_cohorts, in_population, in_outcomes){
  
  
  # in_data <- output_df
  # in_population <- "males"
  # in_cohorts <- 10
  # in_outcomes <- c('inc_num_bl_ihd', 'inc_num_sc_ihd')
  
  age_cohort_list <- list()
  td <- in_data
  aggr <- list()
  l_age <-  min(td$age_cohort)
  for (i in 1:in_cohorts){
    if (l_age <= 100){
      ld <- dplyr::filter(td, age_cohort == l_age)
      
      if (in_population != "total")
        ld <- filter(ld, sex == in_population)
      if (length(in_outcomes) > 0)
        ld <- dplyr::select(ld, age, sex, in_outcomes)
      if (i == 1){
        aggr <- append(aggr, as.list(ld))
        aggr <- as.data.frame(aggr)
        names(aggr) <- paste(names(aggr), l_age, in_population, sep = "_" )
      }
      else {
        n_rows <-  nrow(aggr) - nrow(ld)
        ld[(nrow(ld) + 1):(nrow(ld) + n_rows),] <- NA
        names(ld) <- paste(names(ld), l_age, in_population, sep = "_" )
        aggr <- cbind(aggr, ld)
      }
      
      l_age <- l_age + 5
    }
  }
  
  for (i in 1:length(in_outcomes)){
    aggr[[paste0("total_",in_outcomes[i])]] <- dplyr::select(aggr, starts_with(in_outcomes[i])) %>% rowSums(na.rm = T)
    
  }
  
  aggr
}

# ---- GridArrangSharedLegend ----
# Function to general combined labels for multiple plots in a page
#' @export
GridArrangSharedLegend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right"), mainTitle = "", mainLeft = "", mainBottom = "") {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow, top = mainTitle, left = mainLeft, bottom = mainBottom)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  
  #grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

#' @export
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# ---- GetQualifiedDiseaseName ---- (Belen: check if we need this function)
# Function to get qualified names diseases
#' @export
GetQualifiedDiseaseName <- function (disease){
  if (disease == 'ihd')
    return ('Ischaemic Heart Disease')
  else if (disease == 'bc')
    return ('Breast Cancer')
  else if (disease == 'dm')
    return ('Diabetes')
  else if (disease == 'cc')
    return ('Colon cancer')
  else if (disease == 'is')
    return ('Ischemic stroke')
}

# ---- PlotGBD (may need to update) ----
# Function to generate GBD graphs to compare data national to local (USED in GBD COMPARE############################
#' @export
PlotGBD <- function(in_data1, in_data2, in_sex, in_cause, in_measure) {
  
  # in_data1 <- GBDEngland
  # in_data2 <- GBDGL
  # in_sex <- "male"
  # in_cause <- "all causes"
  # in_measure <- "deaths"
  
  
  data1 <- filter(in_data1, sex == in_sex, cause == in_cause & measure == in_measure) %>% dplyr::select(measure, location, sex, age, metric, cause, one_rate, age_cat)     
  
  data2 <- filter(in_data2, sex == in_sex, cause == in_cause & measure == in_measure) %>% dplyr::select(measure, location, sex, age, metric, cause, one_rate, age_cat)     
  
  
  p <- ggplot(data = data1, aes(age_cat,one_rate)) +
    geom_line(aes(color = "England"))+
    geom_line(data = data2, aes(color = "Greater London"))+
    labs(colour="Locations",x="Age",y= paste(in_cause, in_measure, sep = " "))+
    labs (title = paste("Compare", in_cause, in_measure, in_sex, sep = " "), size=14) + 
    theme_classic()
  print(p)
}

