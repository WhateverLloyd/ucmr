#American Community Survey Data Download

#This script downloads, processes, and organizes 5-year american community 
#survey (ACS) estimates for all US census tracts for the year 2020 via the 
#'tidycensus' R package. ACS variables were selected to represent a wide range 
#of social, demographic, and economic factors in the US population.

#load packages
library(sf)
library(tidycensus)
library(tidyverse)

#set working directory
setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - PFAS + Cancer/UCMR_EA")

#import list of ACS variables to download from tidycensus
acsvars <- read.csv(file = "./data/variables_acs5.csv")

#import list of states with census regions and divisions
stateinfo <- read.csv(file = "./data/StatesFIPS_CensusRegions.csv")
##fix FIPS codes with leading 0's
stateinfo$FIPS <- as.character(str_pad(stateinfo$FIPS, width = 2, side = "left", pad = "0"))
##remove Virgin Islands and Puerto Rico  
stateinfo <- stateinfo %>% filter(!Postal %in% c("VI","PR"))
##get list of states
states <- unique(stateinfo$FIPS)

#tidycensus's 'get_acs' tool only lets you download tract-level data one state
#at a time. Use 'for loop' to download variables from all states and combine
#into single data frame.

##create blank dataframe to store loop outputs
dat2 <- as.data.frame(NULL)

for(i in 1:length(states)){
  
  dat <- NULL
  
  #select state
  st <- states[i]
  
  #download tract level estimates for year 2019 for each state for list of variables
  dat <- get_acs(geography = "tract", state = st, variables = acsvars$acs_code, year = 2020)
  
  #replace variable codes with more descriptive variable labels
  dat <- left_join(dat, acsvars, by = c("variable"="acs_code"))
  
  dat <- dat %>% select(GEOID, variable_label, estimate)
  
  dat2 <- rbind(dat, dat2)
  
}

#reshape data from long to wide giving each variable its own column populated
#with ACS estimates
dat3 <- dat2 %>%
        pivot_wider(names_from = variable_label, values_from = estimate)

#some variables need to be converted to percent by dividing total counts by the
#total population. other variables (eg, no high school education) need to be 
#calculated using estimates from multiple variables.  
acsdat <- dat3 %>%
          mutate(UNDER5_prct = (UNDER5_Mpop + UNDER5_Fpop)/TOTAL_POP*100,
            WHITE_prct = WHITE_total/TOTAL_POP*100,
            BLACK_prct = BLACK_total/TOTAL_POP*100,
            AMIND_prct = AMIND_total/TOTAL_POP*100,
            ASIAN_prct = ASIAN_total/TOTAL_POP*100,
            PACISL_prct = PACISL_total/TOTAL_POP*100,
            OTHER_prct = OTHER_total/TOTAL_POP*100,
            MIXED_prct = MIXED_total/TOTAL_POP*100,
            HISPLAT_prct = HISPLAT_total/TOTAL_POP*100,
            POC_prct = (ALLRACE_pop - WHITE_NH_pop)/ALLRACE_pop*100,
            NOHS_prct = (SCHOOL0_Mpop + SCHOOL1_Mpop + SCHOOL2_Mpop + 
              SCHOOL3_Mpop + SCHOOL4_Mpop + SCHOOL5_Mpop + SCHOOL6_Mpop + 
              SCHOOL7_Mpop + SCHOOL0_Fpop + SCHOOL1_Fpop + SCHOOL2_Fpop + 
              SCHOOL3_Fpop + SCHOOL4_Fpop + SCHOOL5_Fpop + SCHOOL6_Fpop +
              SCHOOL7_Fpop)/TOTALPOP_25over*100,
            LNGISO_prct = (FORLAN1_HH + FORLAN2_HH + FORLAN3_HH + FORLAN4_HH)/TOTAL_HH*100,
            POVERTY_prct = POVERTY_pop/POVERTY_SamplePop*100,
            OWNING_prct = OWNING_total/TOTAL_POP*100,
            RENTING_prct = RENTING_total/TOTAL_POP*100,
            HOME_AGE = ifelse(MED_HOUSE_YEAR < 1900, NA, 2020 - MED_HOUSE_YEAR),
            VETERAN_prct = VETERAN_total/TOTAL_POP*100)

#select and organize variables of interest
acsdat <- acsdat %>% 
          rename("TRACT_ID" = "GEOID") %>%
          select(TRACT_ID, TOTAL_POP, MEDIAN_AGE, UNDER5_prct, OVER18_prct, OVER65_prct, 
            OVER85_prct, MALE_prct, FEMALE_prct, WHITE_prct, BLACK_prct, 
            AMIND_prct, ASIAN_prct, PACISL_prct, OTHER_prct, MIXED_prct, 
            HISPLAT_prct, POC_prct, FOREIGN_BORN_prct, LNGISO_prct, 
            MED_INCOME_usd, GINI, DISABILITY_prct, SNAP_prct, POVERTY_prct, 
            AVG_HOUSEHOLD_SIZE, MED_RENT_usd, MED_HOME_VAL_usd, OWNING_prct, 
            RENTING_prct, HOME_AGE, NOHS_prct, BACHELORS_prct, UNEMPLOYMENT_prct,
            INDUSTRY_AGRICU_prct, INDUSTRY_CONSTR_prct, INDUSTRY_MANUFA_prct,
            INDUSTRY_WHOTRA_prct, INDUSTRY_RETTRA_prct, INDUSTRY_TRANSP_prct,
            INDUSTRY_INFORM_prct, INDUSTRY_FINANC_prct, INDUSTRY_PROSCI_prct,
            INDUSTRY_EDUHEA_prct, INDUSTRY_ARTENT_prct, INDUSTRY_OTHERS_prct,
            INDUSTRY_PUBADM_prct, INSURANCE_PRIVATE_prct, INSURANCE_PUBLIC_prct,
            INSURANCE_NONE_prct, INSURANCE_MEDICARE_prct, INSURANCE_MEDICAID_prct,
            INSURANCE_VA_prct, VETERAN_prct)



#Add Urban-Rural Classification

##import 2020 census tract population centroids
trxy <- read_delim(file = "./data/CenPop2020_Mean_TR.txt", delim = ",")

trxy <- trxy %>%
  mutate(TRACT_ID = as.character(paste(
    str_pad(STATEFP, width = 2, side = "left", pad = "0"),
    str_pad(COUNTYFP, width = 3, side = "left", pad = "0"),
    str_pad(TRACTCE, width = 6, side = "left", pad = "0"),
    sep = ""))) %>%
  select(TRACT_ID, LONGITUDE, LATITUDE)

trxy <- st_as_sf(trxy, coords = c("LONGITUDE","LATITUDE"), crs = 4326)

#import Census's urban area boundary shapefile
urbarea <- st_read("./data/tl_2020_us_uac20/tl_2020_us_uac20.shp")
urbarea <- st_transform(urbarea, crs = 4326)

#intersect with census tract population centroids
trurb <- st_join(trxy, urbarea, join = st_intersects) 

#create urban-rural indicator and join to data
urbrur <- st_drop_geometry(trurb) %>%
  mutate(URBAN_RURAL = ifelse(is.na(UACE20), "Rural", "Urban")) %>%
  select(TRACT_ID, URBAN_RURAL)

acsdat <- left_join(acsdat, urbrur, by = "TRACT_ID")


#export results
write.csv(acsdat, file = "./outputs/ACS5_UStracts2020.csv", row.names = F)
