#PREPARE WORKSPACE ----
#load packages
library(tidyverse)
library(arsenal)

#set working directory
setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - PFAS + Cancer/UCMR_EA")

#IMPORT DATASETS ----

##UCMR5 Tract-level Results
uc5dat <- read.csv(file = "./outputs/UCMR5_by_TRACT_031225.csv")

##Guide to states, regions, divisions, etc.
regs <- read.csv(file = "./data/StatesFIPS_CensusRegions.csv")

##Indices

###5-year American Community Survey
acs <- read.csv("./outputs/ACS5_UStracts2020wPOP.csv")

###NHGIS
nhg <- read.csv(file = "./data/nhgis0002_csv/nhgis0002_ds258_2020_tract.csv")

###CDC Places Data
cdc <- read.csv(file = "./outputs/CensusTract_CDCplaces_healthoutcomes.csv")

#PROCESS DATASETS ----

#ucmr5 data
uc5dat$TRACT_ID <- as.character(str_pad(uc5dat$TRACT_ID, width = 11, side = "left", pad = "0"))

##ACS
acs$TRACT_ID <- as.character(str_pad(acs$TRACT_ID, width = 11, side = "left", pad = "0"))
acs <-  acs %>% 
  dplyr::select(TRACT_ID, TOTAL_POP, MEDIAN_AGE, UNDER5_prct, 
                OVER65_prct, OVER85_prct, MALE_prct, FEMALE_prct, 
                WHITE_prct, BLACK_prct, AMIND_prct, ASIAN_prct, PACISL_prct,
                OTHER_prct, MIXED_prct, HISPLAT_prct, POC_prct,
                FOREIGN_BORN_prct, LNGISO_prct, MED_INCOME_usd, GINI, 
                DISABILITY_prct, SNAP_prct, POVERTY_prct, AVG_HOUSEHOLD_SIZE,
                MED_RENT_usd, MED_HOME_VAL_usd, OWNING_prct, RENTING_prct,
                HOME_AGE, NOHS_prct, BACHELORS_prct, UNEMPLOYMENT_prct, 
                INSURANCE_PRIVATE_prct, INSURANCE_PUBLIC_prct, 
                INSURANCE_NONE_prct, URBAN_RURAL)
#fix poverty percentage (from population not in poverty to population in poverty)
acs$POVERTY_prct <- 100 - acs$POVERTY_prct

#NHGIS
nhg$TRACT_ID <- as.character(str_pad(nhg$GEOCODE, width = 11, side = "left", pad = "0"))

nhg <- nhg %>% 
  rename("POP_TOTAL" = U7I001, 
         "POP_URBAN" = U7I002, 
         "POP_RURAL" = U7I003) %>%
  mutate(AREALAND_KM = AREALAND/1000000) %>%
  mutate(POP_DENS = POP_TOTAL/AREALAND_KM, 
         POP_URBANprct = POP_URBAN/POP_TOTAL*100,
         POP_RURALprct = POP_RURAL/POP_TOTAL*100) %>%
  select(TRACT_ID, AREALAND, POP_TOTAL, POP_DENS, POP_URBAN, 
         POP_URBANprct, POP_RURAL, POP_RURALprct)

nhg[is.na(nhg)] <- 0


#cdc
cdc$TRACT_ID <- as.character(str_pad(cdc$TRACT_ID, width = 11, side = "left",pad = "0"))
cdc <-  cdc %>% 
  select(TRACT_ID, ARTHRITIS_prv, BPHIGH_prv, CANCER_prv, CASTHMA_prv, 
         CHD_prv, COPD_prv, DIABETES_prv, HIGHCHOL_prv, STROKE_prv)

#combine all indices into single dataframe by tract
indxs <-  full_join(acs, nhg, by = "TRACT_ID") %>%
  full_join(., cdc, by = "TRACT_ID") 

#add region/division information
indxs$FIPS <- substr(indxs$TRACT_ID, start = 1, stop = 2)
regs$FIPS <- as.character(str_pad(regs$FIPS, width = 2, side = "left", pad = "0"))
indxs <- left_join(indxs, regs, by = "FIPS")
indxs$EPA_REGION <- paste("EPA_", str_pad(indxs$EPA, width = 2, side = "left", pad = "0"), sep = "")
indxs$RegionDivision <- paste(indxs$Region, indxs$Division, sep = "-")

#join with pfas data
final <- full_join(uc5dat, indxs, by = "TRACT_ID")

#create variable indicating if census tract was in or out of PWS
final$inPWS_yn <- ifelse(is.na(final$N_PWStotal), 0, 1)

#create UCMR status variables (Tested/Not Tested/Detected, 
#and Not Tested (outside PWS)/Not Tested (in PWS)/Detected)
final <-  final %>%
          mutate(UCMRSTATUS1 = case_when(
              is.na(PFAS_yn) == 0 ~ "Not Tested",
              PFAS_yn == 0 ~ "PFAS Not Detected",
              PFAS_yn == 1 ~ "PFAS Detected"),
            UCMRSTATUS2 = case_when(
              is.na(PFAS_yn) & inPWS_yn == 0 ~ "Not Tested Outside PWS",
              is.na(PFAS_yn) & inPWS_yn == 1 ~ "Not Tested Inside PWS",
              PFAS_yn == 0 ~ "PFAS Not Detected",
              PFAS_yn == 1 ~ "PFAS Detected"))

#remove states and territories outside contiguous US
final$FIPS <- as.character(substr(final$TRACT_ID, 1, 2))
final <- final %>% filter(!FIPS %in% c("02","15", "60","66","69","72","78"))

#BUILD TABLE 1 (Not Tested/Not Detected/Detected) ----

#set 'tableby' settings
mycontrols <- tableby.control(numeric.stats = "meansd", numeric.simplify = T, digits = 2, digits.n = NA)

t1v1 <- summary(tableby(formula = UCMRSTATUS1 ~ AMIND_prct+ASIAN_prct+BLACK_prct+WHITE_prct+PACISL_prct+OTHER_prct+MIXED_prct+FOREIGN_BORN_prct+HISPLAT_prct+AVG_HOUSEHOLD_SIZE+HOME_AGE+MED_HOME_VAL_usd+MED_RENT_usd+RENTING_prct+POP_DENS+POP_URBANprct+GINI+MED_INCOME_usd+POVERTY_prct+NOHS_prct+UNEMPLOYMENT_prct+INSURANCE_NONE_prct+ARTHRITIS_prv+BPHIGH_prv+CANCER_prv+CASTHMA_prv+CHD_prv+COPD_prv+DIABETES_prv+HIGHCHOL_prv+STROKE_prv, 
                        data = final, control = mycontrols), text = T) %>% as.data.frame()

t1ur <- summary(tableby(formula = UCMRSTATUS1 ~ URBAN_RURAL+Region, 
                        data = final, control = mycontrols), text = T) %>% as.data.frame()

write.csv(t1v1, "./outputs/NewTable1_part1v1.csv",row.names = F)
write.csv(t1ur, "./outputs/NewTable1_part2v1.csv",row.names = F)


#now summarize variables by ucmr status and population counts

acs2 <- read.csv("./outputs/ACS5_UStracts2020wPOP.csv")
acs2$TRACT_ID <- as.character(str_pad(acs2$TRACT_ID, width = 11, side = "left", pad = "0"))
acs2 <- acs2 %>% 
  dplyr::select(TRACT_ID, URBAN_RURAL, TOTAL_POP, WHITE_pop, BLACK_pop, 
                AMIND_pop, ASIAN_pop, PACISL_pop, OTHER_pop, MIXED_pop, HISPLAT_pop,
                FOREIGN_BORN_pop, RENTING_pop, NOHS_pop, UNEMPLOYMENT_pop, 
                INSURANCE_NONE_pop, POVERTY_pop)

#fix poverty population (from population not in poverty to population in poverty)
acs2$POVERTY_pop <- acs2$TOTAL_POP - acs2$POVERTY_pop

cdc2 <- read.csv(file = "./outputs/CensusTract_CDCplaces_healthoutcomes.csv")
cdc2$TRACT_ID <- as.character(str_pad(cdc2$TRACT_ID, width = 11, side = "left",pad = "0"))
cdc2 <- cdc2 %>% 
  select(TRACT_ID, ARTHRITIS_pop, BPHIGH_pop, CANCER_pop, CASTHMA_pop, 
         CHD_pop, COPD_pop, DIABETES_pop, HIGHCHOL_pop, STROKE_pop)

pfas <- uc5dat %>% 
        select(TRACT_ID, PFAS_yn, N_PWStotal)


final2 <- full_join(pfas, acs2, by = "TRACT_ID") %>%
          full_join(., cdc2, by = "TRACT_ID")

final2$FIPS <- substr(final2$TRACT_ID, 1, 2)
regs2 <- regs %>% select(FIPS, Region)
final2 <- left_join(final2, regs2, by = "FIPS")

final2 <- final2 %>% filter(!FIPS %in% c("02","15", "60","66","69","72","78")) 

final2$inPWS_yn <- ifelse(is.na(final2$N_PWStotal), 0, 1)

final2 <- final2 %>%
  mutate(UCMRSTATUS = case_when(
    is.na(PFAS_yn) ~ "NotTested",
    PFAS_yn == 0 ~ "PFASNotDetected",
    PFAS_yn == 1 ~ "PFASDetected"))

#calculate sums of populations by ucmr status
summs1a <-  final2 %>% 
  select(!c(TRACT_ID, FIPS, PFAS_yn, N_PWStotal, inPWS_yn, Region,  URBAN_RURAL)) %>%
  group_by(UCMRSTATUS) %>%
  summarize_all(.funs = sum, na.rm = T) 

#reshape data from wide to long - and then from long to wide resulting in a 
#column of variables, and seperate columns for each ucmr status filled with
#population totals
summs1b <-  summs1a %>%
  ungroup() %>%
  pivot_longer(!UCMRSTATUS, names_to = "VAR", values_to = "POP") %>%
  pivot_wider(names_from = "UCMRSTATUS", values_from = "POP")

#create similar datasets but with urban rural status, region, epa regions as variables
summsUR <-  final2 %>%
  group_by(UCMRSTATUS, URBAN_RURAL) %>%
  summarize(POP = sum(TOTAL_POP, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = "UCMRSTATUS", values_from = "POP") %>%
  rename("VAR" = URBAN_RURAL)

summsCR <-  final2 %>%
  group_by(UCMRSTATUS, Region) %>%
  summarize(POP = sum(TOTAL_POP, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = "UCMRSTATUS", values_from = "POP") %>%
  rename("VAR" = Region)


#also calculate number of tracts by ucmr status
summsT <- final2 %>% 
  group_by(UCMRSTATUS) %>%
  summarize(N_TRACTS = n()) %>%
  ungroup() %>%
  pivot_longer(!UCMRSTATUS, names_to = "VAR", values_to = "POP") %>%
  pivot_wider(names_from = "UCMRSTATUS", values_from = "POP")

#join all together
summsall <- rbind(summsT, summs1b, summsUR, summsCR)

#calculate overall counts across columns
summsall <- summsall %>%
  mutate(Overall = PFASDetected+PFASNotDetected+NotTested)

write.csv(summsall, "./outputs/NewTable1_part3v1.csv",row.names = F)

#BUILD TABLE 1v2 (Not Tested Outside PWS/Not Tested Inside PWS/Not Detected/Detected) ----

t1v1 <- summary(tableby(formula = UCMRSTATUS2 ~ AMIND_prct+ASIAN_prct+BLACK_prct+WHITE_prct+PACISL_prct+OTHER_prct+MIXED_prct+FOREIGN_BORN_prct+HISPLAT_prct+AVG_HOUSEHOLD_SIZE+HOME_AGE+MED_HOME_VAL_usd+MED_RENT_usd+RENTING_prct+POP_DENS+POP_URBANprct+GINI+MED_INCOME_usd+POVERTY_prct+NOHS_prct+UNEMPLOYMENT_prct+INSURANCE_NONE_prct+ARTHRITIS_prv+BPHIGH_prv+CANCER_prv+CASTHMA_prv+CHD_prv+COPD_prv+DIABETES_prv+HIGHCHOL_prv+STROKE_prv, 
                        data = final, control = mycontrols), text = T) %>% as.data.frame()

t1ur <- summary(tableby(formula = UCMRSTATUS2 ~ URBAN_RURAL+Region, 
                        data = final, control = mycontrols), text = T) %>% as.data.frame()

write.csv(t1v1, "./outputs/NewTable1_part1v2.csv",row.names = F)
write.csv(t1ur, "./outputs/NewTable1_part2v2.csv",row.names = F)


#now summarize variables by ucmr status and population counts

acs2 <- read.csv("./outputs/ACS5_UStracts2020wPOP.csv")
acs2$TRACT_ID <- as.character(str_pad(acs2$TRACT_ID, width = 11, side = "left", pad = "0"))
acs2 <- acs2 %>% 
  dplyr::select(TRACT_ID, URBAN_RURAL, TOTAL_POP, WHITE_pop, BLACK_pop, 
                AMIND_pop, ASIAN_pop, PACISL_pop, OTHER_pop, MIXED_pop, HISPLAT_pop,
                FOREIGN_BORN_pop, RENTING_pop, NOHS_pop, UNEMPLOYMENT_pop, 
                INSURANCE_NONE_pop, POVERTY_pop)
acs2$POVERTY_pop <- acs2$TOTAL_POP - acs2$POVERTY_pop

cdc2 <- read.csv(file = "./outputs/CensusTract_CDCplaces_healthoutcomes.csv")
cdc2$TRACT_ID <- as.character(str_pad(cdc2$TRACT_ID, width = 11, side = "left",pad = "0"))
cdc2 <- cdc2 %>% 
  select(TRACT_ID, ARTHRITIS_pop, BPHIGH_pop, CANCER_pop, CASTHMA_pop, 
         CHD_pop, COPD_pop, DIABETES_pop, HIGHCHOL_pop, STROKE_pop)

pfas <- uc5dat %>% 
  select(TRACT_ID, PFAS_yn, N_PWStotal)


final2 <- full_join(pfas, acs2, by = "TRACT_ID") %>%
  full_join(., cdc2, by = "TRACT_ID")

final2$FIPS <- substr(final2$TRACT_ID, 1, 2)
regs2 <- regs %>% select(FIPS, Region)
final2 <- left_join(final2, regs2, by = "FIPS")

final2 <- final2 %>% filter(!FIPS %in% c("02","15", "60","66","69","72","78")) 

final2$inPWS_yn <- ifelse(is.na(final2$N_PWStotal), 0, 1)

final2 <- final2 %>%
  mutate(UCMRSTATUS = case_when(
    is.na(PFAS_yn) & inPWS_yn == 0 ~ "NotTestedOutsidePWS",
    is.na(PFAS_yn) & inPWS_yn == 1 ~ "NotTestedInsidePWS",
    PFAS_yn == 0 ~ "PFASNotDetected",
    PFAS_yn == 1 ~ "PFASDetected"))

#calculate sums of populations by ucmr status
summs1a <-  final2 %>% 
  select(!c(TRACT_ID, FIPS, PFAS_yn, N_PWStotal, inPWS_yn, Region,  URBAN_RURAL)) %>%
  group_by(UCMRSTATUS) %>%
  summarize_all(.funs = sum, na.rm = T) 

#reshape data from wide to long - and then from long to wide resulting in a 
#column of variables, and seperate columns for each ucmr status filled with
#population totals
summs1b <-  summs1a %>%
  ungroup() %>%
  pivot_longer(!UCMRSTATUS, names_to = "VAR", values_to = "POP") %>%
  pivot_wider(names_from = "UCMRSTATUS", values_from = "POP")

#create similar datasets but with urban rural status, region, epa regions as variables
summsUR <-  final2 %>%
  group_by(UCMRSTATUS, URBAN_RURAL) %>%
  summarize(POP = sum(TOTAL_POP, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = "UCMRSTATUS", values_from = "POP") %>%
  rename("VAR" = URBAN_RURAL)

summsCR <-  final2 %>%
  group_by(UCMRSTATUS, Region) %>%
  summarize(POP = sum(TOTAL_POP, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = "UCMRSTATUS", values_from = "POP") %>%
  rename("VAR" = Region)


#also calculate number of tracts by ucmr status
summsT <- final2 %>% 
  group_by(UCMRSTATUS) %>%
  summarize(N_TRACTS = n()) %>%
  ungroup() %>%
  pivot_longer(!UCMRSTATUS, names_to = "VAR", values_to = "POP") %>%
  pivot_wider(names_from = "UCMRSTATUS", values_from = "POP")

#join all together
summsall <- rbind(summsT, summs1b, summsUR, summsCR)

#calculate overall counts across columns
summsall <- summsall %>%
  mutate(Overall = PFASDetected+PFASNotDetected+NotTestedInsidePWS+NotTestedOutsidePWS)

write.csv(summsall, "./outputs/NewTable1_part3v2.csv",row.names = F)
