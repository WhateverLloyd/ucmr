#load packages
library(tidyverse)

#set working directory
setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - PFAS + Cancer/UCMR_EA")

#import places data
cdc <- read.csv(file = "./data/PLACES__Local_Data_for_Better_Health__Census_Tract_Data_2024_release_20250314.csv")

#fix census tract 11-digit code
cdc$TRACT_ID <- as.character(str_pad(cdc$LocationName, width = 11, side = "left", pad = "0"))

#CDC places data are published for the most recent year (2024); however the rest 
#of the project uses 2020 census data. All tract IDs are the same between 2024 and 2020
#except for connecticut which added new counties in 2022 resulting in new tract IDs.
#I downloaded a crosswalk file from (https://github.com/CT-Data-Collaborative) to fix these
#tracts

xw <- read.csv(file = "./data/2022tractcrosswalk.csv")
xw$TRACT_ID_20 <- as.character(str_pad(xw$tract_fips_2020, width = 11, side = "left", pad = "0"))
xw$TRACT_ID_22 <- as.character(str_pad(xw$Tract_fips_2022, width = 11, side = "left", pad = "0"))
xw <- xw %>% select(TRACT_ID_20, TRACT_ID_22) %>% distinct()

#replace newer connecticut tracts with older 2020 census tract IDs

cdc_ct1 <-  cdc %>% 
            filter(TRACT_ID %in% xw$TRACT_ID_22) %>%
            rename("TRACT_ID_22" = TRACT_ID)

cdc_ct2 <- left_join(cdc_ct1, xw, by = "TRACT_ID_22")

cdc_ct3 <- cdc_ct2 %>% mutate(TRACT_ID = TRACT_ID_20) %>% select(!c(TRACT_ID_20, TRACT_ID_22))

#join these data back to the rest of the cdc data
cdc1 <- cdc %>% filter(!TRACT_ID %in% xw$TRACT_ID_22)

cdc2 <- rbind(cdc1, cdc_ct3)

#select health outcomes of interest and reshape data from long to wide giving each
#outcome its own column and each tract its own row
cdc3 <- cdc2 %>% 
        filter(MeasureId %in% c("HIGHCHOL","BPHIGH","ARTHRITIS",
          "DIABETES","CASTHMA","CANCER","COPD","STROKE","CHD")) %>%
        select(Year, TRACT_ID, MeasureId, Data_Value) %>%
        pivot_wider(id_cols = c("Year", "TRACT_ID"), 
          names_from = "MeasureId", 
          values_from = "Data_Value", 
          names_glue = "{MeasureId}_prv",
          values_fill = NA) %>%
        select(!Year) %>%
        group_by(TRACT_ID) %>%
        summarize_all(.funs = mean, na.rm = T)

#CDC places data provide health outcome data by prevalence. convert these percentages
#to population totals (note: prevalences calculated for adult population >18)

#get tract level population estimates for 2020
pop <-  read.csv("./outputs/ACS5_UStracts2020.csv") %>% 
        select(TRACT_ID, TOTAL_POP, OVER18_prct)
pop$TRACT_ID <- as.character(str_pad(pop$TRACT_ID, width = 11, side = "left", pad = "0"))

#and join to cdc places data
cdc4 <- left_join(cdc3, pop, by = "TRACT_ID")

#get population estimates from prevalence %
##first get population over 18 (population used in prevalence calculations)
cdc4$OVER18_prct <- ifelse(is.na(cdc4$OVER18_prct),0, cdc4$OVER18_prct)
cdc4$OVER18_pop <- round((cdc4$OVER18_prct/100)*cdc4$TOTAL_POP, digits = 0)

cdc4$STROKE_pop <- round((cdc4$STROKE_prv/100)*cdc4$OVER18_pop, digits = 0)
cdc4$ARTHRITIS_pop <- round((cdc4$ARTHRITIS_prv/100)*cdc4$OVER18_pop, digits = 0)
cdc4$DIABETES_pop <- round((cdc4$DIABETES_prv/100)*cdc4$OVER18_pop, digits = 0)
cdc4$CASTHMA_pop <- round((cdc4$CASTHMA_prv/100)*cdc4$OVER18_pop, digits = 0)
cdc4$COPD_pop <- round((cdc4$COPD_prv/100)*cdc4$OVER18_pop, digits = 0)
cdc4$CANCER_pop <- round((cdc4$CANCER_prv/100)*cdc4$OVER18_pop, digits = 0)
cdc4$HIGHCHOL_pop <- round((cdc4$HIGHCHOL_prv/100)*cdc4$OVER18_pop, digits = 0)
cdc4$BPHIGH_pop <- round((cdc4$BPHIGH_prv/100)*cdc4$OVER18_pop, digits = 0)
cdc4$CHD_pop <- round((cdc4$CHD_prv/100)*cdc4$OVER18_pop, digits = 0)
#replace NaN's with NA's
cdc4[is.na(cdc4)] <- NA


#export data
write.csv(cdc4, "./outputs/CensusTract_CDCplaces_healthoutcomes.csv", row.names = F)
