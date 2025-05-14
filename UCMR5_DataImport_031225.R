#UCMR5 Data Import and Processing

#prepare workspace ----

#load packages
library(readxl)
library(sf)
library(tidyverse)

#set working directory
setwd("C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - PFAS + Cancer/UCMR_EA")

#process UCMR5 data ----

#read in most recent release of UCMR-5 data (Jan 2025)
ucmr5 <-  read_delim(file = "../data_sources/ucmr5_jan2025/UCMR5_All.txt", delim = "\t") 

#calculate max and average PFAS values at each PWSID (regardless of PFAS species)
ucmr5a <- ucmr5 %>%
        #remove lithium samples (only non-PFAS contaminant in UCMR5)
          filter(Contaminant != "lithium") %>%
        #convert results from ug/L to ng/L and change non-detects (NAs) to 0
          mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
            AnalyticalResultValue*1000)) %>%
          group_by(PWSID) %>%
          summarize(MAX_PFASngL = max(RESULT_ngL),
            AVG_PFASngL = mean(RESULT_ngL)) %>%
          mutate(PFAS_yn = ifelse(MAX_PFASngL>0,1,0))

#repeat, but calculate max/avg of individual PFAS species at PWSID 
ucmr5b <- ucmr5 %>%
          filter(Contaminant != "lithium") %>%
          mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
            AnalyticalResultValue*1000)) %>%
          group_by(PWSID, Contaminant) %>%
          summarize(max = max(RESULT_ngL)) %>%
          pivot_wider(names_from = Contaminant, 
            values_from = max, values_fill = NA)

#relabel columns where PFAS name has special characters
names(ucmr5b)[names(ucmr5b) == "11Cl-PF3OUdS"] <- "11Cl_PF3OUdS"
names(ucmr5b)[names(ucmr5b) == "4:2 FTS"] <- "4_2FTS"
names(ucmr5b)[names(ucmr5b) == "6:2 FTS"] <- "6_2FTS"
names(ucmr5b)[names(ucmr5b) == "8:2 FTS"] <- "8_2FTS"
names(ucmr5b)[names(ucmr5b) == "9Cl-PF3ONS"] <- "9Cl_PF3ONS"
names(ucmr5b)[names(ucmr5b) == "HFPO-DA"] <- "HFPO_DA"

#rerun PWSID summaries by PFAS groupings (eg, short vs. long)
#import PFAS taxonomy document
chem <- read_xlsx(path = "../data_sources/UCMR5_PFAS_Taxonomy.xlsx")

ucmr5grp <- left_join(ucmr5, chem, by = c("Contaminant" = "Abbreviation"))

#PFAS categories
ucmr5c <- ucmr5grp %>%
          filter(Contaminant != "lithium") %>%
          mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
            AnalyticalResultValue*1000)) %>%
          group_by(PWSID, PFAS_CAT) %>%
          summarize(max = max(RESULT_ngL)) %>%
          pivot_wider(names_from = PFAS_CAT, 
            names_glue = "PFASgroup_{PFAS_CAT}",
            values_from = max, values_fill = NA)

#PFAS chain length
ucmr5d <- ucmr5grp %>%
          filter(Contaminant != "lithium") %>%
          mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
            AnalyticalResultValue*1000)) %>%
          group_by(PWSID, Length) %>%
          summarize(max = max(RESULT_ngL)) %>%
          pivot_wider(names_from = Length,
            names_glue = "ChainLength_{Length}",
            values_from = max, values_fill = NA)

#calculate EPA's hazard index (HI) from UCMR5 data
#HI = (HFPO-DA/10) + (PFBS/2000) + (PFNA/10) + (PFHxS/10)
haz <-  ucmr5 %>%
      #keep only the 4 contaminants used in hazard index
        filter(Contaminant %in% c("HFPO-DA","PFBS","PFNA","PFHxS")) %>%
      #convert non-detects to 0's and results from ug/L (ppm) to ng/L (ppt)
        mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
          AnalyticalResultValue*1000)) %>%
      #find the maximum concentration for each compound by PWSID
        group_by(PWSID, Contaminant) %>%
        summarize(MAX_ngL = max(RESULT_ngL, na.rm = T)) %>%
      #adjust concentrations following hazard index formula
        mutate(MAX_PFAS_adj = case_when(
          Contaminant == "HFPO-DA" ~ MAX_ngL/10,
          Contaminant == "PFBS" ~ MAX_ngL/2000,
          Contaminant == "PFNA" ~ MAX_ngL/10,
          Contaminant == "PFHxS" ~ MAX_ngL/10)) %>%
      #sum adjusted values to get HI and return number of HI compounds sampled
        group_by(PWSID) %>%
        summarize(PFAS_HI = sum(MAX_PFAS_adj, na.rm = T), n_spp = n()) %>%
      #return NA if HI is <1 and not all HI PFAS species were sampled 
        mutate(PFAS_HI = ifelse(n_spp == 4 | PFAS_HI >= 1, PFAS_HI, NA)) %>%
        select(-n_spp)

#Maximum Contaminant Levels (MCL)

#determine if PFAS are above or below MCL for each PWSID
mcl <-  ucmr5 %>%
      #filter results to just contaminants with MCL
        filter(Contaminant %in% c("PFOA", "PFOS", "PFHxS", "HFPO-DA", "PFNA")) %>%
      #convert non-detects to 0's and results from ug/L (ppm) to ng/L (ppt)
        mutate(RESULT_ngL = ifelse(is.na(AnalyticalResultValue), 0, 
          AnalyticalResultValue*1000)) %>%
      #calculate maximum result for each compound by PWSID
        group_by(PWSID, Contaminant) %>%
        summarize(MAX_ngL = max(RESULT_ngL, na.rm = T)) %>%
      #create variable indicating if concentration is above or below MCL
        mutate(PFAS_MCL = case_when(
          Contaminant %in% c("PFOA","PFOS") & MAX_ngL >= 4 ~ 1,
          Contaminant %in% c("PFHxS","HFPO-DA","PFNA") & MAX_ngL >= 10 ~ 1,
          .default = 0)) %>%
      #indicate if PWSID has PFAS above its MCL. Also test to see if all 5 
      #species were tested and return NA if species were missing, unless one 
      #species was above MCL. 
        group_by(PWSID) %>%
        summarize(PFAS_MCL = max(PFAS_MCL), n_spp = n()) %>%
        mutate(PFAS_MCL = ifelse(n_spp == 5 | PFAS_MCL == 1, PFAS_MCL, NA)) %>%
        select(-n_spp)

#join all processed ucmr5 datasets together
ucmr5final <- left_join(ucmr5a, ucmr5b, by = "PWSID") %>%
              left_join(., ucmr5c, by = "PWSID") %>%
              left_join(., ucmr5d, by = "PWSID") %>%
              left_join(., haz, by = "PWSID") %>%
              left_join(., mcl, by = "PWSID")

#add binary variable indicating if PWSID was out of compliance with EPA rules
#either HI >= 1 or any of the 5 regulated compounds exceed their MCL.
ucmr5final <- ucmr5final %>%
              mutate(PFAS_EPA = case_when(is.na(PFAS_MCL) | is.na(PFAS_HI) ~ NA,
                PFAS_MCL == 1 | PFAS_HI >= 1 ~ 1, .default = 0))


#add PWS additional information to UCMR5 pfas data ----

#import UCMR5's 'additional data element' dataset with PWS information 
pwsinfo <- read_delim(file = "../data_sources/ucmr5_jan2025/UCMR5_AddtlDataElem.txt", delim = "\t")

#PWS-reported potential Industry Sources
#any known PFAS industry source impacting PWS
pwssrcA <-  pwsinfo %>% 
            filter(AdditionalDataElement == "PotentialPFASSources") %>%
            mutate(KnownSource = ifelse(Response == "Yes", 1, 0))  %>%
            group_by(PWSID) %>%
            summarize(KnownSource = max(KnownSource, na.rm = T))

#if so, what type of industry 
pwssrcB <-  pwsinfo %>% 
            filter(AdditionalDataElement == "PotentialPFASSourcesDetail")  %>%
            mutate(FLAG = 1) %>%
            group_by(PWSID, Response) %>%
            summarize(FLAG = max(FLAG, na.rm = T))

pwssrcC <-  pwssrcB %>% 
            pivot_wider(id_cols = PWSID, 
              names_from = Response,
              names_glue = "PFASSourceType_{Response}",
              values_from = FLAG, 
              values_fill = 0)

pwssrcD <- full_join(pwssrcA, pwssrcC, by = "PWSID")
pwssrcD[is.na(pwssrcD)] <- 0

#PFAS Treatment 
#(unclear about NMT (not modified after testing) treatment type - does it mean no treatment?)
pwstrtA <-  pwsinfo %>%
            filter(AdditionalDataElement == "PFASTreatment") %>%
            filter(Response != "NMT") %>%
            mutate(FLAG = 1) %>%
            group_by(PWSID, Response) %>%
            summarize(FLAG = max(FLAG, na.rm = T)) %>%
            pivot_wider(id_cols = PWSID, 
              names_from = Response,
              names_glue = "TreatmentType_{Response}",
              values_from = FLAG) %>%
            mutate(PFAStreatment = 1)

pwsinfoA <- full_join(pwstrtA, pwssrcD, by = "PWSID")

pwsinfoA[is.na(pwsinfoA)]<-0
#add PWS info to UCMR5 results
ucmr5final <- full_join(ucmr5final, pwsinfoA, by = "PWSID")

#MATCH PWSIDS to CENSUS TRACTS ----
  
#match PWSID-level UCMR5 results to census tracts by using shapefile of public
#water system service area boundaries and census tract population centroids

#import public water systems shapefile
pwssf <- st_read("./data/Water_System_Boundaries_012825/Final.shp")
#set coordinate system
pwssf <- st_transform(pwssf, crs = 4326)

#import 2020 census tract population centroids
trxy <- read_delim(file = "./data/CenPop2020_Mean_TR.txt", delim = ",")

trxy <- trxy %>%
        mutate(TRACT_ID = as.character(paste(
          str_pad(STATEFP, width = 2, side = "left", pad = "0"),
          str_pad(COUNTYFP, width = 3, side = "left", pad = "0"),
          str_pad(TRACTCE, width = 6, side = "left", pad = "0"),
          sep = ""))) %>%
        select(TRACT_ID, LONGITUDE, LATITUDE)

trxy <- st_as_sf(trxy, coords = c("LONGITUDE","LATITUDE"), crs = 4326)

#match census tract to PWS with "st_join" function returning which PWS each 
#tract's centroid intersects
sf_use_s2(F)

pwstr1 <- st_join(trxy, pwssf, join = st_intersects)

pwstr1 <- st_drop_geometry(pwstr1) %>%
          select(TRACT_ID, PWSID)

#export intersect results for crosswalk (used for other projects)
pwstract_xw <- pwstr1 %>% drop_na() %>% distinct()
write.csv(pwstract_xw, file = "./outputs/TractPopCentr_to_PWSID_xwalk.csv",row.names = F)


#check how many PWSIDs intersected census tract population centroids

##total PWSIDS in shapefile (44,504)
length(unique(pwssf$PWSID))
##PWSIDS which intersect census tract population centroid (10,211)
length(unique(pwstr1$PWSID))

#PWSIDs without shapefile (477)
ucmr5nosf <- ucmr5final %>% filter(!PWSID %in% pwssf$PWSID)
nrow(ucmr5nosf)

#PWSIDs with shapefile, but don't intersect census tract centroid (1,841)
ucmr5notr <- ucmr5final %>% filter(PWSID %in% pwssf$PWSID, !PWSID %in% pwstr1$PWSID) 
nrow(ucmr5notr)

#for those PWSIDS without a shapefile we can use the PWSID to ZIP code crosswalk
#provided with the UCMR5 data release. We determined this crosswalk
#had some errors, so we will reserve the crosswalk only for PWSIDs missing
#shapefiles.

#subset the data by shapefile availability
ucmr5_sf <- ucmr5final %>% filter(PWSID %in% pwssf$PWSID)
ucmr5_nosf <- ucmr5final %>% filter(!PWSID %in% pwssf$PWSID)

#download and prepare crosswalk files

#EPA's UCMR5 PWS-to-ZIP Crosswalk
pwszip <- read_delim(file = "../data_sources/ucmr5_jan2025/UCMR5_ZIPCodes.txt", delim = "\t") 
pwszip$ZIPCODE <- as.character(str_pad(pwszip$ZIPCODE, width = 5, side = "left", pad = "0"))

#ZIP to ZCTA crosswalk
zipzcta <- read.csv(file = "C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - CRB Exposure/Census Geography/zip_to_zcta_2018.csv")
zipzcta$ZCTA <- as.character(str_pad(zipzcta$ZCTA, width = 5, side = "left", pad = "0"))
zipzcta$ZIPCODE <- as.character(str_pad(zipzcta$ZIP_CODE, width = 5, side = "left", pad = "0"))
zipzcta <- zipzcta %>% select(ZIPCODE, ZCTA) 

#Census Tract to ZCTA crosswalk (census 2020)
tractzcta <- read.delim(file = "C:/Users/khill01/OneDrive - Environmental Protection Agency (EPA)/General - CRB Exposure/Census Geography/tab20_zcta520_tract20_natl.txt", sep = "|")
tractzcta$ZCTA <- as.character(str_pad(tractzcta$GEOID_ZCTA5_20, width = 5, side = "left", pad = "0"))
tractzcta$TRACT_ID <- as.character(str_pad(tractzcta$GEOID_TRACT_20, width = 11, side = "left", pad = "0"))
tractzcta <- tractzcta %>% select(ZCTA, TRACT_ID) %>% distinct() %>% na.omit()

#convert UCMR5 PWSID-ZIP Code crosswalk to PWSID-ZCTA
pwszcta <-  left_join(pwszip, zipzcta, by = "ZIPCODE") %>%
            select(PWSID, ZCTA) %>%
            distinct()

#then crosswalk to census tract
pwstract <- left_join(pwszcta, tractzcta, by = "ZCTA", relationship = "many-to-many") %>%
            select(PWSID, TRACT_ID) %>%
            distinct() %>%
            na.omit()


ucmr5_nosfnoxw <- ucmr5_nosf %>% filter(!PWSID %in% pwstract$PWSID)
length(unique(ucmr5_nosfnoxw$PWSID))

#Note, 7 PWSIDS are missing from both the shapefile and UCMR provided ZIP code
#crosswalk. CA1900679, GA0510083, MS0250034, NH0446020, NJ1106319 ,PA2450542,
#TN0001068. I manually searched these PWSIDS using EPA's envirofacts SDWIS 
#search tool and found all to be small water systems servicing businesses, schools,
#or hospitals (ie, not residential). For this reason, I'm okay with excluding 
#them for this analysis

#Summarize PFAS Data by Census Tract ----

#Join PWSID level data to census tracts with crosswalks
#shapefile crosswalk (use "full_join" to include PWSIDS not tested for PFAS)
ucmr5_sf <- full_join(ucmr5_sf, pwstract_xw, by = "PWSID")
#EPA's ZIP code crosswalk
ucmr5_nosf <- left_join(ucmr5_nosf, pwstract, by = "PWSID")
#and recombine datasets
ucmr5_all <- rbind(ucmr5_sf, ucmr5_nosf)

#some census tracts will match to multiple PWSIDs in which case take the maximum
#observed PFAS in that tract regardless of which PWSID it came from

#count PWSIDs by testing/detection status by tract
ucmrtract1 <- ucmr5_all %>%
              filter(!is.na(TRACT_ID)) %>%
              group_by(TRACT_ID) %>%
              summarize(N_PWStotal = n_distinct(PWSID, na.rm = T),
                N_PWStested = n_distinct(PWSID[!is.na(PFAS_yn)], na.rm = T),
                N_PWSdetected = n_distinct(PWSID[PFAS_yn == 1], na.rm = T))

#calculate the sum of additional PWS info variables (treated/source data) by tract
ucmrtract2 <- ucmr5_all %>%
              filter(!is.na(TRACT_ID)) %>%
              group_by(TRACT_ID) %>%
              summarize_at(vars(TreatmentType_OTH:PFASSourceType_PP), sum, na.rm = T)

#calculate maximum values of PFAS measurements by tract
ucmrtract3 <- ucmr5_all %>%
              filter(!is.na(PFAS_yn)) %>%
              group_by(TRACT_ID) %>%
              summarize_at(vars(MAX_PFASngL:PFAS_EPA), max, na.rm = T)
ucmrtract3[ucmrtract3==-Inf] <- NA

#join tract-level summaries into single data frame
final <-  left_join(ucmrtract1, ucmrtract2, by = "TRACT_ID") %>%
          left_join(ucmrtract3, by = "TRACT_ID")

#Incorporate Private Well Usage by Tract
##these data will estimate the population living within census tract who are not
##drinking from a public water system

#read in block group private well estimates
prvw <- read.csv(file = "./data/PrivateWellUse_BlockGroups_2020.csv")

#convert from block groups to tracts
prvw$BLOCKGRP_ID <- as.character(str_pad(prvw$BLOCKGRP_ID, width = 12, side = "left", pad = "0"))
prvw$TRACT_ID <- substr(prvw$BLOCKGRP_ID, start = 1, stop = 11)

prvw1 <-  prvw %>% 
          group_by(TRACT_ID) %>% 
          summarize(TotalPop = sum(Population_2020, na.rm = T),
            PrivWellPop = sum(WellPop_2020, na.rm = T)) %>%
          mutate(PrivWellPrct = ifelse(TotalPop == 0, 0, PrivWellPop/TotalPop*100)) %>%
          select(TRACT_ID, PrivWellPrct)

#join with UCMR data
final <- left_join(final, prvw1, by = "TRACT_ID")

#replace tract level estimates of private well use with NA with 0's (tracts with
#no population drinking from private wells were not included in private well use
#database
final$PrivWellPrct[is.na(final$PrivWellPrct)] <- 0


#Export Data ----
#export final dataset
write.csv(final, "./outputs/UCMR5_by_TRACT_031225.csv", row.names = F)
